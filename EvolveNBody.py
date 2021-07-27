#from matplotlib import pyplot
#from matplotlib.animation as animation
import time
import pickle
import os

from amuse.units.quantities import AdaptingVectorQuantity
from amuse.lab import *
from amuse.community.gadget2.interface import Gadget2
from amuse.community.fi.interface import Fi
from amuse.community.hermite0.interface import  Hermite
from amuse.units import units , nbody_system
from amuse.units.units import *
from amuse.lab import Particles
from amuse.io import read_set_from_file, write_set_to_file

from amuse.couple.bridge import Bridge, CalculateFieldForParticles, CalculateFieldForCodesUsingReinitialize
from amuse.community.mi6.interface import MI6
from amuse.community.huayno.interface import Huayno
from amuse.ext.sink import new_sink_particles

import StarModels
import BinaryCalculations
#import TCEPlotting

def DynamicsForBinarySystem(dynamicsCode, semmiMajor, binary, outputDirectory="/vol/sci/astro/home/glanz"):
    unitConverter = nbody_system.nbody_to_si(binary.total_mass(), semmiMajor)
    try:
        os.mkdir(outputDirectory)
    except(OSError):
        pass
    system = dynamicsCode(unitConverter, redirection="file", redirect_file=outputDirectory + "/dynamics_code_out{0}.log"
                     .format(str(time.localtime().tm_year) + "-" +
                            str(time.localtime().tm_mon) + "-" + str(time.localtime().tm_mday) + "-" +
                            str(time.localtime().tm_hour) + ":" + str(time.localtime().tm_min) + ":" +
                            str(time.localtime().tm_sec)))
    system.parameters.epsilon_squared = 0 | units.m**2
    system.parameters.inttype_parameter = system.inttypes.SHARED10
    system.parameters.timestep_parameter = 0.02
    system.particles.add_particles(binary)
    return system

def FindSmallestCell(gas):
        smallestRadius = gas.total_radius()
        for gasParticle in gas:
            if gasParticle.radius < smallestRadius:
                smallestRadius = gasParticle.radius
        return smallestRadius

def FindLowestNumberOfNeighbours(gas):
        numberOfNeighbours = len(gas)
        for gasParticle in gas:
            if gasParticle.num_neighbours < numberOfNeighbours:
                numberOfNeighbours = gasParticle.num_neighbours
        return numberOfNeighbours

def CheckMerger(pointParticles):
    for i in xrange(len(pointParticles)):
        for j in xrange(i+1,len(pointParticles)):
            if BinaryCalculations.CalculateVectorSize(
                    BinaryCalculations.CalculateSeparation(pointParticles[i],pointParticles[j])) <= \
                    max(pointParticles[i].radius, pointParticles[j].radius):
                print "merger between particle", i, " and particle ",j, "!"
                return True
    return False


def HydroSystem(sphCode, envelope, core, t_end, n_steps, beginTime, core_radius, numberOfWorkers = 1, outputDirectory=""):
    unitConverter = nbody_system.nbody_to_si(envelope.total_mass() + core.mass, core_radius*100)
    print "preparing the system with ",numberOfWorkers, " workers"
    if outputDirectory == "":
        outputDirectory = "code_output"
    try:
        os.makedirs(outputDirectory)
    except(OSError):
        pass
    system = sphCode(unitConverter, mode="adaptivegravity", redirection="file", redirect_file= outputDirectory + "/sph_code_out{0}.log"
                     .format(str(time.localtime().tm_year) + "-" +
                            str(time.localtime().tm_mon) + "-" + str(time.localtime().tm_mday) + "-" +
                            str(time.localtime().tm_hour) + ":" + str(time.localtime().tm_min) + ":" +
                            str(time.localtime().tm_sec)), number_of_workers=numberOfWorkers)
    if sphCode.__name__ == "Fi":
        system.parameters.timestep = t_end / n_steps
        system.parameters.eps_is_h_flag = True
    else:
        system.parameters.time_step = t_end / n_steps
        system.parameters.gadget_output_directory = outputDirectory
    #system.parameters.begin_time = beginTime
    system.parameters.time_limit_cpu = 7200000000 | units.s
    print "core radius:",core.radius.as_string_in(units.RSun), core.radius
    print "current timestep_accuracy= ", system.parameters.timestep_accuracy_parameter
    system.parameters.timestep_accuracy_parameter = 0.05
    print "current time max= ", system.parameters.time_max
    system.parameters.time_max = t_end * 1.5
    system.dm_particles.add_particle(core)
    print "core added to hydro"
    system.gas_particles.add_particles(envelope)
    print "envelope added to hydro"
    print system.dm_particles
    #print core
    print system.parameters.epsilon_squared
    print system.parameters.gas_epsilon
    print system.parameters.timestep
    print system.parameters.begin_time
    print system.parameters.time_max
    print "output directory: ", system.parameters.gadget_output_directory
    return system

def CoupledSystem(hydroSystem, binarySystem, t_end, n_steps, beginTime, relax = False, numberOfWorkers = 8, outputDirectory=""):
    unitConverter = nbody_system.nbody_to_si(binarySystem.particles.total_mass(), t_end)
    try:
        os.makedirs(outputDirectory)
    except(OSError):
        pass
    kickerCode = MI6(unitConverter,number_of_workers= numberOfWorkers, redirection='file', redirect_file=outputDirectory + '/kicker_code_mi6_out{0}.log'.format(str(time.localtime().tm_year) + "-" +
                            str(time.localtime().tm_mon) + "-" + str(time.localtime().tm_mday) + "-" +
                            str(time.localtime().tm_hour) + ":" + str(time.localtime().tm_min) + ":" +
                            str(time.localtime().tm_sec)))
    print "kicker code intialized"
    epsilonSquared = (hydroSystem.dm_particles.radius[0]/ 2.8)**2
    kickerCode.parameters.epsilon_squared = epsilonSquared
    print epsilonSquared
    kickFromBinary = CalculateFieldForCodesUsingReinitialize(kickerCode, (binarySystem,))
    print "creating bridge"
    coupledSystem = Bridge(timestep=(t_end / (2 * n_steps)), verbose=False, use_threading= not relax)
    if not relax:
        kick_from_hydro = CalculateFieldForParticles(particles=hydroSystem.particles, gravity_constant=constants.G)
        kick_from_hydro.smoothing_length_squared = epsilonSquared
        coupledSystem.add_system(binarySystem, (kick_from_hydro,), False)
    coupledSystem.add_system(hydroSystem, (kickFromBinary,), False)
    return coupledSystem

def PrintEnergies(coupledSystem):
    print "potential energy: ", coupledSystem.potential_energy
    print "cores potential: ", coupledSystem.dm_particles.potential_energy()
    try:
        for sys in coupledSystem.systems:
            if sys.__class__.__name__ == Huayno.__name__:
                print "potential on companions: ", coupledSystem.partners[sys]
    except:
        pass
    print "kinetic energy: ", coupledSystem.kinetic_energy
    print "thermal energy: ", coupledSystem.thermal_energy
    print "total energy: ", coupledSystem.potential_energy + coupledSystem.kinetic_energy + coupledSystem.thermal_energy

def RunSystem(system=None, endTime=10000 | units.yr, timeSteps=3,
        savedVersionPath="", saveAfterMinute=1, step=-1, relax=False,initialCOM=None,
                    initialCOMV=None):


    '''
    Now check if there is a saved state
    '''
    if relax:
        adding = "relaxation"
    else:
        adding = "evolution"

    try:
        os.makedirs(savedVersionPath + "/" + adding)
    except(OSError):
        pass

    try:
        os.makedirs(savedVersionPath + '/pics/')
    except(OSError):
        pass

    timeStep = endTime / timeSteps
    currentTime = 0.0 | units.Myr
    currentSimulationTime = currentTime

    if step != -1:
        currentTime = step * timeStep

    coupledSystem = system
    try:
        dm = coupledSystem.dm_particles.copy()
        gas = coupledSystem.gas_particles.copy()
    except(Exception):
        coupledSystem.dm_particles = coupledSystem.particles
        dm = coupledSystem.particles.copy()
        gas = None

    if initialCOM is None:
        if relax and step==-1:
            initialCOM = coupledSystem.particles.center_of_mass()
        else:
            initialCOM =  [0.0, 0.0, 0.0] | units.m
    if initialCOMV is None:
        if relax and step==-1:
            initialCOMV = coupledSystem.particles.center_of_mass_velocity()
        else:
            initialCOMV = [0.0, 0.0, 0.0] | units.m / units.s

    print "initial com: ", initialCOM
    print "initial com v: ", initialCOMV

    # if not relax:
    #    sinks = new_sink_particles(coupledSystem.codes[0].particles, sink_radius= stars.radius[0]*2) #sink radius is the particle radius * 2

    currentSecond = time.time()
    # coupledSystem.time = currentTime
    print "starting SPH " + adding
    print coupledSystem.dm_particles
    print "evolving from step ", step + 1
    #print "beggining time: ", coupledSystem.time + currentTime
    if step == -1:
        try:
            StarModels.SaveGas(savedVersionPath + "/" + adding + "/gas_00.amuse", gas)
            StarModels.SaveDm(savedVersionPath + "/" + adding + "/dm_00.amuse", dm)
        except(Exception):
            print "didnt save gas"
            StarModels.SaveDm(savedVersionPath + "/" + adding + "/dm_00.amuse", coupledSystem.particles)
        print "pre state saved - {0}".format(savedVersionPath) + "/" + adding
        PrintEnergies(coupledSystem)

    while currentTime < endTime:
        step += 1
        particles = coupledSystem.particles
        if relax:
            print "com: ", particles.center_of_mass()
            if (particles.center_of_mass() != initialCOM).all():
                particles.position += (initialCOM - particles.center_of_mass())
            print "com: ", particles.center_of_mass()
            print "com v: ", particles.center_of_mass_velocity()
            relaxingVFactor = (step * 1.0 / timeSteps)
            particles.velocity = relaxingVFactor * (
                        particles.velocity - particles.center_of_mass_velocity()) + initialCOMV
            print "com v: ", particles.center_of_mass_velocity()
        else:
            # check if there is a merger - don't continue
            if CheckMerger(coupledSystem.dm_particles):
                coupledSystem.stop()
                return coupledSystem.gas_particles, coupledSystem.dm_particles
            #    sinks.accrete(coupledSystem.gas_particles)

        coupledSystem.evolve_model(currentSimulationTime + timeStep)
        print "   Evolved to:", (currentTime + timeStep).as_quantity_in(units.day)
        # print "time step is - ", coupledSystem.get_time_step()
        currentTime += timeStep
        currentSimulationTime += timeStep
        if (time.time() - currentSecond) > saveAfterMinute * 60:
            if savedVersionPath != "":
                try:
                    StarModels.SaveGas(savedVersionPath + "/" + adding + "/gas_{0}.amuse".format(step),
                                       coupledSystem.gas_particles)
                    StarModels.SaveDm(savedVersionPath + "/" + adding + "/dm_{0}.amuse".format(step),
                                  coupledSystem.dm_particles)
                except(Exception):
                    print "didnt save gas"
                    StarModels.SaveDm(savedVersionPath + "/" + adding + "/dm_{0}.amuse".format(step),
                                  coupledSystem.particles)
                print "state saved - {0}".format(savedVersionPath) + "/" + adding
                print coupledSystem.dm_particles
                PrintEnergies(coupledSystem)
                currentSecond = time.time()
        #dm = coupledSystem.dm_particles.copy()
        #gas = coupledSystem.gas_particles.copy()

        # if not relax:
        #    print "masses: ", sinks.mass.as_quantity_in(units.MSun)
    coupledSystem.stop()

    return gas, dm


def Run(totalMass, semmiMajor, sphEnvelope, sphCore, stars, endTime= 10000 | units.yr, timeSteps = 3 ,
        savedVersionPath = "", saveAfterMinute = 1, step = -1, relax = False, sphCode = Gadget2, dynamicsCode = Huayno,
         numberOfWorkers = 1, takeCompanionInRelaxation = True, system=None, dmToSave=None, gasToSave=None,initialCOM=None,
                    initialCOMV=None):
    '''

    Args:
        totalMassrelaxation tie: the triple summarized mass
        semmiMajor:
        gasParticles: all the gas particles in the system
        dmParticles: all the dark matter particles in the system
        endTime: when should end the evolution
        timeSteps:
        savedVersionPath:
        saveAfterMinute:elocity
    starCore.velocity += giant.velocity

    return giant.mass, starEnvelope,
        step: the begining step of the simulation
        relax: if it is a relaxation simulation or a normal evolution
        hydroCode: which sph code to use (default Gadget2)elocity
    starCore.velocity += giant.velocity

    return giant.mass, starEnvelope,
        dynamicCode: which dynamic code to use (default Huayno)

    Returns:

    '''


    '''
    Now check if there is a saved state
    '''
    if relax:
        adding = "relaxation"
    else:
        adding = "evolution"


    try:
        os.makedirs(savedVersionPath + "/" + adding)
    except(OSError):
        pass

    try:
        os.makedirs(savedVersionPath + '/pics/')
    except(OSError):
        pass

    timeStep = endTime / timeSteps
    currentTime = 0.0 | units.Myr
    currentSimulationTime = currentTime

    if step!= -1:
        currentTime = step * timeStep

    if system is None:
        outputDirectory = savedVersionPath + "/codes_output_{0}".format(str(time.localtime().tm_year) + "-" +
                            str(time.localtime().tm_mon) + "-" + str(time.localtime().tm_mday) + "-" +
                            str(time.localtime().tm_hour) + ":" + str(time.localtime().tm_min) + ":" +
                            str(time.localtime().tm_sec))
        os.makedirs(outputDirectory)
        if relax and step==-1:
            coreParticleRadius = sphCore.radius * 20.0 * (250.0 * 1000.0 / len(sphEnvelope)) # will be 20 for 250K, 10 for 500K and less for better resolution
        else:
            try:
                coreParticleRadius = sphCore.epsilon
            except:
                coreParticleRadius = sphCore.radius
        hydroSystem = HydroSystem(sphCode, sphEnvelope, sphCore, endTime, timeSteps, currentTime, coreParticleRadius, numberOfWorkers, outputDirectory=outputDirectory + "/hydro")

        #hydroSystem.time = currentTime
        if not relax or takeCompanionInRelaxation:
            print "\nSetting up {0} to simulate triple system".format(dynamicsCode.__name__)
            binarySystem = DynamicsForBinarySystem(dynamicsCode, semmiMajor, stars, outputDirectory=outputDirectory + "/dynamics")
            binarySystem.time = currentTime

            print "\nSetting up Bridge to simulate triple system"
            coupledSystem = CoupledSystem(hydroSystem, binarySystem, endTime, timeSteps, currentTime, relax=relax, numberOfWorkers=numberOfWorkers, outputDirectory=outputDirectory + "/coupled")
        else:
            coupledSystem = hydroSystem
    else: # got it from the outside
        coupledSystem = system

    return RunSystem(system=coupledSystem, endTime=endTime, timeSteps=timeSteps,
        savedVersionPath=savedVersionPath, saveAfterMinute=saveAfterMinute, step=step, relax=relax,initialCOM=initialCOM,
                    initialCOMV=initialCOMV)


def EvolveBinary(totalMass, semmiMajor, sphEnvelope, sphCore, stars, endTime= 10000 | units.yr, timeSteps = 3 ,
        savedVersionPath = "", saveAfterMinute = 0, step = -1, relax = False, sphCode = Gadget2, dynamicsCode = Huayno,
         numberOfWorkers = 1, takeCompanionInRelaxation = True, system=None,initialCOM=None,
                    initialCOMV=None):
    '''

    Args:
        totalMass: total sph star mass
        semmiMajor: semmimajor of the binary, there is no use here
        sphEnvelope: gas
        sphCore: core
        stars: the binary
        endTime: end of simulation
        timeSteps:
        savedVersionPath:
        saveAfterMinute:
        step:
        relax:
        sphCode: default gadget2
        dynamicsCode: default Huayno
        numberOfWorkers: not in use
        takeCompanionInRelaxation: should take the companioin count during relaxation?
        system: if want to prepare the evolving system from the outside and not use the default one

    Returns:

    '''
    #Now check if there is a saved state

    if relax:
        adding = "relaxation"
    else:
        adding = "evolution"


    try:
        os.makedirs(savedVersionPath + "/" + adding)
    except(OSError):
        pass

    try:
        os.makedirs(savedVersionPath + '/pics/')
    except(OSError):
        pass

    timeStep = endTime / timeSteps
    currentTime = 0.0 | units.Myr
    currentSimulationTime = currentTime

    if step != -1:
        currentTime = step * timeStep
        
    if system is  None:
        outputDirectory = savedVersionPath + "/codes_output_{0}".format(str(time.localtime().tm_year) + "-" +
                            str(time.localtime().tm_mon) + "-" + str(time.localtime().tm_mday) + "-" +
                            str(time.localtime().tm_hour) + ":" + str(time.localtime().tm_min) + ":" +
                            str(time.localtime().tm_sec))
        os.makedirs(outputDirectory)
        if relax and step==-1:
            coreParticleRadius = sphCore.radius * 20.0 * (250.0 * 1000.0 / len(sphEnvelope)) # will be 10 for 250K, 5 for 500K and less for better resolution
        else:
            coreParticleRadius = sphCore.epsilon
        hydroSystem = HydroSystem(sphCode, sphEnvelope, sphCore, endTime, timeSteps, currentTime, coreParticleRadius, numberOfWorkers, outputDirectory=outputDirectory + "/hydro")
        if not relax:
            hydroSystem.dm_particles.add_particle(stars[-1])
            coupledSystem = hydroSystem
        elif takeCompanionInRelaxation:
            companionField = CalculateFieldForParticles(
                Particles(particles=[stars[-1]]), gravity_constant=constants.G)
            coupledSystem = Bridge(timestep=timeStep, verbose=False, use_threading=False)
            coupledSystem.add_system(hydroSystem, (companionField,), False, h_smooth_is_eps=False)
        else:
            coupledSystem = hydroSystem
    else:
        coupledSystem = system

    return RunSystem(system=coupledSystem, endTime=endTime, timeSteps=timeSteps,
        savedVersionPath=savedVersionPath, saveAfterMinute=saveAfterMinute, step=step, relax=relax,initialCOM=initialCOM,
                    initialCOMV=initialCOMV)

