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
#import TCEPlotting

def DynamicsForBinarySystem(dynamicsCode, semmiMajor, binary):

    unitConverter = nbody_system.nbody_to_si(binary.total_mass(), semmiMajor)
    system = dynamicsCode(unitConverter, redirection="file", redirect_file="dynamics_code_out.log",info_file="info2.txt")
    system.parameters.epsilon_squared = 0 | units.m**2
    system.parameters.inttype_parameter = system.inttypes.SHARED10
    system.parameters.timestep_parameter = 0.2
    system.particles.add_particles(binary)
    return system

def HydroSystem(sphCode, envelope, core, t_end, n_steps, beginTime, core_radius, numberOfWorkers = 1):
    if sphCode.__name__ =="Gadget2":
        unitConverter = nbody_system.nbody_to_si(envelope.total_mass() + core.mass, core_radius*1000*2)
    else:
        unitConverter = nbody_system.nbody_to_si(envelope.total_mass() + core.mass, core_radius*1000)
    system = sphCode(unitConverter, redirection="file", redirect_file="sph_code_out{0}.log"
                     .format(str(time.localtime().tm_year) + "-" +
                            str(time.localtime().tm_mon) + "-" + str(time.localtime().tm_mday) + "-" +
                            str(time.localtime().tm_hour) + ":" + str(time.localtime().tm_min) + ":" +
                            str(time.localtime().tm_sec)), number_of_workers=numberOfWorkers)
    if sphCode.__name__ == "Fi":
        system.parameters.timestep = t_end / n_steps
        system.parameters.eps_is_h_flag = True
    system.parameters.begin_time = beginTime
    system.parameters.time_limit_cpu = 7200000000 | units.s
    print "core radius:",core.radius.as_string_in(units.RSun), core.radius
    system.dm_particles.add_particle(core)
    system.gas_particles.add_particles(envelope)
    system.parameters.timestep_accuracy_parameter = 0.05
    system.parameters.time_max = t_end * 1.5
    system.parameters.cpu_file = "cpu_code_out_{0}.txt".format(str(time.localtime().tm_year) + "-" +
                            str(time.localtime().tm_mon) + "-" + str(time.localtime().tm_mday) + "-" +
                            str(time.localtime().tm_hour) + ":" + str(time.localtime().tm_min) + ":" +
                            str(time.localtime().tm_sec))
    system.parameters.energy_file = "energy_out_{0}.txt".format(str(time.localtime().tm_year) + "-" +
                            str(time.localtime().tm_mon) + "-" + str(time.localtime().tm_mday) + "-" +
                            str(time.localtime().tm_hour) + ":" + str(time.localtime().tm_min) + ":" +
                            str(time.localtime().tm_sec))
    system.parameters.info_file = "info_out_{0}.txt".format(str(time.localtime().tm_year) + "-" +
                            str(time.localtime().tm_mon) + "-" + str(time.localtime().tm_mday) + "-" +
                            str(time.localtime().tm_hour) + ":" + str(time.localtime().tm_min) + ":" +
                            str(time.localtime().tm_sec))
    print system.dm_particles
    print system.parameters.epsilon_squared
    print system.parameters.gas_epsilon
    print system.parameters.timestep
    print "info file: ", system.parameters.info_file, "output: ", system.parameters.gadget_output_directory
    return system

def CoupledSystem(hydroSystem, binarySystem, t_end, n_steps, beginTime, relax = False):
    unitConverter = nbody_system.nbody_to_si(binarySystem.particles.total_mass(), t_end)
    kickerCode = MI6(unitConverter,number_of_workers= 8, redirection='file', redirect_file='kicker_code_mi6_out.log')
    epsilonSquared = (hydroSystem.dm_particles.radius[0]/ 2.8)**2
    kickerCode.parameters.epsilon_squared = epsilonSquared
    kickFromBinary = CalculateFieldForCodesUsingReinitialize(kickerCode, (binarySystem,))
    coupledSystem = Bridge(timestep=(t_end / (2 * n_steps)), verbose=False, use_threading= not relax)
    if not relax:
        kick_from_hydro = CalculateFieldForParticles(particles=hydroSystem.particles, gravity_constant=constants.G)
        kick_from_hydro.smoothing_length_squared = epsilonSquared
        coupledSystem.add_system(binarySystem, (kick_from_hydro,), False)
    coupledSystem.add_system(hydroSystem, (kickFromBinary,), False)
    return coupledSystem



def Run(totalMass, semmiMajor, sphEnvelope, sphCore, stars, endTime= 10000 | units.yr, timeSteps = 3 ,
        savedVersionPath = "", saveAfterMinute = 1, step = -1, relax = False, sphCode = Gadget2, dynamicsCode = Huayno,
         numberOfWorkers = 1, takeCompanionInRelaxation = True, system=None, dmToSave=None, gasToSave=None):
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

    if step!= -1:
        currentTime = step * timeStep

    if system is None:
        if step == -1:# the radius of the core is now 10 times the real one becuase of epsilon = 10r_c
            hydroSystem = HydroSystem(sphCode, sphEnvelope, sphCore, endTime, timeSteps, currentTime, sphCore.radius, numberOfWorkers)
        else:# if its not the first step we shouldn't multiply by 20 again...
            hydroSystem = HydroSystem(sphCode, sphEnvelope, sphCore, endTime, timeSteps, currentTime, sphCore.radius/(10*2), numberOfWorkers)
        if not relax or takeCompanionInRelaxation:
            print "\nSetting up {0} to simulate triple system".format(dynamicsCode.__name__)
            binarySystem = DynamicsForBinarySystem(dynamicsCode, semmiMajor, stars.stars)

            print "\nSetting up Bridge to simulate triple system"
            coupledSystem = CoupledSystem(hydroSystem, binarySystem, endTime, timeSteps, currentTime, relax=relax)
        else:
            coupledSystem = hydroSystem
    else: # got it from the outside
        coupledSystem = system
    if dmToSave != None:
        dm = dmToSave
        gas = gasToSave
        systemParticles = Particles(particles=dm)
        systemParticles.add_particle(gas)
    else:
        dm = coupledSystem.dm_particles.copy()
        gas = coupledSystem.gas_particles.copy()
        systemParticles = coupledSystem.particles

    centerOfMassRadius = systemParticles.center_of_mass()
    centerOfMassV = systemParticles.center_of_mass_velocity()

    #if not relax:
    #    sinks = new_sink_particles(coupledSystem.codes[0].particles, sink_radius= stars.radius[0]*2) #sink radius is the particle radius * 2

    currentSecond = time.time()

    print "starting SPH " + adding
    print "evolving from step ", step + 1
    if step ==-1:
        StarModels.SaveGas(savedVersionPath + "/" + adding + "/gas_00.amuse", gasToSave)
        StarModels.SaveDm(savedVersionPath + "/" + adding + "/dm_00.amuse".format(step+1), dmToSave)
        print "pre state saved - {0}".format(savedVersionPath) + "/" + adding
    
    while currentTime < endTime:
        step += 1
        particles = systemParticles
        if relax:
            particles.position += (centerOfMassRadius - particles.center_of_mass())
            relaxingVFactor = (step * 1.0 / timeSteps)
            particles.velocity = relaxingVFactor * (particles.velocity - particles.center_of_mass_velocity()) + centerOfMassV
        #else:
        #    sinks.accrete(coupledSystem.gas_particles)
        coupledSystem.evolve_model(currentTime)
        print "   Evolved to:", currentTime.as_quantity_in(units.day)

        currentTime += timeStep
        if (time.time() - currentSecond) > saveAfterMinute * 60:
            if savedVersionPath != "":
                StarModels.SaveGas(savedVersionPath + "/" + adding + "/gas_{0}.amuse".format(step), systemParticles.gas_particles)
                StarModels.SaveDm(savedVersionPath + "/" + adding + "/dm_{0}.amuse".format(step), systemParticles.dm_particles)
                print "state saved - {0}".format(savedVersionPath) + "/" + adding
                print systemParticles.dm_particles
                print len(systemParticles.gas_particles)
                currentSecond = time.time()
        dm = systemParticles.dm_particles.copy()
        gas = systemParticles.gas_particles.copy()
        #if not relax:
        #    print "masses: ", sinks.mass.as_quantity_in(units.MSun)
    coupledSystem.stop()

    return gas, dm

def EvolveBinary(totalMass, semmiMajor, sphEnvelope, sphCore, stars, endTime= 10000 | units.yr, timeSteps = 3 ,
        savedVersionPath = "", saveAfterMinute = 0, step = -1, relax = False, sphCode = Gadget2, dynamicsCode = Huayno,
         numberOfWorkers = 1, takeCompanionInRelaxation = True, system=None):
    '''

    Args:
        totalMass: total sph star mass
        semmiMajor: semmimajor of the binary
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

    if step != -1:
        currentTime = step * timeStep
        
    if system is  None:
        if step == -1:
            hydroSystem = HydroSystem(sphCode, sphEnvelope, sphCore, endTime, timeSteps, currentTime, sphCore.radius, numberOfWorkers)
        else: # if its not the first step we shouldn't multiply by 20 again...
            hydroSystem = HydroSystem(sphCode, sphEnvelope, sphCore, endTime, timeSteps, currentTime, sphCore.radius/(10*2), numberOfWorkers)
        if not relax or takeCompanionInRelaxation:
            hydroSystem.dm_particles.add_particle(stars.stars[-1])
            coupledSystem = hydroSystem
        else:
            coupledSystem = hydroSystem
    else:
        coupledSystem = system

    dm = coupledSystem.dm_particles.copy()
    gas = coupledSystem.gas_particles.copy()

    centerOfMassRadius = coupledSystem.particles.center_of_mass()
    centerOfMassV = coupledSystem.particles.center_of_mass_velocity()

    #if not relax:#sinks = new_sink_particles(coupledSystem.codes[0].particles, sink_radius= stars.radius[-1]*2)
    #    sinks = new_sink_particles(coupledSystem.dm_particles[-1:], sink_radius= stars.radius[-1]*2)

    currentSecond = time.time()

    print "starting SPH " + adding
    print "evolving from step ", step + 1

    particles = coupledSystem.particles.copy()
    if step == -1:
        StarModels.SaveGas(savedVersionPath + "/" + adding + "/gas_00.amuse",coupledSystem.gas_particles)
        StarModels.SaveDm(savedVersionPath + "/" + adding + "/dm_00.amuse",coupledSystem.dm_particles)
        print "pre - state saved"
    
    while currentTime < endTime:
        step += 1
        particles = coupledSystem.particles
        if relax:
            particles.position = particles.position + (centerOfMassRadius - particles.center_of_mass())
            relaxingVFactor = (step / timeSteps)
            particles.velocity = relaxingVFactor * (particles.velocity - particles.center_of_mass_velocity()) + centerOfMassV
        #else:
        #    sinks.accrete(coupledSystem.gas_particles)

        coupledSystem.evolve_model(currentTime)
        print "   Evolved to:", currentTime.as_quantity_in(units.day)

        currentTime += timeStep
        if (time.time() - currentSecond) > saveAfterMinute * 60:
            if savedVersionPath != "":
                StarModels.SaveGas(savedVersionPath + "/" + adding + "/gas_{0}.amuse".format(step+1), coupledSystem.gas_particles)
                StarModels.SaveDm(savedVersionPath + "/" + adding + "/dm_{0}.amuse".format(step+1), coupledSystem.dm_particles)
                print coupledSystem.dm_particles
                #pickle.dump(StarModels.SphMetaData(coupledSystem.gas_particles),open(savedVersionPath+"/metaData_{0}.p".format(step), 'wb'), pickle.HIGHEST_PROTOCOL)
                print "state saved - {0}".format(savedVersionPath) + "/" + adding
                currentSecond = time.time()
        print len(coupledSystem.gas_particles)
        dm = coupledSystem.dm_particles.copy()
        gas = coupledSystem.gas_particles.copy()
        #if not relax:
        #    print "masses: ", sinks.mass.as_quantity_in(units.MSun)
    coupledSystem.stop()
    return gas, dm
