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
import TCEPlotting

def DynamicsForBinarySystem(dynamicsCode, semmiMajor, binary):

    unitConverter = nbody_system.nbody_to_si(binary.total_mass(), semmiMajor)
    system = dynamicsCode(unitConverter, redirection="file", redirect_file="dynamics_code_out.log")
    system.parameters.epsilon_squared = 0 | units.m**2
    system.parameters.inttype_parameter = system.inttypes.SHARED10
    system.parameters.timestep_parameter = 0.2
    system.particles.add_particles(binary)
    return system

def HydroSystem(sphCode, envelope, core, t_end, n_steps, beginTime, core_radius, numberOfWorkers = 1):
    unitConverter = nbody_system.nbody_to_si(envelope.total_mass() + core.mass, core_radius*1000)
    system = sphCode(unitConverter, redirection="file", redirect_file="sph_code_out.log",cpu_file="cpu_code_out.txt")
    if sphCode.__name__ == "Fi":
        system.parameters.timestep = t_end / n_steps
        system.parameters.eps_is_h_flag = True
    system.parameters.begin_time = beginTime
    #if sphCode.__name__ =="Gadget2":
        #system.parameters.number_of_workers = numberOfWorkers
    system.parameters.time_limit_cpu = 7200000000 | units.s
    print "core radius before: ", core.radius
    if sphCode.__name__ == "Gadget2":
        core.radius = core_radius * 2 
        #core.radius = core_radius
    else:
        core.radius = core_radius
    print "core radius:",core.radius.as_string_in(units.RSun)
    #print "core: ", core
    system.dm_particles.add_particle(core)
    system.gas_particles.add_particles(envelope)
    #print system.parameters.gas_epsilon
    #print system.parameters.radius
    print system.dm_particles
    #print system.particles
    #print system.parameters.epsilon_squared
    #system.parameters.epsilon_squared = core_radius ** 2#TODO:check this!
    print system.parameters.epsilon_squared
    print system.parameters.gas_epsilon
    print system.parameters.timestep
    print "info file: ", system.parameters.info_file
    return system

def CoupledSystem(hydroSystem, binarySystem, t_end, n_steps, beginTime, relax = False):
    unitConverter = nbody_system.nbody_to_si(binarySystem.particles.total_mass(), t_end)
    kickerCode = MI6(unitConverter,number_of_workers= 8, redirection='file', redirect_file='kicker_code_mi6_out.log')
    epsilonSquared = hydroSystem.parameters.epsilon_squared
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
         numberOfWorkers = 1):
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
    hydroSystem = HydroSystem(sphCode, sphEnvelope, sphCore, endTime, timeSteps, currentTime, sphCore.radius, numberOfWorkers)

    print "\nSetting up {0} to simulate triple system".format(dynamicsCode.__name__)
    binarySystem = DynamicsForBinarySystem(dynamicsCode, semmiMajor, stars.stars)

    print "\nSetting up Bridge to simulate triple system"
    coupledSystem = CoupledSystem(hydroSystem, binarySystem, endTime, timeSteps, currentTime, relax=relax)

    dm = coupledSystem.dm_particles.copy()
    gas = coupledSystem.gas_particles.copy()

    centerOfMassRadius = coupledSystem.particles.center_of_mass()
    centerOfMassV = coupledSystem.particles.center_of_mass_velocity()

    if not relax:
        sinks = new_sink_particles(coupledSystem.codes[0].particles, sink_radius= stars.radius[0]*2)

    #if step!= 0:
    #    x, y, z = pickle.load(open(savedVersionPath+"xyz.p", 'rb'))
    currentSecond = time.time()

    potential_energies = hydroSystem.potential_energy.as_vector_with_length(1).as_quantity_in(units.J)
    kinetic_energies = hydroSystem.kinetic_energy.as_vector_with_length(1).as_quantity_in(units.J)
    thermal_energies = coupledSystem.thermal_energy.as_vector_with_length(1).as_quantity_in(units.J)

    print "starting SPH " + adding
    print "evolving from step ", step + 1
    StarModels.SaveGas(savedVersionPath + "/" + adding + "/gas_0.amuse", coupledSystem.gas_particles)
    StarModels.SaveDm(savedVersionPath + "/" + adding + "/dm_0.amuse".format(step+1), coupledSystem.dm_particles)
    print "pre state saved - {0}".format(savedVersionPath) + "/" + adding
    
    while currentTime < endTime:
        step += 1
        particles = coupledSystem.particles
        if relax:
            particles.position += (centerOfMassRadius - particles.center_of_mass())
            relaxingVFactor = (step * 1.0 / timeSteps)
            particles.velocity = relaxingVFactor * (particles.velocity - particles.center_of_mass_velocity()) + centerOfMassV
        else:
            sinks.accrete(coupledSystem.gas_particles)
        coupledSystem.evolve_model(currentTime)
        print "   Evolved to:", currentTime.as_quantity_in(units.day)
        potential_energies.append(coupledSystem.potential_energy)
        kinetic_energies.append(coupledSystem.kinetic_energy)
        thermal_energies.append(coupledSystem.thermal_energy)
        print "   Energies calculated"

        currentTime += timeStep
        if (time.time() - currentSecond) > saveAfterMinute * 60:
            if savedVersionPath != "":
                StarModels.SaveGas(savedVersionPath + "/" + adding + "/gas_{0}.amuse".format(step), coupledSystem.gas_particles)
                StarModels.SaveDm(savedVersionPath + "/" + adding + "/dm_{0}.amuse".format(step), coupledSystem.dm_particles)
                print "state saved - {0}".format(savedVersionPath) + "/" + adding
                #TCEPlotting.PlotDensity(hydroSystem.gas_particles,hydroSystem.dm_particles, binarySystem.particles,
                #                        step=step, plottingPath=savedVersionPath + '/pics/' + adding )
                currentSecond = time.time()
        dm = coupledSystem.dm_particles.copy()
        gas = coupledSystem.gas_particles.copy()
        if not relax:
            print "masses: ", sinks.mass.as_quantity_in(units.MSun)
    coupledSystem.stop()
    return gas, dm

def EvolveBinary(totalMass, semmiMajor, sphEnvelope, sphCore, stars, endTime= 10000 | units.yr, timeSteps = 3 ,
        savedVersionPath = "", saveAfterMinute = 0, step = -1, relax = False, sphCode = Gadget2, dynamicsCode = Huayno,
         numberOfWorkers = 1, takeCompanionInRelaxation = True):
    '''

    Args:
        totalMass:
        semmiMajor:
        sphEnvelope:
        sphCore:
        stars:the binary stars
        endTime:
        timeSteps:
        savedVersionPath:
        saveAfterMinute:
        step:
        relax:
        sphCode:
        dynamicsCode:
        numberOfWorkers:

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
    #hydroSystem = HydroSystem(sphCode, sphEnvelope, sphCore, endTime, timeSteps, currentTime, 20.1731812713 | units.RSun, numberOfWorkers)
    hydroSystem = HydroSystem(sphCode, sphEnvelope, sphCore, endTime, timeSteps, currentTime, sphCore.radius, numberOfWorkers)
    if not relax or takeCompanionInRelaxation:
        #print "\nSetting up {0} to simulate binary system".format(dynamicsCode.__name__)
        #nbody = nbody_system.nbody_to_si(stars.stars.total_mass(), endTime)
        #binarySystem = ph4(nbody)
        #binarySystem.particles.add_particle(stars.stars[-1])

        #print "\nSetting up Bridge to simulate binary system"
        #coupledSystem = CoupledSystem(hydroSystem, binarySystem, endTime, timeSteps, currentTime, relax=relax)
        if step == -1:
            hydroSystem.dm_particles.add_particle(stars.stars[-1])
        coupledSystem = hydroSystem
    else:
        coupledSystem = hydroSystem

    dm = coupledSystem.dm_particles.copy()
    gas = coupledSystem.gas_particles.copy()

    centerOfMassRadius = coupledSystem.particles.center_of_mass()
    centerOfMassV = coupledSystem.particles.center_of_mass_velocity()

    if not relax:#sinks = new_sink_particles(coupledSystem.codes[0].particles, sink_radius= stars.radius[-1]*2)
        sinks = new_sink_particles(coupledSystem.dm_particles[-1:], sink_radius= stars.radius[-1]*2)

    currentSecond = time.time()

    potential_energies = hydroSystem.potential_energy.as_vector_with_length(1).as_quantity_in(units.J)
    kinetic_energies = hydroSystem.kinetic_energy.as_vector_with_length(1).as_quantity_in(units.J)
    thermal_energies = coupledSystem.thermal_energy.as_vector_with_length(1).as_quantity_in(units.J)

    print "starting SPH " + adding
    print "evolving from step ", step + 1

    particles = coupledSystem.particles.copy()
    
    StarModels.SaveGas(savedVersionPath + "/" + adding + "/gas_00.amuse",coupledSystem.gas_particles)
    StarModels.SaveDm(savedVersionPath + "/" + adding + "/dm_00.amuse",coupledSystem.dm_particles)
    print "pre - state saved"
    
    while currentTime < endTime:
        step += 1
        particles = coupledSystem.particles
        if relax:
            particles.position = particles.position + (centerOfMassRadius - particles.center_of_mass())
            relaxingVFactor = (step / timeSteps)
            particles.velocity = relaxingVFactor * (particles.velocity - gas.center_of_mass_velocity()) + centerOfMassV
        else:
            sinks.accrete(coupledSystem.gas_particles)

        coupledSystem.evolve_model(currentTime)
        print "   Evolved to:", currentTime.as_quantity_in(units.day)
        potential_energies.append(coupledSystem.potential_energy)
        kinetic_energies.append(coupledSystem.kinetic_energy)
        thermal_energies.append(coupledSystem.thermal_energy)

        print "   Energies calculated"

        currentTime += timeStep
        if (time.time() - currentSecond) > saveAfterMinute * 60:
            if savedVersionPath != "":
                StarModels.SaveGas(savedVersionPath + "/" + adding + "/gas_{0}.amuse".format(step+1), coupledSystem.gas_particles)
                StarModels.SaveDm(savedVersionPath + "/" + adding + "/dm_{0}.amuse".format(step+1), coupledSystem.dm_particles)
                #pickle.dump(StarModels.SphMetaData(coupledSystem.gas_particles),open(savedVersionPath+"/metaData_{0}.p".format(step), 'wb'), pickle.HIGHEST_PROTOCOL)
                print "state saved - {0}".format(savedVersionPath) + "/" + adding
                currentSecond = time.time()
        print len(coupledSystem.gas_particles)
        dm = coupledSystem.dm_particles.copy()
        gas = coupledSystem.gas_particles.copy()
        if not relax:
            print "masses: ", sinks.mass.as_quantity_in(units.MSun)
    coupledSystem.stop()
    return gas, dm
