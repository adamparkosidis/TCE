from matplotlib import pyplot
#from matplotlib.animation as animation
import time
import pickle
import os

from amuse.units.quantities import AdaptingVectorQuantity
from amuse.community.gadget2.interface import Gadget2
from amuse.community.fi.interface import Fi
from amuse.community.hermite0.interface import  Hermite
from amuse.units import units , nbody_system
from amuse.units.units import *
from amuse.plot import plot as aplot, native_plot, sph_particles_plot
from amuse.lab import Particles
from amuse.io import read_set_from_file, write_set_to_file

from amuse.couple.bridge import Bridge, CalculateFieldForCodesUsingReinitialize
from amuse.community.mi6.interface import MI6
from amuse.community.huayno.interface import Huayno
from amuse.ext.sink import new_sink_particles

import StarModels
import TCEPlotting

def DynamicsForBinarySystem(dynamicsCode, totalMass, semmiMajor, binary):

    unitConverter = nbody_system.nbody_to_si(totalMass, semmiMajor)
    system = dynamicsCode(unitConverter, redirection="file", redirect_file="dynamics_code_out.log")
    system.parameters.epsilon_squared = 0 | units.m**2
    system.parameters.inttype_parameter = system.inttypes.SHARED10
    system.parameters.timestep_parameter = 0.2
    system.particles.add_particles(binary)
    return system

def HydroSystem(sphCode, envelope, core, t_end, n_steps, beginTime, core_radius):
    envelopeMass = StarModels.CalculateTotalMass(envelope)
    unitConverter = nbody_system.nbody_to_si(envelopeMass + core.mass, t_end)
    system = sphCode(unitConverter, redirection="file", redirect_file="sph_code_out.log")
    system.parameters.timestep = t_end / n_steps
    system.parameters.eps_is_h_flag = True
    system.parameters.begin_time = beginTime
    if sphCode.__name__()=="Gasget2":
        system.parameters.number_of_workers = 7
    core.radius = core_radius * 2
    system.dm_particles.add_particle(core)
    system.gas_particles.add_particles(envelope)
    return system

def CoupledSystem(hydroSystem, binary, t_end, n_steps, beginTime):
    unitConverter = nbody_system.nbody_to_si(binary.particles.total_mass(), t_end)
    kickerCode = MI6(unitConverter,number_of_workers= 8, redirection='file', redirect_file='kicker_code_mi6_out.log')
    kickerCode.parameters.epsilon_squared = 1.0 | units.RSun**2
    kickFromBinary = CalculateFieldForCodesUsingReinitialize(kickerCode, (binary,))

    coupledSystem = Bridge(timestep=(t_end / (2 * n_steps)), verbose=False, use_threading=False)
    coupledSystem.add_system(hydroSystem, (kickFromBinary,), False)

    return coupledSystem



def Run(totalMass, semmiMajor, sphEnvelope,sphCore, stars, endTime= 10000 | units.yr, timeSteps = 3 ,
        savedVersionPath = "", saveAfterMinute = 15, step = 0, relax = False, sphCode = Gadget2, dynamicsCode = Huayno,  sinkRadius = 4.0| units.RSun):
    '''

    Args:
        totalMass: the triple summarized mass
        semmiMajor:
        gasParticles: all the gas particles in the system
        dmParticles: all the dark matter particles in the system
        endTime: when should end the evolution
        timeSteps:
        savedVersionPath:
        saveAfterMinute:
        step: the begining step of the simulation
        relax: if it is a relaxation simulation or a normal evolution
        hydroCode: which sph code to use (default Gadget2)
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

    if step!= 0:
        sphEnvelope= StarModels.LoadGas(savedVersionPath + "/" + adding + "_gas_{0}.amuse".format(step))
        sphCore= StarModels.LoadDm(savedVersionPath + "/" + adding + "_dm_{0}.amuse".format(step))
        currentTime = step * timeStep

    hydroSystem = HydroSystem(sphCode, sphEnvelope, sphCore, endTime, timeSteps, currentTime, sphCore.radius)
    print "evolving from step ", step

    native_plot.figure(figsize=(20, 20), dpi=60)

    print "\nSetting up {0} to simulate triple system".format(dynamicsCode.__name__)
    binarySystem = DynamicsForBinarySystem(dynamicsCode, stars, currentTime)
    print "\nSetting up Bridge to simulate triple system"
    coupledSystem = CoupledSystem(hydroSystem, binarySystem, endTime, timeSteps, currentTime)

    dm = hydroSystem.dm_particles.copy()
    gas = hydroSystem.gas_particles.copy()

    sph_particles_plot(hydroSystem.gas_particles)
    #native_plot.show()
    native_plot.savefig(savedVersionPath + '/pics/' + adding + '_0.jpg')
    centerOfMassRadius = coupledSystem.particles.center_of_mass()
    centerOfMassV = coupledSystem.particles.center_of_mass_velocity()

    if not relax:
        sinks = new_sink_particles(coupledSystem.codes[0].particles, sink_radius= sinkRadius)

    #if step!= 0:
    #    x, y, z = pickle.load(open(savedVersionPath+"xyz.p", 'rb'))
    currentSecond = time.time()
    timeToSave = saveAfterMinute * 60

    potential_energies = hydroSystem.potential_energy.as_vector_with_length(1).as_quantity_in(units.J)
    kinetic_energies = hydroSystem.kinetic_energy.as_vector_with_length(1).as_quantity_in(units.J)
    thermal_energies = coupledSystem.thermal_energy.as_vector_with_length(1).as_quantity_in(units.J)

    print "starting SPH " + adding

    while currentTime < endTime:
        step += 1
        if relax:
            particles.position += (centerOfMassRadius - particles.center_of_mass())
            relaxingVFactor = (step / timeSteps)
            particles.velocity = relaxingVFactor * (particles.velocity - gas.center_of_mass_velocity()) + centerOfMassV
        else:
            sinks.accrete(coupledSystem.gas_particles)

        coupledSystem.evolve_model(time)
        print "   Evolved to:", time
        potential_energies.append(coupledSystem.potential_energy)
        kinetic_energies.append(coupledSystem.kinetic_energy)
        thermal_energies.append(coupledSystem.thermal_energy)
        print "   Energies calculated"

        print "current time = ", coupledSystem.model_time.as_quantity_in(units.yr)
        currentTime += timeStep
        gas = hydroSystem.gas_particles.copy()
        sph_particles_plot(gas)
        native_plot.savefig(savedVersionPath + "/pics/" + adding + "_{0}".format(step))
        print "pic {0} saved".format(step)
        particles = coupledSystem.particles

        coupledSystem.evolve_model(currentTime)
        if (time.time() - currentSecond) % timeToSave :
            if savedVersionPath != "":
                StarModels.SaveGas(savedVersionPath + "/" + adding + "/gas_{0}.amuse".format(step), hydroSystem.gas_particles)
                StarModels.SaveDm(savedVersionPath + "/" + adding + "/dm_{0}.amuse".format(step), hydroSystem.dm_particles)
                print "state saved - {0}".format(savedVersionPath) + "/" + adding
                TCEPlotting.PlotDensity(hydroSystem.gas_particles,hydroSystem.dm_particles, binarySystem.particles)
    #native_plot.show()

    coupledSystem.stop()
    return hydroSystem.gas_particles, hydroSystem.dm_particles

def EvolveBinary(totalMass, semmiMajor, binary, endTime= 10000 | units.yr, timeSteps = 10,
                 orbitPlotPath = 'Binary_dynamics.eps'):
    '''

    :param totalMass: the triple summarized mass
    :param semmiMajor:
    :param gasParticles: all the gas particles in the system
    :param dmParticles: all the dark matter particles in the system
    :param endTime: when should end the evolution
    :param timeSteps: in how many steps you want to evolve
    :return: None
    '''

    nbody = nbody_system.nbody_to_si(totalMass, semmiMajor)
    # evolve
    evolutionCode = Hermite(nbody)
    evolutionCode.parameters.epsilon_squared = 0.0 | units.AU**2
    evolutionCode.particles.add_particles(binary)
    print "starting system evolution"
    x =  AdaptingVectorQuantity()
    y =  AdaptingVectorQuantity()
    z =  AdaptingVectorQuantity()
    x.append(evolutionCode.particles.x)
    y.append(evolutionCode.particles.y)
    z.append(evolutionCode.particles.z)
    timeStep = endTime / timeSteps
    currentTime = 0.0 | units.Myr
    while currentTime < endTime:
        evolutionCode.evolve_model(currentTime)
        print "current time = ", evolutionCode.model_time.as_quantity_in(units.yr)
        currentTime += timeStep
        x.append(evolutionCode.particles.x)
        y.append(evolutionCode.particles.y)
        z.append(evolutionCode.particles.z)
    x= x.value_in(units.AU)
    y= y.value_in(units.AU)
    pyplot.figure(figsize= (20, 20), dpi= 80)
    pyplot.plot(x[:, 0], y[:, 0], 'r.', ms= 10.0, )
    pyplot.plot(x[:, 1], y[:, 1], 'g.')
    pyplot.xlim(-1, 1)
    pyplot.ylim(-1, 1)
    pyplot.xlabel('AU')
    pyplot.savefig(orbitPlotPath)
    evolutionCode.stop()

