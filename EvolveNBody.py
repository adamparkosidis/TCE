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


def Run(totalMass, semmiMajor, gasParticles, dmParticles, endTime= 10000 | units.yr, timeSteps = 3 ,
        savedVersionPath = "", saveAfterMinute = 15, step = 0, relax = False):
    '''

    :param totalMass: the triple summarized mass
    :param semmiMajor:
    :param gasParticles: all the gas particles in the system
    :param dmParticles: all the dark matter particles in the system
    :param endTime: when should end the evolution
    :param timeSteps
    :param step: the begining step of the simulation
    :param relax: if it is a relaxation simulation or a normal evolution
    evolutionCode.parameters.code_time_unit =  units.yr: in how many steps you want to evolve
    :return: None
    '''

    # creating the NBody system with the 3
    nbody = nbody_system.nbody_to_si(totalMass, semmiMajor)

    # evolve
    #evolutionCode = Fi(nbody)
    evolutionCode = Gadget2(nbody, number_of_workers=7)
    evolutionCode.parameters.time_max = 1000. | units.yr
    evolutionCode.parameters.time_limit_cpu = 1000000 | units.s
    evolutionCode.parameters.eps_is_h_flag = True
    timeStep = endTime / timeSteps
    evolutionCode.parameters.timestep = timeStep
    currentTime = 0.0 | units.Myr
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
    if step!= 0:
        evolutionCode.gas_particles.add_particles(read_set_from_file(savedVersionPath + "/" + adding + "_gas_{0}.hdf5".format(step),
                                                                     'amuse', close_file= True))
        evolutionCode.dm_particles.add_particle(read_set_from_file(savedVersionPath + "/" + adding + "_dm_{0}.hdf5".format(step),
                                                                   'amuse', close_file= True)[0])
        evolutionCode.parameters.begin_time = (step * timeStep)
        currentTime = step * timeStep
        print "evolving from step " , step
    else:
        for gasParticle in gasParticles:
            evolutionCode.gas_particles.add_particles(gasParticle)
        for dmParticle in dmParticles:
            evolutionCode.dm_particles.add_particle(dmParticle)

    native_plot.figure(figsize=(20, 20), dpi=60)

    dm = evolutionCode.dm_particles.copy()
    gas = evolutionCode.gas_particles.copy()

    sph_particles_plot(evolutionCode.gas_particles)
    #native_plot.show()
    native_plot.savefig(savedVersionPath + '/pics/' + adding + '_0.jpg')
    centerOfMassRadius = evolutionCode.gas_particles.center_of_mass()
    centerOfMassV = evolutionCode.gas_particles.center_of_mass_velocity()
    x =  AdaptingVectorQuantity()
    y =  AdaptingVectorQuantity()
    z =  AdaptingVectorQuantity()
    x.append(evolutionCode.particles.x)
    y.append(evolutionCode.particles.y)
    z.append(evolutionCode.particles.z)
    #if step!= 0:
    #    x, y, z = pickle.load(open(savedVersionPath+"xyz.p", 'rb'))
    currentSecond = time.time()
    timeToSave = saveAfterMinute * 60

    print "starting SPH " + adding
    while currentTime < endTime:
        step += 1
        evolutionCode.evolve_model(currentTime)
        if (time.time() - currentSecond) % timeToSave :
            if savedVersionPath != "":
                write_set_to_file(evolutionCode.gas_particles, savedVersionPath + "/" + adding + "/gas_{0}.hdf5".format(step),
                                  'amuse' , append_to_file= False)
                write_set_to_file(Particles(particles = evolutionCode.dm_particles),
                                  savedVersionPath + "/" + adding + "/dm_{0}.hdf5".format(step), 'amuse', append_to_file= False)
                pickle.dump([x,y,z], open(savedVersionPath + "/" + adding +"xyz.p", 'wb'), pickle.HIGHEST_PROTOCOL)
                print "state saved - {0}".format(savedVersionPath) + "/" + adding

        print "current time = ", evolutionCode.model_time.as_quantity_in(units.yr)
        currentTime += timeStep
        gas = evolutionCode.gas_particles.copy()
        sph_particles_plot(gas)
        native_plot.savefig(savedVersionPath + "/pics/" + adding + "_{0}".format(step))
        print "pic {0} saved".format(step)
        if relax:
            gas.add_particle(evolutionCode.dm_particles)
            evolutionCode.gas_particles.position += (centerOfMassRadius - gas.center_of_mass())
            evolutionCode.dm_particles.position += (centerOfMassRadius - gas.center_of_mass())
            relaxingVFactor = (step / timeSteps)
            evolutionCode.gas_particles.velocity = relaxingVFactor * (evolutionCode.gas_particles.velocity -
                                                                                gas.center_of_mass_velocity()) + centerOfMassV
            evolutionCode.dm_particles.velocity = relaxingVFactor * (evolutionCode.dm_particles.velocity -
                                                                               gas.center_of_mass_velocity()) + centerOfMassV

        gas = evolutionCode.gas_particles.copy()
        dm = evolutionCode.dm_particles.copy()
        x.append(evolutionCode.particles.x)
        y.append(evolutionCode.particles.y)
        z.append(evolutionCode.particles.z)
    #native_plot.show()
    x=x.value_in(units.AU)
    y=y.value_in(units.AU)
    pyplot.figure(figsize= (20, 20), dpi= 80)
    pyplot.plot(x[:, 0], y[:, 0], 'r.', ms= 20.0, )
    pyplot.plot(x[:, 1], y[:, 1], 'g.')
    pyplot.plot(x[:, 2], y[:, 2], 'b.')
    pyplot.plot(x[:, 0], y[:, 0], 'r.', ms= 10.0, )
    pyplot.plot(x[:, 1], y[:, 1], 'g.')
    pyplot.plot(x[:, 2], y[:, 2], 'b.')
    pyplot.xlim(-20, 20)
    pyplot.ylim(-20, 20)
    pyplot.xlabel('AU')
    pyplot.savefig('savings/pics/{0}_dynamics.jpg'.format(adding))
    evolutionCode.stop()
    return gas, dm

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
    pyplot.plot(x[:, 0], y[:, 0], 'r.', ms= 10.0, )
    pyplot.plot(x[:, 1], y[:, 1], 'g.')
    pyplot.xlim(-1, 1)
    pyplot.ylim(-1, 1)
    pyplot.xlabel('AU')
    pyplot.savefig(orbitPlotPath)
    evolutionCode.stop()

