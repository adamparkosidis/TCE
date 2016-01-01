from matplotlib import pyplot
#from matplotlib.animation as animation
import time

from amuse.units.quantities import AdaptingVectorQuantity
from amuse.community.gadget2.interface import Gadget2
from amuse.community.hermite0.interface import  Hermite
from amuse.units import units , nbody_system
from amuse.units.units import *
from amuse.plot import plot as aplot, native_plot, sph_particles_plot
from amuse.lab import Particles
from amuse.io import read_set_from_file, write_set_to_file


def Run(totalMass, semmiMajor, gasParticles, dmParticles, endTime= 10000 | units.yr, timeSteps = 3 ,
        savedVersionPath = "", saveAfterMinute = 15, step = 1):
    '''

    :param totalMass: the triple summarized mass
    :param semmiMajor:
    :param gasParticles: all the gas particles in the system
    :param dmParticles: all the dark matter particles in the system
    :param endTime: when should end the evolution
    :param timeSteps
    evolutionCode.parameters.code_time_unit =  units.yr: in how many steps you want to evolve
    :return: None
    '''
    # creating the NBody system with the 3
    nbody = nbody_system.nbody_to_si(totalMass, semmiMajor)

    # evolve
    evolutionCode = Gadget2(nbody, number_of_workers=7)
    evolutionCode.parameters.time_max = 1000. | units.yr
    #evolutionCode.parameters.timestep = 1.0 | units.yr
    timeStep = endTime / timeSteps
    currentTime = 0.0 | units.Myr
    '''
    Now check if there is a saved state
    '''
    if step!= 1:
        evolutionCode.gas_particles.add_particles(read_set_from_file(savedVersionPath+"_gas_{0}.hdf5".format(step),
                                                                     'amuse', close_file= True))
        evolutionCode.dm_particles.add_particle(read_set_from_file(savedVersionPath+"_dm_{0}.hdf5".format(step),
                                                                   'amuse', close_file= True)[0])
        evolutionCode.parameters.begin_time = ((step - 1) * timeStep)
        currentTime = (step - 1) * timeStep
    else:
        for gasParticle in gasParticles:
            evolutionCode.gas_particles.add_particles(gasParticle)
        for dmParticle in dmParticles:
            evolutionCode.dm_particles.add_particle(dmParticle)

    print "starting system evolution"
    native_plot.figure(figsize=(20, 20), dpi=100)
    parts = evolutionCode.gas_particles.copy()
    sph_particles_plot(parts)
    #native_plot.show()
    native_plot.savefig('savings/pics/evolution_0.jpg')
    x =  AdaptingVectorQuantity()
    y =  AdaptingVectorQuantity()
    z =  AdaptingVectorQuantity()
    x.append(evolutionCode.particles.x)
    y.append(evolutionCode.particles.y)
    z.append(evolutionCode.particles.z)
    currentSecond = time.time()
    timeToSave = saveAfterMinute * 60
    while currentTime < endTime:
        evolutionCode.evolve_model(currentTime)
        if (time.time() - currentSecond) % timeToSave :
            if savedVersionPath != "":
                write_set_to_file(evolutionCode.gas_particles, savedVersionPath+"_gas_{0}.hdf5".format(step), 'amuse' ,
                                  append_to_file= False)
                write_set_to_file(Particles(particles = evolutionCode.dm_particles),
                                  savedVersionPath+"_dm_{0}.hdf5".format(step), 'amuse', append_to_file= False)
                print "state saved - {0}".format(savedVersionPath)

        print "current time = ", evolutionCode.model_time.as_quantity_in(units.yr)
        currentTime += timeStep
        step += 1
        parts = evolutionCode.gas_particles.copy()
        sph_particles_plot(parts)
        native_plot.savefig(savedVersionPath + "savings/pics/evolution_{0}.jpg".format(step))
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
    pyplot.savefig('savings/pics/TCE_dynamics.jpg')
    evolutionCode.stop()

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
