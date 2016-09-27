import os
import time
from amuse.units import units , nbody_system
from amuse.units.units import *
from amuse.plot import native_plot, sph_particles_plot
from amuse.ext.star_to_sph import pickle_stellar_model
from amuse.datamodel import Particles
from amuse.io import read_set_from_file
from amuse.community.gadget2.interface import Gadget2
import StarModels
import EvolveNBody

def Run(configurationFile, mesaPath = ""):
    '''
    creating the binary
    :return:main star's mass, the envelope particles, the core particles, the binary stars and the binary semmimajor
    '''
    giant = StarModels.CreatePointStar(configurationFile, configurationSection="Star")

    sphStar = StarModels.SphStar(giant,configurationFile,configurationSection="Star",
                                savedMesaStarPath = mesaPath, takeSavedMesa=True)

    print "Now having the sph star and the binaries, ready for relaxing"
    starEnvelope, dmStars = Relax(sphStar.gas_particles, sphStar.core_particle, endTime= sphStar.relaxationTime, timeSteps=sphStar.relaxationTimeSteps ,
        savedVersionPath = "mesaPath", saveAfterMinute = 1, step = -1, sphCode = Gadget2,
          numberOfWorkers = sphStar.numberOfWorkers)
    starCore = dmStars[-1]
    starCore.radius = sphStar.core_particle.radius
    sphMetaData = StarModels.SphMetaData(sphStar)

    #saved state
    StarModels.SaveState(mesaPath, giant.mass, starEnvelope, dmStars,0, sphMetaData)

def HydroSystem(sphCode, envelope, core, t_end, n_steps, beginTime, core_radius, numberOfWorkers = 1):
    unitConverter = nbody_system.nbody_to_si(envelope.total_mass() + core.mass, t_end)
    system = Gadget2(unitConverter, redirection="file", redirect_file="sph_code_out.log", number_of_workers=numberOfWorkers)
    print "hiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii"
    system.parameters.begin_time = beginTime
    #if sphCode.__name__ =="Gadget2":
        #system.parameters.number_of_workers = numberOfWorkers
    print "core radius before: ", core.radius
    if sphCode.__name__ == "Gadget2":
        core.radius = core_radius * 2
        #core.radius = core_radius
    else:
        core.radius = core_radius
    print "core radius:",core.radius.as_string_in(units.RSun)
    system.dm_particles.add_particle(core)
    system.gas_particles.add_particles(envelope)
    return system

def Relax(sphEnvelope, sphCore, endTime= 10000 | units.yr, timeSteps = 3 ,
        savedVersionPath = "", saveAfterMinute = 1, step = -1, sphCode = Gadget2,
          numberOfWorkers = 1):

    adding = "relaxation"
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

    dm = hydroSystem.dm_particles.copy()
    gas = hydroSystem.gas_particles.copy()

    centerOfMassRadius = hydroSystem.particles.center_of_mass()
    centerOfMassV = hydroSystem.particles.center_of_mass_velocity()

    #    x, y, z = pickle.load(open(savedVersionPath+"xyz.p", 'rb'))
    currentSecond = time.time()

    print "starting SPH " + adding
    print "evolving from step ", step + 1

    while currentTime < endTime:
        step += 1
        particles = hydroSystem.particles
        particles.position += (centerOfMassRadius - particles.center_of_mass())
        relaxingVFactor = (step * 1.0 / timeSteps)
        particles.velocity = relaxingVFactor * (particles.velocity - particles.center_of_mass_velocity()) + centerOfMassV
        hydroSystem.evolve_model(currentTime)
        print "   Evolved to:", currentTime.as_quantity_in(units.day)
        print "   Energies calculated"

        currentTime += timeStep
        if (time.time() - currentSecond) > saveAfterMinute * 60:
            if savedVersionPath != "":
                StarModels.SaveGas(savedVersionPath + "/" + adding + "/gas_{0}.amuse".format(step), hydroSystem.gas_particles)
                StarModels.SaveDm(savedVersionPath + "/" + adding + "/dm_{0}.amuse".format(step), hydroSystem.dm_particles)
                print "state saved - {0}".format(savedVersionPath) + "/" + adding
                currentSecond = time.time()
        dm = hydroSystem.dm_particles.copy()
        gas = hydroSystem.gas_particles.copy()

    hydroSystem.stop()
    return gas, dm

if __name__ == "__main__":
    Run(configurationFile="configuration.ini", mesaPath = "../amuse-master/")
