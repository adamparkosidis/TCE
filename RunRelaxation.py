import os
import time
import ConfigParser
import sys
import glob
from amuse.lab import *
from amuse.units import units , nbody_system
from amuse.datamodel import Particle
from amuse.units.units import *
from amuse.plot import native_plot, sph_particles_plot
from amuse.ext.star_to_sph import pickle_stellar_model
from amuse.datamodel import Particles
from amuse.io import read_set_from_file
from amuse.community.gadget2.interface import Gadget2
import StarModels
import EvolveNBody

def Run(configurationFile, mesaPath = "", withCoreParticle=False, coreMass = 0|units.MSun):
    '''
    creating the binary
    :return:main star's mass, the envelope particles, the core particles, the binary stars and the binary semmimajor
    '''
    parser = ConfigParser.ConfigParser()
    parser.read(configurationFile)
    sphParticles = float(parser.get("Star", "sphParticles"))

    internal_structure= CreateMesaDictionaryFromFiles(mesaPath)
    mesa=MESA()
    mesa.new_particle_from_model(internal_structure)

    if withCoreParticle:
        sphStar = convert_stellar_model_to_SPH(mesa, sphParticles,
                                               with_core_particle = withCoreParticle, target_core_mass  = coreMass ,
                                                           do_store_composition = False,base_grid_options=dict(type="fcc"))
    else:
        sphStar = convert_stellar_model_to_SPH(mesa, sphParticles,
                                                       do_store_composition = False,base_grid_options=dict(type="fcc"))
    print "Now having the sph star and the binaries, ready for relaxing"
    starEnvelope, dmStars = Relax(sphStar.gas_particles, sphStar.core_particle, endTime= sphStar.relaxationTime, timeSteps=sphStar.relaxationTimeSteps ,
        savedVersionPath = "mesaPath", saveAfterMinute = 1, step = -1, sphCode = Gadget2,
          numberOfWorkers = sphStar.numberOfWorkers)
    starCore = dmStars[-1]
    starCore.radius = sphStar.core_particle.radius
    sphMetaData = StarModels.SphMetaData(sphStar)

    #saved state
    StarModels.SaveState(mesaPath, sphStar.total_mass(), starEnvelope, dmStars,0, sphMetaData)

def HydroSystem(sphCode, envelope, core, t_end, n_steps, beginTime, core_radius, numberOfWorkers = 1):
    unitConverter = nbody_system.nbody_to_si(envelope.total_mass() + core.mass, t_end)
    system = Gadget2(unitConverter, redirection="file", redirect_file="sph_code_out.log", number_of_workers=numberOfWorkers)
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

def CreateArrayFromFile(filePath):
    file = open(filePath,"r")
    array = file.readlines()

    return array

def CreateMesaDictionaryFromFiles(fileDirectory):
    internal_structure = dict()
    files = os.listdir(fileDirectory)
    for file in files:
        if os.path.isfile(fileDirectory+"/"+file):
            internal_structure[str(file)]=CreateArrayFromFile(fileDirectory+"/"+file)

    return internal_structure

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
    Run("AGBConfiguration.ini", mesaPath = "../../../BIGDATA/yossef/WDRelaxation/AGB")
