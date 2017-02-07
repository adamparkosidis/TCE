import os
import time
import ConfigParser
import sys
import glob
from amuse.lab import *
from amuse.units import units , nbody_system
from amuse.datamodel import Particle
from amuse.units.units import *
from amuse.ext.sph_to_star import *
from amuse.plot import native_plot, sph_particles_plot
from amuse.ext.star_to_sph import pickle_stellar_model
from amuse.ext.sph_to_star import *
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
    internal_structure = CreateMesaDictionaryFromFiles(mesaPath)
    internal_structure = AddUnits(internal_structure)
    #stellarModel = derive_stellar_structure(internal_structure)
    mesa= MESA()
    mesa.initialize_code()
    #mesa.parameters.stabilize_new_stellar_model_flag = False
    
    mesaParticle =  mesa.new_particle_from_model(internal_structure, 0.0 | units.Myr)
    print mesaParticle
    if withCoreParticle:
        sphStar = convert_stellar_model_to_SPH(mesa, sphParticles,
                                               with_core_particle = withCoreParticle, target_core_mass  = coreMass ,
                                                           do_store_composition = True,base_grid_options=dict(type="fcc"))
    else:
        sphStar = convert_stellar_model_to_SPH(mesa, sphParticles,
                                                       do_store_composition = True,base_grid_options=dict(type="fcc"))
    print "Now having the sph star and the binaries, ready for relaxing"
    starEnvelope, dmStars = Relax(sphStar.gas_particles, sphStar.core_particle, endTime= sphStar.relaxationTime, timeSteps=sphStar.relaxationTimeSteps,
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
    print "core radius:", core.radius.as_string_in(units.RSun)
    system.dm_particles.add_particle(core)
    system.gas_particles.add_particles(envelope)
    return system

def derive_stellar_structure(internal_structure):
        stellar_model = Grid()
        stellar_model.dmass = internal_structure['dmass']
        stellar_model.mass = stellar_model.dmass.accumulate()
        stellar_model.rho = internal_structure['rho']
        stellar_model.radius = internal_structure['radius']
        stellar_model.temperature = internal_structure['temperature']
        stellar_model.luminosity = internal_structure['luminosity']
        setattr(stellar_model, 'X_H', internal_structure['X_H'])
        setattr(stellar_model, 'X_He', internal_structure['X_He'])
        setattr(stellar_model, 'X_C', internal_structure['X_C'])
        setattr(stellar_model, 'X_N', internal_structure['X_N'])
        setattr(stellar_model, 'X_O', internal_structure['X_O'])
        setattr(stellar_model, 'X_Ne', internal_structure['X_Ne'])
        setattr(stellar_model, 'X_Mg', internal_structure['X_Mg'])
        setattr(stellar_model, 'X_Si', internal_structure['X_Si'])
        setattr(stellar_model, 'X_Fe', numpy.zeros(len(stellar_model.dmass)))
        return stellar_model
        return stellar_model

def CreateArrayFromFile(filePath):
    file = open(filePath,"r")
    array = file.readlines()
    newArray = []
    for element in array:
        element = element[:-1]
        newArray.append(float(element))
    return newArray

    return array

def CreateMesaDictionaryFromFiles(fileDirectory):
    internal_structure = dict()
    files = os.listdir(fileDirectory)
    onlyFiles =  [file for file in files if os.path.isfile(fileDirectory+"/"+file)]
    for file in onlyFiles:
        internal_structure[str(file)] = CreateArrayFromFile(fileDirectory+"/"+file)

    return internal_structure

def ConvertUnits(listOfElements, factor):
    return [float(element) * factor for element in listOfElements]


def AddUnits(internal_structure):
    internal_structure['dmass'] = internal_structure['dmass'] | units.MSun
    #internal_structure['radius'] = internal_structure['radius'] | units.m
    internal_structure['radius'] = ConvertUnits(internal_structure['radius'], 6.957 * 10**10) | units.cm
    internal_structure['rho'] = internal_structure['rho'] | units.g/units.cm **3
    internal_structure['temperature'] = internal_structure['temperature'] | units.K
    internal_structure['X_H'] = internal_structure['X_H']
    internal_structure['luminosity'] = ConvertUnits(internal_structure['luminosity'], 3.826 * 10**33) | units.erg/ units.s
    #internal_structure['luminosity'] = internal_structure['luminosity'] | units.LSun
    internal_structure['X_He'] = internal_structure['X_He']
    internal_structure['X_C'] = internal_structure['X_C']
    internal_structure['X_N'] = internal_structure['X_N']
    internal_structure['X_O'] = internal_structure['X_O']
    internal_structure['X_Ne'] = internal_structure['X_Ne']
    internal_structure['X_Mg'] = internal_structure['X_Mg']
    internal_structure['X_Si'] = internal_structure['X_Si']
    internal_structure['X_Fe'] = [0.0 for element in internal_structure['X_Si']]

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
    Run("COConfiguration.ini", mesaPath = "/BIGDATA/yossef/WDRelaxation/CO")
