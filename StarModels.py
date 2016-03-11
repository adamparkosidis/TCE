import ConfigParser
import numpy
import pickle
import math

from amuse.lab import *
from amuse.units import *
from amuse.datamodel import Particle
from amuse.ext import orbital_elements
from amuse.plot import plot, native_plot, sph_particles_plot

import EvolveNBody

def get_relative_velocity_at_apastron(total_mass, semimajor_axis, ecc):
    return (constants.G * total_mass * ((1.0 - ecc)/(1.0 + ecc)) / semimajor_axis).sqrt()

def CreatePointStar(configurationFile="", configurationSection=""):
    star = Particle()
    print 'parsing configurations'
    parser = ConfigParser.ConfigParser()
    parser.read(configurationFile)
    star.mass = float(parser.get(configurationSection, "mass")) | units.MSun
    star.metalicity = float(parser.get(configurationSection, "metalicity"))
    star.radius = float(parser.get(configurationSection, "radius")) | units.AU
    return star

class SphStar:
    def __init__(self, pointStar, configurationFile="", configurationSection="", savedMesaStarPath = "",savedGas="", savedDm=""):
        print 'parsing configurations'
        parser = ConfigParser.ConfigParser()
        parser.read(configurationFile)
        self.pointStar = pointStar
        self.sphParticles = float(parser.get(configurationSection, "sphParticles"))
        self.coreMass = float(parser.get(configurationSection, "coreMass")) | units.MSun
        self.relaxationTime = float(parser.get(configurationSection, "relaxationTime")) | units.yr
        self.relaxationTimeSteps = float(parser.get(configurationSection, "relaxationTimeSteps"))

        # Convert the star to SPH model ###
        if savedGas != "" and savedDm != "":
            self.sphStar.gas_particles = LoadGas(savedGas)
            self.sphStar.core_particle = LoadDm(savedDm)
        else:
            if savedMesaStarPath == "":
                mesaStar = self.EvolveStarWithStellarCode(MESA)

                self.sphStar = convert_stellar_model_to_SPH(mesaStar, self.sphParticles, do_relax = False, with_core_particle=True,
                                                    target_core_mass = self.coreMass)
            else:
                self.sphStar = convert_stellar_model_to_SPH(None, self.sphParticles, pickle_file = savedMesaStarPath,
                                                       with_core_particle = True, target_core_mass  = self.coreMass ,
                                                       do_store_composition = False,base_grid_options=dict(type="fcc"))




    def  EvolveStarWithStellarCode(self, code = MESA):
        '''
        evolve with (default) MESA or other
        :return: the star after has been created with MESA
        '''
        evolutionType = code()
        print "evolving with MESA"
        mainStar = evolutionType.particles.add_particle(self.pointStar)
        print "particle added, current radius = ", mainStar.radius, "target radius = ", self.pointStar.radius
        while mainStar.radius < self.pointStar.radius:
            mainStar.evolve_one_step()
        return mainStar

def LoadGas(savedGas):
    gas = read_set_from_file(savedGas, format='amuse')
    return gas

def LoadDm(savedDm):
    dms = read_set_from_file(savedDm, format='amuse')[-1]
    return dms

def SaveGas(savingPath,gas):
     write_set_to_file(gas, savingPath, format='amuse')


def SaveDm(savingPath,dms):
    write_set_to_file(dms, savingPath, format='amuse')


class Star:

    def __init__(self, configurationFile="", configurationSection="", savedMesaStarPath = "", takeSavedMesa = False):
        print 'Creating Particle'
        self.star = Particle()
        if configurationSection == "":
            self.star.mass = 2.0 | units.MSun
            self.coreMass = 0.1 | units.MSun
            self.star.metalicity = 0.02
            self.radius = 1 | units.AU
            self.star.radius = 1 | units.AU
            self.sphParticles = 100.0
            self.envelopeRadius = 1 | units.AU
            self.relaxationTime = 10
            self.relaxationTimeSteps = 1
        else:
            print 'parsing configurations'
            parser = ConfigParser.ConfigParser()
            parser.read(configurationFile)
            self.star.mass = float(parser.get(configurationSection, "mass")) | units.MSun
            self.coreMass = float(parser.get(configurationSection, "coreMass")) | units.MSun
            self.star.metalicity = float(parser.get(configurationSection, "metalicity"))
            self.radius = float(parser.get(configurationSection, "radius")) | units.AU
            self.star.radius = float(parser.get(configurationSection, "radius")) | units.AU
            self.sphParticles = float(parser.get(configurationSection, "sphParticles"))
            self.envelopeRadius = float(parser.get(configurationSection, "envelopeRadius")) | units.AU
            self.relaxationTime = float(parser.get(configurationSection, "relaxationTime"))
            self.relaxationTimeSteps = float(parser.get(configurationSection, "relaxationTimeSteps"))
            self.semiMajorAxis = float(parser.get(configurationSection, "semiMajorAxis"))
            self.eccentricity = float(parser.get(configurationSection, "eccentricity"))
        self.star.position = self.semiMajorAxis * (1 + self.eccentricity) * ([1, 0, 0] | units.none)
        self.savedPath = savedMesaStarPath
        gas_particles, core_particles = self.GetRelaxedSphModel(takeSavedMesa)
        native_plot.figure(figsize=(60, 60), dpi=100)
        sph_particles_plot(gas_particles)
        #native_plot.show()

        self.envelope = gas_particles
        self.core = core_particles


    def CreateSphModel(self):
        '''
        :return: the sph model and the maximum radius of this star
        '''
        time = 0 | units.yr
        mesaStar = self.EvolveStarWithMesa()
        # Convert the star to SPH model ###
        sphModel = convert_stellar_model_to_SPH(mesaStar, self.sphParticles, do_relax = False, with_core_particle=True,
                                                target_core_mass = self.coreMass)
        return sphModel, mesaStar.radius

    def GetRelaxedSphModel(self,takeSavedMesa):
        '''

        :return: sph star after relaxation
        '''
        if takeSavedMesa:
            totalMass,radius, starEnvelopeRadius= pickle.load(open(self.savedPath + "/mass.p", 'rb'))
            gasParticles = read_set_from_file(self.savedPath + "/gasParticles.hdf5",'amuse', close_file= True)
            coreParticles = read_set_from_file(self.savedPath + "/coreParticles.hdf5",'amuse', close_file= True)[0]
        else:
            sphStar, radius = self.CreateSphModel()
            #save everything
            write_set_to_file(sphStar.gas_particles, self.savedPath + "/gasParticles.hdf5", 'amuse' , append_to_file= False)
            write_set_to_file(Particles(particles = [sphStar.core_particle]), self.savedPath + "/coreParticles.hdf5", 'amuse', append_to_file= False)
            pickle.dump([self.star.mass, radius, self.envelopeRadius], open(self.savedPath + "/mass.p", 'wb'), pickle.HIGHEST_PROTOCOL)
            totalMass = self.star.mass + self.coreMass
            starEnvelopeRadius = self.envelopeRadius
            gasParticles = sphStar.gas_particles
            coreParticles = sphStar.core_particle
        starVolume = 4.0*numpy.pi*(radius**3)/3.0
        starAverageDensity = self.star.mass / starVolume
        relaxationTime = 2.0 / (constants.G*starAverageDensity).sqrt() # dynamical time
        return EvolveNBody.Run(totalMass, starEnvelopeRadius, [gasParticles], [coreParticles], relaxationTime.as_quantity_in(units.yr),
                               self.relaxationTimeSteps, savedVersionPath= self.savedPath, relax= True)
        #return sphStar

def GetRelativeVelocityAtApastron(total_mass, semimajor_axis, ecc):
    return (constants.G * total_mass * ((1.0 - ecc)/(1.0 + ecc)) / semimajor_axis).sqrt()

def CalculateTotalMass(particles):
    mass = 0.0 | units.MSun
    for particle in particles:
        mass += particle.mass
    return mass

class Binary:
    '''
    creating a binary star
    :param configurationFile: where to take the binary's attributes
    :param configurationSection: in which section in the configuration file
    :param binaryMasses: what are the two masses?
    :param binarySemimajorAxis: the semmimajor
    :return: binary star model
    '''
    def __init__(self, configurationFile="", configurationSection=""):
        print "Initializing ", configurationSection

        parser = ConfigParser.ConfigParser()
        parser.read(configurationFile)
        masses = [float(parser.get(configurationSection, "mass1")) | units.MSun,
                  float(parser.get(configurationSection, "mass2")) | units.MSun]
        self.semimajorAxis = float(parser.get(configurationSection, "semmimajor")) | units.AU
        self.inclination = float(parser.get(configurationSection, "inclination"))
        self.eccentricity = float(parser.get(configurationSection, "eccentricity"))

        stars =  Particles(2)
        stars.mass = masses
        stars.position = [0.0, 0.0, 0.0] | units.AU
        stars.velocity = [0.0, 0.0, 0.0] | units.km / units.s
        stars[0].y = self.semimajorAxis
        stars[0].vx = -math.cos(self.inclination)*get_relative_velocity_at_apastron(
            stars.total_mass(), self.semimajorAxis, 0)
        stars[0].vz = math.sin(self.inclination)*get_relative_velocity_at_apastron(
            stars.total_mass(), self.semimajorAxis, 0)
        stars.move_to_center()
        self.stars = stars
