import ConfigParser
import numpy

from amuse.lab import *
from amuse.units import *
from amuse.datamodel import Particle
from amuse.ext import orbital_elements
from amuse.plot import plot, native_plot, sph_particles_plot

import RelaxModel

class Star:

    def __init__(self, configurationFile="", configurationSection=""):
        print 'Creating Particle'
        self.star = Particle()
        if configurationSection == "":
            self.star.mass = 2.0 | units.MSun
            self.coreMass = 0.1 | units.MSun
            self.star.metalicity = 0.02
            self.radius = 1 | units.AU
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
            self.sphParticles = float(parser.get(configurationSection, "sphParticles"))
            self.envelopeRadius = float(parser.get(configurationSection, "envelopeRadius")) | units.AU
            self.relaxationTime = float(parser.get(configurationSection, "relaxationTime"))
            self.relaxationTimeSteps = float(parser.get(configurationSection, "relaxationTimeSteps"))
        relaxedModel = self.GetRelaxedSphModel()
        print relaxedModel.gas_particles
        native_plot.figure(figsize=(60, 60), dpi=100)
        sph_particles_plot(relaxedModel.gas_particles)
        #native_plot.show()

        self.envelope = relaxedModel.gas_particles
        self.core = relaxedModel.core_particle

    def  EvolveStarWithMesa(self):
        '''
        evolve with MESA
        :return: the star after has been created with MESA
        '''
        evolutionType = MESA()
        print "evolving with MESA"
        mainStar = evolutionType.particles.add_particle(self.star)
        print "particle added, current radius = ", mainStar.radius, "target radius = ", self.radius
        while mainStar.radius < self.radius:
            mainStar.evolve_one_step()
        return mainStar

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

    def GetRelaxedSphModel(self):
        '''

        :return: sph star after relaxation
        '''
        sphStar, radius = self.CreateSphModel()
        starVolume = 4.0*numpy.pi*(radius**3)/3.0
        starAverageDensity = self.star.mass / starVolume
        relaxationTime = 2.0 / (constants.G*starAverageDensity).sqrt() # dynamical time
        sphStar.core_particle
        write_set_to_file(sphStar, "sph.hdf5", 'amuse' , append_to_file= False)
        write_set_to_file(sphStar.core_particle, "coreParticle.hdf5", 'amuse' , append_to_file= False)
        write_set_to_file(Particles(particles = sphStar.gas_particles), "gasParticles.hdf5", 'amuse', append_to_file= False)
        return RelaxModel.RelaxedModel(self.star.mass + self.coreMass, self.envelopeRadius, [sphStar.gas_particles],
                                [sphStar.core_particle], relaxationTime.as_quantity_in(units.yr), self.relaxationTimeSteps)
        #return sphStar


def CreateBinary(configurationFile="", configurationSection="", binaryMasses = [1.0 | units.MSun, 1.0 | units.MSun],
                 binarySemimajorAxis = 0.1 | units.AU):
    '''
    creating a binary stard
    :param configurationFile: where to take the binary's attributes
    :param configurationSection: in which section in the configuration file
    :param binaryMasses: what are the two masses?
    :param binarySemimajorAxis: the semmimajor
    :return: binary star model
    '''
    if configurationFile == "":
        masses = binaryMasses
        semimajor_axis = binarySemimajorAxis
    else:
        parser = ConfigParser.ConfigParser()
        parser.read(configurationFile)
        masses = [float(parser.get(configurationSection, "mass1")) | units.MSun,
                  float(parser.get(configurationSection, "mass2")) | units.MSun]
        semimajor_axis = float(parser.get(configurationSection, "semmimajor")) | units.AU
    binary = orbital_elements.new_binary_from_orbital_elements(
        masses[0],
        masses[1],
        semimajor_axis, G=constants.G)
    binary.move_to_center()
    return binary
