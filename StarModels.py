import ConfigParser
import numpy
import pickle
import math
import os
from amuse.lab import *
from amuse.units import *
from amuse.datamodel import Particle
from amuse.ext import orbital_elements
from amuse.plot import plot, native_plot, sph_particles_plot

import EvolveNBody


def CreatePointStar(configurationFile="", configurationSection=""):
    star = Particle()
    print 'parsing configurations'
    parser = ConfigParser.ConfigParser()
    parser.read(configurationFile)
    star.mass = float(parser.get(configurationSection, "mass")) | units.MSun
    star.metalicity = float(parser.get(configurationSection, "metalicity"))
    star.radius = float(parser.get(configurationSection, "radius")) | units.AU
    star.position = [0.0, 0.0, 0.0] | units.AU
    star.velocity = [0.0, 0.0, 0.0] | units.m/units.s
    return star
'''
def CreatePointStar(mass = 1.0 | units.MSun, radius = 1.0 | units.RSun):
    star = Particle()
    star.mass = mass
    star.radius = radius
    return star
'''

class SphStar:
    def __init__(self, pointStar, configurationFile="", configurationSection="", savedMesaStarPath = "", takeSavedMesa = False,savedGas="", savedDm=""):
        print 'parsing configurations'
        parser = ConfigParser.ConfigParser()
        parser.read(configurationFile)
        self.pointStar = pointStar
        self.sphParticles = float(parser.get(configurationSection, "sphParticles"))
        self.coreMass = float(parser.get(configurationSection, "coreMass")) | units.MSun
        self.relaxationTime = float(parser.get(configurationSection, "relaxationTime")) | units.day
        self.relaxationTimeSteps = float(parser.get(configurationSection, "relaxationTimeSteps"))
        self.evolutionTime = float(parser.get(configurationSection, "evolutionTime")) | units.day
        self.evolutionTimeSteps = float(parser.get(configurationSection, "evolutionTimeSteps"))
        self.numberOfWorkers = float(parser.get(configurationSection, "numberOfWorkers"))

        # Convert the star to SPH model ###
        if savedGas != "" and savedDm != "":
            self.gas_particles = LoadGas(savedGas)
            self.core_particle = LoadDm(savedDm)
        else:
            if takeSavedMesa:
                print "taking save state from: ", savedMesaStarPath + "/" + MESA.__name__
                self.sphStar = convert_stellar_model_to_SPH(None, self.sphParticles, pickle_file = savedMesaStarPath + "/" + MESA.__name__,
                                                       with_core_particle = True, target_core_mass  = self.coreMass ,
                                                       do_store_composition = False,base_grid_options=dict(type="fcc"))
            else:
                mesaStar = self.EvolveStarWithStellarCode(MESA, savedMesaStarPath)

                self.sphStar = convert_stellar_model_to_SPH(mesaStar, self.sphParticles, do_relax = False, with_core_particle=True,
                                                    target_core_mass = self.coreMass, base_grid_options=dict(type="fcc"))
            self.gas_particles = self.sphStar.gas_particles
            self.core_particle = self.sphStar.core_particle
            self.gas_particles.position += pointStar.position
            self.gas_particles.velocity += pointStar.velocity
            self.core_particle.position += pointStar.position
            self.core_particle.velocity += pointStar.velocity



    def  EvolveStarWithStellarCode(self, code = MESA, savingPath = ""):
        '''
        evolve with (default) MESA or other
        :return: the star after has been created with MESA
        '''
        evolutionType = code()
        print "evolving with MESA"
        mainStar = evolutionType.particles.add_particle(self.pointStar)
        print "particle added, current radius = ", mainStar.radius.as_quantity_in(units.AU), "target radius = ", self.pointStar.radius
        while mainStar.radius < self.pointStar.radius:
            mainStar.evolve_one_step()
        try:
            os.makedirs(savingPath)
        except(OSError):
            pass
        pickle_stellar_model(mainStar, savingPath + "/" + code.__name__)
        print "star saved to: ", savingPath + "/" + code.__name__ , "mass: ",mainStar.mass, "stellar type:", mainStar.stellar_type
        return mainStar

class SphMetaData:
    def __init__(self,sphStar):
        self.relaxationTime = sphStar.relaxationTime
        self.relaxationTimeSteps = sphStar.relaxationTimeSteps
        self.evolutionTime = sphStar.evolutionTime
        self.evolutionTimeSteps = sphStar.evolutionTimeSteps
        self.numberOfWorkers = sphStar.numberOfWorkers

def LoadGas(savedGas):
    gas = read_set_from_file(savedGas, format='amuse')
    return gas

def LoadDm(savedDm):
    dms = read_set_from_file(savedDm, format='amuse')
    return dms

def SaveGas(savingPath,gas):
     write_set_to_file(gas, savingPath, format='amuse')


def SaveDm(savingPath,dms):
    write_set_to_file(dms,savingPath, 'amuse')


def TakeTripleSavedState(savedVersionPath, configurationFile, step = -1 ):
    '''
    :param savedVersionPath: the path to where you have your saved state
    :return: the saved system
    '''
    print "using saved state file - {0}".format(savedVersionPath)
    if step > -1:
        starEnvelope= LoadGas(savedVersionPath + "/gas_{0}.amuse".format(step))
        load= LoadDm(savedVersionPath + "/dm_{0}.amuse".format(step))
        starCore=load[-1]
        innerBinary = Binary(particles=Particles(2, particles=[load[0], load[1]]))
    else:
        starEnvelope = LoadGas(savedVersionPath+"/envelope.amuse")
        load = LoadDm(savedVersionPath + "/dm.amuse")
        starCore=load[-1]
        innerBinary = Binary(configurationFile, configurationSection="InnerBinary")
        outerBinary = Binary(configurationFile, configurationSection="OuterBinary")
        #fix the position
        giant = Particle()
        giant.mass = starEnvelope.total_mass() + starCore.mass
        giant.position = outerBinary.semimajorAxis * (1 + outerBinary.eccentricity) * ([1, 0, 0] | units.none)
        giant.velocity = GetRelativeVelocityAtApastron(
            giant.mass + innerBinary.stars.total_mass(),
            outerBinary.semimajorAxis, outerBinary.eccentricity) * ([0, 1, 0] | units.none)
        triple = innerBinary.stars
        giantInSet = triple.add_particle(giant)

        triple.move_to_center()
        innerBinary.stars = triple - giantInSet
        tripleSemmimajor = outerBinary.semimajorAxis

    starMass = starEnvelope.total_mass() + starCore.mass
    if step > -1:
        tripleVelocityDifference = BinaryCalculations.CalculateVelocityDifference(innerBinary, giant)
        tripleSeparation = BinaryCalculations.CalculateSeparation(innerBinary, giant)
        tripleSemmimajor = BinaryCalculations.CalculateSemiMajor(tripleVelocityDifference, tripleSeparation,
                                                             starMass + innerBinary.stars.total_mass())

    sphMetaData = pickle.load(open(savedVersionPath + "/metaData.p", "rb"))

    return starMass, starEnvelope, starCore, innerBinary, tripleSemmimajor, sphMetaData

def TakeBinarySavedState(savedVersionPath, configurationFile, step = -1 ):
    '''
    :param savedVersionPath: the path to where you have your saved state
    :return: the saved system
    '''
    print "using saved state file - {0}".format(savedVersionPath)
    if step > -1:
        starEnvelope = LoadGas(savedVersionPath+"/gas_{0}.amuse".format(step))
        load = LoadDm(savedVersionPath + "/dm_{0}.amuse".format(step))
        starCore=load[-1]
        starMass = starEnvelope.total_mass() + starCore.mass
        binary = Binary(particles=Particles(2, particles=[load[0], load[1]]))
    else:
        starEnvelope = LoadGas(savedVersionPath+"/envelope.amuse")
        load = LoadDm(savedVersionPath + "/dm.amuse")
        binary = Binary(configurationFile, configurationSection="Binary")
        binary.stars.radius = binary.radius
        starCore=load[-1]
        starMass = starEnvelope.total_mass() + starCore.mass
        #changing the mass to the one after relaxation
        binary.stars[0].mass = starMass
	vx, vy, vz = starEnvelope.center_of_mass_velocity()
        starEnvelopeV = (vx, vy, vz)
        print starEnvelopeV * starEnvelope.total_mass() / starMass + (starCore.vx,starCore.vy,starCore.vz)*starCore.mass / starMass
        binary.stars[0].velocity = (starEnvelopeV * starEnvelope.total_mass() +
                          (starCore.vx, starCore.vy, starCore.vz) * starCore.mass) / starMass
        print "(giant, star): ", binary.stars    
    sphMetaData = pickle.load(open(savedVersionPath + "/metaData.p", "rb"))
    return starEnvelope, starCore, binary, binary.semimajorAxis, sphMetaData

def SaveState(savedVersionPath, starMass, starEnvelope, dms, tripleSemmimajor, sphMetaData):
    '''
    :param savedVersionPath:  the path to where you want to save the state after creating the system
    :param starMass:
    :param starEnvelope: sphParticles
    :param starCore: dm particles after sph
    :param binary: binary star
    :param tripleSemmimajor: semmimajor of the triple system
    :return: None
    '''

    try:
        os.makedirs(savedVersionPath)
    except(OSError):
        pass

    pickle.dump(sphMetaData,open(savedVersionPath+"/metaData.p", 'wb'), pickle.HIGHEST_PROTOCOL)
    SaveDm(savedVersionPath+"/dm.amuse", dms)
    SaveGas(savedVersionPath+"/envelope.amuse", starEnvelope)
    print "state saved - {0}".format(savedVersionPath)
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
            self.relaxationTime = 10| units.day
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
            self.relaxationTime = float(parser.get(configurationSection, "relaxationTime")) | units.day
            self.relaxationTimeSteps = float(parser.get(configurationSection, "relaxationTimeSteps"))
            self.numberOfWorkers = float(parser.get(configurationSection, "numberOfWorkers"))
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
    def __init__(self, configurationFile="", configurationSection="", particles=None):
        if particles is None:
            print "Initializing ", configurationSection
            parser = ConfigParser.ConfigParser()
            parser.read(configurationFile)
            masses = [float(parser.get(configurationSection, "mass1")) | units.MSun,
                      float(parser.get(configurationSection, "mass2")) | units.MSun]
            self.semimajorAxis = float(parser.get(configurationSection, "semmimajor")) | units.AU
            self.inclination = math.radians(float(parser.get(configurationSection, "inclination")))
            self.eccentricity = float(parser.get(configurationSection, "eccentricity"))
            self.radius = [float(parser.get(configurationSection, "radius1")) | units.AU,
                      float(parser.get(configurationSection, "radius2")) | units.AU]

            stars = Particles(2)
            stars.mass = masses

            stars.position = [0.0, 0.0, 0.0] | units.AU
            stars.velocity = [0.0, 0.0, 0.0] | units.km / units.s
            stars[0].y = self.semimajorAxis
            stars[0].vx = -math.cos(self.inclination)*GetRelativeVelocityAtApastron(
                stars.total_mass(), self.semimajorAxis, self.eccentricity)
            stars[0].vz = math.sin(self.inclination)*GetRelativeVelocityAtApastron(
                stars.total_mass(), self.semimajorAxis, self.eccentricity)
            stars.move_to_center()

            self.stars = stars
        else:
            self.LoadBinary(particles)
    def LoadBinary(self, particles):

        print "loading binary "
        self.stars = particles
        masses = [particles[0].mass,particles[1].mass]
        velocityDifference = BinaryCalculations.CalculateVelocityDifference(self.stars[0], self.stars[1])
        separation = BinaryCalculations.CalculateSeparation(self.stars[0], self.stars[1])
        mass = particles.total_mass()

        self.semimajorAxis = BinaryCalculations.CalculateSemiMajor(velocityDifference,separation,mass)
        self.inclination = math.radians(BinaryCalculations.CalculateInclination((0,0,0) | (units.m / units.s),(0,0,0)| units.m,
                                                                   velocityDifference, separation))
        self.eccentricity =BinaryCalculations.CalculateEccentricity(self.stars[0],self.stars[1],self.semimajorAxis)

    def CalculateEccentricity(self):
        V= self.CalculateVelocityDifference()
        R = self.CalculateSeparation()
        h = CalculateSpecificMomentum(V,R)
        hSize = CalculateVectorSize(h)
        miu = constants.G*(self.stars.total_mass())
        element2 = (R[0].value_in(units.m),R[1].value_in(units.m),R[2].value_in(units.m))/(CalculateVectorSize(R).value_in(units.m))
        vxh = VectorCross(V,h)
        element10 = (vxh[0].value_in(units.m**3*units.s**-2),vxh[1].value_in(units.m**3*units.s**-2),vxh[2].value_in(units.m**3*units.s**-2))
        element11 = miu.value_in(units.m**3*units.s**-2)
        element1 = element10/element11
        eAttitude1 = CalculateVectorSize(element1 - element2)
        eAttitude2 = (1 + (2 * CalculateSpecificEnergy(V,R,self.stars[0],self.stars[1])*(hSize**2))/(constants.G*(self.stars.total_mass()))**2)**0.5
        return eAttitude1


