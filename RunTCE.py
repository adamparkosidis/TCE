import pickle
import os

from amuse.units import units
from amuse.units.units import *
from amuse.lab import Particles
from amuse.io import read_set_from_file, write_set_to_file
from amuse.plot import native_plot, sph_particles_plot

import StarModels
import EvolveNBody



def CreateTripleSystem(configurationFile, savedPath = "", takeSavedSPH = False, takeSavedMesa = False):
    '''
    creating the TCE
    :return:main star's mass, the envelope particles, the core particles, the binary stars and the triple semmimajor
    '''
    giant = StarModels.CreatePointStar(configurationFile,configurationSection="MainStar")
    innerBinary = StarModels.Binary(configurationFile, configurationSection="InnerBinary")
    outerBinary = StarModels.Binary(configurationFile, configurationSection="OuterBinary")

    giant.position = outerBinary.semimajorAxis * (1 + outerBinary.eccentricity) * ([1, 0, 0] | units.none);
    giant.velocity = StarModels.GetRelativeVelocityAtApastron(
        giant.mass + innerBinary.stars.total_mass(),
        outerBinary.semimajorAxis, outerBinary.eccentricity) * ([0, 1, 0] | units.none)
    triple = innerBinary.stars
    giantInSet = triple.add_particle(giant)

    triple.move_to_center()

    if takeSavedSPH:
        print "using saved ph saved file - {0}".format(savedPath)
        starEnvelope = StarModels.LoadGas(savedPath +"_envelope.amuse")
        starCore = StarModels.LoadDm(savedPath + "_core.amuse")
        native_plot.figure(figsize=(30, 30), dpi=60)
        sph_particles_plot(starEnvelope)
        #native_plot.show()
    else:
        if takeSavedMesa:
            sphStar = StarModels.SphStar(giantInSet,configurationFile,configurationSection="MainStar", savedMesaStarPath = savedPath)
        else:
            sphStar = StarModels.SphStar(giantInSet,configurationFile,configurationSection="MainStar")

        starEnvelope = sphStar.sphStar.gas_particles
        starCore = sphStar.sphStar.core_particle

    starEnvelope, relaxedDM = EvolveNBody.Run(totalMass= giant.mass + innerBinary.stars[0].mass + innerBinary.stars[1].mass,
                    semmiMajor= outerBinary.semimajorAxis, sphEnvelope= [starEnvelope], sphCore=starCore,stars=innerBinary.stars,
                    endTime= sphStar.relaxationTime, timeSteps= sphStar.relaxationTimeSteps, relax=True)
    starCore = relaxedDM[0]

    # fixing positions
    starEnvelope.position += giant.position
    starCore.position += giant.position

    # fixing velocities
    starEnvelope.velocity += giant.velocity
    starCore.velocity += giant.velocity

    return giant.mass, starEnvelope, starCore, innerBinary, outerBinary.semmimajorAxis

def TakeSavedState(savedVersionPath, configurationFile):
    '''
    :param savedVersionPath: the path to where you have your saved state
    :return: the saved system
    '''
    print "using saved state file - {0}".format(savedVersionPath)
    starMass, binary, tripleSemmimajor = pickle.load(open(savedVersionPath+".p", 'rb'))
    starEnvelope = StarModels.LoadGas(savedVersionPath+"_envelope.amuse")
    starCore = StarModels.LoadDm(savedVersionPath+"_core.amuse")

    # create the binary
    newBinary = StarModels.Binary(configurationFile, "BinaryStar").stars
    newBinary.position += 1.0 | units.AU

    native_plot.figure(figsize=(30, 30), dpi=60)
    sph_particles_plot(starEnvelope)
    #native_plot.show()

    return starMass, starEnvelope, starCore, newBinary, tripleSemmimajor

def SaveState(savedVersionPath, starMass, starEnvelope, starCore, binary, tripleSemmimajor):
    '''

    :param savedVersionPath:  the path to where you want to save the state after creating the system
    :param starMass:
    :param starEnvelope: sphParticles
    :param starCore: dm particles after sph
    :param binary: binary star
    :param tripleSemmimajor: semmimajor of the triple system
    :return: None
    '''
    StarModels.SaveGas(savedVersionPath+"_envelope.amuse", starEnvelope)
    StarModels.SaveDm(savedVersionPath+"_core.amuse", starCore)
    print "state saved - {0}".format(savedVersionPath)



def Start(savedVersionPath = "savings/TCETry", takeSavedState = "False", configurationFile = "TCEConfiguration.ini"):
    '''
    This is the main function of our simulation
    :param savedVersionPath: path to the saved state
    :param takeSavedState: do you have a saved state you want to use? True or False if it is all saved ,
                            Mesa if its only the Mesa Star, Sph - if its only the sph star (after relaxation)
    :return: None
    '''
    try:
        os.makedirs(savedVersionPath + "/pics")
    except(OSError):
        pass
    # creating the triple system
    if takeSavedState == "True":
        starMass, starEnvelope, starCore, binary, tripleSemmimajor = TakeSavedState(savedVersionPath, configurationFile)
    else:
        if takeSavedState == "Sph":
            starMass, starEnvelope, starCore, binary, tripleSemmimajor = CreateTripleSystem(configurationFile, savedVersionPath, takeSavedSPH= True)
        elif takeSavedState == "Relax":
            starMass, starEnvelope, starCore, binary, tripleSemmimajor = CreateTripleSystem(configurationFile, savedVersionPath, takeSavedSPH= True)
        elif takeSavedState == "Mesa":
            starMass, starEnvelope, starCore, binary, tripleSemmimajor = CreateTripleSystem(configurationFile, savedVersionPath, takeSavedMesa= True)
        else:
            starMass, starEnvelope, starCore, binary, tripleSemmimajor = CreateTripleSystem(configurationFile, savedVersionPath)
        SaveState(savedVersionPath, starMass, starEnvelope, starCore, binary, tripleSemmimajor)


    #EvolveNBody.Run(totalMass= starMass, semmiMajor= tripleSemmimajor, gasParticles= [starEnvelope],
    #               dmParticles= [starCore], endTime= 1000. | units.yr, timeSteps= 12 ,
    #               savedVersionPath= savedVersionPath, step= 0)

    #EvolveNBody.Run(totalMass= starMass + binary[0].mass,
    #                semmiMajor= tripleSemmimajor, gasParticles= [starEnvelope], dmParticles= [starCore , binary[0]],
    #                endTime= 10. | units.yr, timeSteps= 5, savedVersionPath= savedVersionPath)


    #EvolveNBody.EvolveBinary(totalMass= binary[0].mass + binary[1].mass,
    #                semmiMajor= 0.15 | units.AU, binary= binary , endTime= 100 | units.yr, timeSteps = 2)

    # creating the NBody system with the 3 and evolving
    EvolveNBody.Run(totalMass= starMass + binary[0].mass + binary[1].mass,
                    semmiMajor= tripleSemmimajor, gasParticles= [starEnvelope], dmParticles= [starCore , binary],
                    endTime= 100.0 | units.yr, timeSteps= 20)
    print "****************** Simulation Completed ******************"
if __name__ == "__main__":
    Start()

def MakeAMovieFromSavedState(savedVersionPath= "savings/TCE500000" , steps = []):
    #TODO: do something
    print "blabla"
