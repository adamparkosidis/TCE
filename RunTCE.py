import pickle
import os.path

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
    if takeSavedSPH:
        print "using saved ph saved file - {0}".format(savedPath)
        starMass , tempBinary, starEnvelopeRadius= pickle.load(open(savedPath +".p", 'rb'))
        starEnvelope = read_set_from_file(savedPath +"_envelope.hdf5",'amuse', close_file= True)
        starCore = read_set_from_file(savedPath + "_core.hdf5",'amuse', close_file= True)[0]
        native_plot.figure(figsize=(30, 30), dpi=60)
        sph_particles_plot(starEnvelope)
        #native_plot.show()
    else:
        if takeSavedMesa:
            star = StarModels.Star(configurationFile, "MainStar", savedMesaStarPath = savedPath, takeSavedMesa= True)
        else:
            star = StarModels.Star(configurationFile, "MainStar", savedMesaStarPath = savedPath, takeSavedMesa= False)
        starMass = star.star.mass
        starEnvelope = star.envelope
        starEnvelopeRadius = star.envelopeRadius
        starCore = star.core
    # create the binary
    binary = StarModels.CreateBinary(configurationFile, "BinaryStar")

    # create the binary with the main star
    tripleMass = [starMass, binary[0].mass + binary[1].mass]
    tripleSemmimajor = starEnvelopeRadius
    triple = StarModels.CreateBinary(binaryMasses= tripleMass, binarySemimajorAxis= tripleSemmimajor)

    # fixing positions
    starEnvelope.position += triple.position[0]
    starCore.position += triple.position[0]
    binary.position += triple.position[1]

    # fixing velocities
    starEnvelope.velocity += triple.velocity[0]
    starCore.velocity += triple.velocity[0]
    binary.velocity += triple.velocity[1]

    return starMass, starEnvelope, starCore, binary, tripleSemmimajor

def TakeSavedState(savedVersionPath, configurationFile):
    '''
    :param savedVersionPath: the path to where you have your saved state
    :return: the saved system
    '''
    print "using saved state file - {0}".format(savedVersionPath)
    starMass, binary, tripleSemmimajor = pickle.load(open(savedVersionPath+".p", 'rb'))
    starEnvelope = read_set_from_file(savedVersionPath+"_envelope.hdf5",'amuse', close_file= True)
    starCore = read_set_from_file(savedVersionPath+"_core.hdf5",'amuse', close_file= True)[0]
    #TODO: check this...
    starEnvelope.move_to_center()
    starCore.position = [0.0, 0.0, 0.0] | units.m



    # create the binary
    newBinary = StarModels.CreateBinary(configurationFile, "BinaryStar")
    newBinary.position += 1.0 | units.AU

    native_plot.figure(figsize=(30, 30), dpi=60)
    sph_particles_plot(starEnvelope)
    native_plot.show()

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
    pickle.dump([starMass, binary, tripleSemmimajor], open(savedVersionPath+".p", 'wb'), pickle.HIGHEST_PROTOCOL)
    write_set_to_file(starEnvelope, savedVersionPath+"_envelope.hdf5", 'amuse' , append_to_file= False)
    write_set_to_file(Particles(particles = [starCore]), savedVersionPath+"_core.hdf5", 'amuse', append_to_file= False)
    print "state saved - {0}".format(savedVersionPath)


def Start(savedVersionPath = "savings/Passy500000", takeSavedState = "False", configurationFile = "PassyConfiguration.ini"):
    '''
    This is the main function of our simulation
    :param savedVersionPath: path to the saved state
    :param takeSavedState: do you have a saved state you want to use? True or False if it is all saved ,
                            Mesa if its only the Mesa Star, Sph - if its only the sph star (after relaxation)
    :return: None
    '''

    # creating the triple system
    if takeSavedState == "True":
        starMass, starEnvelope, starCore, binary, tripleSemmimajor = TakeSavedState(savedVersionPath, configurationFile)
    else:
        if takeSavedState == "Sph":
            starMass, starEnvelope, starCore, binary, tripleSemmimajor = CreateTripleSystem(configurationFile, savedVersionPath, takeSavedSPH= True)
        elif takeSavedState == "Mesa":
            starMass, starEnvelope, starCore, binary, tripleSemmimajor = CreateTripleSystem(configurationFile, savedVersionPath, takeSavedMesa= True)
        else:
            starMass, starEnvelope, starCore, binary, tripleSemmimajor = CreateTripleSystem(configurationFile, savedVersionPath)
        SaveState(savedVersionPath, starMass, starEnvelope, starCore, binary, tripleSemmimajor)


    #EvolveNBody.Run(totalMass= starMass, semmiMajor= tripleSemmimajor, gasParticles= [starEnvelope],
    #               dmParticles= [starCore], endTime= 1000. | units.yr, timeSteps= 12 ,
    #               savedVersionPath= savedVersionPath, step= 1)

    EvolveNBody.Run(totalMass= starMass + binary[0].mass,
                    semmiMajor= tripleSemmimajor, gasParticles= [starEnvelope], dmParticles= [starCore , binary[0]],
                    endTime= 10. | units.yr, timeSteps= 5, savedVersionPath= savedVersionPath)


    #EvolveNBody.EvolveBinary(totalMass= binary[0].mass + binary[1].mass,
    #                semmiMajor= 0.15 | units.AU, binary= binary , endTime= 100 | units.yr, timeSteps = 2)

    # creating the NBody system with the 3 and evolving
    #EvolveNBody.Run(totalMass= starMass + binary[0].mass + binary[1].mass,
    #                semmiMajor= tripleSemmimajor, gasParticles= [starEnvelope], dmParticles= [starCore , binary],
    #                endTime= 100.0 | units.yr, timeSteps= 20)
    print "****************** Simulation Completed ******************"

def MakeAMovieFromSavedState(savedVersionPath= "savings/TCE500000" , steps = []):
    #TODO: do something
    print "blabla"