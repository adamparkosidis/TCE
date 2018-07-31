import pickle
import os
import sys

from amuse.units import units
from amuse.units.units import *
from amuse.plot import native_plot, sph_particles_plot
from amuse.ext.star_to_sph import pickle_stellar_model
from amuse.datamodel import Particles
from amuse.io import read_set_from_file
import StarModels
import EvolveNBody



def CreateBinarySystem(configurationFile, savedPath = "", takeSavedSPH = False, takeSavedMesa = False):
    '''
    creating the binary
    :return:main star's mass, the envelope particles, the core particles, the binary stars and the binary semmimajor
    '''
    binary = StarModels.Binary(configurationFile, configurationSection="Binary")
    binary.stars.radius = binary.radius
    giant = binary.stars[0]
    print "giant: ", giant
    #create the sph giant
    sphStar = StarModels.SphStar(giant, configurationFile,configurationSection="MainStar",
                                savedMesaStarPath = savedPath, takeSavedMesa=takeSavedMesa)
    print "Now having the sph star and the binaries, ready for relaxing"
    starEnvelope, dmStars = EvolveNBody.EvolveBinary(totalMass= binary.stars.total_mass(),
                    semmiMajor= binary.semimajorAxis, sphEnvelope= sphStar.gas_particles, sphCore=sphStar.core_particle,
                                             stars=binary, endTime= sphStar.relaxationTime,
                                             timeSteps= sphStar.relaxationTimeSteps, relax=True,
                                              numberOfWorkers= sphStar.numberOfWorkers, savedVersionPath=savedPath, saveAfterMinute=5, takeCompanionInRelaxation= False)
    starCore = dmStars[0]
    starCore.radius = sphStar.core_particle.radius

    sphMetaData = StarModels.SphMetaData(sphStar)
    #save state after relaxation
    StarModels.SaveState(savedPath, starEnvelope.total_mass() + starCore.mass, starEnvelope, dmStars, binary.semimajorAxis, sphMetaData)

    #moving the main star back to the center
    diffPosition = starCore.position - giant.position
    diffVelocity = (starCore.velocity*starCore.mass + starEnvelope.center_of_mass_velocity() * starEnvelope.total_mass())/ giant.mass
    starEnvelope.position -= diffPosition
    starCore.position -= diffPosition
    starEnvelope.velocity -= diffVelocity
    starCore.velocity -= diffVelocity

    return starEnvelope, starCore, binary, binary.semimajorAxis, sphMetaData

def CreateTwoSPHBinarySystem(configurationFile, savedPath = "", takeSavedSPH = False, takeSavedMesa = False):
    '''
    creating the TCE
    :return:main star's mass, the envelope particles, the core particles, the binary stars and the binary semmimajor
    '''
    binary = StarModels.Binary(configurationFile, configurationSection="Binary")
    binary.stars.radius = binary.radius
    print binary

    sphStar1 = StarModels.SphStar(binary[0],configurationFile,configurationSection="SphStar1",
                                savedMesaStarPath = savedPath, takeSavedMesa=takeSavedMesa)

    print "Now having the first sph star , ready for relaxing"
    star1Envelope, dmStars1 = EvolveNBody.Run(totalMass= binary.stars.total_mass(),
                    semmiMajor= binary.semimajorAxis, sphEnvelope= sphStar1.gas_particles, sphCore=sphStar1.core_particle,
                                             stars=binary, endTime= sphStar1.relaxationTime,
                                             timeSteps= sphStar1.relaxationTimeSteps, relax=True,
                                              numberOfWorkers= sphStar1.numberOfWorkers, savedVersionPath=savedPath, saveAfterMinute=15)
    star1Core = dmStars1[-1]
    sph1MetaData = StarModels.SphMetaData(sphStar1)
    #saved state
    StarModels.SaveState(savedPath + "/sph1", binary.mass, star1Envelope, dmStars1, binary.semimajorAxis, sph1MetaData)

    print "first sph star is relaxed"

    sphStar2 = StarModels.SphStar(binary[1],configurationFile,configurationSection="SphStar2",
                                savedMesaStarPath = savedPath, takeSavedMesa=takeSavedMesa)

    print "Now having the second sph star , ready for relaxing"
    star2Envelope, dmStars2 = EvolveNBody.Run(totalMass= binary.stars.total_mass(),
                    semmiMajor= binary.semimajorAxis, sphEnvelope= sphStar2.gas_particles, sphCore=sphStar2.core_particle,
                                             stars=binary, endTime= sphStar2.relaxationTime,
                                             timeSteps= sphStar2.relaxationTimeSteps, relax=True,
                                              numberOfWorkers= sphStar2.numberOfWorkers, savedVersionPath=savedPath, saveAfterMinute=15)
    star2Core = dmStars1[-1]
    sph2MetaData = StarModels.SphMetaData(sphStar2)
    #saved state
    StarModels.SaveState(savedPath + "/sph2", binary.mass, star1Envelope, dmStars2, binary.semimajorAxis, sph1MetaData)
    print "second sph star is relaxed and saved"

    return [star1Envelope, star2Envelope], [star1Core,star2Core] , binary, binary.semimajorAxis, sph1MetaData


def Start(savedVersionPath = "/vol/sci/astro/bigdata/code/amuse-10.0/Glanz/savings/Passy/50000", takeSavedState = "False", step = -1, configurationFile = "/vol/sci/astro/bigdata/code/amuse-10.0/Glanz/savings/Passy/50000/PassyConfiguration.ini"):
    '''
    This is the main function of our simulation
    :param savedVersionPath: path to the saved state
    :param takeSavedState: do you have a saved state you want to use? True or False if it is all saved right before the evolution,
			    Relax if its in the middle of the relaxation, Evolve if its in the evolutionProcess,
                            Mesa if its only the Mesa Star
    :return: None
    '''
    try:
        os.makedirs(savedVersionPath + "/pics")
    except(OSError):
        pass
    # creating the triple system
    if takeSavedState == "True":
        starEnvelope, starCore, binary, semmimajor, sphMetaData = \
            StarModels.TakeBinarySavedState(savedVersionPath, configurationFile, step= -1)
    elif takeSavedState == "Evolve":
        starEnvelope, starCore, binary, semmimajor,sphMetaData = \
            StarModels.TakeBinarySavedState(savedVersionPath + "/evolution", configurationFile, step) 
        print starCore.radius, starEnvelope.velocity, binary.velocity
    else:
        if takeSavedState == "Mesa":
            starEnvelope, starCore, binary, semmimajor, sphMetaData = CreateBinarySystem(configurationFile, savedVersionPath, takeSavedMesa= True)
            print starCore
        else:
            starEnvelope, starCore, binary, semmimajor, sphMetaData = CreateBinarySystem(configurationFile, savedVersionPath)

    # creating the NBody system with the 3 and evolving
    EvolveNBody.EvolveBinary(totalMass= binary.stars.total_mass(),
                    semmiMajor= semmimajor, sphEnvelope= starEnvelope,
                    sphCore=starCore, stars=binary,
                    endTime= sphMetaData.evolutionTime, timeSteps= sphMetaData.evolutionTimeSteps, numberOfWorkers= sphMetaData.numberOfWorkers, step= step,
                    savedVersionPath=savedVersionPath,relax= False)

    print "****************** Simulation Completed ******************"

if __name__ == "__main__":
    args = sys.argv
    if len(args) > 1:
        Start(savedVersionPath=args[1],takeSavedState=args[2], step=int(args[3]), configurationFile=args[1] + "/PassyConfiguration.ini")
    else:
        Start(takeSavedState="Evolve", step=1625)
