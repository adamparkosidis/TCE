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
from StarModels import GiantSPHCenterOfMassPosition, GiantSPHCenterOfMassVelocity


def CreateTripleSystem(configurationFile, savedPath = "", takeSavedSPH = False, takeSavedMesa = False):
    '''
    creating the TCE
    :return:main star's mass, the envelope particles, the core particles, the binary stars and the triple semmimajor
    '''
    print configurationFile
    giant = StarModels.CreatePointStar(configurationFile,configurationSection="MainStar")
    innerBinary = StarModels.Binary(configurationFile, configurationSection="InnerBinary")

    #now setting up the giant (want it to be relaxed and spinning)
    outerBinary = StarModels.Binary(configurationFile, configurationSection="OuterBinary")
    #notice that the giant is the binary.stars[0], the companions are the next

    #the inner binary's center of mass is the second star of the outer binary. so move the center of mass to that place.
    innerBinary.stars.position += outerBinary.stars[1].position
    innerBinary.stars.velocity += outerBinary.stars[1].velocity

    giant.position = outerBinary.stars[0].position
    giant.velocity = outerBinary.stars[0].velocity

    triple = innerBinary.stars
    giantInSet = triple.add_particle(giant)
    innerBinary.stars = triple - giantInSet

    triple.position -= giantInSet.position
    triple.velocity -= giantInSet.velocity

    print triple

    sphStar = StarModels.SphStar(giantInSet,configurationFile,configurationSection="MainStar",
                                savedMesaStarPath = savedPath, takeSavedMesa=takeSavedMesa)
    print sphStar.core_particle

    print "Now having the sph star and the binaries, ready for relaxing"
    starEnvelope, dmStars = EvolveNBody.Run(totalMass= giantInSet.mass + innerBinary.stars.total_mass(),
                    semmiMajor= outerBinary.semimajorAxis, sphEnvelope= sphStar.gas_particles, sphCore=sphStar.core_particle,
                                             stars=innerBinary, endTime= sphStar.relaxationTime,
                                             timeSteps= sphStar.relaxationTimeSteps, relax=True,
                                              numberOfWorkers= 2, savedVersionPath=savedPath, saveAfterMinute=10)
    starCore = dmStars[-1]
    #starCore.radius = sphStar.core_particle.radius
    print starCore
    print "moving the main star back to the center"
    diffPosition = GiantSPHCenterOfMassPosition(starEnvelope, starCore) - giantInSet.position
    #diffVelocity = GiantSPHCenterOfMassVelocity(starEnvelope, starCore) -giantInSet.velocity
    starEnvelope.position -= diffPosition
    starCore.position -= diffPosition
    '''starEnvelope.velocity -= diffVelocity
    starCore.velocity -= diffVelocity'''
    starEnvelope.velocity = giantInSet.velocity
    starCore.velocity = giantInSet.velocity
    dmStars[-1].position = starCore.position
    dmStars[-1].velocity = starCore.velocity
    print diffPosition
    print starCore
    sphMetaData = StarModels.SphMetaData(sphStar)
    
    #saved state
    StarModels.SaveState(savedPath, starEnvelope.total_mass() + starCore.mass, starEnvelope, dmStars, outerBinary.semimajorAxis, sphMetaData)



    return giant.mass, starEnvelope, starCore, innerBinary, outerBinary.semimajorAxis, sphMetaData



def Start(savedVersionPath = "/home/hilaglanz/Documents/80265", takeSavedState = "False", step = -1, configurationFile = "/home/hilaglanz/Documents/80265/TCEConfiguration.ini"):
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

    continueEvolutionSimulation = False
    
    # creating the triple system
    if takeSavedState == "True":
        starMass, starEnvelope, starCore, binary, tripleSemmimajor, sphMetaData = \
            StarModels.TakeTripleSavedState(savedVersionPath, configurationFile, step= -1)
    elif takeSavedState == "Evolve":
        starMass, starEnvelope, starCore, binary, tripleSemmimajor,sphMetaData = \
            StarModels.TakeTripleSavedState(savedVersionPath + "/evolution", configurationFile, step)
    else:
        if takeSavedState == "Mesa":
            starMass, starEnvelope, starCore, binary, tripleSemmimajor, sphMetaData = CreateTripleSystem(configurationFile, savedVersionPath, takeSavedMesa= True)            
        else:
            starMass, starEnvelope, starCore, binary, tripleSemmimajor, sphMetaData = CreateTripleSystem(configurationFile, savedVersionPath)
        step=-1

    # creating the NBody system with the 3 and evolving
    EvolveNBody.Run(totalMass= starMass + binary.stars.total_mass(),
                    semmiMajor= tripleSemmimajor, sphEnvelope= starEnvelope,
                    sphCore=starCore, stars=binary,
                    endTime= sphMetaData.evolutionTime, timeSteps= sphMetaData.evolutionTimeSteps, numberOfWorkers= sphMetaData.numberOfWorkers , step= step,
                    savedVersionPath=savedVersionPath, saveAfterMinute= 0)

    print "****************** Simulation Completed ******************"
if __name__ == "__main__":
    args = sys.argv
    if len(args) > 1:
        Start(savedVersionPath=args[1],takeSavedState=args[2], step=int(args[3]), configurationFile=args[1] + "/TCEConfiguration.ini")
    else:
        Start(takeSavedState="No", step=-1)

