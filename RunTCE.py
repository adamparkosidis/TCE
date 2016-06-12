import pickle
import os

from amuse.units import units
from amuse.units.units import *
from amuse.plot import native_plot, sph_particles_plot
from amuse.ext.star_to_sph import pickle_stellar_model
from amuse.datamodel import Particles
from amuse.io import read_set_from_file
import StarModels
import EvolveNBody

class SphMetaData:
    def __init__(self,sphStar):
        self.relaxationTime = sphStar.relaxationTime
        self.relaxationTimeSteps = sphStar.relaxationTimeSteps
        self.evolutionTime = sphStar.evolutionTime
        self.evolutionTimeSteps = sphStar.evolutionTimeSteps
        self.numberOfWorkers = sphStar.numberOfWorkers

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
    innerBinary.stars = triple - giantInSet

    sphStar = StarModels.SphStar(giantInSet,configurationFile,configurationSection="MainStar",
                                savedMesaStarPath = savedPath, takeSavedMesa=takeSavedMesa)
    print "Now having the sph star and the binaries, ready for relaxing"
    starEnvelope, dmStars = EvolveNBody.Run(totalMass= giant.mass + innerBinary.stars[0].mass +
                                                         innerBinary.stars[1].mass,
                    semmiMajor= outerBinary.semimajorAxis, sphEnvelope= sphStar.gas_particles, sphCore=sphStar.core_particle,
                                             stars=innerBinary, endTime= sphStar.relaxationTime,
                                             timeSteps= sphStar.relaxationTimeSteps, relax=True,
                                              numberOfWorkers= sphStar.numberOfWorkers, savedVersionPath=savedPath, saveAfterMinute=15)
    starCore = dmStars[-1]
    sphMetaData = SphMetaData(sphStar)

    return giant.mass, starEnvelope, starCore, innerBinary, outerBinary.semimajorAxis, sphMetaData



def Start(savedVersionPath = "Glanz/savings/TCETry", takeSavedState = "False", step = -1, configurationFile = "TCEConfiguration.ini"):
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
        starMass, starEnvelope, starCore, binary, tripleSemmimajor, sphMetaData = StarModels.TakeSavedState(savedVersionPath, configurationFile, step= -1)
    elif takeSavedState == "Evolve":
        starMass, starEnvelope, starCore, binart, tripleSemmimajor,sphMetaData = StarModels.TakeSavedState(savedVersionPath + "/evolution", configurationFile, step)
    else:
        if takeSavedState == "Mesa":
            starMass, starEnvelope, starCore, binary, tripleSemmimajor, sphMetaData = CreateTripleSystem(configurationFile, savedVersionPath, takeSavedMesa= True)
        else:
            starMass, starEnvelope, starCore, binary, tripleSemmimajor, sphMetaData = CreateTripleSystem(configurationFile, savedVersionPath)


    #TODO: change positions of binary stars

    # creating the NBody system with the 3 and evolving
    EvolveNBody.Run(totalMass= starMass + binary.stars[0].mass + binary.stars[1].mass,
                    semmiMajor= tripleSemmimajor, sphEnvelope= starEnvelope,
                    sphCore=starCore, stars=binary,
                    endTime= sphMetaData.evolutionTime, timeSteps= sphMetaData.evolutionTimeSteps, numberOfWorkers= sphMetaData.numberOfWorkers, step= step,
                    savedVersionPath=savedVersionPath)

    print "****************** Simulation Completed ******************"
if __name__ == "__main__":
    Start(takeSavedState="True")

def MakeAMovieFromSavedState(savedVersionPath= "savings/TCE500000" , steps = []):
    #TODO: do something
    print "blabla"
