import pickle
import os

from amuse.units import units
from amuse.community.gadget2.interface import Gadget2
from amuse.units.units import *
from amuse.plot import native_plot, sph_particles_plot
from amuse.ext.star_to_sph import pickle_stellar_model
from amuse.datamodel import Particles
from amuse.io import read_set_from_file
import StarModels
import EvolveNBody



def CreateTripleSystem(configurationFile, savedPath = "", takeSavedSPH = False, takeSavedMesa = False):
    '''
    creating the TCE
    the inner binary is made by the giant and one MS star, the outer binary is made of the center of mass of the inner binary with another far MS star.
    :return:main star's mass, the envelope particles, the core particles, the binary stars and the triple semmimajor

    '''
    giant = StarModels.CreatePointStar(configurationFile,configurationSection="MainStar")
    innerBinary = StarModels.Binary(configurationFile, configurationSection="InnerBinary")

    #now setting up the giant (want it to be relaxed and spinning)
    outerBinary = StarModels.Binary(configurationFile, configurationSection="OuterBinary")

    outerBinary.stars.position += innerBinary.stars.center_of_mass()
    outerBinary.stars.velocity += innerBinary.stars.center_of_mass_velocity()
    giant.position = innerBinary.stars[0].position
    giant.velocity = innerBinary.stars[0].velocity

    sphStar = StarModels.SphStar(giant,configurationFile,configurationSection="MainStar",
                                savedMesaStarPath = savedPath, takeSavedMesa=takeSavedMesa)

    print "Now having the sph star and the binaries, ready for relaxing"

    hydroSystem = EvolveNBody.HydroSystem(Gadget2, sphStar.gas_particles, sphStar.core_particle, sphStar.evolutionTime,
                                          sphStar.evolutionTimeSteps, 0.0 | units.Myr, sphStar.core_particle.radius,
                                          sphStar.numberOfWorkers)
    #adding the companions
    hydroSystem.dm_particles.add_particle(innerBinary.stars[-1])
    hydroSystem.dm_particles.add_particle(outerBinary.stars[-1])

    print hydroSystem.dm_particles
    starEnvelope, dmStars = EvolveNBody.Run(totalMass= outerBinary.stars.total_mass(),
                    semmiMajor= outerBinary.semimajorAxis, sphEnvelope= sphStar.gas_particles, sphCore=sphStar.core_particle,
                                             stars=None, endTime= sphStar.relaxationTime,
                                             timeSteps= sphStar.relaxationTimeSteps, relax=True,
                                              numberOfWorkers= sphStar.numberOfWorkers, savedVersionPath=savedPath,
                                            saveAfterMinute=10, system=hydroSystem)
    starCore = dmStars[0]
    starCore.radius = sphStar.core_particle.radius
    sphMetaData = StarModels.SphMetaData(sphStar)

    #saved state
    StarModels.SaveState(savedPath, starEnvelope.total_mass() + starCore.mass, starEnvelope, dmStars, outerBinary.semimajorAxis, sphMetaData)

    return giant.mass, starEnvelope, starCore, innerBinary, outerBinary, sphMetaData



def Start(savedVersionPath = "Glanz/savings/TCEBecomming/300000", takeSavedState = "False", step = -1, configurationFile = "Glanz/savings/TCEBecomming/300000/Configuration.ini"):
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
        starMass, starEnvelope, starCore,innerBinary, outerBinary, sphMetaData = \
            StarModels.TakeTripleSavedState(savedVersionPath, configurationFile, step= -1, opposite=True)
    elif takeSavedState == "Evolve":
        starMass, starEnvelope, starCore,innerBinary, outerBinary, sphMetaData = \
            StarModels.TakeTripleSavedState(savedVersionPath + "/evolution", configurationFile, step, opposite=True)
    else:
        if takeSavedState == "Mesa":
            starMass, starEnvelope, starCore, innerBinary, outerBinary, sphMetaData = CreateTripleSystem(configurationFile, savedVersionPath, takeSavedMesa= True)
        else:
            starMass, starEnvelope, starCore, innerBinary, outerBinary, sphMetaData = CreateTripleSystem(configurationFile, savedVersionPath)

    # creating the NBody system with the 3 and evolving

    hydroSystem = EvolveNBody.HydroSystem(Gadget2, starEnvelope, starCore, sphMetaData.evolutionTime,
                                          sphMetaData.evolutionTimeSteps, 0.0 | units.Myr, starCore.radius,
                                          sphMetaData.numberOfWorkers)
    #adding the companions
    hydroSystem.dm_particles.add_particle(innerBinary.stars[-1])
    hydroSystem.dm_particles.add_particle(outerBinary.stars[-1])

    EvolveNBody.Run(totalMass= starMass + innerBinary.stars.mass[-1] + outerBinary.stars.mass[-1],
                    semmiMajor= outerBinary.semimajorAxis, sphEnvelope= starEnvelope,
                    sphCore=starCore, stars=innerBinary,
                    endTime= sphMetaData.evolutionTime, timeSteps= sphMetaData.evolutionTimeSteps, numberOfWorkers= sphMetaData.numberOfWorkers, step= step,
                    savedVersionPath=savedVersionPath, saveAfterMinute= 0)

    print "****************** Simulation Completed ******************"
if __name__ == "__main__":
    Start(takeSavedState="True")

