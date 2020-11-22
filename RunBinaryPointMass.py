import pickle
import os, time
import sys

from amuse.community.gadget2.interface import Gadget2
from amuse.community.huayno.interface import Huayno
from amuse.units import units
from amuse.units.units import *
from amuse.plot import native_plot, sph_particles_plot
from amuse.ext.star_to_sph import pickle_stellar_model, unpickle_stellar_structure
from amuse.datamodel import Particles
from amuse.io import read_set_from_file
import StarModels
from StarModels import GiantSPHCenterOfMassVelocity, GiantSPHCenterOfMassPosition
import EvolveNBody



def CreateBinarySystem(configurationFile, savedPath = "", takeSavedSPH = False, takeSavedMesa = False):
    '''
    creating the binary
    :return:main star's mass, the envelope particles, the core particles, the binary stars and the binary semmimajor
    '''
    binary = StarModels.Binary(configurationFile, configurationSection="Binary")
    binary.stars.radius = binary.radius
    giant = binary.stars[0]

    #put the giant in the center
    binary.stars.position -= giant.position
    binary.stars.velocity -= giant.velocity


    sphStar = StarModels.SphStar(giant, configurationFile, configurationSection="MainStar",
                                 savedMesaStarPath=savedPath, takeSavedMesa=takeSavedMesa)
    metaData = StarModels.SphMetaData(sphStar)
    binary.stars[0].mass = sphStar.particles.total_mass()
    binary.stars[0].epsilon = 10 * sphStar.core_particle.radius

    print "binary: ", binary.stars

    return binary, metaData

def Start(savedVersionPath = "/vol/sci/astro/bigdata/code/amuse-10.0/Glanz/savings/Passy/500000", takeSavedState = "False", step = -1, configurationFile = "/vol/sci/astro/bigdata/code/amuse-10.0/Glanz/savings/Passy/500000/PassyConfiguration.ini"):
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
        starEnvelope, starCore, binary, semmimajor, metaData = \
            StarModels.TakeBinarySavedState(savedVersionPath, configurationFile, step= -1)
    elif takeSavedState == "Evolve":
        starEnvelope, starCore, binary, semmimajor,metaData = \
            StarModels.TakeBinarySavedState(savedVersionPath + "/evolution", configurationFile, step)
    else:
        if takeSavedState == "Mesa":
            binary, metaData = CreateBinarySystem(configurationFile, savedVersionPath, takeSavedMesa= True)
        else:
            binary, metaData = CreateBinarySystem(configurationFile, savedVersionPath)

        step = -1

    system = EvolveNBody.DynamicsForBinarySystem(Huayno,binary.semimajorAxis,binary,savedVersionPath)

    # creating the NBody system with the binary and evolving
    EvolveNBody.Run(totalMass=binary.stars.total_mass(),
                    semmiMajor=binary.semimajorAxis, sphEnvelope=None,
                    sphCore=binary[0], stars=None,
                    endTime=metaData.evolutionTime, timeSteps=metaData.evolutionTimeSteps,
                    numberOfWorkers=metaData.numberOfWorkers, step=step,
                    savedVersionPath=savedVersionPath, saveAfterMinute=0, system=system)

    print "****************** Simulation Completed ******************"

if __name__ == "__main__":
    args = sys.argv
    if len(args) > 1:
        Start(savedVersionPath=args[1],takeSavedState=args[2], step=int(args[3]), configurationFile=args[1] + "/BinaryConfiguration.ini")
    else:
        Start(takeSavedState="Evolve", step=1625)
