import pickle
import os
import sys, time
from amuse.units import units
from amuse.units.units import *
from amuse.units import units
from amuse.community.gadget2.interface import Gadget2
from amuse.datamodel import Particles
import StarModels
import EvolveNBody


def Start(savedVersionPath="/home/hilaglanz/Documents/80265", step=-1,
          configurationFile="/home/hilaglanz/Documents/80265/TCEConfiguration.ini", mergedParticles=0o1, opposite=False):
    '''
    This is the main function of our simulation
    :param savedVersionPath: path to the saved state
    :param mergedParticles: which particles should be combined where 0 is the core:
            01: 0 + 1
            02: 0 + 2
            12: 1 + 2
    :return: None
    '''
    try:
        os.makedirs(savedVersionPath + "/pics")
    except(OSError):
        pass

    starMass, starEnvelope, starCore, binary, tripleSemmimajor, sphMetaData = \
        StarModels.TakeTripleSavedState(savedVersionPath+ "/evolution", configurationFile, step=step, opposite= opposite)

    # creating the new binary system
    if mergedParticles == "01":
        if not opposite:
            starCore = StarModels.MergeParticles(Particles(particles=[starCore,binary.stars[0]]))
            companion = binary.stars[1]
        else:
            starCore = StarModels.MergeParticles(Particles(particles=[starCore,binary[0]]))
            companion = binary[1]

    elif mergedParticles == "02":
        if not opposite:
            starCore = StarModels.MergeParticles(Particles(particles=[starCore, binary.stars[1]]))
            companion = binary.stars[0]
        else:
            starCore = StarModels.MergeParticles(Particles(particles=[starCore, binary[1]]))
            companion = binary[0]
    else:
        if not opposite:
            companion = StarModels.MergeParticles(Particles(particles=[binary.stars[0], binary.stars[1]]))
        else:
            companion = StarModels.MergeParticles(Particles(particles=[binary[0], binary[1]]))

    print("core: ", starCore)
    print("companion: ", companion)

    outputDirectory = savedVersionPath + "/codes_output_{0}".format(str(time.localtime().tm_year) + "-" +
                            str(time.localtime().tm_mon) + "-" + str(time.localtime().tm_mday) + "-" +
                            str(time.localtime().tm_hour) + ":" + str(time.localtime().tm_min) + ":" +
                            str(time.localtime().tm_sec))

    hydroSystem = EvolveNBody.HydroSystem(Gadget2, starEnvelope, starCore, sphMetaData.evolutionTime,
                                          sphMetaData.evolutionTimeSteps, 0.0 | units.Myr, starCore.radius,
                                          sphMetaData.numberOfWorkers, outputDirectory)

    hydroSystem.dm_particles.add_particle(companion)
    coupledSystem = hydroSystem

    EvolveNBody.EvolveBinary(totalMass=starEnvelope.total_mass() + starCore.mass + companion.mass,
                    semmiMajor=tripleSemmimajor, sphEnvelope=starEnvelope,
                    sphCore=starCore, stars=None,
                    endTime= sphMetaData.evolutionTime, timeSteps= sphMetaData.evolutionTimeSteps, numberOfWorkers= sphMetaData.numberOfWorkers, step= step,
                    savedVersionPath=savedVersionPath, saveAfterMinute= 0, system=coupledSystem)

    print("****************** Simulation Completed ******************")


if __name__ == "__main__":
    args = sys.argv
    if len(args) > 1:
        if len(args) < 5:
            Start(savedVersionPath=args[1], step=int(args[2]),
              configurationFile=args[1] + "/TCEConfiguration.ini", mergedParticles=args[3])
        else:
            Start(savedVersionPath=args[1], step=int(args[2]),
              configurationFile=args[1] + "/TCEConfiguration.ini", mergedParticles=args[3], opposite=bin(int(args[4])))
    else:
        print("can do nothing without arguments")

