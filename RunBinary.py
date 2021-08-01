import pickle
import os, time
import sys

from amuse.community.gadget2.interface import Gadget2
from amuse.units import units
from amuse.units.units import *
from amuse.plot import native_plot, sph_particles_plot
from amuse.ext.star_to_sph import pickle_stellar_model
from amuse.datamodel import Particles
from amuse.io import read_set_from_file
import StarModels
from StarModels import GiantSPHCenterOfMassVelocity, GiantSPHCenterOfMassPosition
import EvolveNBody



def CreateBinarySystem(configurationFile, savedPath = "", takeSavedSPH = 0, takeSavedMesa = False,
                       doubleSPH=False, step=-1):
    '''
    creating the binary
    :return:main star's mass, the envelope particles, the core particles, the binary stars and the binary semmimajor
    '''
    #check if 2SPH and if yes return CreateTwoSPHBinarySystem
    if doubleSPH:
        return CreateTwoSPHBinarySystem(configurationFile,savedPath,takeSavedSPH,takeSavedMesa, step)

    binary = StarModels.Binary(configurationFile, configurationSection="Binary")
    binary.stars.radius = binary.radius
    giant = binary.stars[0]
    print "binary: ", binary.stars

    #put the giant in the center
    binary.stars.position -= giant.position
    binary.stars.velocity -= giant.velocity

    #create the sph giant
    sphStar = StarModels.SphStar(giant, configurationFile,configurationSection="MainStar",
                                savedMesaStarPath = savedPath, takeSavedMesa=takeSavedMesa)
    binary.stars[0].mass = sphStar.particles.total_mass()
    binary.UpdateWithMassChange()


    print "Now having the sph star and the binaries, ready for relaxing"
    starEnvelope, dmStars = EvolveNBody.EvolveBinary(totalMass= binary.stars.total_mass(),
                    semmiMajor= binary.semimajorAxis, sphEnvelope= sphStar.gas_particles, sphCore=sphStar.core_particle,
                                             stars=binary.stars, endTime= sphStar.relaxationTime,
                                             timeSteps= sphStar.relaxationTimeSteps, relax=True,
                                              numberOfWorkers= sphStar.numberOfWorkers, savedVersionPath=savedPath, saveAfterMinute=5, takeCompanionInRelaxation= False)
    '''
    binary.stars.move_to_center()
    giant = binary.stars[0]
    '''
    starCore = dmStars[0]
    starCore.radius = sphStar.core_particle.radius

    sphMetaData = StarModels.SphMetaData(sphStar)


    #moving the main star back to the center
    diffPosition = GiantSPHCenterOfMassPosition(starEnvelope, starCore) - giant.position
    print "diff position: ", diffPosition
    starEnvelope.position -= diffPosition
    starCore.position -= diffPosition
    starEnvelope.velocity = giant.velocity
    starCore.velocity = giant.velocity

    #save state after relaxation
    StarModels.SaveState(savedPath, starEnvelope.total_mass() + starCore.mass, starEnvelope, dmStars, binary.semimajorAxis, sphMetaData)


    return starEnvelope, starCore, binary, binary.semimajorAxis, sphMetaData

def CreateTwoSPHBinarySystem(configurationFile, savedPath = "", takeSavedSPH = 0, takeSavedMesa = False, step=-1):
    '''
    creating the TCE
    :return:main star's mass, the envelope particles, the core particles, the binary stars and the binary semmimajor
    '''

    print("creating a double SPH system")
    binary = StarModels.Binary(configurationFile, configurationSection="Binary")
    binary.stars.radius = binary.radius
    print binary


    #put the giant in the center
    binary.stars.position -= binary.stars[0].position
    binary.stars.velocity -= binary.stars[0].velocity

    sphStar1 = StarModels.SphStar(binary.stars[0],configurationFile,configurationSection="SphStar1",
                                savedMesaStarPath = savedPath + "/sph1", takeSavedMesa=(takeSavedMesa or takeSavedSPH!=0))

    print "Now having the first sph star , ready for relaxing"


    binary.stars.move_to_center()
    binary.stars.position -= binary.stars[1].position
    binary.stars.velocity -= binary.stars[1].velocity

    sphStar2 = StarModels.SphStar(binary.stars[1],configurationFile,configurationSection="SphStar2",
                                savedMesaStarPath = savedPath + "/sph2", takeSavedMesa=(takeSavedMesa or takeSavedSPH==2))

    print "Now having the second sph star , ready for relaxing"

    core_radius = max(sphStar1.core_particle.radius, sphStar2.core_particle.radius)
    sphStar1.core_particle.radius = core_radius
    sphStar2.core_particle.radius = core_radius

    if takeSavedSPH == 0:
        star1Envelope, dmStars1 = EvolveNBody.Run(totalMass= binary.stars.total_mass(),
                        semmiMajor= binary.semimajorAxis, sphEnvelope= sphStar1.gas_particles, sphCore=sphStar1.core_particle,
                                                 stars=binary.stars, endTime= sphStar1.relaxationTime,
                                                 timeSteps= sphStar1.relaxationTimeSteps, relax=True, takeCompanionInRelaxation=False,
                                                  numberOfWorkers= sphStar1.numberOfWorkers, savedVersionPath=savedPath + "/sph1", saveAfterMinute=15)

    else:
        if takeSavedSPH == 2:
            star1Envelope, dmStars1 = StarModels.TakeSPHSavedState(savedPath + "/sph1", step=-1)
        else:
            star1Envelope, dmStars1 = StarModels.TakeSPHSavedState(savedPath + "/sph1/relaxation", step=step)


    star1Core = dmStars1[-1]
    star1Core.radius = sphStar1.core_particle.radius
    sph1MetaData = StarModels.SphMetaData(sphStar1)

    # saved state
    StarModels.SaveState(savedPath + "/sph1", star1Envelope.total_mass() + star1Core.mass, star1Envelope, dmStars1,
                         binary.semimajorAxis, sph1MetaData)

    print("first sph star is relaxed")

    if takeSavedSPH == 2:
        star2Envelope, dmStars2 = StarModels.TakeSPHSavedState(savedPath + "/sph2/relaxation", step)
    else:
        star2Envelope, dmStars2 = EvolveNBody.Run(totalMass= binary.stars.total_mass(),
                        semmiMajor= binary.semimajorAxis, sphEnvelope= sphStar2.gas_particles, sphCore=sphStar2.core_particle,
                                                 stars=binary.stars, endTime= sphStar2.relaxationTime,
                                                 timeSteps= sphStar2.relaxationTimeSteps, relax=True, takeCompanionInRelaxation=False,
                                                  numberOfWorkers= sphStar2.numberOfWorkers, savedVersionPath=savedPath + "/sph2", saveAfterMinute=15)
    star2Core = dmStars2[-1]
    star2Core.radius = sphStar2.core_particle.radius
    sph2MetaData = StarModels.SphMetaData(sphStar2)
    #saved state
    StarModels.SaveState(savedPath + "/sph2", star2Envelope.total_mass() + star2Core.mass, star2Envelope, dmStars2, binary.semimajorAxis, sph2MetaData)
    print("second sph star is relaxed and saved")

    binary.stars.move_to_center()
    binary.stars[0].mass = sphStar1.particles.total_mass()
    binary.UpdateWithMassChange()

    binary.stars[1].mass = sphStar2.particles.total_mass()
    binary.UpdateWithMassChange()

    diffPos1 = sphStar1.particles.center_of_mass() - binary.stars[0].position
    diffPos2 = sphStar2.particles.center_of_mass() - binary.stars[1].position
    diffVel1 = sphStar1.particles.center_of_mass_velocity() - binary.stars[0].velocity
    diffVel2 = sphStar2.particles.center_of_mass_velocity() - binary.stars[1].velocity
    star1Envelope.position -= diffPos1
    star1Core.position -= diffPos1
    star2Envelope.position -= diffPos2
    star2Core.position -= diffPos2
    star1Envelope.velocity-= diffVel1
    star1Core.velocity-= diffVel1
    star2Envelope.velocity-= diffVel2
    star2Core.velocity-= diffVel2

    totalEnvelope = Particles()
    totalEnvelope.add_particles(star1Envelope)
    totalEnvelope.add_particles(star2Envelope)
    totalCores = Particles()
    totalCores.add_particle(star1Core)
    totalCores.add_particle(star2Core)
    StarModels.SaveState(savedPath, totalEnvelope.total_mass() + totalCores.mass, totalEnvelope, totalCores, binary.semimajorAxis, sph1MetaData)
    #now a trick so that the usuall code will include all particles in the hydro- all envelopes and core1 will be added from binary0 and core2 as the companion.
    binary.stars[-1] = star2Core
    binary.stars[0].mass += star2Envelope
    return totalEnvelope, star1Core , binary, binary.semimajorAxis, sph1MetaData


def Start(savedVersionPath = "/vol/sci/astro/bigdata/code/amuse-10.0/Glanz/savings/Passy/500000",
          takeSavedState = "False", step = -1,
          configurationFile = "/vol/sci/astro/bigdata/code/amuse-10.0/Glanz/savings/Passy/500000/PassyConfiguration.ini",
          doubleSPH=False):
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
    relax = False
    simulationTime = None
    simulationTimeSteps = None
    initialCOM = None
    initialCOMV = None
    if takeSavedState == "Single":#continue the relaxation but without forcing the com to stay in place
        loadingStep = -1
        savedModelPath = savedVersionPath

        if step > -1:
            savedModelPath=savedVersionPath+"/evolution"
            loadingStep = step

        starEnvelope, starCore, binary, semmimajor, sphMetaData = \
            StarModels.TakeBinarySavedState(savedModelPath, configurationFile, step= loadingStep)
        outputDirectory = savedVersionPath + "/codes_output_{0}".format(str(time.localtime().tm_year) + "-" +
                                                                        str(time.localtime().tm_mon) + "-" + str(
            time.localtime().tm_mday) + "-" +
                                                                        str(time.localtime().tm_hour) + ":" + str(
            time.localtime().tm_min) + ":" +
                                                                        str(time.localtime().tm_sec))
        os.makedirs(outputDirectory)
        try:
            coreParticleRadius = starCore.epsilon
        except:
            coreParticleRadius = starCore.radius

        currentTime = 0.0 | units.Myr
        system = EvolveNBody.HydroSystem(Gadget2, starEnvelope, starCore, sphMetaData.relaxationTime,
                                         sphMetaData.relaxationTimeSteps, currentTime, coreParticleRadius,
                                  sphMetaData.numberOfWorkers, outputDirectory=outputDirectory + "/hydro")


        EvolveNBody.RunSystem(system,sphMetaData.relaxationTime,sphMetaData.relaxationTimeSteps,savedVersionPath,0,step,False)

        print "****************** Simulation Completed ******************"
        return

    # creating the triple system
    if takeSavedState == "True":
        starEnvelope, starCore, binary, semmimajor, sphMetaData = \
            StarModels.TakeBinarySavedState(savedVersionPath, configurationFile, step= -1, doubleSPH=doubleSPH)
    elif takeSavedState == "Evolve":
        starEnvelope, starCore, binary, semmimajor,sphMetaData = \
            StarModels.TakeBinarySavedState(savedVersionPath + "/evolution", configurationFile, step, doubleSPH=doubleSPH)
    elif takeSavedState == "Relax": # this option is currently supported only for the circumstellar case, for the other need to form the companions
        if doubleSPH:
            starEnvelope, starCore, binary, semmimajor, sphMetaData = CreateBinarySystem(configurationFile,
                                                                                         savedVersionPath, takeSavedSPH=1,
                                                                                         doubleSPH=doubleSPH, step=step)
        else:
            starEnvelope, starCore, binary, semmimajor,sphMetaData = \
                StarModels.TakeBinarySavedState(savedVersionPath + "/relaxation", configurationFile, step=step)
        relax=True
        simulationTime = sphMetaData.relaxationTime
        simulationTimeSteps = sphMetaData.relaxationTimeSteps
        try:
            initialCOM = sphMetaData.initialCOM
            initialCOMV = sphMetaData.initialCOMV
        except:
            print("couldn't rertrieve initial com")
    elif takeSavedState == "Relax2":
        if not doubleSPH:
            print("no implementation for relax2 without doubleSPH")
            return
        starEnvelope, starCore, binary, semmimajor, sphMetaData = CreateBinarySystem(configurationFile,
                                                                                     savedVersionPath,
                                                                                     doubleSPH=doubleSPH,
                                                                                     takeSavedSPH=2, step=step)

    else:
        if takeSavedState == "Mesa":
            starEnvelope, starCore, binary, semmimajor, sphMetaData = CreateBinarySystem(configurationFile,
                                                                                         savedVersionPath,
                                                                                         takeSavedMesa= True,
                                                                                         doubleSPH=doubleSPH)
            print starCore
        else:
            starEnvelope, starCore, binary, semmimajor, sphMetaData = CreateBinarySystem(configurationFile,
                                                                                         savedVersionPath,
                                                                                         doubleSPH=doubleSPH)

        step = -1
    if simulationTime is None:
        simulationTime = sphMetaData.evolutionTime
        simulationTimeSteps= sphMetaData.evolutionTimeSteps

    # creating the NBody system with the 3 and evolving
    EvolveNBody.EvolveBinary(totalMass= binary.stars.total_mass(),
                    semmiMajor= semmimajor, sphEnvelope= starEnvelope,
                    sphCore=starCore, stars=binary.stars,
                    endTime= simulationTime, timeSteps= simulationTimeSteps, numberOfWorkers= sphMetaData.numberOfWorkers, step= step,
                    savedVersionPath=savedVersionPath,relax= relax, initialCOM=initialCOM,
                    initialCOMV=initialCOMV)

    print "****************** Simulation Completed ******************"

def IsTrue(str):
    if str.upper =="TRUE" or str == "1":
        return True
    return False

if __name__ == "__main__":
    args = sys.argv
    doubleSPH = "False"
    if len(args)>4:
        doubleSPH = args[4]

    if len(args) > 1:
        Start(savedVersionPath=args[1],takeSavedState=args[2], step=int(args[3]),
              configurationFile=args[1] + "/BinaryConfiguration.ini", doubleSPH=IsTrue(doubleSPH))
    else:
        Start(takeSavedState="Evolve", step=1625)
