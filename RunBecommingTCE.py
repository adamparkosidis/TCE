import pickle
import os, sys, time

from amuse.units import units
from amuse.community.gadget2.interface import Gadget2
from amuse.community.huayno.interface import Huayno
from amuse.community.hermite0.interface import Hermite
from amuse.community.mi6.interface import MI6
from amuse.couple.bridge import Bridge, CalculateFieldForParticles, CalculateFieldForCodesUsingReinitialize
from amuse.units import units , nbody_system
from amuse.ext import bridge
from amuse.units.units import *
from amuse.plot import native_plot, sph_particles_plot
from amuse.ext.star_to_sph import pickle_stellar_model
from amuse.datamodel import Particles
from amuse.io import read_set_from_file
from amuse.units.nbody_system import nbody_to_si
import StarModels
import EvolveNBody


def CreateTripleSystem(configurationFile, savedPath = "", takeSavedSPH = False, takeSavedMesa = False):
    '''
    creating the TCE
    the inner binary is made by the giant and one compact object, the outer binary is made of another compact object
    and the center of mass of the inner binary.
    :return:main star's mass, the envelope particles, the core particles, the binary stars and the triple semmimajor

    '''
    giant = StarModels.CreatePointStar(configurationFile,configurationSection="MainStar")
    innerBinary = StarModels.Binary(configurationFile, configurationSection="InnerBinary")

    #now setting up the giant (want it to be relaxed and spinning)
    outerBinary = StarModels.Binary(configurationFile, configurationSection="OuterBinary")


    sphStar = StarModels.SphStar(giant,configurationFile,configurationSection="MainStar",
                                savedMesaStarPath = savedPath, takeSavedMesa=takeSavedMesa)

    innerBinary.stars[0].mass = sphStar.particles.total_mass()
    innerBinary.UpdateWithMassChange()

    # the inner binary's center of mass is the second star of the outer binary. so move the center of mass to that place.
    innerBinary.stars.position += outerBinary.stars[1].position
    innerBinary.stars.velocity += outerBinary.stars[1].velocity

    # we now move the system so the giant will be in the middle
    giantPossitionDiff = innerBinary.stars[0].position
    giantVelocityDiff = innerBinary.stars[0].velocity
    innerBinary.stars.position -= giantPossitionDiff
    innerBinary.stars.velocity -= giantVelocityDiff
    outerBinary.stars.position -= giantPossitionDiff
    outerBinary.stars.velocity -= giantVelocityDiff

    giant.position = innerBinary.stars[0].position
    giant.velocity = innerBinary.stars[0].velocity

    print "Now having the sph star and the binaries, ready for relaxing"
    outputDirectory = savedPath + "/codes_output_{0}".format(str(time.localtime().tm_year) + "-" +
                            str(time.localtime().tm_mon) + "-" + str(time.localtime().tm_mday) + "-" +
                            str(time.localtime().tm_hour) + ":" + str(time.localtime().tm_min) + ":" +
                            str(time.localtime().tm_sec))
    os.makedirs(outputDirectory)
    hydroSystem = EvolveNBody.HydroSystem(Gadget2, sphStar.gas_particles, sphStar.core_particle, sphStar.relaxationTime,
                                          sphStar.relaxationTimeSteps, 0.0 | units.Myr,
                                          sphStar.core_particle.radius * 20.0 * (250.0 * 1000.0 / len(sphStar.gas_particles)),
                                          sphStar.numberOfWorkers, outputDirectory)
    companionField = CalculateFieldForParticles(Particles(particles=[innerBinary.stars[1],outerBinary.stars[0]]), gravity_constant=constants.G)
    coupledSystem = Bridge(timestep=(sphStar.relaxationTime / (2 * sphStar.relaxationTimeSteps)), verbose=False, use_threading= False)
    coupledSystem.add_system(hydroSystem,(companionField,),False,h_smooth_is_eps=True)
    '''
    #kickFromBinary = CalculateFieldForParticles(particles=Particles(particles=[innerBinary.stars[1],outerBinary.stars[1]]), gravity_constant=constants.G)
    nbodyConverter = nbody_to_si(innerBinary.stars[1].mass, sphStar.relaxationTime)
    companion1 = Hermite(nbodyConverter)
    companion1.particles.add_particle(innerBinary.stars[1])
    nbodyConverter = nbody_to_si(outerBinary.stars[1].mass, sphStar.relaxationTime)
    companion2 = Hermite(nbodyConverter)
    companion2.particles.add_particle(outerBinary.stars[1])
    epsilonSquared = (hydroSystem.dm_particles.radius[0] / 2.8)**2
    #kickFromBinary.smoothing_length_squared = epsilonSquared
    coupledSystem.add_system(hydroSystem, (companion1, companion2,), False)

    coupledSystem= hydroSystem'''
    starEnvelope, dmStars = EvolveNBody.Run(totalMass= outerBinary.stars.total_mass(),
                    semmiMajor= outerBinary.semimajorAxis, sphEnvelope= sphStar.gas_particles, sphCore=sphStar.core_particle,
                                             stars=None, endTime= sphStar.relaxationTime,
                                             timeSteps= sphStar.relaxationTimeSteps, relax=True,
                                              numberOfWorkers= sphStar.numberOfWorkers, savedVersionPath=savedPath, saveAfterMinute=10,system=coupledSystem)

    starCore = dmStars[0]
    starCore.radius = sphStar.core_particle.radius
    #moving the com back to the center
    giant.position += giantPossitionDiff
    giant.velocity += giantVelocityDiff
    innerBinary.stars.position += giantPossitionDiff
    innerBinary.stars.velocity += giantVelocityDiff
    outerBinary.stars.position += giantPossitionDiff
    outerBinary.stars.velocity += giantVelocityDiff


    sphMetaData = StarModels.SphMetaData(sphStar)


    #saved state
    StarModels.SaveState(savedPath, starEnvelope.total_mass() + starCore.mass, starEnvelope, dmStars, outerBinary.semimajorAxis, sphMetaData)

    #moving the main star back to the center
    diffPosition = starCore.position - giant.position
    #diffVelocity = (starCore.velocity*starCore.mass + starEnvelope.center_of_mass_velocity() * starEnvelope.total_mass())/ giant.mass
    starEnvelope.position -= diffPosition
    starCore.position = giant.position
    starEnvelope.velocity = giant.velocity
    starCore.velocity = giant.velocity
    companions = Particles()
    companions.add_particle(innerBinary.stars[1])
    companions.add_particle(outerBinary.stars[0])

    return giant.mass, starEnvelope, starCore, companions, outerBinary.semimajorAxis, sphMetaData



def Start(savedVersionPath = "Glanz/savings/TCEBecomming/500000/nbody", takeSavedState = "False", step = -1, configurationFile = "Glanz/savings/TCEBecomming/500000/Configuration.ini"):
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
        starMass, starEnvelope, starCore,companions, outerSemmimajor, sphMetaData = \
            StarModels.TakeTripleSavedState(savedVersionPath, configurationFile, step= -1, opposite=True)
    elif takeSavedState == "Evolve":
        starMass, starEnvelope, starCore,companions, outerSemmimajor, sphMetaData = \
            StarModels.TakeTripleSavedState(savedVersionPath + "/evolution", configurationFile, step, opposite=True)
    else:
        if takeSavedState == "Mesa":
            starMass, starEnvelope, starCore,companions, outerSemmimajor, sphMetaData = CreateTripleSystem(configurationFile, savedVersionPath, takeSavedMesa= True)
        else:
            starMass, starEnvelope, starCore,companions, outerSemmimajor, sphMetaData = CreateTripleSystem(configurationFile, savedVersionPath)

    # creating the NBody system with the 3 and evolving
    outputDirectory = savedVersionPath + "/codes_output_{0}".format(str(time.localtime().tm_year) + "-" +
                            str(time.localtime().tm_mon) + "-" + str(time.localtime().tm_mday) + "-" +
                            str(time.localtime().tm_hour) + ":" + str(time.localtime().tm_min) + ":" +
                            str(time.localtime().tm_sec))
    try:
        coreRadius = starCore.epsilon
    except:
        coreRadius = starCore.radius
    hydroSystem = EvolveNBody.HydroSystem(Gadget2, starEnvelope, starCore, sphMetaData.evolutionTime,
                                          sphMetaData.evolutionTimeSteps, 0.0 | units.Myr, coreRadius,
                                          sphMetaData.numberOfWorkers, outputDirectory)

    hydroSystem.dm_particles.add_particle(companions[1])
    hydroSystem.dm_particles.add_particle(companions[0])
    coupledSystem = hydroSystem

    '''unitConverter = nbody_system.nbody_to_si(outerBinary.stars[1].mass + innerBinary.stars[1].mass, sphMetaData.evolutionTime)
    binarySystem = Huayno(unitConverter)
    binarySystem.particles.add_particle(innerBinary.stars[1])
    binarySystem.particles.add_particle(starCore)
    coupledSystem = Bridge(timestep=(sphMetaData.evolutionTime / (2 * sphMetaData.evolutionTimeSteps)), verbose=False, use_threading= False)
    coupledSystem.add_system(binarySystem)
    print coupledSystem.particles
    coupledSystem.add_system(hydroSystem)
    print "bridging between ", hydroSystem.dm_particles[::-1][1:]
    coupledSystem.channels.add_channel(binarySystem.particles.new_channel_to(hydroSystem.dm_particles[::-1][1:]))
    starsToSave = Particles(particles=[innerBinary.stars[1], starCore])
    binarySystem.particles.new_channel_to(starsToSave)
    hydroSystem.gas_particles.new_channel_to(starEnvelope)'''
    EvolveNBody.Run(totalMass= companions.total_mass()+starMass,
                    semmiMajor= outerSemmimajor, sphEnvelope= starEnvelope,
                    sphCore=starCore, stars=None,
                    endTime= sphMetaData.evolutionTime, timeSteps= sphMetaData.evolutionTimeSteps, numberOfWorkers= sphMetaData.numberOfWorkers, step= step,
                    savedVersionPath=savedVersionPath, saveAfterMinute= 0, system=coupledSystem)

    print "****************** Simulation Completed ******************"
if __name__ == "__main__":
    args = sys.argv
    if len(args) > 1:
        Start(savedVersionPath=args[1],takeSavedState=args[2], step=int(args[3]), configurationFile=args[1] + "/Configuration.ini")
    else:
        Start()

