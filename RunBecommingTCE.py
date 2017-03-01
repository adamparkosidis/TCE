import pickle
import os

from amuse.units import units
from amuse.community.gadget2.interface import Gadget2
from amuse.community.huayno.interface import Huayno
from amuse.community.mi6.interface import MI6
from amuse.couple.bridge import Bridge, CalculateFieldForParticles, CalculateFieldForCodesUsingReinitialize
from amuse.units import units , nbody_system
from amuse.ext import bridge
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


    #hydroSystem.dm_particles.add_particles(innerBinary.stars[1])

    unitConverter = nbody_system.nbody_to_si(outerBinary.stars.total_mass(), sphStar.relaxationTime)
    kickerCode = MI6(unitConverter,number_of_workers= 8, redirection='file', redirect_file='kicker_code_mi6_out.log')
    binarySystem = Huayno(unitConverter)
    binarySystem.particles.add_particles(outerBinary.stars)

    epsilonSquared = (hydroSystem.dm_particles.radius[0]/ 2.8)**2
    kickerCode.parameters.epsilon_squared = 0
    coupledSystem = Bridge(timestep=(sphStar.relaxationTime / (2 * sphStar.relaxationTimeSteps)), verbose=False, use_threading= False)
    kick_from_hydro = CalculateFieldForParticles(particles=hydroSystem.particles, gravity_constant=constants.G)
    kick_from_hydro.smoothing_length_squared = epsilonSquared
    coupledSystem.add_system(binarySystem, (kick_from_hydro,), False)
    coupledSystem.add_system(hydroSystem, (kickerCode,), False)
    coupledSystem.channels.add_channel(binarySystem.particles[0].new_channel_to(hydroSystem.particles))

    print hydroSystem.dm_particles
    print coupledSystem.particles
    starEnvelope, dmStars = EvolveNBody.Run(totalMass= outerBinary.stars.total_mass(),
                    semmiMajor= outerBinary.semimajorAxis, sphEnvelope= sphStar.gas_particles, sphCore=sphStar.core_particle,
                                             stars=None, endTime= sphStar.relaxationTime,
                                             timeSteps= sphStar.relaxationTimeSteps, relax=True,
                                              numberOfWorkers= sphStar.numberOfWorkers, savedVersionPath=savedPath,
                                            saveAfterMinute=10, system=coupledSystem)
    starCore = dmStars[0]
    starCore.radius = sphStar.core_particle.radius

    sphMetaData = StarModels.SphMetaData(sphStar)

    #saved state
    StarModels.SaveState(savedPath, starEnvelope.total_mass() + starCore.mass, starEnvelope, dmStars, outerBinary.semimajorAxis, sphMetaData)

    #moving the main star back to the center
    diffPosition = starCore.position - giant.position
    diffVelocity = (starCore.velocity*starCore.mass + starEnvelope.center_of_mass_velocity() * starEnvelope.total_mass())/ giant.mass
    starEnvelope.position -= diffPosition
    starCore.position -= diffPosition
    starEnvelope.velocity -= diffVelocity
    starCore.velocity -= diffVelocity



    return giant.mass, starEnvelope, starCore, innerBinary, outerBinary, sphMetaData



def Start(savedVersionPath = "Glanz/savings/TCEBecomming/300000/3AU", takeSavedState = "False", step = -1, configurationFile = "Glanz/savings/TCEBecomming/300000/3AU/Configuration.ini"):
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
    NBodySystem = EvolveNBody.DynamicsForBinarySystem(Huayno, innerBinary.semimajorAxis, innerBinary.stars)
    NBodySystem.particles.add_particles(outerBinary.stars[-1])

    unitConverter = nbody_system.nbody_to_si(outerBinary.particles.total_mass(), sphMetaData.evolutionTime)
    #kickerCode = MI6(unitConverter,number_of_workers= 8, redirection='file', redirect_file='kicker_code_mi6_out.log')
    #kickFromCompanions = CalculateFieldForCodesUsingReinitialize(kickerCode, (innerBinary.stars[-1], outerBinary.stars[-1]))
    #kick_from_hydro = CalculateFieldForParticles(particles=hydroSystem.gas_particles, gravity_constant=constants.G)

    kickerCodeOuter = CalculateFieldForParticles(particles=outerBinary.stars[1], gravity_constant=constants.G)
    kick_from_hydro = CalculateFieldForParticles(particles=hydroSystem.gas_particles, gravity_constant=constants.G)

    hydroSystem.dm_particles.add_particles(innerBinary.stars[-1])
    #hydroSystem.dm_particles.add_particles(outerBinary.stars[-1])
    coupledSystem = Bridge()
    coupledSystem.add_system(hydroSystem, (kickerCodeOuter,), False)
    coupledSystem.add_system(kickerCodeOuter, (kick_from_hydro, ), False)

    print innerBinary.stars
    print outerBinary.stars
    print coupledSystem.particles
    EvolveNBody.Run(totalMass= starMass + innerBinary.stars.mass[-1] + outerBinary.stars.mass[-1],
                    semmiMajor= outerBinary.semimajorAxis, sphEnvelope= starEnvelope,
                    sphCore=starCore, stars=innerBinary,
                    endTime= sphMetaData.evolutionTime, timeSteps= sphMetaData.evolutionTimeSteps, numberOfWorkers= sphMetaData.numberOfWorkers, step= step,
                    savedVersionPath=savedVersionPath, saveAfterMinute= 0, system=coupledSystem)

    print "****************** Simulation Completed ******************"
if __name__ == "__main__":
    Start(takeSavedState="Mesa")

