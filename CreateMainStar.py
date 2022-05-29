from amuse.units import units
from amuse.units.units import *
from amuse.community.mesa.interface import MESA
from amuse.ext.star_to_sph import convert_stellar_model_to_SPH

import RelaxModel

### evolve with MESA ###
def  EvolveStarWithMesa(star,radius):
    evolutionType = MESA()
    print("evolving with MESA")
    mainStar = evolutionType.particles.add_particle(star)
    print("particle added, current radius = ",mainStar.radius, "target radius = ", radius)
    while mainStar.radius < radius:
        mainStar.evolve_one_step()
        #print "radius = " ,mainStar.radius
    return mainStar

### returns the sph model and the maximum radius of this star ###
def CreateSphModel(star):
    time = 0 | units.yr
    mesaStar = EvolveStarWithMesa(star.star,star.radius)
    # Convert the star to SPH model ###
    sphModel = convert_stellar_model_to_SPH(mesaStar, star.sphParticles, do_relax = True, with_core_particle=True,
                                            target_core_mass = star.coreMass)
    return sphModel, mesaStar.radius


def GetRelaxedSphModel(star):
    sphStar, radius = CreateSphModel(star)
    starVolume = 4.0*numpy.pi*(radius**3)/3.0
    starAverageDensity = star.star.mass / starVolume
    relaxationTime = 1.0 / (constants.G*starAverageDensity).sqrt() # dynamical time
    #return RelaxModel.Relax(star.star.mass + star.coreMass, star.envelopeRadius, [sphStar.gas_particles], [sphStar.core_particle], relaxationTime.as_quantity_in(yr), star.relaxationTimeSteps)
    return sphStar



