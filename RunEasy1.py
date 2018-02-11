import ConfigParser
import numpy
import pickle
import math
import os
from amuse.lab import *
from amuse.units import *
from amuse.datamodel import Particle
from amuse.ext import orbital_elements
from amuse.plot import plot, native_plot, sph_particles_plot
import time
from amuse.units.quantities import AdaptingVectorQuantity
from amuse.lab import *
from amuse.community.gadget2.interface import Gadget2
from amuse.community.fi.interface import Fi
from amuse.community.hermite0.interface import  Hermite
from amuse.units import units , nbody_system
from amuse.units.units import *
from amuse.lab import Particles
from amuse.io import read_set_from_file, write_set_to_file

from amuse.couple.bridge import Bridge, CalculateFieldForParticles, CalculateFieldForCodesUsingReinitialize
from amuse.community.mi6.interface import MI6
from amuse.community.huayno.interface import Huayno
from amuse.ext.sink import new_sink_particles

if __name__ == "__main__":
    part = Particle()
    evolutionCode = MESA()
    mainStar = evolutionCode.particles.add_particle(part)
    mainStar.evolve_one_step()
    sphStar = convert_stellar_model_to_SPH(mainStar, 500000, do_relax = False, with_core_particle=True,
                                                    target_core_mass = 2.4, base_grid_options=dict(type="fcc"))
    parts = Particles(2)
    parts.mass = [1.0 | units.MSun , 1.0 | units.MSun]
    unitConverter = nbody_system.nbody_to_si(10.0 | units.MSun, 0.511 | units.AU)
    dynamicSystem = Huayno(unitConverter, redirection="file")
    dynamicSystem.particles.add_particles(parts)
    unitConverter = nbody_system.nbody_to_si(8.0 | units.MSun, 0.01*1000*2 | units.RSun)
    hydroSystem = Gadget2(unitConverter, redirection="file",  number_of_workers=12)
    hydroSystem.dm_particles.add_particle(sphStar.core_particle)
    hydroSystem.gas_particles.add_particles(sphStar.gas_particles)

    unitConverter = nbody_system.nbody_to_si(2.0 | units.MSun, 1400 | units.day)
    kickerCode = MI6(unitConverter,number_of_workers= 8, redirection='file', redirect_file='kicker_code_mi6_out.log')
    epsilonSquared = (hydroSystem.dm_particles.radius[0]/ 2.8)**2
    kickerCode.parameters.epsilon_squared = epsilonSquared
    kickFromBinary = CalculateFieldForCodesUsingReinitialize(kickerCode, (dynamicSystem,))
    coupledSystem = Bridge(timestep=(1400.0 / (2 * 7000.0)), verbose=False, use_threading= False)

    kick_from_hydro = CalculateFieldForParticles(particles=hydroSystem.particles, gravity_constant=constants.G)
    kick_from_hydro.smoothing_length_squared = epsilonSquared
    coupledSystem.add_system(dynamicSystem, (kick_from_hydro,), False)
    coupledSystem.add_system(hydroSystem, (kickFromBinary,), False)

    coupledSystem.evolve_model(0.2 | units.day)