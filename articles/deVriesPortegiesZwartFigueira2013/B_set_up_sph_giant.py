import os
import os.path
import shutil
import math

from amuse.units import units, constants, nbody_system
from amuse.datamodel import Particles, Particle
from amuse.io import write_set_to_file
from amuse.ext.star_to_sph import convert_stellar_model_to_SPH
from amuse.couple.bridge import Bridge, CalculateFieldForCodesUsingReinitialize

from amuse.community.fi.interface import Fi
from amuse.community.huayno.interface import Huayno
from amuse.community.mi6.interface import MI6

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot
from amuse.plot import scatter, xlabel, ylabel, plot, pynbody_column_density_plot, HAS_PYNBODY

from xiTau_parameters import triple_parameters


def new_working_directory():
    i = 0
    current_directory = os.getcwd()
    while os.path.exists(os.path.join(current_directory, "run_{0:=03}".format(i))):
        i += 1
    new_directory = os.path.join(current_directory, "run_{0:=03}".format(i))
    os.mkdir(new_directory)
    print "Created new directory for output:", new_directory
    os.mkdir(os.path.join(new_directory, "plots"))
    shutil.copy(__file__, new_directory)
    os.chdir(new_directory)

def get_relative_velocity_at_apastron(total_mass, semimajor_axis, ecc):
    return (constants.G * total_mass * ((1.0 - ecc)/(1.0 + ecc)) / semimajor_axis).sqrt()

def set_up_initial_conditions(relative_inclination):
    stars = set_up_inner_binary(relative_inclination)
    giant = set_up_outer_star(stars.total_mass())
    giant_in_set = stars.add_particle(giant)
    stars.move_to_center()
    return giant_in_set, stars - giant_in_set

def set_up_inner_binary(inclination):
    semimajor_axis = triple_parameters["a_in"]
    masses = triple_parameters["masses_in"]
    print "   Initializing inner binary"
    stars =  Particles(2)
    stars.mass = masses
    stars.position = [0.0, 0.0, 0.0] | units.AU
    stars.velocity = [0.0, 0.0, 0.0] | units.km / units.s
    stars[0].y = semimajor_axis
    stars[0].vx = -math.cos(inclination)*get_relative_velocity_at_apastron(
        stars.total_mass(), semimajor_axis, 0)
    stars[0].vz = math.sin(inclination)*get_relative_velocity_at_apastron(
        stars.total_mass(), semimajor_axis, 0)
    stars.move_to_center()
    return stars

def set_up_outer_star(inner_binary_mass):
    semimajor_axis = triple_parameters["a_out"]
    eccentricity = triple_parameters["eccentricity_out"]
    print "   Initializing outer star"
    giant = Particle()
    giant.mass = triple_parameters["mass_out"]
    giant.position = semimajor_axis * (1 + eccentricity) * ([1, 0, 0] | units.none)
    giant.velocity = get_relative_velocity_at_apastron(
        giant.mass + inner_binary_mass,
        semimajor_axis, eccentricity) * ([0, 1, 0] | units.none)
    return giant

def set_up_sph_giant(giant, stellar_structure_file, number_of_sph_particles):
    model = convert_stellar_model_to_SPH(
        None, 
        number_of_sph_particles, 
        pickle_file = stellar_structure_file, 
        with_core_particle = True,
        target_core_mass  = 2.0 | units.MSun,
        do_store_composition = False,
        base_grid_options=dict(type="fcc")
    )
    sph_giant = model.gas_particles
    core = model.core_particle
    sph_giant.position += giant.position
    sph_giant.velocity += giant.velocity
    core.position += giant.position
    core.velocity += giant.velocity
    print "Core radius:", model.core_radius.as_string_in(units.RSun)
    return sph_giant, core, model.core_radius

def new_dynamics_for_binary(dynamics_code, binary_particles):
    unit_converter = nbody_system.nbody_to_si(binary_particles.total_mass(), 1.0|units.AU)
    system = dynamics_code(unit_converter, redirection="file", redirect_file="dynamics_code_out.log")
    system.parameters.epsilon_squared = 0 | units.m**2
    system.parameters.inttype_parameter = system.inttypes.SHARED10
    system.parameters.timestep_parameter = 0.2
    system.particles.add_particles(binary_particles)
    return system

def new_hydro(sph_code, sph_giant, core, t_end, n_steps, core_radius):
    unit_converter = nbody_system.nbody_to_si(sph_giant.total_mass() + core.mass, t_end)
    system = sph_code(unit_converter, redirection="file", redirect_file="sph_code_out.log")
    system.parameters.timestep = t_end / n_steps
    system.parameters.eps_is_h_flag = True
    core.radius = core_radius * 2
    system.dm_particles.add_particle(core)
    system.gas_particles.add_particles(sph_giant)
    return system

def new_coupled_system(hydro, binary_system, t_end, n_steps):
    unit_converter = nbody_system.nbody_to_si(binary_system.particles.total_mass(), t_end)
    kicker_code = MI6(unit_converter, redirection='file', redirect_file='kicker_code_mi6_out.log')
    kicker_code.parameters.epsilon_squared = 100 | units.RSun**2
    kick_from_binary = CalculateFieldForCodesUsingReinitialize(kicker_code, (binary_system,))
    
    coupled_system = Bridge(timestep=(t_end / (2 * n_steps)), verbose=False, use_threading=False)
    coupled_system.add_system(hydro, (kick_from_binary,), False)
    return coupled_system

def evolve_system(coupled_system, t_end, n_steps):
    times = (t_end * range(1, n_steps+1) / n_steps).as_quantity_in(units.day)
    hydro = coupled_system.codes[0].code # only calculate potential energy for the giant (SPH particles)
    potential_energies = hydro.potential_energy.as_vector_with_length(1).as_quantity_in(units.J)
    kinetic_energies = hydro.kinetic_energy.as_vector_with_length(1).as_quantity_in(units.J)
    thermal_energies = hydro.thermal_energy.as_vector_with_length(1).as_quantity_in(units.J)
    
    giant = coupled_system.particles
    giant_initial_position = giant.center_of_mass()
    giant_initial_velocity = giant.center_of_mass_velocity()
    for i_step, time in enumerate(times):
        giant.position += giant_initial_position - giant.center_of_mass()
        giant.velocity = (giant.velocity - giant.center_of_mass_velocity()) * (i_step * 1.0 / n_steps) + giant_initial_velocity
        
        print "   Evolving...",
        coupled_system.evolve_model(time)
        print "   Evolved to:", time,
        
        potential_energies.append(hydro.potential_energy)
        kinetic_energies.append(hydro.kinetic_energy)
        thermal_energies.append(hydro.thermal_energy)
        print "   Energies calculated"
        density_plot(coupled_system, i_step)
        
    snapshotfile = "hydro_giant_gas.amuse"
    write_set_to_file(coupled_system.gas_particles, snapshotfile, format='amuse')
    snapshotfile = "hydro_giant_dm.amuse"
    write_set_to_file(coupled_system.dm_particles, snapshotfile, format='amuse')
    
    coupled_system.stop()
    
    energy_evolution_plot(times[:len(kinetic_energies)-1], kinetic_energies, 
            potential_energies, thermal_energies)

def density_plot(coupled_system, i_step):
    if not HAS_PYNBODY:
        return
    figname = os.path.join("plots", "hydro_giant{0:=04}.png".format(i_step))
    print "  -   Hydroplot saved to: ", figname
    pynbody_column_density_plot(coupled_system.gas_particles, width=3|units.AU, vmin=26, vmax=33)
    scatter(coupled_system.dm_particles.x, coupled_system.dm_particles.y, c="w")
    pyplot.savefig(figname)
    pyplot.close()

def energy_evolution_plot(time, kinetic, potential, thermal, figname = "energy_evolution.png"):
    time.prepend(0.0 | units.day)
    pyplot.figure(figsize = (5, 5))
    plot(time, kinetic, label='K')
    plot(time, potential, label='U')
    plot(time, thermal, label='Q')
    plot(time, kinetic + potential + thermal, label='E')
    xlabel('Time')
    ylabel('Energy')
    pyplot.legend(prop={'size':"x-small"}, loc=4)
    pyplot.savefig(figname)
    pyplot.close()


if __name__ == "__main__":
    dynamics_code = Huayno
    sph_code = Fi
    number_of_sph_particles = 100000
    stellar_structure_file = os.path.join(os.getcwd(), "giant_models", "model_0100_100.2")
    
    relative_inclination = math.radians(0.0)
    
    t_end = 100.0 | units.day
    n_steps = 1000
    
    new_working_directory()
    
    print "Initializing triple"
    giant, binary = set_up_initial_conditions(relative_inclination)
    print "\nInitialization done:\n", giant + binary
    
    print "Converting giant from file {0} to {1} SPH particles".format(stellar_structure_file, number_of_sph_particles)
    sph_giant, core, core_radius = set_up_sph_giant(giant, stellar_structure_file, number_of_sph_particles)


    print "\nSetting up {0} to simulate triple system".format(sph_code.__name__)
    hydro = new_hydro(sph_code, sph_giant, core, t_end, n_steps, core_radius)
    print "\nSetting up {0} to simulate triple system".format(dynamics_code.__name__)
    binary_system = new_dynamics_for_binary(dynamics_code, binary)
    print "\nSetting up Bridge to simulate triple system"
    coupled_system = new_coupled_system(hydro, binary_system, t_end, n_steps)
    
    print "\nEvolving to {0}".format(t_end)
    evolve_system(coupled_system, t_end, n_steps)
    
    print "Done"
