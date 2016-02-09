import os
import os.path
import shutil
import math

from amuse.units import units, constants, nbody_system
from amuse.units.quantities import AdaptingVectorQuantity
from amuse.datamodel import Particles, Particle
from amuse.io import write_set_to_file, read_set_from_file
from amuse.ext.star_to_sph import convert_stellar_model_to_SPH
from amuse.ext.sink import new_sink_particles
from amuse.couple.bridge import Bridge, CalculateFieldForParticles, CalculateFieldForCodesUsingReinitialize

from amuse.community.fi.interface import Fi
from amuse.community.gadget2.interface import Gadget2
from amuse.community.huayno.interface import Huayno
from amuse.community.mi6.interface import MI6


import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot
from amuse.plot import scatter, xlabel, ylabel, plot, native_plot, sph_particles_plot


from B_set_up_sph_giant import set_up_initial_conditions, new_dynamics_for_binary, new_hydro

def new_working_directory():
    i = 0
    current_directory = os.getcwd()
    while os.path.exists(os.path.join(current_directory, "run_{0:=03}".format(i))):
        i += 1
    new_directory = os.path.join(current_directory, "run_{0:=03}".format(i))
    os.mkdir(new_directory)
    print "Created new directory for output:", new_directory
    os.mkdir(os.path.join(new_directory, "plots"))
    os.mkdir(os.path.join(new_directory, "snapshots"))
    shutil.copy(__file__, new_directory)
    os.chdir(new_directory)

def load_sph_giant(gas_particles_file, dm_particles_file):
    sph_giant = read_set_from_file(gas_particles_file, format='amuse')
    core = read_set_from_file(dm_particles_file, format='amuse')[-1]
    return sph_giant, core

def new_coupled_system(hydro, binary_system, t_end, n_steps):
    kick_from_hydro = CalculateFieldForParticles(particles=hydro.particles, gravity_constant=constants.G)
    kick_from_hydro.smoothing_length_squared = 1.0 | units.RSun**2
    
    unit_converter = nbody_system.nbody_to_si(binary_system.particles.total_mass(), t_end)
    kicker_code = MI6(unit_converter, redirection='file', redirect_file='kicker_code_mi6_out.log')
    kicker_code.parameters.epsilon_squared = 1.0 | units.AU**2
    kick_from_binary = CalculateFieldForCodesUsingReinitialize(kicker_code, (binary_system,))
    
    coupled_system = Bridge(timestep=(t_end / (2 * n_steps)), verbose=False, use_threading=True)
    coupled_system.add_system(binary_system, (kick_from_hydro,), False)
    coupled_system.add_system(hydro, (kick_from_binary,), False)
    return coupled_system

def evolve_system(coupled_system, t_end, n_steps):
    times = (t_end * range(1, n_steps+1) / n_steps).as_quantity_in(units.day)
    
    sinks = new_sink_particles(coupled_system.codes[0].particles, sink_radius=1.0|units.RSun)

    hydro = coupled_system.codes[0].code # only calculate potential energy for the giant (SPH particles)
    potential_energies = hydro.potential_energy.as_vector_with_length(1).as_quantity_in(units.J)
    kinetic_energies = hydro.kinetic_energy.as_vector_with_length(1).as_quantity_in(units.J)
    thermal_energies = coupled_system.thermal_energy.as_vector_with_length(1).as_quantity_in(units.J)

    x =  AdaptingVectorQuantity()
    y =  AdaptingVectorQuantity()
    z =  AdaptingVectorQuantity()
    x.append(hydro.particles.x)
    y.append(hydro.particles.y)
    z.append(hydro.particles.z)

    native_plot.figure(figsize=(20, 20), dpi=60)
    sph_particles_plot(coupled_system.gas_particles)
    native_plot.savefig('plots/0.jpg')

    for i_step, time in enumerate(times):
        sinks.accrete(coupled_system.gas_particles)
        coupled_system.evolve_model(time)
        print "   Evolved to:", time
        potential_energies.append(coupled_system.potential_energy)
        kinetic_energies.append(coupled_system.kinetic_energy)
        thermal_energies.append(coupled_system.thermal_energy)
        print "   Energies calculated"
        #density_plot(coupled_system, i_step)
        if i_step % 1 == 1:
            snapshotfile = os.path.join("snapshots", "hydro_triple_{0:=06}_gas.amuse".format(i_step))
            write_set_to_file(coupled_system.gas_particles, snapshotfile, format='amuse')
            snapshotfile = os.path.join("snapshots", "hydro_triple_{0:=06}_dm.amuse".format(i_step))
            write_set_to_file(coupled_system.dm_particles, snapshotfile, format='amuse')
            sph_particles_plot(coupled_system.gas_particles)
            native_plot.savefig('plots/{0}.jpg'.format(i_step))
        x.append(hydro.particles.x)
        y.append(hydro.particles.y)
        z.append(hydro.particles.z)
    #energy_evolution_plot(times[:len(kinetic_energies)-1], kinetic_energies,
    #        potential_energies, thermal_energies)

    x=x.value_in(units.AU)
    y=y.value_in(units.AU)
    pyplot.figure(figsize= (20, 20), dpi= 80)
    pyplot.plot(x[:, 0], y[:, 0], 'r.', ms= 20.0, )
    pyplot.plot(x[:, 1], y[:, 1], 'g.')
    pyplot.plot(x[:, 2], y[:, 2], 'b.')
    pyplot.plot(x[:, 0], y[:, 0], 'r.', ms= 10.0, )
    pyplot.plot(x[:, 1], y[:, 1], 'g.')
    pyplot.plot(x[:, 2], y[:, 2], 'b.')
    pyplot.xlim(-20, 20)
    pyplot.ylim(-20, 20)
    pyplot.xlabel('AU')
    pyplot.savefig('{0}_dynamics.jpg'.format("evolution"))

    coupled_system.stop()

'''
def density_plot(coupled_system, i_step):
    if not HAS_PYNBODY:
        print "problem plotting"
        #return
    figname = os.path.join("plots", "hydro_giant{0:=04}.png".format(i_step))
    print "  -   Hydroplot saved to: ", figname
    pynbody_column_density_plot(coupled_system.gas_particles, width=5|units.AU, vmin=20, vmax=100)
    scatter(coupled_system.dm_particles.x, coupled_system.dm_particles.y, c="w")
    pyplot.savefig(figname)
    pyplot.close()
'''
def energy_evolution_plot(time, kinetic, potential, thermal, figname = "energy_evolution1.png"):
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
    gas_particles_file = os.path.join(os.getcwd(), "run_001", "hydro_giant_gas.amuse")
    dm_particles_file = os.path.join(os.getcwd(), "run_001", "hydro_giant_dm.amuse")
    
    # Output from set_up_sph_giant in B_set_up_sph_giant.py:
    core_radius = 0.7047075092| units.RSun
    
    relative_inclination = math.radians(9.0)
    
    t_end = 1400.0 | units.day
    n_steps = 7000
    
    new_working_directory()
    print "Initializing triple"
    giant, binary = set_up_initial_conditions(relative_inclination)
    print "\nInitialization done:\n", giant + binary
    
    sph_giant, core = load_sph_giant(gas_particles_file, dm_particles_file)
    
    print "\nSetting up {0} to simulate triple system".format(sph_code.__name__)
    hydro = new_hydro(sph_code, sph_giant, core, t_end, n_steps, core_radius)
    print "\nSetting up {0} to simulate triple system".format(dynamics_code.__name__)
    binary_system = new_dynamics_for_binary(dynamics_code, binary)
    print "\nSetting up Bridge to simulate triple system"
    coupled_system = new_coupled_system(hydro, binary_system, t_end, n_steps)
    
    print "\nEvolving to {0}".format(t_end)
    evolve_system(coupled_system, t_end, n_steps)
    print "Done"
