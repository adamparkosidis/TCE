import os
import os.path

from amuse.units import units
from amuse.datamodel import Particle

from amuse.ext.star_to_sph import pickle_stellar_model
from amuse.community.mesa.interface import MESA as stellar_evolution_code

from xiTau_parameters import triple_parameters


def evolve_giant(giant, stop_radius):
    stellar_evolution = stellar_evolution_code()
    giant_in_code = stellar_evolution.particles.add_particle(giant)
    while (giant_in_code.radius < 0.7 | units.AU):
        giant_in_code.evolve_one_step()

    
    print "Giant starts to ascend the giant branch, now saving model every step..."
    print giant_in_code.as_set()
    
    i = 0
    while (giant_in_code.radius < stop_radius):
        giant_in_code.evolve_one_step()
        print giant_in_code.radius, giant_in_code.age
        pickle_file_name = "./model_{0:=04}_".format(i) + "%0.1f"%(giant_in_code.radius.value_in(units.AU))
        pickle_stellar_model(giant_in_code, pickle_file_name)
        i += 1
    

if __name__ == "__main__":
    model_directory = os.path.join("../../../../../BIGDATA/code/amuse-10.0", "giant_models")
    if not os.path.exists(model_directory):
        os.mkdir(model_directory)
    os.chdir(model_directory)
    
    giant = Particle(mass = triple_parameters["mass_out"])
    
    print "\nEvolving with", stellar_evolution_code.__name__
    evolve_giant(giant, 1.0 | units.AU)
    print "Done"
