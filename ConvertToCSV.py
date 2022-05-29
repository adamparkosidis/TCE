import sys
from amuse.io import write_set_to_file
from amuse.datamodel import Particles
from DensityPlotting import SphGiant


def main(args):
    if len(args) <= 2: #dont have the path
        print("no path received")
    else:
        if len(args) > 2:#opposite case
            giant = SphGiant(args[1], args[2],opposite=True)
        else:
            giant = SphGiant(args[1], args[2])
        particles= giant.gasParticles
        particles.add_particle(giant.core)
        print("writing csv file to ", args[1].replace("gas","giant",1))
        write_set_to_file(particles, args[1].replace("gas", "giant", 1))

if __name__ == "__main__":
    main(sys.argv)
