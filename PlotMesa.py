import os
import os.path
import sys
import threading
import multiprocessing
import shutil
import math
import pickle
import gc
import h5py
from glob import glob
import matplotlib
matplotlib.use('Agg')
from amuse.units import units, constants, nbody_system
from amuse.units import *

from amuse.lab import *
from amuse.units.quantities import AdaptingVectorQuantity
from amuse.datamodel import Particles, Particle
from amuse.ext import sph_to_star, star_to_sph
from amuse.io import write_set_to_file, read_set_from_file

from amuse.plot import scatter, xlabel, ylabel, plot, pynbody_column_density_plot, HAS_PYNBODY, _smart_length_units_for_pynbody_data, convert_particles_to_pynbody_data, UnitlessArgs, semilogx, semilogy, loglog, xlabel, ylabel

from matplotlib import pyplot
import numpy as np
import pynbody
import pynbody.plot.sph as pynbody_sph
from amuse.plot import scatter, xlabel, ylabel, plot, native_plot, sph_particles_plot

class Star:
    def __init__(self, pickleFile):
        self.pickle_file = pickleFile
        self.unpickle_stellar_structure()

    def unpickle_stellar_structure(self):
        if os.path.isfile(self.pickle_file):
            infile = open(self.pickle_file, 'rb')
        else:
            raise BaseException("Input pickle file '{0}' does not exist".format(self.pickle_file))
        structure = pickle.load(infile)
        self.mass   = structure['mass']
        self.radius = structure['radius']
        self.number_of_zones     = structure['number_of_zones']
        self.number_of_species   = structure['number_of_species']
        self.species_names       = structure['species_names']
        self.density_profile     = structure['density_profile']
        self.radius_profile      = structure['radius_profile']
        self.mu_profile          = structure['mu_profile']
        self.composition_profile = structure['composition_profile']
        self.specific_internal_energy_profile = structure['specific_internal_energy_profile']
        self.midpoints_profile   = structure['midpoints_profile']
        self.temperature =  self.specific_internal_energy_profile * self.mu_profile / (1.5 * constants.kB)
        self.pressure = (2.0/3)*self.specific_internal_energy_profile * self.density_profile
        self.sound_speed = (((5.0/3.0) * constants.kB * self.temperature /
                              self.mu_profile)**0.5).as_quantity_in(units.m / units.s)
        print len(self.specific_internal_energy_profile), len(self.temperature), len(self.sound_speed)

def mu(X = None, Y = 0.25, Z = 0.02, x_ion = 0.1):
    """
    Compute the mean molecular weight in kg (the average weight of particles in a gas)
    X, Y, and Z are the mass fractions of Hydrogen, of Helium, and of metals, respectively.
    x_ion is the ionisation fraction (0 < x_ion < 1), 1 means fully ionised
    """
    if X is None:
        X = 1.0 - Y - Z
    elif abs(X + Y + Z - 1.0) > 1e-6:
        print "Error in calculating mu: mass fractions do not sum to 1.0"
    return constants.proton_mass / (X*(1.0+x_ion) + Y*(1.0+2.0*x_ion)/4.0 + Z*x_ion/2.0)



def structure_from_star(star):
    radius_profile = star.radius_profile
    density_profile = star.density_profile
    #print radius_profile
    #print star.temperature
    #TODO: check if the mass profile is correct...
    return dict(
        radius = radius_profile.as_quantity_in(units.RSun),
        density = density_profile,
        mass = star.mu_profile * star.mass,
        temperature = star.temperature,
        specific_energy= star.specific_internal_energy_profile,
        pressure = star.pressure,
        sound_speed = star.sound_speed
    )

def temperature_density_plot(outputDir, pickleFile):
    if not HAS_PYNBODY:
        print "problem plotting"
        return
    width = 5.0 | units.AU
    length_unit, pynbody_unit = _smart_length_units_for_pynbody_data(width)

    star = Star(outputDir + "/" + pickleFile)
    print star.number_of_zones
    data = structure_from_star(star)
    adding = pickleFile.split("MESA")[-1]
    figure = pyplot.figure(figsize = (8, 10))
    pyplot.subplot(1, 1, 1)
    ax = pyplot.gca()
    plotT = semilogy(data["radius"], data["temperature"].as_quantity_in(units.K), 'r-', label = r'$T(r)$')
    xlabel('Radius')
    ylabel('Temperature')
    ax.twinx()
    plotrho = semilogy(data["radius"], data["density"], 'g-', label = r'$\rho(r)$')
    plots = plotT + plotrho
    labels = [one_plot.get_label() for one_plot in plots]
    ax.legend(plots, labels, loc=3)
    ylabel('Density')

    #plot to file
    textFile = open(outputDir + '/radial_profile/temperature_' + adding + '.txt', 'w')
    textFile.write(', '.join([str(y) for y in data["temperature"]]))
    textFile.close()
    textFile = open(outputDir + '/radial_profile/density_' + adding + '.txt', 'w')
    textFile.write(', '.join([str(y) for y in data["density"]]))
    textFile.close()
    textFile = open(outputDir + '/radial_profile/radius_' + adding +  '.txt', 'w')
    textFile.write(', '.join([str(y) for y in data["radius"]]))
    textFile.close()
    textFile = open(outputDir + '/radial_profile/sound_speed_' + adding +  '.txt', 'w')
    textFile.write(', '.join([str(y) for y in data["sound_speed"]]))
    textFile.close()
    textFile = open(outputDir + '/radial_profile/pressure_' + adding +  '.txt', 'w')
    textFile.write(', '.join([str(y) for y in data["pressure"]]))
    textFile.close()
    textFile = open(outputDir + '/radial_profile/specific_energy_' + adding +  '.txt', 'w')
    textFile.write(', '.join([str(y) for y in data["specific_energy"]]))
    textFile.close()
    #print "saved"
    pyplot.legend()
    pyplot.suptitle('Structure of a ' + adding + ' mass star')
    pyplot.savefig(outputDir + "/radial_profile/temperature_radial_profile_" + adding + ".jpg")
    pyplot.close()

    native_plot.figure()
    r, th = np.meshgrid(data["radius"][:-1000].value_in(units.RSun),np.linspace(0,2*np.pi,len(data["radius"][:-1000])))
    native_plot.subplot(projection="polar")
    native_plot.xlabel("radius")
    native_plot.ylabel("radius")
    density2D, r2 = np.meshgrid(data["density"][:-1000].value_in(units.kg / units.m**3) , data["radius"][:-1000].value_in(units.RSun))
    #print density2D
    #native_plot.pcolormesh(density2D)
    native_plot.pcolormesh(th, r, density2D)
    native_plot.plot(np.linspace(0,2*np.pi,len(data["radius"][:-1000])), r, color='k', ls='none')
    native_plot.colorbar()
    #native_plot.imshow(density, cmap="hot", interpolation="nearest", origin='lower')
    native_plot.savefig(outputDir + "/radial_profile/density_" + adding + ".jpg")
    native_plot.close()

def Plot1Axe(x, fileName, outputDir, timeStep= 1400.0/7000.0, beginTime = 0):
    if len(x) == 0:
        return
    timeLine = [beginTime + time * timeStep for time in xrange(len(x))] | units.day
    native_plot.figure(figsize= (20, 20), dpi= 80)
    plot(timeLine,x)
    xlabel('time[days]')
    native_plot.savefig(outputDir + '/' + fileName + '.jpg')
    textFile = open(outputDir + '/' + fileName + '.txt', 'w')
    textFile.write(', '.join([str(y) for y in x]))
    textFile.close()


def GetArgs(args):
    if len(args) > 1:
        directory=args[1]
    else:
        directory = args[0]

    return directory

def InitializeSnapshots(savingDir):
    '''
    taking the snapshots directory of past run
    Returns: sorted mesa pickle files

    '''
    #insideDirectories = glob(savingDir + "/*/")
    snapshots = os.listdir(os.path.join(os.getcwd(),savingDir))
    picklerFiles = []
    for snapshotFile in snapshots:
        if 'MESA' in snapshotFile: #if the word dm is in the filename
            picklerFiles.append(snapshotFile)
    picklerFiles.sort()
    return picklerFiles



def main(args= ["/BIGDATA/code/amuse-10.0/Glanz/savings/MesaModels"]):
    savingDir = GetArgs(args)
    print "plotting pics to " +  savingDir + "/radial_profile" +  " from " +  os.path.join(os.getcwd(),savingDir)
    try:
        os.makedirs(savingDir + "/radial_profile")
    except(OSError):
        print ""
    pickleMesaFiles = InitializeSnapshots(savingDir)
    print "found ",len(pickleMesaFiles) ," pickled files"
    for pickleFile in pickleMesaFiles:
        temperature_density_plot(savingDir, pickleFile)


if __name__ == "__main__":
    main(sys.argv)

