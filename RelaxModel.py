from amuse.community.gadget2.interface import Gadget2
from amuse.units import *
from amuse.units.units import *
from amuse.units import units , nbody_system
from amuse.plot import plot, native_plot, sph_particles_plot

class RelaxedModel:
    def __init__(self, totalMass, semmiMajor, gasParticles, dmParticles, endTime=10000 | units.yr, timeSteps=7):
        # creating the NBody system with the 3
        nbody = nbody_system.nbody_to_si(totalMass, semmiMajor)

        # evolve
        evolutionCode = Gadget2(nbody, number_of_workers=2)
        for gasParticle in gasParticles:
            evolutionCode.gas_particles.add_particles(gasParticle)
        for dmParticle in dmParticles:
            evolutionCode.dm_particles.add_particle(dmParticle)

        parts = evolutionCode.gas_particles.copy()
        centerOfMassRadius = parts.center_of_mass()
        centerOfMassV = parts.center_of_mass_velocity()

        print "starting SPH relaxation"
        native_plot.figure(figsize=(10, 10), dpi=60)
        timeStep = endTime / timeSteps
        currentTime = 0.0 | units.Myr
        while currentTime < endTime:
            evolutionCode.evolve_model(currentTime)
            print "current time = ", evolutionCode.model_time.as_quantity_in(units.yr)
            print "radius = ", evolutionCode.particles.radius[:-1]
            currentTime += timeStep
            parts = evolutionCode.gas_particles.copy()
            sph_particles_plot(parts)
            # native_plot.show()
            #native_plot._imsave(string.format("relax_{0}", currentTime))

            f = ((currentTime/timeStep) * 1.0)/timeSteps
            parts.add_particle(evolutionCode.dm_particles)
            evolutionCode.gas_particles.position += centerOfMassRadius - parts.center_of_mass()
            evolutionCode.dm_particles.position += centerOfMassRadius - parts.center_of_mass()
            evolutionCode.gas_particles.velocity = f * (evolutionCode.gas_particles.velocity - parts.center_of_mass_velocity()) + centerOfMassV
            evolutionCode.dm_particles.velocity = f * (evolutionCode.dm_particles.velocity - parts.center_of_mass_velocity()) + centerOfMassV

        # TODO: should I change the center of mass?

        evolutionCode.stop()
        self.gas_particles = evolutionCode.gas_particles.copy()
        self.dm_particles = evolutionCode.dm_particles.copy()
        native_plot.figure(figsize=(10, 10), dpi=60)
        sph_particles_plot(self.gas_particles)
        native_plot.show()