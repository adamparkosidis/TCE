from amuse.community.gadget2.interface import Gadget2
from amuse.units import *
from amuse.units.units import *
from amuse.units import units , nbody_system
from amuse.plot import plot, native_plot, sph_particles_plot

class RelaxedModel:
    def __init__(self, totalMass, semmiMajor, gasParticles, dmParticles, endTime=10000 | units.yr, timeSteps=7, savedVersionPath = ""):
        # creating the NBody system with the 3
        nbody = nbody_system.nbody_to_si(totalMass, semmiMajor)

        # evolve
        evolutionCode = Gadget2(nbody, number_of_workers=7)
        evolutionCode.parameters.time_limit_cpu = 1000000 | units.s
        for gasParticle in gasParticles:
            evolutionCode.gas_particles.add_particles(gasParticle)
        for dmParticle in dmParticles:
            evolutionCode.dm_particles.add_particle(dmParticle)

        gas = evolutionCode.gas_particles.copy()
        dm = evolutionCode.dm_particles.copy()
        centerOfMassRadius = gas.center_of_mass()
        centerOfMassV = gas.center_of_mass_velocity()

        print "starting SPH relaxation"
        native_plot.figure(figsize=(30, 30), dpi=60)
        timeStep = endTime / timeSteps
        currentTime = 0.0 | units.Myr
        currentStep = 0
        while currentTime < endTime:
            evolutionCode.evolve_model(currentTime)
            print "current time = ", evolutionCode.model_time.as_quantity_in(units.yr)
            currentTime += timeStep
            gas = evolutionCode.gas_particles.copy()
            sph_particles_plot(gas)
            native_plot.savefig(savedVersionPath + "/pics/relax_{0}".format(currentTime))
            gas.add_particle(evolutionCode.dm_particles)
            evolutionCode.gas_particles.position += (centerOfMassRadius - gas.center_of_mass())
            evolutionCode.dm_particles.position += (centerOfMassRadius - gas.center_of_mass())
            relaxingVFactor = (currentStep / timeSteps)
            evolutionCode.gas_particles.velocity = relaxingVFactor * (evolutionCode.gas_particles.velocity -
                                                                                gas.center_of_mass_velocity()) + centerOfMassV
            evolutionCode.dm_particles.velocity = relaxingVFactor * (evolutionCode.dm_particles.velocity -
                                                                               gas.center_of_mass_velocity()) + centerOfMassV
            currentStep += 1
            gas = evolutionCode.gas_particles.copy()
            dm = evolutionCode.dm_particles.copy()

        evolutionCode.stop()
        self.gas_particles = gas
        self.core_particle = dm