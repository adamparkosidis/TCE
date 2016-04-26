import os
import os.path
import shutil
import math

from amuse.units import units, constants, nbody_system
from amuse.units import *

from amuse.lab import *
from amuse.units.quantities import AdaptingVectorQuantity
from amuse.datamodel import Particles, Particle
from amuse.io import write_set_to_file, read_set_from_file

from amuse.plot import scatter, xlabel, ylabel, plot, pynbody_column_density_plot, HAS_PYNBODY, _smart_length_units_for_pynbody_data, convert_particles_to_pynbody_data, UnitlessArgs

from matplotlib import pyplot
import pynbody.plot.sph as pynbody_sph
from amuse.plot import scatter, xlabel, ylabel, plot, native_plot, sph_particles_plot


class Star:
    def __init__(self, particle1,particle2):
        if particle1 != None and particle2 != None:
            self.Star(particle1,particle2)
        else:
            self.position = (0.0, 0.0, 0.0) | units.m
            self.vx, self.vy, self.vz = (0.0 , 0.0, 0.0 ) | units.m / units.s
            self.v = 0.0 | units.m / units.s
            self.mass = 0.0 | units.kg

    def Star(self,particle1,particle2):
        self.vx = (particle1.vx* particle1.mass + particle2.vx * particle2.mass)/(particle1.mass + particle2.mass)
        self.vy = (particle1.vy* particle1.mass + particle2.vy * particle2.mass)/(particle1.mass + particle2.mass)
        self.vz = (particle1.vz* particle1.mass + particle2.vz * particle2.mass)/(particle1.mass + particle2.mass)
        self.v = CalculateVectorSize((self.vx,self.vy,self.vz))
        self.x = (particle1.x* particle1.mass + particle2.x * particle2.mass)/(particle1.mass + particle2.mass)
        self.y = (particle1.y * particle1.mass + particle2.y * particle2.mass)/(particle1.mass + particle2.mass)
        self.z = (particle1.z * particle1.mass + particle2.z * particle2.mass)/(particle1.mass + particle2.mass)
        self.position =(self.x, self.y, self.z)
        self.mass  = particle1.mass + particle2.mass
        self.velocityDifference = CalculateVelocityDifference(particle1,particle2)
        self.separation = CalculateSeparation(particle1,particle2)
        self.specificEnergy = CalculateSpecificEnergy(self.velocityDifference,self.separation,particle1,particle2)

class SphGiant:
    def __init__(self, gas_particles_file, dm_particles_file):
        self.gasParticles = read_set_from_file(gas_particles_file, format='amuse')
        self.core = read_set_from_file(dm_particles_file, format='amuse')[-1]

        self.gas = Star(None, None)
        self.gas.mass = self.gasParticles.total_mass()
        self.gas.position = self.gasParticles.center_of_mass()
        self.gas.x , self.gas.y, self.gas.z = self.gasParticles.center_of_mass()
        self.gas.vx, self.gas.vy, self.gas.vz = self.gasParticles.center_of_mass_velocity()
        self.gas.v = (self.gas.vx, self.gas.vy, self.gas.vz)

        self.mass = self.gas.mass + self.core.mass
        self.x , self.y, self.z = (self.gas.position * self.gas.mass + self.core.position * self.core.mass) / self.mass
        self.vx, self.vy, self.vz = (self.gas.v * self.gas.mass + (self.core.vx, self.core.vy, self.core.vz) * self.core.mass) / self.mass
        self.v = (self.vx, self.vy, self.vz)


    def CalculateInnerSPH(self, relativeParticle):
        self.innerGas = Star(None, None)
        radius = CalculateVectorSize(CalculateSeparation(relativeParticle, self.core))
        self.innerGas.mass = self.CalculateTotalGasMassInsideRadius(radius)
        self.innerGas.position = self.core.position
        self.innerGas.x , self.innerGas.y, self.innerGas.z = self.innerGas.position
        self.CalculateSphVelocityInsideRadius(radius)

    def CalculateTotalGasMassInsideRadius(self, radius):
        innerMass = self.core.mass
        i = 0
        for particle in self.gasParticles:
            separation = CalculateVectorSize(CalculateSeparation(particle, self.core))
            if separation < radius:
                innerMass += particle.mass
        return innerMass
    '''
    innerMass = self.core.mass
    i = 0
    for particleRadius in self.gasParticles.radius:
        #separation = CalculateVectorSize(CalculateSeparation(particle, self.core))
        print radius.as_quantity_in(units.AU), particleRadius.as_quantity_in(units.AU)
        if particleRadius < radius:
            innerMass += self.gasParticles.mass[i]
        else:
            print particleRadius
        i += 1
    return innerMass
    '''

    def CalculateSphVelocityInsideRadius(self,radius):
        self.innerGas.vxTot , self.innerGas.vyTot , self.innerGas.vzTot = ( 0.0 , 0.0, 0.0 )| units.m * units.s**-1
        particles = 0
        for particle in self.gasParticles:
            separation = CalculateVectorSize(CalculateSeparation(particle, self.core))
            if separation < radius:
                self.innerGas.vxTot += particle.vx
                self.innerGas.vyTot += particle.vy
                self.innerGas.vzTot += particle.vz
                particles += 1
        if particles > 0:
            self.innerGas.vxTot /= particles
            self.innerGas.vyTot /= particles
            self.innerGas.vzTot /= particles
        self.innerGas.v = (self.innerGas.vxTot, self.innerGas.vyTot, self.innerGas.vzTot)

    def CountLeavingParticlesInsideRadius(self):
        leavingParticles = 0
        for particle in self.gasParticles:
            specificEnergy = CalculateSpecificEnergy(CalculateVelocityDifference(particle,self.gas),CalculateSeparation(particle, self.gas), particle, self.gas)
            if specificEnergy > 0:
                leavingParticles += 1
        return leavingParticles

def LoadBinaries(file):
    load = read_set_from_file(file, format='amuse')
    stars = Particles(2, particles= [load[0], load[1]])
    return stars

def get_relative_velocity_at_apastron(total_mass, semimajor_axis, ecc):
    return (constants.G * total_mass * ((1.0 - ecc)/(1.0 + ecc)) / semimajor_axis).sqrt()

def VectorCross(v1,v2):
    return (v1[1]*v2[2]-v1[2]*v2[1],v1[2]*v2[0]-v1[0]*v2[2],v1[0]*v2[1]-v1[1]*v2[0])

def SqalarMul(v1,v2):
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]

def CalculateSemiMajor(V,R,M):#as in 2.134 in Solar System Dynamics bd C.D Murray and S.F Dermott
    a = 1/((2/CalculateVectorSize(R))-(CalculateVectorSize(V)**2)/(constants.G*(particle1.mass+particle2.mass)))
    return abs(a)

def CalculateSpecificMomentum(V,R):
    return VectorCross(R,V)

def CalculateEccentricity(particle1,particle2,a):
    #h = ((particle1.mass * GetVelocitySize(particle1)) + (particle2.mass * GetVelocitySize(particle2)))/ miu
    V= CalculateVelocityDifference(particle1,particle2)
    R = CalculateSeparation(particle1,particle2)
    h = CalculateSpecificMomentum(V,R)
    hSize = CalculateVectorSize(h)
    '''
    if a < 0.0 |units.m:
        e =  (1 + ((hSize**2)/(constants.G*(particle1.mass + particle2.mass)*a)))**0.5
    else:
        e =  (1 - ((hSize**2)/(constants.G*(particle1.mass + particle2.mass)*a)))**0.5
    '''
    miu = constants.G*(particle1.mass + particle2.mass)
    element2 = (R[0].value_in(units.m),R[1].value_in(units.m),R[2].value_in(units.m))/(CalculateVectorSize(R).value_in(units.m))
    vxh = VectorCross(V,h)
    element10 = (vxh[0].value_in(units.m**3*units.s**-2),vxh[1].value_in(units.m**3*units.s**-2),vxh[2].value_in(units.m**3*units.s**-2))
    element11 = miu.value_in(units.m**3*units.s**-2)
    element1 = element10/element11
    eAttitude1 = CalculateVectorSize(element1 - element2)
    eAttitude2 = (1 + (2 * CalculateSpecificEnergy(V,R,particle1,particle2)*(hSize**2))/(constants.G*(particle1.mass+particle2.mass))**2)**0.5
    return eAttitude1

def CalculateInclination(V12,R12,V23,R23):
    h12 = CalculateSpecificMomentum(V12,R12)
    h23 = CalculateSpecificMomentum(V23,R23)
    return math.degrees(math.acos(SqalarMul(h12,h23)/(CalculateVectorSize(h12)*CalculateVectorSize(h23))))
        

def CalculateSpecificEnergy(V,R,particle1,particle2):
    return -constants.G*(particle1.mass+particle2.mass)/(CalculateVectorSize(R)) + 0.5*(CalculateVectorSize(V)**2)


def CalculateVectorSize(vector):
    returnValue = (vector[0]**2 + vector[1]**2 + vector[2]**2)**0.5
    return returnValue

def GetVelocitySize(particle):
    return  CalculateVectorSize((particle.vx,particle.vy,particle.vz))

def CalculateSeparation(particle1, particle2):
    x = particle1.x - particle2.x
    y = particle1.y - particle2.y
    z = particle1.z - particle2.z
    return (x,y,z)

def CalculateVelocityDifference(particle1, particle2):
    vx = particle1.vx - particle2.vx
    vy = particle1.vy - particle2.vy
    vz = particle1.vz - particle2.vz
    return (vx,vy,vz)


def GetPositionSize(particle):
    return CalculateVectorSize((particle.x ,particle.y,  particle.z))

def PlotDensity(sphGiant,core,binary,i):
    if not HAS_PYNBODY:
        print "problem plotting"
        return
    pynbody_column_density_plot(sphGiant, width=5|units.AU, vmin= 1e16, vmax= 1e34)
    scatter(core.x, core.y, c="r")
    scatter(binary.x, binary.y, c="w")
    pyplot.savefig("tempPics/plotting_{0}.jpg".format(i))
    pyplot.close()

def PlotVelocity(sphGiant,core,binary,i):
    if not HAS_PYNBODY:
        print "problem plotting"
        return
    width = 2.0 * sphGiant.position.lengths_squared().amax().sqrt()
    length_unit, pynbody_unit = _smart_length_units_for_pynbody_data(width)
    pyndata = convert_particles_to_pynbody_data(sphGiant, length_unit, pynbody_unit)
    UnitlessArgs.strip([1]|length_unit, [1]|length_unit)
    units = 'm_p cm^-2'
    pynbody_sph.velocity_image(pyndata, width=width.value_in(length_unit), units=units,vmin= 1e16, vmax= 1e34)
    UnitlessArgs.current_plot = native_plot.gca()
    scatter(core.x, core.y, c="r")
    scatter(binary.x, binary.y, c="w")
    pyplot.savefig("tempPics/velocity/velocity_plotting_{0}.jpg".format(i))
    pyplot.close()

def PlotOrbitalParameters(x,y,plotFile ='dynamics.jpg' ):

    pyplot.figure(figsize= (20, 20), dpi= 80)
    pyplot.plot(x[:,0], y[:,0], 'r.')
    pyplot.plot(x[:,1], y[:,1], 'g.')
    pyplot.plot(x[:,2], y[:, 2], 'b.')
    pyplot.xlim(-20, 20)
    pyplot.ylim(-20, 20)
    pyplot.xlabel('RSun')
    pyplot.savefig(plotFile)
    pyplot.show()

def Plot1Axe(x, fileName):
    timeStep = 1400.0/7000.0
    timeLine = [time*timeStep*3600*24 for time in xrange(len(x))] | units.day
    native_plot.figure(figsize= (20, 20), dpi= 80)
    #PlotOrbitalParameters(x,steps,plotFile=fileName)
    plot(timeLine,x)
    xlabel('time[s]')
    native_plot.savefig('tempPics/' + fileName + '.jpg')

def PlotSemiMajorAxis(semimajors):
    for a in semimajors:
        if a[0]:
            Plot1Axe(a[0], a[1])

def PlotEccentricity(eccentricities):
    for e in eccentricities:
        if e[0] != []:
            Plot1Axe(e[0], e[1])

def PlotBinaryDistance(distances):
    Plot1Axe(distances, "distances")


if __name__ == "__main__":
    beginStep = 0
    snapshots = os.listdir(os.path.join(os.getcwd(), "run_008/snapshots"))
    dmFiles = []
    gasFiles = []
    for snapshotFile in snapshots:
        if 'dm' in snapshotFile: #if the word dm is in the filename
            dmFiles.append(snapshotFile)
        if 'gas' in snapshotFile:
            gasFiles.append(snapshotFile)
    dmFiles.sort()
    gasFiles.sort()
    particle1x =  []
    particle2x = []
    corex = []
    particle1y = []
    particle2y = []
    corey = []
    particle1z = []
    particle2z = []
    corez = []
    binaryDistances = AdaptingVectorQuantity()
    aInners = AdaptingVectorQuantity()
    aOuters = AdaptingVectorQuantity()
    aOuters1 = AdaptingVectorQuantity() # for the couple particle1 + giant
    aOuters2 = AdaptingVectorQuantity() # for the couple particle2 + giant

    eInners = []
    eOuters = []
    eOuters1 = [] # for the couple particle1 + giant
    eOuters2 = [] # for the couple particle2 + giant
    giantMasses = []
    particle1Masses = []
    particle2Masses = []
    inclinations = []
    print len(dmFiles)
    for i in xrange(len(dmFiles)):
        gas_particles_file = os.path.join(os.getcwd(), "run_008/snapshots", gasFiles[i])
        dm_particles_file = os.path.join(os.getcwd(), "run_008/snapshots", dmFiles[i])
        relative_inclination = math.radians(9.0)

        binary = LoadBinaries(dm_particles_file)
        particle1 , particle2 = binary[0] , binary[1]

        innerBinary = Star(particle1,particle2)

        sphGiant = SphGiant(gas_particles_file, dm_particles_file)

        triple1 = Star(particle1, sphGiant)
        triple2 = Star(particle2, sphGiant)

        aInner = CalculateSemiMajor(innerBinary.velocityDifference,innerBinary.separation, innerBinary.mass)
        eInner = CalculateEccentricity(particle1,particle2, aInner)

        if CalculateVectorSize(innerBinary.separation) < 13.0*(10**8) | units.m: #if its closer than 2 solar radiuses
            print "merger between the inner binary!" , innerBinary.separation

        if CalculateVectorSize(triple1.separation) < 13.0*(10**8) | units.m: #if its closer than 2 solar raiuses
            print "merger between particle1 and the giant!" , triple1.separation
            break

        if CalculateVectorSize(triple2.separation) < 13.0*(10**8) | units.m: #if its closer than 2 solar raiuses
            print "merger between particle 2 and the giant!" , triple2.separation
            break

        if innerBinary.specificEnergy > 0 | (units.m **2 / units.s **2):
            print "binary is breaking up", innerBinary.specificEnergy

            if triple1.specificEnergy > 0 | (units.m **2 / units.s **2):
                print "triple1 is breaking up", triple1.specificEnergy
                sphGiant.CalculateInnerSPH(particle2)
                print "innerGasMas: ", sphGiant.innerGas.mass.value_in(units.MSun)
                newBinaryVelocityDifference = CalculateVelocityDifference(particle2, sphGiant.innerGas)
                newBinarySeparation = CalculateSeparation(particle2, sphGiant.innerGas)
                newBinaryMass = particle2.mass + sphGiant.innerGas.mass

                aOuter1 = CalculateSemiMajor(newBinaryVelocityDifference, newBinarySeparation, newBinaryMass)
                eOuter1 = CalculateEccentricity(particle2, sphGiant.innerGas, aOuter1)
                aOuters1.append(aOuter1)
                eOuters1.append(eOuter1)

            if triple2.specificEnergy > 0 | (units.m **2 / units.s **2):
                print "triple2 is breaking up", triple2.specificEnergy
                sphGiant.CalculateInnerSPH(particle1)
                print "innerGasMas: ", sphGiant.innerGas.mass.value_in(units.MSun)
                newBinaryVelocityDifference = CalculateVelocityDifference(particle1, sphGiant.innerGas)
                newBinarySeparation = CalculateSeparation(particle1, sphGiant.innerGas)
                newBinaryMass = particle1.mass + sphGiant.innerGas.mass

                aOuter2 = CalculateSemiMajor(newBinaryVelocityDifference, newBinarySeparation, newBinaryMass)
                eOuter2 = CalculateEccentricity(particle1, sphGiant.innerGas, aOuter2)
                eOuter2 = CalculateEccentricity(particle1, sphGiant.innerGas, aOuter2)
                aOuters2.append(aOuter2)
                eOuters2.append(eOuter2)

        else:
            #tripleVCOM = (innerBinary.mass * innerBinary.v + sphGiant.gas.mass*sphGiant.gas.v)/(innerBinary.mass + sphGiant.gas.mass)
            #tripleCOM = (innerBinary.mass * innerBinary.position + sphGiant.gas.mass * sphGiant.gas.position) / (innerBinary.mass + sphGiant.gas.mass)
            #print "System's Center Of Mass: ", tripleCOM[0].as_quantity_in(units.AU), tripleCOM[1].as_quantity_in(units.AU), tripleCOM[2].as_quantity_in(units.AU)
            tripleMass = innerBinary.mass + sphGiant.mass


            tripleVelocityDifference = CalculateVelocityDifference(innerBinary,sphGiant.gas)
            tripleSeparation = CalculateSeparation(innerBinary,sphGiant.gas)

            aOuter = CalculateSemiMajor(tripleVelocityDifference, tripleSeparation, tripleMass)
            eOuter = CalculateEccentricity(innerBinary,sphGiant.gas, aOuter)

            inclination = CalculateInclination(tripleVelocityDifference, tripleSeparation, innerBinary.velocityDifference, innerBinary.separation)
            tripleSpecificEnergy = CalculateSpecificEnergy(innerBinary.velocityDifference, innerBinary.separation, particle1, particle2)

            binaryDistances.append(CalculateVectorSize(innerBinary.separation))

            aInners.append(aInner)
            aOuters.append(aOuter)
            eInners.append(eInner)
            eOuters.append(eOuter)
            inclinations.append(inclination)

        particle1x.append(particle1.x| units.RSun)
        particle2x.append(particle2.x| units.RSun)
        corex.append(sphGiant.core.x| units.RSun)
        particle1y.append(particle1.y| units.RSun)
        particle2y.append(particle2.y| units.RSun)
        corey.append(sphGiant.core.y| units.RSun)
        particle1z.append(particle1.z| units.RSun)
        particle2z.append(particle2.z| units.RSun)
        corez.append(sphGiant.core.z| units.RSun)


        #PlotDensity(sphGiant.gasParticles,sphGiant.core,binary,i + beginStep)
        #PlotVelocity(sphGiant.gasParticles,sphGiant.core,binary,i + beginStep)

        #print  aOuter / aInner
    PlotBinaryDistance(binaryDistances)
    PlotSemiMajorAxis([(aInners,"aInners"),(aOuters, "aOuters"), (aOuters1, "aOuters1"), (aOuters2, "aOuters2")])
    PlotEccentricity([(eInners, "eInners"), (eOuters, "eOuters"), (eOuters1, "eOuters1"), (eOuters2, "eOuters2")])

    Plot1Axe(inclinations,"inclinations")
    #PlotOrbitalParameters(particle1x,particle1y, "tempPics/orbitalsParticle1.jpg")


