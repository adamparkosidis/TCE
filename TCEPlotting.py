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

from amuse.plot import scatter, xlabel, ylabel, plot, pynbody_column_density_plot, HAS_PYNBODY


from matplotlib import pyplot
from amuse.plot import scatter, xlabel, ylabel, plot, native_plot, sph_particles_plot


def load_sph_giant(gas_particles_file, dm_particles_file):
    sph_giant = read_set_from_file(gas_particles_file, format='amuse')
    core = read_set_from_file(dm_particles_file, format='amuse')[-1]

    return sph_giant, core


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
    return a

def CalculateSpecificMomentum(V,R):
    return VectorCross(R,V)

def CalculateEccentricity(V,R,a):
    #h = ((particle1.mass * GetVelocitySize(particle1)) + (particle2.mass * GetVelocitySize(particle2)))/ miu
    h = CalculateSpecificMomentum(V,R)
    hSize = CalculateVectorSize(h)
    e =  (1 - ((hSize**2)/(constants.G*(particle1.mass + particle2.mass)*a)))**0.5
    return e

def CalculateInclination(V12,R12,V23,R23):
    h12 = CalculateSpecificMomentum(V12,R12)
    h23 = CalculateSpecificMomentum(V23,R23)
    return math.degrees(math.acos(SqalarMul(h12,h23)/(CalculateVectorSize(h12)*CalculateVectorSize(h23))))
        

def CalculateSpecificEnergy(V,R,particle1,particle2):
    return -constants.G*(particle1.mass+particle2.mass)/(CalculateVectorSize(R)) + 0.5*(CalculateVectorSize(V)**2)

def PlotDensity(sph_giant,core,binary,step ,plottingPath = "/tempPics"):
    if not HAS_PYNBODY:
        print "problem plotting"
        #return
    pynbody_column_density_plot(sph_giant, width=5|units.AU)
    scatter(core.x, core.y, c="r")
    scatter(binary.x, binary.y, c="w")
    pyplot.savefig(plottingPath + "/DensityPlot_{0}.jpg".format(step))
    pyplot.close()

def GetSphMass(sphStar):
    totalMass = 0.0 | units.kg
    for mass in sphStar.mass:
        totalMass += mass
    return totalMass

def GetSphVelocity(sphStar):
    vxTot , vyTot , vzTot = ( 0.0 , 0.0, 0.0 )| units.m * units.s**-1
    for vx in sphStar.vx:
        vxTot += vx

    for vy in sphStar.vy:
        vyTot += vy

    for vz in sphStar.vz:
        vzTot += vz

    return vx, vy, vz

def GetSphPosition(sphStar):
    xTot , yTot , zTot = ( 0.0 , 0.0, 0.0 )| units.m
    for x in sphStar.x:
        xTot += x

    for y in sphStar.y:
        yTot += y

    for z in sphStar.z:
        zTot += z

    return x, y, z

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
    vx =  particle1.vx - particle2.vx
    vy = particle1.vy - particle2.vy
    vz = particle1.vz - particle2.vz
    return (vx,vy,vz)

class Star:
    def __init__(self, particle1,particle2):
        if particle1 != None and particle2 != None:
            self.Star(particle1,particle2)
        else:
            self.position = 0.0 | units.m
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
        self.position = CalculateVectorSize((self.x,self.y,self.z))
        self.mass  = particle1.mass + particle2.mass
        self.velocityDifference = CalculateVelocityDifference(particle1,particle2)
        self.separation = CalculateSeparation(particle1,particle2)
        self.specificEnergy = CalculateSpecificEnergy(self.velocityDifference,self.separation,particle1,particle2)

def GetPositionSize(particle):
    return CalculateVectorSize((particle.x ,particle.y,  particle.z))

def CheckLeavingParticle(sphStar,core,mass):
    for particle in xrange(len(sphStar.vx)):
        energy = -constants.G*(sphStar.mass[particle]+mass)/(CalculateVectorSize((sphStar.x[particle],sphStar.y[particle],sphStar.z[particle])))\
                 + 0.5*(CalculateVectorSize((sphStar.vx[particle],sphStar.vy[particle],sphStar.vz[particle]))**2)
        if energy > 0 | (units.m **2 / units.s **2):
            print "particle ", particle, " has left the system"

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
    timeLine = [time*1400.0*3600*24/7000.0 for time in xrange(len(x))]
    native_plot.figure(figsize= (20, 20), dpi= 80)
    #PlotOrbitalParameters(x,steps,plotFile=fileName)
    plot(timeLine,x)
    xlabel('time[s]')
    native_plot.savefig('tempPics/' + fileName + '.jpg')

def PlotSemiMajorAxis(aInners,aOuters):
    Plot1Axe(aInners, "aInners")
    Plot1Axe(aOuters, "aOuters")

def PlotEccentricity(eInners,eOuters):
    Plot1Axe(eInners, "eInners")
    Plot1Axe(eOuters, "eOuters")

def PlotBinaryDistance(distances):
    Plot1Axe(distances,"distances")

if __name__ == "__main__":
    snapshots = os.listdir(os.path.join(os.getcwd(), "run_006/snapshots/"))
    dmFiles = []
    gasFiles = []
    print "giant: ", giant.position, "sph: ", starEnvelope
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
    eInners = []
    eOuters = []
    inclinations = []
    print len(dmFiles)
    for i in xrange(len(dmFiles)):
        gas_particles_file = os.path.join(os.getcwd(), "run_006/snapshots/", gasFiles[i])
        dm_particles_file = os.path.join(os.getcwd(), "run_006/snapshots/", dmFiles[i])
        relative_inclination = math.radians(9.0)

        sph_giant, core = load_sph_giant(gas_particles_file, dm_particles_file)
        giant = Star(None, None)
        giant.mass = GetSphMass(sph_giant) + core.mass
        giant.vx ,giant.vy, giant.vz =  GetSphVelocity(sph_giant)
        sphV = (giant.vx **2 + giant.vy**2 + giant.vz**2)**0.5
        giant.v = (sphV * GetSphMass(sph_giant) + (GetVelocitySize(core) * core.mass)) / giant.mass
        giant.x ,giant.y, giant.z  =  GetSphPosition(sph_giant)
        sphPos = (giant.x **2 + giant.y**2 + giant.z**2)**0.5


        giant.position = (sphPos * GetSphMass(sph_giant) + GetPositionSize(core) * core.mass) / giant.mass

        binary = LoadBinaries(dm_particles_file)
        particle1 , particle2 = binary[0] , binary[1]
        binaryStar = Star(particle1,particle2)

        triple1 = Star(particle1,giant)
        triple2 = Star(particle2,giant)

        aInner = CalculateSemiMajor(binaryStar.velocityDifference,binaryStar.separation, binaryStar.mass)
        eInner = CalculateEccentricity(binaryStar.velocityDifference,binaryStar.separation, aInner)

        if CalculateVectorSize(binaryStar.separation) < 13.0*(10**8) | units.m: #if its closer than 2 solar raiuses
            print "merger!" , binaryStar.separation

        if binaryStar.specificEnergy > 0 | (units.m **2 / units.s **2):
            print "binary is breaking up", binaryStar.specificEnergy

            if triple1.specificEnergy > 0 | (units.m **2 / units.s **2):
                        print "triple1 is breaking up", triple1.specificEnergy

            if triple2.specificEnergy > 0 | (units.m **2 / units.s **2):
                        print "triple2 is breaking up", triple2.specificEnergy
            continue


        tripleVCOM = (binaryStar.mass * binaryStar.v + giant.mass*giant.v)/(binaryStar.mass + giant.mass)
        tripleCOM = (binaryStar.mass * binaryStar.position + giant.mass * giant.position) / (binaryStar.mass + giant.mass)
        tripleMass = binaryStar.mass + giant.mass


        tripleVelocityDifference = CalculateVelocityDifference(binaryStar,giant)
        tripleSeparation = CalculateSeparation(binaryStar,giant)

        aOuter = CalculateSemiMajor(tripleVelocityDifference, tripleSeparation, tripleMass)
        eOuter = CalculateEccentricity((binaryStar.vx,binaryStar.vy,binaryStar.vz), tripleSeparation, aOuter)
        print "eOuter: ", eOuter

        tripleSpecificEnergy = CalculateSpecificEnergy(binaryStar.velocityDifference,binaryStar.separation,particle1,particle2)

        particle1x.append(particle1.x| units.RSun)
        particle2x.append(particle2.x| units.RSun)
        corex.append(core.x| units.RSun)
        particle1y.append(particle1.y| units.RSun)
        particle2y.append(particle2.y| units.RSun)
        corey.append(core.y| units.RSun)
        particle1z.append(particle1.z| units.RSun)
        particle2z.append(particle2.z| units.RSun)
        corez.append(core.z| units.RSun)

        binaryDistances.append(CalculateVectorSize(binaryStar.separation))

        PlotDensity(sph_giant,core,binary, step)
        aInners.append(abs(aInner))
        aOuters.append(abs(aOuter))
        if eInner < 0 or  eInner > 1:
                eInner = -1
        eInners.append(eInner)
        if eOuter < 0 or  eOuter > 1:
                eOuter = -1
        eOuters.append(eOuter)
        
        inclination = CalculateInclination(tripleVelocityDifference,tripleSeparation,binaryStar.velocityDifference,binaryStar.separation)
        print inclination
        inclinations.append(CalculateInclination(tripleVelocityDifference,tripleSeparation,binaryStar.velocityDifference,binaryStar.separation))
        #print  aOuter / aInner
    PlotBinaryDistance(binaryDistances)
    PlotSemiMajorAxis(aInners,aOuters)
    PlotEccentricity(eInners,eOuters)
    Plot1Axe(inclinations,"inclinations")
    #PlotOrbitalParameters(particle1x,particle1y, "tempPics/orbitalsParticle1.jpg")


