import os
import os.path
import sys
import threading
import shutil
import math
import pickle
import gc
import h5py

import matplotlib
matplotlib.use('Agg')
from amuse.units import units, constants, nbody_system
from amuse.units import *

from amuse.lab import *
from amuse.units.quantities import AdaptingVectorQuantity
from amuse.datamodel import Particles, Particle
from amuse.io import write_set_to_file, read_set_from_file

from amuse.plot import scatter, xlabel, ylabel, plot, pynbody_column_density_plot, HAS_PYNBODY, _smart_length_units_for_pynbody_data, convert_particles_to_pynbody_data, UnitlessArgs

from matplotlib import pyplot
import pynbody
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
        self.core = read_set_from_file(dm_particles_file, format='amuse')[0]
        self.gas = Star(None, None)
        self.gas.mass = self.gasParticles.total_mass()
        self.gas.position = self.gasParticles.center_of_mass()
        self.gas.x , self.gas.y, self.gas.z = self.gasParticles.center_of_mass()
        self.gas.vx, self.gas.vy, self.gas.vz = self.gasParticles.center_of_mass_velocity()
        self.gas.v = (self.gas.vx, self.gas.vy, self.gas.vz)
        print "gas position: ", self.gas.position, " core position: ", self.core.position
        self.mass = self.gas.mass + self.core.mass
        print "core mass: ",self.core.mass.as_quantity_in(units.MSun)," gas mass: ", self.gas.mass.as_quantity_in(units.MSun), " total star mass: ", self.mass.as_quantity_in(units.MSun)
        self.x , self.y, self.z = (self.gas.position * self.gas.mass + self.core.position * self.core.mass) / self.mass
        #self.x , self.y, self.z = self.core.position
        self.vx, self.vy, self.vz = (self.gas.v * self.gas.mass + (self.core.vx, self.core.vy, self.core.vz) * self.core.mass) / self.mass
        #self.vx, self.vy, self.vz =  self.core.vx, self.core.vy, self.core.vz
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
    print load
    stars = Particles(2, particles= [load[0], load[1]])
    return stars

def get_relative_velocity_at_apastron(total_mass, semimajor_axis, ecc):
    return (constants.G * total_mass * ((1.0 - ecc)/(1.0 + ecc)) / semimajor_axis).sqrt()

def VectorCross(v1,v2):
    return (v1[1]*v2[2]-v1[2]*v2[1],v1[2]*v2[0]-v1[0]*v2[2],v1[0]*v2[1]-v1[1]*v2[0])

def SqalarMul(v1,v2):
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]

def CalculateSemiMajor(V,R,M):#as in 2.134 in Solar System Dynamics bd C.D Murray and S.F Dermott
    print "R: ", R ," V: ",V," M: ",M
    a = 1/((2/CalculateVectorSize(R))-(CalculateVectorSize(V)**2)/(constants.G*M))
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
    #TODO: should maybe take into acount only the inner mass?
    return -constants.G*(particle1.mass+particle2.mass)/(CalculateVectorSize(R)) + 0.5*(CalculateVectorSize(V)**2)


def CalculateVectorSize(vector):
    returnValue = (vector[0]**2 + vector[1]**2 + vector[2]**2)**0.5
    return returnValue

def GetVelocitySize(particle):
    return  CalculateVectorSize((particle.vx,particle.vy,particle.vz))

def CalculateSeparation(particle1, particle2):
    x = (particle1.x - particle2.x).as_quantity_in(units.AU)
    y = (particle1.y - particle2.y).as_quantity_in(units.AU)
    z = (particle1.z - particle2.z).as_quantity_in(units.AU)
    return (x,y,z)

def CalculateVelocityDifference(particle1, particle2):
    vx = particle1.vx - particle2.vx
    vy = particle1.vy - particle2.vy
    vz = particle1.vz - particle2.vz
    return (vx,vy,vz)

def GetPositionSize(particle):
    return CalculateVectorSize((particle.x ,particle.y,  particle.z))

def PlotDensity(sphGiant,core,binary,i, outputDir, vmin, vmax):
    if not HAS_PYNBODY:
        print "problem plotting"
        return
    pynbody_column_density_plot(sphGiant ,resolution=1000, width=2.0|units.AU,vmin= vmin, vmax= vmax,cmap= "hot")
    scatter(core.x, core.y, c="r")
    scatter(binary.x, binary.y, c="w")
    pyplot.savefig(outputDir + "/plotting_{0}.jpg".format(i))
    pyplot.close()

def PlotVelocity(sphGiant,core,binary,step, outputDir, vmin, vmax):
    if not HAS_PYNBODY:
        print HAS_PYNBODY
        print "problem plotting"
        return
    
    width = 1.7 * sphGiant.position.lengths_squared().amax().sqrt()
    length_unit, pynbody_unit = _smart_length_units_for_pynbody_data(width)
    pyndata = convert_particles_to_pynbody_data(sphGiant, length_unit, pynbody_unit)
    UnitlessArgs.strip([1]|length_unit, [1]|length_unit)
    units = 'm_p cm^-2'
    pynbody_sph.velocity_image(pyndata, width=width.value_in(length_unit), units=units,vmin= vmin, vmax= vmax)
    UnitlessArgs.current_plot = native_plot.gca()
    scatter(core.x, core.y, c="r")
    scatter(binary.x, binary.y, c="w")
    pyplot.savefig(outputDir + "/velocity/velocity_plotting_{0}.jpg".format(step))
    pyplot.close()

'''
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
'''
def Plot1Axe(x, fileName, outputDir, beginTime = 0):
    if len(x) == 0:
        return
    timeStep = 1400.0/7000.0
    timeLine = [beginTime + time * timeStep for time in xrange(len(x))] | units.day
    native_plot.figure(figsize= (20, 20), dpi= 80)
    plot(timeLine,x)
    xlabel('time[days]')
    native_plot.savefig(outputDir + '/' + fileName + '.jpg')

def PlotSemiMajorAxis(semimajors, outputDir, beginTime = 0):
    for a in semimajors:
        if a[0]:
            Plot1Axe(a[0], a[1], outputDir)

def PlotEccentricity(eccentricities, outputDir, beginTime = 0):
    for e in eccentricities:
        if e[0] != []:
            Plot1Axe(e[0], e[1], outputDir)

def PlotBinaryDistance(distances, outputDir, beginTime = 0):
    for d in distances:
        if d[0]:
            Plot1Axe(d[0], d[1], outputDir)

def AnalyzeBinaryChunk(savingDir,gasFiles,dmFiles,outputDir,chunk, vmin, vmax):
    for i in chunk:
        print "step #",i
        gas_particles_file = os.path.join(os.getcwd(), savingDir,gasFiles[i])
        dm_particles_file = os.path.join(os.getcwd(),savingDir, dmFiles[i])
        binaryDistances = AdaptingVectorQuantity()
        semmimajors = AdaptingVectorQuantity()
        eccentricities = []
        sphGiant = SphGiant(gas_particles_file, dm_particles_file)
        try:
            binary = LoadBinaries(dm_particles_file)
            companion = binary[1]
        except: #no binary
            binary = []
            companion = sphGiant
        #print binary
        if len(binary) > 1:
            isBinary= True
            binary = Star(companion, sphGiant)
        else:
            isBinary=False
            binary = Star(sphGiant, sphGiant)

        if CalculateVectorSize(CalculateSeparation(sphGiant.core,companion)) < sphGiant.core.radius * 2:
            print "merger between companion and the giant!"
            #break

        if isBinary:
            #check if the companion is inside, take into account only the inner mass of the companion's orbit
            sphGiant.CalculateInnerSPH(companion)
            print "innerGasMass: ", sphGiant.innerGas.mass.value_in(units.MSun)

            #check if the binary is breaking up
            if binary.specificEnergy > 0 | (units.m **2 / units.s **2):
                print "binary is breaking up", binary.specificEnergy

            newBinaryVelocityDifference = CalculateVelocityDifference(companion, sphGiant.innerGas)
            newBinarySeparation = CalculateSeparation(companion, sphGiant.innerGas)
            newBinaryMass = companion.mass + sphGiant.innerGas.mass
            newBinarySpecificEnergy = CalculateSpecificEnergy(newBinaryVelocityDifference,newBinarySeparation,sphGiant.innerGas,companion)
            semmimajor = CalculateSemiMajor(newBinaryVelocityDifference, newBinarySeparation, newBinaryMass).as_quantity_in(units.AU)
            eccentricity = CalculateEccentricity(companion, sphGiant.innerGas, semmimajor)
            semmimajors.append(semmimajor)
            eccentricities.append(eccentricity)
            binaryDistances.append(CalculateVectorSize(newBinarySeparation).as_quantity_in(units.RSun))

            #check if the binary is breaking up
            if newBinarySpecificEnergy > 0 | (units.m **2 / units.s **2):
                print "binary is breaking up", binary.specificEnergy
                #break

            print "semmimajor: ",semmimajor," eccentricity: ",eccentricity, "separation: ",CalculateVectorSize(newBinarySeparation).as_quantity_in(units.RSun)

            print newBinarySeparation
            PlotDensity(sphGiant.gasParticles,sphGiant.core,companion, i , outputDir, vmin, vmax)
            PlotVelocity(sphGiant.gasParticles,sphGiant.core,companion, i, outputDir, vmin, vmax)

            binaryDistancesFile = open(savingDir+"InnerBinaryDistances{0}".format(chunk[0]),"w")
            binaryDistancesFile.writelines(binaryDistances)
            binaryDistancesFile.close()
            semmimajorsFile = open(savingDir+"aInners{0}".format(chunk[0]),"w")
            semmimajorsFile.writelines(semmimajors)
            semmimajorsFile.close()
            eccentricitiesFile = open(savingDir+"eInners{0}".format(chunk[0]),"w")
            eccentricitiesFile.writelines(eccentricities)
            eccentricitiesFile.close()

            #close opened handles
            for f in [obj for obj in gc.get_objects() if isinstance(obj,h5py.File)]:
                try:
                    f.close()
                except:
                    pass

def AnalyzeBinary(beginStep, lastStep, dmFiles, gasFiles, savingDir, outputDir, vmin, vmax):
    binaryDistances = AdaptingVectorQuantity()
    semmimajors = AdaptingVectorQuantity()
    eccentricities = []
    inclinations = []
    separationTime = 0
    if lastStep == 0 : # no boundary on last step
        lastStep = len(dmFiles)
    print lastStep
    chunkSize= 5
    chunks = [xrange(i,i+chunkSize) for i in xrange(beginStep,lastStep,chunkSize)]
    chunks[-1]= xrange(int((lastStep-beginStep)/chunkSize)*chunkSize,int((lastStep-beginStep)/chunkSize)*chunkSize +
                       lastStep-int((lastStep-beginStep)/chunkSize)*chunkSize)
    threads = []
    for chunk in chunks[1:]:
        threads.append(threading._start_new_thread(AnalyzeBinaryChunk,[savingDir,gasFiles,outputDir,chunk, vmin, vmax]))
    AnalyzeBinary(savingDir,gasFiles,dmFiles,outputDir,chunks[0], vmin, vmax)

    for chunk in chunks: #combine them
        binaryDistancesFile = open(savingDir + "/InnerBinaryDistances{0}".format(chunk[0]),"r")
        binaryDistances.append(binaryDistancesFile.readlines())
        binaryDistancesFile.close()
        os.remove(binaryDistancesFile)
        semmimajorsFile = open(savingDir + "/aInners{0}".format(chunk[0]),"r")
        semmimajors.append(semmimajorsFile.readlines())
        semmimajorsFile.close()
        os.remove(semmimajorsFile)
        eccentricitiesFile = open(savingDir + "/eInners{0}".format(chunk[0]),"r")
        eccentricities.append(eccentricitiesFile.readlines())
        eccentricitiesFile.close()
        os.remove(eccentricitiesFile)

    PlotBinaryDistance([(binaryDistances, "InnerBinaryDistances")], outputDir + "/graphs")
    PlotSemiMajorAxis([(semmimajors,"aInners")], outputDir+"/graphs")
    PlotEccentricity([(eccentricities, "eInners")], outputDir + "/graphs")



def AnalyzeTriple(beginStep, lastStep, dmFiles, gasFiles, savingDir, outputDir, vmin, vmax ):
    particle2x = []
    particle1x =  []
    corex = []
    particle1y = []
    particle2y = []
    corey = []
    particle1z = []
    particle2z = []
    corez = []
    binaryDistances = AdaptingVectorQuantity()
    triple1Distances = AdaptingVectorQuantity()
    triple2Distances = AdaptingVectorQuantity()
    aInners = AdaptingVectorQuantity()
    aOuters = AdaptingVectorQuantity()
    aOuters1 = AdaptingVectorQuantity() # for the couple particle1 + giant
    aOuters2 = AdaptingVectorQuantity() # for the couple particle2 + giant

    eInners = []
    eOuters = []
    eOuters1 = [] # for the couple particle1 + giant
    eOuters2 = [] # for the couple particle2 + giant

    inclinations = []
    separationTime = 0
    if lastStep == 0 : # no boundary on last step
        lastStep = len(dmFiles)
    print lastStep
    for i in xrange(beginStep,lastStep,5):
        print "step #",i
        gas_particles_file = os.path.join(os.getcwd(), savingDir,gasFiles[i])
        dm_particles_file = os.path.join(os.getcwd(),savingDir, dmFiles[i])

        sphGiant = SphGiant(gas_particles_file, dm_particles_file)

        #binary = Particles(2,pickle.load(open(os.path.join(os.getcwd(),savingDir,"binary.p"),"rb")))
        binary = LoadBinaries(dm_particles_file)
        #print binary

        particle1 , particle2 = binary[1] , binary[2]

        innerBinary = Star(particle1,particle2)
        triple1 = Star(particle1, sphGiant)
        triple2 = Star(particle2, sphGiant)

        aInner = CalculateSemiMajor(innerBinary.velocityDifference,innerBinary.separation, innerBinary.mass)
        eInner = CalculateEccentricity(particle1,particle2, aInner)

        if CalculateVectorSize(innerBinary.separation) < 13.0*(10**8) | units.m: #if its closer than 2 solar radiuses
            print "merger between the inner binary!" , innerBinary.separation.as_quantity_in(units.RSun)

        if CalculateVectorSize(CalculateSeparation(sphGiant.core,particle1)) < sphGiant.core.radius * 2:
            print "merger between particle1 and the giant!"
            #break

        if CalculateVectorSize(CalculateSeparation(sphGiant.core, particle2)) < sphGiant.core.radius * 2:
            print "merger between particle 2 and the giant!"
            #break

        #check if the binry is breaking up
        if innerBinary.specificEnergy > 0 | (units.m **2 / units.s **2):
            print "binary is breaking up", innerBinary.specificEnergy

            #check if the couple particle1 + giant are breaking up
            if triple1.specificEnergy > 0 | (units.m **2 / units.s **2):
                print "triple1 is breaking up", triple1.specificEnergy
                if separationTime == 0:
                    separationTime = i * 1400/7000
                #check if the couple particle2 + giant are also breaking up
                if triple2.specificEnergy > 0 | (units.m **2 / units.s **2):
                    print "triple2 is also breaking up", triple2.specificEnergy
                    break
                else:
                    sphGiant.CalculateInnerSPH(particle2)
                    print "innerGasMass: ", sphGiant.innerGas.mass.value_in(units.MSun)
                    newBinaryVelocityDifference = CalculateVelocityDifference(particle2, sphGiant.innerGas)
                    newBinarySeparation = CalculateSeparation(particle2, sphGiant.innerGas)
                    newBinaryMass = particle2.mass + sphGiant.innerGas.mass

                    aOuter1 = CalculateSemiMajor(newBinaryVelocityDifference, newBinarySeparation, newBinaryMass)
                    eOuter1 = CalculateEccentricity(particle2, sphGiant.innerGas, aOuter1)
                    aOuters1.append(aOuter1)
                    eOuters1.append(eOuter1)
                    triple2Distances.append(CalculateVectorSize(newBinarySeparation))
                    print newBinarySeparation

            #check if the couple particle2 + giant are breaking up
            if triple2.specificEnergy > 0 | (units.m **2 / units.s **2):
                print "triple2 is breaking up", triple2.specificEnergy
                if separationTime == 0:
                    separationTime = i * 1400/7000
                #check if the couple particle1 + giant are also breaking up
                if triple1.specificEnergy > 0 | (units.m **2 / units.s **2):
                    print "triple1 is also breaking up", triple1.specificEnergy
                    break

                else:
                    sphGiant.CalculateInnerSPH(particle1)
                    print "innerGasMass: ", sphGiant.innerGas.mass.value_in(units.MSun)
                    newBinaryVelocityDifference = CalculateVelocityDifference(particle1, sphGiant.innerGas)
                    newBinarySeparation = CalculateSeparation(particle1, sphGiant.innerGas)
                    newBinaryMass = particle1.mass + sphGiant.innerGas.mass

                    aOuter2 = CalculateSemiMajor(newBinaryVelocityDifference, newBinarySeparation, newBinaryMass)
                    eOuter2 = CalculateEccentricity(particle1, sphGiant.innerGas, aOuter2)
                    aOuters2.append(aOuter2)
                    eOuters2.append(eOuter2)
                    triple1Distances.append(CalculateVectorSize(newBinarySeparation))

            else:#all the three are connected
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

        PlotDensity(sphGiant.gasParticles,sphGiant.core,binary,i + beginStep, outputDir, vmin, vmax)
        PlotVelocity(sphGiant.gasParticles,sphGiant.core,binary,i + beginStep, outputDir, vmin, vmax)

        #close opened handles
        for f in [obj for obj in gc.get_objects() if isinstance(obj,h5py.File)]:
            try:
                f.close()
            except:
                pass

    PlotBinaryDistance([(binaryDistances, "InnerBinaryDistances"), (triple1Distances, "triple1Distances"), (triple2Distances, "triple2Distances")], outputDir + "/graphs")
    PlotSemiMajorAxis([(aInners,"aInners"),(aOuters, "aOuters")], outputDir+"/graphs")
    PlotSemiMajorAxis([(aOuters1, "aOuters1"), (aOuters2, "aOuters2")], outputDir+ "/graphs", separationTime)
    PlotEccentricity([(eInners, "eInners"), (eOuters, "eOuters")], outputDir + "/graphs")
    PlotEccentricity([(eOuters1, "eOuters1"), (eOuters2, "eOuters2")],outputDir + "/graphs", separationTime)
    Plot1Axe(inclinations,"inclinations", outputDir+"/graphs")
    #PlotOrbitalParameters(particle1x,particle1y, "tempPics/orbitalsParticle1.jpg")

def GetArgs(args):
    if len(args) > 1:
        directory=args[1]
    else:
        directory = args[0]
    if len(args) > 2:
        savingDir = directory + "/" + args[2]
        if args[2] == "snapshots":
            toCompare = False
        else:
            toCompare = True
    else:
        savingDir = directory + "/evolution"
        toCompare = True
    if len(args) > 3:
        beginStep = int(args[3])
    else:
        beginStep = 0
    if (args) > 4:
        lastStep = int(args[4])
    else:
        lastStep = 0
    if len(args) > 5:
        vmin= float(args[5])
    else:
        vmin = 1e16
    if len(args) > 6:
        vmax = float(args[6])
    else:
        vmax= 1e34
    outputDir = savingDir + "/pics"
    return savingDir, toCompare, beginStep, lastStep, vmin, vmax, outputDir

def InitializeSnapshots(savingDir, toCompare=False):
    '''
    taking the snapshots directory of past run
    Returns: sorted dm snapshots and gas snapshots

    '''
    snapshots = os.listdir(os.path.join(os.getcwd(),savingDir))
    numberOfSnapshots = len(snapshots) / 2
    dmFiles = []
    gasFiles = []
    for snapshotFile in snapshots:
        if 'dm' in snapshotFile: #if the word dm is in the filename
            dmFiles.append(snapshotFile)
        if 'gas' in snapshotFile:
            gasFiles.append(snapshotFile)
    if toCompare:
        dmFiles.sort(cmp=compare)
        gasFiles.sort(cmp= compare)
    else:
        dmFiles.sort()
        gasFiles.sort()
    numberOfCompanion = 0
    if len(dmFiles) > 0:
        numberOfCompanion = len(read_set_from_file(os.path.join(os.getcwd(), savingDir,dmFiles[0]), format='amuse'))
    return gasFiles, dmFiles, numberOfCompanion

def compare(st1, st2):
    num1 = int(st1.split("_")[1].split(".")[0])
    num2 = int(st2.split("_")[1].split(".")[0])
    if num1 < num2:
        return -1
    return 1


def main(args= ["../../BIGDATA/code/amuse-10.0/runs200000/run_003","evolution",0,1e16,1e34]):
    savingDir, toCompare, beginStep, lastStep, vmin, vmax, outputDir = GetArgs(args)
    print "plotting pics to " +  outputDir +  " from " +  savingDir +" begin step = " , beginStep , " vmin, vmax = " , vmin, vmax, "special comparing = ", toCompare
    try:
        os.makedirs(outputDir)
    except(OSError):
        pass
    try:
        os.makedirs(outputDir + "/velocity")
    except(OSError):
        pass
    try:
        os.makedirs(outputDir + "/graphs")
    except (OSError):
        pass
    gasFiles, dmFiles, numberOfCompanion = InitializeSnapshots(savingDir, toCompare)

    if numberOfCompanion <= 2: #binary
        print "analyzing binary"
        AnalyzeBinary(beginStep, lastStep, dmFiles, gasFiles, savingDir, outputDir, vmin, vmax)
    elif numberOfCompanion ==3: #triple
        AnalyzeTriple(beginStep, lastStep, dmFiles, gasFiles, savingDir, outputDir, vmin, vmax)

if __name__ == "__main__":
    main(sys.argv)
