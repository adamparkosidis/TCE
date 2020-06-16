import os
import os.path
import shutil
import math
import pickle

from amuse.units import units, constants, nbody_system
from amuse.units import *
from amuse.datamodel import Particles, Particle
def get_relative_velocity_at_apastron(total_mass, semimajor_axis, ecc):
    return (constants.G * total_mass * ((1.0 - ecc)/(1.0 + ecc)) / semimajor_axis).sqrt()

def VectorCross(v1,v2):
     vx = v1[1]*v2[2]-v1[2]*v2[1]
     vy = v1[2]*v2[0]-v1[0]*v2[2]
     vz = v1[0]*v2[1]-v1[1]*v2[0]
     v = v1
     v[0] = vx
     v[1] = vy
     v[2] = vz
     return v

def SqalarMul(v1,v2):
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]

def CalculateSemiMajor(V,R,M):#as in 2.134 in Solar System Dynamics bd C.D Murray and S.F Dermott
    #print "R: ", R ," V: ",V," M: ",M
    a = 1/((2/CalculateVectorSize(R))-(CalculateVectorSize(V)**2)/(constants.G*M))
    return abs(a)

def CalculateSpecificMomentum(V,R):
    return VectorCross(R,V)

def CalculateEccentricity(particle1,particle2):
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
    #print eAttitude1, eAttitude2
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

def CalculateBinaryParameters(particle, sphGiant):
    sphGiant.CalculateInnerSPH(particle)
    innerMass = sphGiant.innerGas.mass.value_in(units.MSun)
    binaryVelocityDifference = CalculateVelocityDifference(particle, sphGiant.innerGas)
    binarySeparation = CalculateSeparation(particle, sphGiant.innerGas)
    binaryMass = particle.mass + sphGiant.innerGas.mass

    aOuter1 = CalculateSemiMajor(binaryVelocityDifference, binarySeparation, binaryMass)
    eOuter1 = CalculateEccentricity(particle, sphGiant.innerGas)
    return innerMass, aOuter1.value_in(units.AU), eOuter1, CalculateVectorSize(binarySeparation).value_in(units.RSun)

def CalculatePotentialEnergy(particle1, particle2):
    return -constants.G * particle1.mass * particle2.mass \
            / CalculateVectorSize(CalculateSeparation(particle1, particle2))


def CalculateOmega(particles):
    r =  ((particles.x*particles.x+particles.y*particles.y+particles.z*particles.z)**0.5).reshape((-1,1))
    return -0.5 * CalculateVectorSize((particles.mass.reshape((-1, 1)) * (particles.position.cross(particles.velocity)/r)**2).sum(axis=0))
