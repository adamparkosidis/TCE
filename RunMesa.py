import os
import time
import argparse
import pickle
import ConfigParser
from amuse.units import units, constants
from amuse.lab import *
from amuse.community.mesa.interface import MESA
import StarModels

types = dict()
types["MS"] = (1 + 1) | units.stellar_type
types["RGB"] = (1 + 9) | units.stellar_type
types["AGB"] = (1 + 9) | units.stellar_type



class SphStar:
    def __init__(self, pointStar, configurationFile="", configurationSection="", savedMesaStarPath = "", takeSavedMesa = False,savedGas="", savedDm=""):
        print 'parsing configurations'
        parser = ConfigParser.ConfigParser()
        parser.read(configurationFile)
        self.pointStar = pointStar
        self.sphParticles = float(parser.get(configurationSection, "sphParticles"))
        self.stellar_type = types[(parser.get(configurationSection, "stellar_type"))]
        self.saving_name = parser.get(configurationSection, "name")
        try:
            self.radius = float(parser.get(configurationSection, "radius")) | units.RSun
        except:
            self.radius = 0 | units.RSun

        if takeSavedMesa:
            self.sphStar = convert_stellar_model_to_SPH(pointStar,self.sphParticles, do_relax = False, with_core_particle=True,
                                                        target_core_mass= float(parser.get(configurationSection, "coreMass")) | units.MSun,
                                            base_grid_options=dict(type="fcc"), pickle_file=savedMesaStarPath)
        else:
            mesaStar = self.EvolveStarWithStellarCode(MESA, savedMesaStarPath + "/" + self.saving_name, stellar_type= self.stellar_type)
            self.sphStar = convert_stellar_model_to_SPH(mesaStar, self.sphParticles, do_relax = False, with_core_particle=True, target_core_mass= mesaStar.core_mass,
                                            base_grid_options=dict(type="fcc"))
        self.gas_particles = self.sphStar.gas_particles
        self.core_particle = self.sphStar.core_particle

    def CalculateBindingEnergy(self, star):
        mass_profile = star.get_mass_profile()
        cumulative_mass_profile = star.get_cumulative_mass_profile()
        thermal_energy_profile = star.get_thermal_energy_profile()
        radius_profile = star.get_radius_profile()
        E = 0.0 |units.erg
        for i in range(len(mass_profile)):
            if (cumulative_mass_profile[i] | units.MSun) <= star.core_mass:
                continue
            E += (thermal_energy_profile[i] -constants.G*(cumulative_mass_profile[i] | units.MSun)/radius_profile[i])*(mass_profile[i] | units.MSun)
            print cumulative_mass_profile[i], mass_profile[i], radius_profile[i]

        return E

    def CheckLimitType(self, starType):
            return starType >= 16 or starType == 7 or starType < self.stellar_type.value_in(units.stellar_type)

    def EvolveStarWithStellarCode(self, code = MESA, savingPath = "", stellar_type= 3 | units.stellar_type ):
        '''
        evolve with (default) MESA or other
        :return: the star after has been created with MESA
        '''
        evolutionType = code()
        #evolutionType2=code()
        print "evolving with MESA"
        radiuses = []
        times = []
        p = Particle(mass= self.pointStar.mass)
        p.metalicity = self.pointStar.metalicity
        p.radius = self.pointStar.radius
        p.position = self.pointStar.position
        p.velocity = self.pointStar.velocity
        #mainStar = evolutionType.pre_ms_stars.add_particle(self.pointStar)
        #mainStar2 = evolutionType2.pre_ms_stars.add_particle(p)
        mainStar = evolutionType.particles.add_particle(self.pointStar)
        print "particle with mass=",mainStar.mass,"  added, current radius = ", mainStar.radius.as_quantity_in(units.RSun)," current type=",mainStar.stellar_type,"  target radius = ", self.radius, " target type = ",stellar_type
        oldStellarType = mainStar.stellar_type.value_in(units.stellar_type)
        maxRadii = mainStar.radius
        oldTime = time.time()
        #print evolutionType.pre_ms_stars
        #print evolutionType2.pre_ms_stars
        try:
            os.makedirs(savingPath)
        except(OSError):
            pass
        mainStar.reset_number_of_backups_in_a_row()
        while mainStar.radius < self.radius or \
                (self.radius.value_in(units.RSun) == 0 and
                     self.CheckLimitType(mainStar.stellar_type.value_in(units.stellar_type))):
            try:
                mainStar.evolve_one_step()
            except:
                try:
                    radiuses.append(mainStar.radius)
                    times.append(mainStar.age)
                    if mainStar.radius > maxRadii:
                        maxRadii = mainStar.radius
                    if mainStar.stellar_type.value_in(units.stellar_type)!= oldStellarType:
                        oldStellarType = mainStar.stellar_type.value_in(units.stellar_type)
                        #print mainStar.stellar_type, mainStar.radius, maxRadii
                        E = self.CalculateBindingEnergy(mainStar)
                        print "E=", E
                        print "lamda=", (-constants.G*mainStar.mass*(mainStar.mass-mainStar.core_mass)/mainStar.radius) / E
                    mainStar.reset_number_of_backups_in_a_row()
                    mainStar.evolve_one_step()
                except Exception as e:
                    print "could not evolve further", e
                    break
            currTime = time.time()
            if (currTime-oldTime) > 120:
                print "*",mainStar.radius
            oldTime = currTime
            #old=mainStar2.age
            #mainStar2.evolve_for(70|units.yr)
            #mainStar.evolve_for(mainStar2.age-old)
            #evolutionType.evolve_model(mainStar.age + (100|units.yr))
            #print "***************"
            #rint mainStar.age, mainStar2.age
            #print mainStar.radius, mainStar2.radius
            #print mainStar.temperature, mainStar2.temperature
            #print mainStar.get_radius_profile()[0], mainStar2.get_radius_profile()[0]
            #cm= mainStar.get_cumulative_mass_profile()
            #print cm[0], cm[-1]
            #u =  mainStar.get_thermal_energy_profile()
            #print u[0], u[-1]
            # =  mainStar.get_luminosity_profile()
            #print l[-1], mainStar2.get_luminosity_profile()[-1]
            #[0] = l[0] * 1000
            #ainStar.set_luminosity_profile(l)
            #print mainStar.get_luminosity_profile()[-1] - l[-1], l[-1]
            #mainStar.luminosity += 100.0 | units.K
            radiuses.append(mainStar.radius)
            times.append(mainStar.age)
            if mainStar.stellar_type.value_in(units.stellar_type) != oldStellarType:
                print mainStar.stellar_type, mainStar.radius , maxRadii
                E = self.CalculateBindingEnergy(mainStar)
                print "E=", E
                print "lamda=", (-constants.G*mainStar.mass*(mainStar.mass-mainStar.core_mass)/mainStar.radius) / E
                #save a pickle of all the mesa properties
                if not os.path.isfile(savingPath + "/" + code.__name__ + "_" + str(mainStar.mass.value_in(units.MSun)) + "_" + str(mainStar.stellar_type.value_in(units.stellar_type))):
                    pickle_stellar_model(mainStar, savingPath + "/" + code.__name__ + "_" + str(mainStar.mass.value_in(units.MSun)) + "_" + str(mainStar.stellar_type.value_in(units.stellar_type)))
                oldStellarType = mainStar.stellar_type.value_in(units.stellar_type)
            else:
                if maxRadii < mainStar.radius:
                    maxRadii = mainStar.radius
        radiuses.append(mainStar.radius)

        print evolutionType
        print mainStar
        if not os.path.isfile(savingPath + "/" + code.__name__ + "_" + str(mainStar.mass.value_in(units.MSun)) + "_" + str(mainStar.stellar_type.value_in(units.stellar_type))):
            pickle_stellar_model(mainStar, savingPath + "/" + code.__name__ + "_" + str(mainStar.mass.value_in(units.MSun)) + "_" + str(mainStar.stellar_type.value_in(units.stellar_type)))

        textFile = open(savingPath + '/radiuses_' + str(mainStar.mass.value_in(units.MSun)) + '.txt', 'w')
        textFile.write(', '.join([str(y) for y in radiuses]))
        textFile.close()

        textFile = open(savingPath + '/times_' + str(mainStar.mass.value_in(units.MSun)) + '.txt', 'w')
        textFile.write(', '.join([str(y) for y in times]))
        textFile.close()
        print "star saved to: ", savingPath + "/" + code.__name__ , "mass: ",mainStar.mass, "stellar type:", mainStar.stellar_type
        return mainStar

def Start(savedVersionPath = "/BIGDATA/code/amuse-10.0/Glanz/savings/MesaModels", takeSavedState = False, step = -1, configurationFile = "/BIGDATA/code/amuse-10.0/Glanz/savings/MesaModels/1RGBConfiguration.ini", savedMesaPath=""):
    #print types["RGB"], types["AGB"]
    giant = StarModels.CreatePointStar(configurationFile,configurationSection="MainStar")
    if takeSavedState:
        sphStar = SphStar(giant, configurationFile, configurationSection="MainStar",
                                savedMesaStarPath= savedMesaPath, takeSavedMesa=True)
    else:
        sphStar = SphStar(giant, configurationFile, configurationSection="MainStar",
                                savedMesaStarPath = savedVersionPath, takeSavedMesa=False)
    #saved state
    unitConverter = nbody_system.nbody_to_si(sphStar.gas_particles.total_mass() + sphStar.core_particle.mass, sphStar.core_particle.radius*1000*2)
    system=Gadget2(unitConverter)
    system.dm_particles.add_particle(sphStar.core_particle)
    system.gas_particles.add_particles(sphStar.gas_particles)
    StarModels.SaveDm(savedVersionPath+"/dm_" + str(sphStar.pointStar.mass.value_in(units.MSun)) + "_" + str(sphStar.stellar_type.value_in(units.stellar_type)) +".amuse", system.dm_particles)
    StarModels.SaveGas(savedVersionPath+"/gas_" + str(sphStar.pointStar.mass.value_in(units.MSun)) + "_" + str(sphStar.stellar_type.value_in(units.stellar_type)) +".amuse", system.gas_particles)
    print "state saved - {0}".format(savedVersionPath)

    print "****************** Simulation Completed ******************"

def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--saving_path', type=str, required=True,  help='path for saving models')
    parser.add_argument('--configuration_file', type=str,  help='path to where to config file is located', default="")

    return parser

if __name__ == "__main__":
    parser = InitParser()
    args = parser.parse_args()
    if args.configuration_file == "":
        configuration_files = [config_file for config_file in os.listdir(args.saving_path) if "Configuration.ini" in config_file]
        if len(configuration_files) > 1:
            print "cant figure out which configuration file to use"
            raise Exception("too many confing files")
        configuration_file = os.path.join(args.saving_path,configuration_files[0])
    else:
        configuration_file = args.configuration_file
    print "taking config file: ", configuration_file
    Start(savedVersionPath = args.saving_path, takeSavedState = False, step = -1, configurationFile = configuration_file)
    #Start(savedVersionPath="/vol/sci/astro/bigdata/glanz/amuse10/savings/MesaModels/HotJupiter", takeSavedState=False,step=-1,configurationFile="/vol/sci/astro/bigdata/glanz/amuse10/savings/MesaModels/HotJupiter/HJConfiguration.ini")

    #Start(savedVersionPath = "/home/hilaglanz/Dropbox/TCE/plots/MesaModels",savedMesaPath = "/home/hilaglanz/Dropbox/TCE/plots/MesaModels/1_0AGB/MESA_0.663429839586_5", takeSavedState = "True", step = -1, configurationFile = "/home/hilaglanz/Dropbox/TCE/plots/MesaModels/1AGBConfiguration.ini")
    '''
    Start(savedVersionPath = "/BIGDATA/code/amuse-10.0/Glanz/savings/MesaModels", takeSavedState = "False", step = -1, configurationFile = "/BIGDATA/code/amuse-10.0/Glanz/savings/MesaModels/1RGBConfiguration.ini")
    Start(savedVersionPath = "/BIGDATA/code/amuse-10.0/Glanz/savings/MesaModels", takeSavedState = "False", step = -1, configurationFile = "/BIGDATA/code/amuse-10.0/Glanz/savings/MesaModels/1.5AGBConfiguration.ini")
    Start(savedVersionPath = "/BIGDATA/code/amuse-10.0/Glanz/savings/MesaModels", takeSavedState = "False", step = -1, configurationFile = "/BIGDATA/code/amuse-10.0/Glanz/savings/MesaModels/2AGBConfiguration.ini")
    Start(savedVersionPath = "/BIGDATA/code/amuse-10.0/Glanz/savings/MesaModels", takeSavedState = "False", step = -1, configurationFile = "/BIGDATA/code/amuse-10.0/Glanz/savings/MesaModels/2.5AGBConfiguration.ini")
    Start(savedVersionPath = "/BIGDATA/code/amuse-10.0/Glanz/savings/MesaModels", takeSavedState = "False", step = -1, configurationFile = "/BIGDATA/code/amuse-10.0/Glanz/savings/MesaModels/3AGBConfiguration.ini")
'''
