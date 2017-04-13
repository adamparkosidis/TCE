import os
import ConfigParser
from amuse.units import units
from amuse.lab import *
from amuse.community.mesa.interface import MESA
import StarModels

types = dict()
types["MS"] = 1 | units.stellar_type
types["RGB"] = 3 | units.stellar_type
types["AGB"] = 4 | units.stellar_type


class SphStar:
    def __init__(self, pointStar, configurationFile="", configurationSection="", savedMesaStarPath = "", takeSavedMesa = False,savedGas="", savedDm=""):
        print 'parsing configurations'
        parser = ConfigParser.ConfigParser()
        parser.read(configurationFile)
        self.pointStar = pointStar
        self.sphParticles = parser.get(configurationSection, "sphParticles")
        self.stellar_type = types[(parser.get(configurationSection, "stellar_type"))]

        mesaStar = self.EvolveStarWithStellarCode(MESA, savedMesaStarPath, stellar_type= self.stellar_type)

        self.sphStar = convert_stellar_model_to_SPH(mesaStar, self.sphParticles, do_relax = False, with_core_particle=False,
                                            base_grid_options=dict(type="fcc"))
        self.gas_particles = self.sphStar.gas_particles
        self.core_particle = self.sphStar.core_particle

    def  EvolveStarWithStellarCode(self, code = MESA, savingPath = "", stellar_type= 3 | units.stellar_type ):
        '''
        evolve with (default) MESA or other
        :return: the star after has been created with MESA
        '''
        evolutionType = code()
        print "evolving with MESA"
        mainStar = evolutionType.particles.add_particle(self.pointStar)
        print "particle added, current radius = ", mainStar.radius.as_quantity_in(units.AU), "target type = ",stellar_type
        while mainStar.stellar_type < stellar_type:
            mainStar.evolve_one_step()
        try:
            os.makedirs(savingPath)
        except(OSError):
            pass
        print evolutionType
        print mainStar
        pickle_stellar_model(mainStar, savingPath + "/" + code.__name__ + "_" + str(mainStar.mass.value_in(units.RSun)) + "_" + str(mainStar.stellar_type.value_in(units.stellar_type)))
        print "star saved to: ", savingPath + "/" + code.__name__ , "mass: ",mainStar.mass, "stellar type:", mainStar.stellar_type
        return mainStar

def Start(savedVersionPath = "/BIGDATA/code/amuse-10.0/Glanz/savings/MesaModels", takeSavedState = "False", step = -1, configurationFile = "/BIGDATA/code/amuse-10.0/Glanz/savings/MesaModels/1RGBConfiguration.ini"):

    giant = StarModels.CreatePointStar(configurationFile,configurationSection="MainStar")
    sphStar = SphStar(giant, configurationFile, configurationSection="MainStar",
                                savedMesaStarPath = "", takeSavedMesa=False)
    #saved state
    StarModels.SaveDm(savedVersionPath+"/dm_" + sphStar.pointStar.mass.value_in(units.RSun) + "_" + str(sphStar.stellar_type.value_in(units.stellar_type)) +".amuse", [sphStar.core_particle])
    StarModels.SaveGas(savedVersionPath+"/envelope_" + sphStar.pointStar.mass.value_in(units.RSun) + "_" + str(sphStar.stellar_type.value_in(units.stellar_type)) +".amuse", sphStar.gas_particles)
    print "state saved - {0}".format(savedVersionPath)

    print "****************** Simulation Completed ******************"

if __name__ == "__main__":
    Start()

