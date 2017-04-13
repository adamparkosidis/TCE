import os
import StarModels
from amuse.units import units
from amuse.lab import *
from amuse.community.mesa.interface import MESA

def  EvolveStarWithStellarCode(self, code = MESA, stellar_type= 3 | units.stellar_type, savingPath = ""):
        '''
        evolve with (default) MESA or other
        :return: the star after has been created with MESA
        '''
        evolutionType = code()
        print "evolving with MESA"
        mainStar = evolutionType.particles.add_particle(self.pointStar)
        print "particle added, current radius = ", mainStar.radius.as_quantity_in(units.AU), "target radius = ", self.pointStar.radius
        while mainStar.stellar_type < stellar_type:
            mainStar.evolve_one_step()
        try:
            os.makedirs(savingPath)
        except(OSError):
            pass
        print evolutionType
        print mainStar
        pickle_stellar_model(mainStar, savingPath + "/" + code.__name__ + "_" + mainStar.mass + "_" + mainStar.stellar_type)
        print "star saved to: ", savingPath + "/" + code.__name__ , "mass: ",mainStar.mass, "stellar type:", mainStar.stellar_type
        return mainStar

def Start(savedVersionPath = "/BIGDATA/code/amuse-10.0/Glanz/savings/MesaModels", takeSavedState = "False", step = -1, configurationFile = "/BIGDATA/code/amuse-10.0/Glanz/savings/MesaModels/1RGBConfiguration.ini"):
    giant = StarModels.CreatePointStar(configurationFile,configurationSection="MainStar")

    sphStar = StarModels.SphStar(giant, configurationFile, configurationSection="MainStar",
                                savedMesaStarPath = "", takeSavedMesa=False)

    sphMetaData = StarModels.SphMetaData(sphStar)

    #saved state
    StarModels.SaveState(savedVersionPath, sphStar.gas_particles.total_mass() + sphStar.core_particle.mass,
                         sphStar.gas_particles, [sphStar.core_particle], 0 | units.AU, sphMetaData)

    print "****************** Simulation Completed ******************"

if __name__ == "__main__":
    Start()

