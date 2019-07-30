import os, io, glob
import gc
import h5py
import argparse
from amuse.io import  read_set_from_file
from amuse.units import units

def GetHeadresFromObject(object, finiteChar=""):
    print help(object)
    print object.__array_interface__
    print object.__dir__.__dict__
    csvHeaders = ""
    for key in object.__dict__:
        csvHeaders += key + finiteChar + ","
    if len(csvHeaders) > 0 :
        return csvHeaders[:-1]
    return csvHeaders

def GetHeadersOfParticle(particle, finiteChar=""):
    return "mass" + finiteChar + ",radius" + finiteChar + ",vx" + finiteChar + ",vy" + finiteChar + ",vz" + finiteChar + \
           ",x" + finiteChar + ",y" + finiteChar + ",z" + finiteChar +",ax" + finiteChar + ",ay" + finiteChar + ",az" + finiteChar
    '''return "ax" + finiteChar + ",ay" + finiteChar + ",az" + finiteChar + ",epsilon" + finiteChar + ",mass" + \
           finiteChar + ",radius" + finiteChar + ",vx" + finiteChar + ",vy" + finiteChar + ",vz" + finiteChar + \
           ",x" + finiteChar + ",y" + finiteChar + ",z" + finiteChar'''

def GetValuesOfObject(object, finiteChar=""):
    csvValues = ""
    for key, value in object.__dict__:
        csvValues += value + finiteChar + ","
    if len(csvValues) > 0 :
        return csvValues[:-1]
    return csvValues

def GetValuesOfParticle(particle):
    result =  str(particle.mass.value_in(units.g)) + "," \
    + str(particle.radius.value_in(units.cm)) + "," + str(particle.vx.value_in(units.cm / units.s)) + "," + \
           str(particle.vy.value_in(units.cm / units.s)) + "," + str(particle.vz.value_in(units.cm / units.s)) + "," + \
    str(particle.x.value_in(units.cm)) + "," + str(particle.y.value_in(units.cm)) + "," + str(particle.z.value_in(units.cm))
    try:
        result += ", " +str(particle.ax.value_in(units.cm / units.s**2)) + "," + str(particle.ay.value_in(units.cm / units.s**2)) + "," \
                  + str(particle.az.value_in(units.cm / units.s**2))
    except:
        result += ",,,"

    return result

def GetHeadersOfMultipleParticles(multipleObjects):
    headers = ""
    for i, obj in enumerate(multipleObjects):
        objHeaders = GetHeadersOfParticle(obj, str(i + 1))
        if len(headers) > 0 and len(objHeaders) > 0:
           headers += ","
        headers += objHeaders
    return headers

def GetValuesOfMultipleParticles(multipleObjects):
    values = ""
    for i, obj in enumerate(multipleObjects):
        objValues = GetValuesOfParticle(obj)
        if len(values) > 0 and len(objValues) > 0 :
            values += ","
        values += objValues
    return values


def GetTimeOfFile(fileNumber, defaultTimeStep = 0.2 | units.day):
    return (fileNumber * defaultTimeStep)

def GetBinaryStateFromFile(directoryPath, fileNumber):
    try:
        return read_set_from_file(os.path.join(directoryPath, "dm_" + fileNumber + ".amuse"), format='amuse')
    except Exception as ex:
        if fileNumber == 0 or fileNumber =='0':
            return read_set_from_file(os.path.join(directoryPath, "dm_00.amuse"), format='amuse')
        else:
            raise Exception(ex)

def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--first', type=int,  help='first file number')
    parser.add_argument('--last', type=int,  help='last file number')
    parser.add_argument('--time_step', type=float,  help='time between files in days')
    parser.add_argument('--source_dir', type=str,  help='path to amuse files directory')
    return parser

if __name__ == "__main__":
    parser = InitParser()
    args = parser.parse_args()
    firstBinary = GetBinaryStateFromFile(args.source_dir, str(args.first))
    csvData = GetHeadersOfMultipleParticles(firstBinary) + ",time" + "\r\n"

    if args.last == 0:
        snapshots = glob.glob(args.source_dir + "/dm*.amuse")
        numberOfSnapshots = len(snapshots) - 1
    else:
        numberOfSnapshots = args.last
    print numberOfSnapshots
    if args.time_step is None:
        args.time_step = 0.2
    args.time_step = args.time_step | units.day
    for n in xrange(args.first, numberOfSnapshots):
        csvData += GetValuesOfMultipleParticles(GetBinaryStateFromFile(args.source_dir, str(n))) + ", " + str(GetTimeOfFile(n, args.time_step).value_in(units.s)) + '\r\n'
        for f in [obj for obj in gc.get_objects() if isinstance(obj,h5py.File)]:
            try:
                f.close()
            except:
                pass

    file = open(args.source_dir + "/resultCSV.csv", mode="w")
    file.write(csvData)
    file.close()
    print "file saved to ", args.source_dir + "/resultCSV.csv"

