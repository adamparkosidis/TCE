import os, io
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
    return "ax" + finiteChar + ",ay" + finiteChar + ",az" + finiteChar + ",epsilon" + finiteChar + ",mass" + \
           finiteChar + ",radius" + finiteChar + ",vx" + finiteChar + ",vy" + finiteChar + ",vz" + finiteChar + \
           ",x" + finiteChar + ",y" + finiteChar + ",z"

def GetValuesOfObject(object, finiteChar=""):
    csvValues = ""
    for key, value in object.__dict__:
        csvValues += value + finiteChar + ","
    if len(csvValues) > 0 :
        return csvValues[:-1]
    return csvValues

def GetValuesOfParticle(particle):
    return particle.ax.as_string_in(units.m / units.s**2) + "," + particle.ay.as_string_in(units.m / units.s**2) + "," + particle.az.as_string_in(units.m / units.s**2) \
           + "," + particle.epsilon.as_string_in(units.m) + "," + particle.mass.as_string_in(units.kg) + "," \
    + particle.radius.as_string_in(units.m) + "," + particle.vx.as_string_in(units.m / units.s) + "," + \
           particle.vy.as_string_in(units.m / units.s) + "," + particle.vz.as_string_in(units.m / units.s) + "," + \
    particle.x.as_string_in(units.m) + "," + particle.y.as_string_in(units.m) + "," + particle.z.as_string_in(units.m)

def GetHeadersOfBinaryObject(binary):
    headersOfFirst = GetHeadresFromObject(binary[0], "1")
    headersOfSecond = GetHeadresFromObject(binary[1], "2")
    headers = headersOfFirst
    if len(headersOfFirst) > 0 and len(headersOfSecond) > 0:
        headers += ","
    headers += headersOfSecond
    return headers

def GetValuesOfBinaryObject(binary):
    valuesOfFirst = GetValuesOfObject(binary[0], "1")
    valuesOfSecond = GetValuesOfObject(binary[1], "2")
    values = valuesOfFirst
    if len(valuesOfFirst) > 0 and len(valuesOfSecond) > 0:
        values += ","
    values += valuesOfSecond
    return values

def GetHeadersOfBinaryParticles(binary):
    headersOfFirst = GetHeadersOfParticle(binary[0], "1")
    headersOfSecond = GetHeadersOfParticle(binary[1], "2")
    headers = headersOfFirst
    if len(headersOfFirst) > 0 and len(headersOfSecond) > 0:
        headers += ","
    headers += headersOfSecond
    return headers

def GetValuesOfBinaryParticle(binary):
    valuesOfFirst = GetValuesOfParticle(binary[0])
    valuesOfSecond = GetValuesOfParticle(binary[1])
    values = valuesOfFirst
    if len(valuesOfFirst) > 0 and len(valuesOfSecond) > 0:
        values += ","
    values += valuesOfSecond
    return values

def GetTimeOfFile(fileNunber, defaultTimeStep = 0.2 | units.day):
    return (fileNunber * defaultTimeStep).as_string_in(units.day)

def GetBinaryStateFromFile(directoryPath, fileNumber):
    return read_set_from_file(os.path.join(directoryPath, "dm_" + fileNumber + ".amuse"), format='amuse')

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
    csvData = GetHeadersOfBinaryParticles(firstBinary) + ",time" + "\r\n"

    if args.last == 0:
        snapshots = os.listdir(args.source_dir)
        numberOfSnapshots = len(snapshots) / 2
    else:
        numberOfSnapshots = args.last
    if args.time_step is None:
        args.time_step = 0.2 | units.day
    for n in xrange(args.first, numberOfSnapshots + 1):
        csvData += GetValuesOfBinaryParticle(GetBinaryStateFromFile(args.source_dir, str(n))) + ", " + GetTimeOfFile(n, args.time_step | units.day)\
                   + "\r\n"
        for f in [obj for obj in gc.get_objects() if isinstance(obj,h5py.File)]:
            try:
                f.close()

    file = open(args.source_dir + "/resultCSV.csv", mode="w")
    file.write(csvData)
    file.close()
    print "file saved to ", args.source_dir + "/resultCSV.csv"

