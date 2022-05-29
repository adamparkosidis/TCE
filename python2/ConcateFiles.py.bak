import os, sys, glob
import argparse

def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--first', type=str, help='first day')
    parser.add_argument('--last', type=str, help='last day')
    parser.add_argument('--name', type=str, help='name of files without times')
    parser.add_argument('--source_dir', type=str, help='path to files directory')
    parser.add_argument('--timestep', type=float,default=0.2, help='timesteps in days')
    return parser

def GetLastDayOfFile(fileName):
    return fileName.split('_to_')[-1].split('days')[0]

def GetListOfFilesToConcat(timesDict, curr_first, last, timestep):
    if curr_first == last:
        return []

    if curr_first not in timesDict:
        print curr_first, " is not the first in any file"
        return None

    for f in timesDict[curr_first]:
        print f
        l = GetListOfFilesToConcat(timesDict, str(float(GetLastDayOfFile(f)) + timestep), last, timestep)
        if l is not None:
            return [f] + l

        l = GetListOfFilesToConcat(timesDict, str(float(GetLastDayOfFile(f))), last, timestep)
        if l is not None:
            return [f] + l
        '''
        for step in range(0,20):
            l = GetListOfFilesToConcat(timesDict, str(float(GetLastDayOfFile(f)) + 0.1 * step), last, timestep)
            if l is not None:
                return [f] + l
        '''

    return None

def ConcateFiles(files):
    data = ""
    first = True
    for f in files:
        file = open(f,mode="r")
        if not first:
            data += ", "
        data += file.readline()
        first = False

    return data

if __name__ == "__main__":
    parser = InitParser()
    args = parser.parse_args()
    files = glob.glob(args.source_dir + "\\" + args.name + "*")
    begintimeDict = dict()
    for f in files:
        beginTime = f.split('time_')[-1].split('_to')[0]
        if beginTime not in begintimeDict:
            begintimeDict[beginTime] = []
        begintimeDict[beginTime].append(f)

    filesToConcat = GetListOfFilesToConcat(begintimeDict, args.first, args.last, args.timestep)
    if filesToConcat is None:
        print "could not concat those files..."
    else:
        newFile = open(args.source_dir + "\\" + args.name + "time_" + args.first + "_to_" + args.last, mode="w")
        newFile.write(ConcateFiles(filesToConcat))
        newFile.close()
        print "new file has been saved to ", args.source_dir + "\\" + args.name + "time_" + args.first + "_to_" + args.last