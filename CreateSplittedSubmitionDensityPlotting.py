import os
import argparse
import subprocess

def main(args):
    chunkSize = (args.LAST_STEP - args.BEGIN_STEP) / args.NODES
    if chunkSize == 0:
        if args.LAST_STEP - args.BEGIN_STEP == 0:
            return
        else:
            chunkSize = 1

    for i in range(args.BEGIN_STEP,args.LAST_STEP,chunkSize):
        jobName = "density" + args.RADII + str(i)
        innerDir = args.INNER_DIR
        if innerDir is None:
            innerDir = ""
        jobOutputs = os.path.join(args.SIMULATIONS_DIR, args.RADII, str(args.GIANT_MASS)+"MSun",str(args.PHASE)+"Phase",
                                  str(args.INNER_SEPARATION)+"RSun",str(args.INCLINATION)+"inclin", innerDir)
        jobOutput = os.path.join(jobOutputs,jobName + ".out.txt")
        jobError = os.path.join(jobOutputs,jobName + ".err.txt")
        submitionBashCommand = "sbatch -n 18 -N 1 --mem-per-cpu=10G -J " + jobName +\
                               " --output=" + jobOutput + " --error=" + jobError + \
                               "--mail-type=ALL --mail-user=hilaglanz@gmail.com" \
                               "mpiexec -n 1 /usr/local/astro/amuse/11.2/./amuse-new.sh " + args.FILETORUN + " " + \
                               jobOutputs + " " + args.SNAPSHOTS_DIR + " " + str(i) + " " + str(i+chunkSize) + " " + \
                               args.V_MIN + " " + args.V_MAX + " " + str(args.PLOT)
        process = subprocess.Popen(submitionBashCommand, stdout= subprocess.PIPE)
        output, error = process.communicate()

def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--RADII', type=str)
    parser.add_argument('--GIANT_MASS', type=int)
    parser.add_argument('--PHASE', type=int)
    parser.add_argument('--INNER_SEPARATION', type=int)
    parser.add_argument('--INCLINATION', type=int)

    parser.add_argument('--SNAPSHOTS_DIR', type=str,  help='path to amuse files directory')
    parser.add_argument('--BEGIN_STEP', type=int,  help='first file number')
    parser.add_argument('--LAST_STEP', type=int,  help='last file number')
    parser.add_argument('--V_MIN', type=str, help= "vmin for plotting density profiles")
    parser.add_argument('--V_MAX', type=str, help ="vmax for plotting density profiles")
    parser.add_argument('--INNER_DIR', type=str, help = "inner simulation directory if exists")
    parser.add_argument('--PLOT', type=int,  help='true if should density profiles')

    parser.add_argument('--SIMULATIONS_DIR', type=str,  help='path to simulations directory')
    parser.add_argument('--NODES', type=int, help="number of nodes to spread to")

    parser.add_argument('--FILETORUN', type=str, help="path to file to run")
    return parser

if __name__ == "__main__":
    parser = InitParser()
    args = parser.parse_args()
    print(args)
    main(args)