"""
clusterfit.py
May 3rd, 2016

Uses mpi4py to leverage the power of any beowulf cluster.  The problem of
    fitting functions to the allele data is the epitome of an
    "embarassingly parallel" problem, so it is logical to use an MPI
    to solve it.
"""
import argparse
import csv
import os
import sys
from sys import argv
from mpi4py import MPI
import time
from hzarloop import Hzarloop as hzar


def clusterFit(inputfile, outputcsv, locinames, pickAlleles=False, iteration=1, modelIII=True, chainlength="1e5"):
    """
    use mpi4py and a message passing interface such as mpich or openmpi
        to make the fits in parallel across a beowulf cluster.
    """
    comm = MPI.COMM_WORLD           # The mpi object used to communicate between cores
    rank = comm.Get_rank()          # The core number to differentiate the processes
    size = comm.Get_size()          # number of cores specified
    if rank == 0:
        start = time.time()
        # print("Rank 0: {}.".format(sys.version))               # Used for debugging
        hz = hzar(inputfile, outputcsv, locinames, pickAlleles, iteration, modelIII, chainlength)  # create Hzarloop object
        with open(hz.inputfile, "r") as f:
            csvstring = f.read()
            csvstring = csvstring.replace(".A", ".a").replace(".B", ".b").replace(".N", ".n")
            csvstring = csvstring.replace("_", "")
            # print(csvstring)
        with open("rfile.R", "w") as f:
            f.write(hz.rawrscript.format(**{"csvtext": csvstring, "chainlength": hz.chainlength}))
        jobs = []
        numloci = len(hz.locinames)
        if numloci >= size-1:
            for i in range(size-1):
                job = []
                for j in range(i, numloci, size-1):
                    job.append(hz.locinames[j])
                jobs.append([job, csvstring])
        else:
            raise Exception("Use at most {} cores.".format(numloci - 1))

        # print([job[0] for job in jobs])                         # check to make sure we are divying up work correctly

        for i, job in enumerate(jobs):
            dest = i + 1
            comm.send(job, dest=dest, tag=1)

        with open(outputcsv, "w") as f:
            writer = csv.DictWriter(f, fieldnames=hz.headers, lineterminator="\n")
            writer.writeheader()
        with open("done.txt", "w") as f:
            pass

        n = 0
        donecount = 0
        totalanalyzed = 0
        doneset = set()
        while donecount < size-1:
            if n not in doneset:
                print("checking source {}.".format(n+1))
                if totalanalyzed != 0:
                    predict = (((time.time() - start) / totalanalyzed) * numloci)/3600
                    print("~~ Predicted time: {} hours or {} minutes.".format(round(predict, 2), round(predict * 60, 2)))
                    print("~~ Iteration: {}.".format(iteration + 1))
                print("~~ Checking core:  {}.".format(n + 1))
                print("~~ Cores finished: {}.".format(donecount))
                print("~~ Alleles fitted: {} of {}.".format(totalanalyzed, numloci))
                print("\n\n\n\n\nWAITING TO RECEIVE {}".format(n+1))
                towrite = comm.recv(source=n+1, tag=2)
                # print(towrite)
                if towrite != "done":
                    totalanalyzed += 1
                    with open(outputcsv, "a") as f:                           # Intentionally opening and closing many times so that data will not be lost in case of error
                        writer = csv.DictWriter(f, fieldnames=hz.headers, lineterminator="\n")
                        print("Writing new row.")
                        writer.writerow(towrite)
                else:
                    donecount += 1
                    doneset.add(n)
            n = (n+1) % (size-1)
        end = time.time()
        print("Finished analyzing {} alleles in {} seconds with {} cores.".format(numloci, end-start, size))

    else:
        # print("Rank {}: {}.".format(rank, sys.version))         # Used for debugging
        # print("Rank {}")
        jobs = comm.recv(source=0, tag=1)
        # print(job[1])
        for job in jobs[0]:
            towrite = hzar.makeFit([job], hzar.rawrscript, jobs[1], rank)
            # print("towrite: {}.".format(towrite))
            comm.send(towrite[0], dest=0, tag=2)
        comm.send("done", dest=0, tag=2)


def main():
    parser = argparse.ArgumentParser(
        description="Use hzar r package to fit cline models to multiple loci",
        )
    parser.add_argument('--input', metavar='FN', help='hzar input')
    parser.add_argument('--locinames', metavar='LOC', help='file containing names of loci', default="None")
    parser.add_argument('--output', metavar='DIR', help='directory for output files', default='None')
    parser.add_argument('--repeat', metavar='REP', help='Number of times to repeat the fitting process.', default="1")
    parser.add_argument('--pickalleles', metavar='PA', help='Whether or not to let the program decide which alleles to fit. Takes in a filename to write to.', default='None')
    parser.add_argument('--modelIII', metavar='MIII', help='Whether or not to fit modelIII.  True or False.', default="True")
    parser.add_argument('--chainlength', metavar="CL", help='Chainlength.  Default 1e5', default="1e5")
    if len(argv) == 1:
        parser.print_usage()
        exit(0)
    args = parser.parse_args()
    try:
        os.mkdir(args.output)
    except:
        print("directory exists")
    for i in range(int(args.repeat)):
        locinames = args.locinames if args.locinames != "None" else None
        if args.output != "None":
            fname = args.output.split("/")[-1]
            outname = "{}/{}_{}.csv".format(args.output, fname, i)
        else:
            outname = None
        pick = args.pickalleles if args.pickalleles != "None" else False
        usemodeliii = True if args.modelIII == "True" else False
        chainlength = args.chainlength
        clusterFit(args.input, outname, args.locinames, pick, i, usemodeliii, chainlength)


def test1000():
    clusterFit("hzar_10loci.csv", "small_output/qsub1000.csv", "alleleList.txt")


def testSingle():
    clusterFit("hzar_10loci.csv", "small_output/singleacisstest.csv", "locinames_short.txt")


def test10():
    clusterFit("hzar_10loci.csv", "small_output/aciss10.csv", "hzar_10loci_names.txt")


def test100():
    clusterFit("hzar_10loci.csv", "small_output/qsubtest.csv", "hzar_100loci_names.txt")


if __name__ == '__main__':
    main()
    # test10()
    # test100()
