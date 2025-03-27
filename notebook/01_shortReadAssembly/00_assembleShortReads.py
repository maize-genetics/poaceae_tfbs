#!/usr/bin/env python3
# Written by Aimee Schulz, adapted by Charlie Hale
# Usage: script.py <read_paths.tsv> <threads> OPTIONAL: <--interleaved> (default is to use separate files for forward and reverse reads) <--irods> (default is to not download from irods)  
# "read_paths.tsv" should be in the following format:
# Sample	FilePath1	FilePath2
# MySample1	path/to/file1	path/to/file2 # This is for a sample with a single set of paired-end reads
# MySample2	path/to/file1a,path/to/file1b	path/to/file2a,path/to/file2b #This is for a sample that needs multiple merged raw read files

import multiprocessing
import subprocess
import logging
from pandas import *
import os
import traceback
import argparse
import math
import glob

# Set root directory as current working directory
rootDir = os.getcwd()

# Create the parser to read command line arguments
parser = argparse.ArgumentParser(description='Process command line arguments')

# Add an argument for the TSV file path
parser.add_argument('tsv_path', help='Path to the TSV file')

# Add argument for number of threads to use
parser.add_argument('threads', help='Number of threads to use')

# Add an optional argument to determine number of raw read files to use
# Defaults to false unless specified
parser.add_argument('--interleaved', action='store_true', help='Should be set to true if paired-end reads are merged into same file')

# Add an optional argument to check whether to download from irods
# defaults to False unless specified
parser.add_argument('--irods', action='store_true', help='Enable irods (default: False)')

# Parse the arguments
args = parser.parse_args()

# Read the TSV file using the provided path
#dataLocation = pandas.read_csv(args.csv_path)
dataLocation = pandas.read_csv(args.tsv_path, delimiter='\t')

#Check the value of irods argument
if args.irods:
    from irods.session import iRODSSession

    # Perform actions if irods is True
    print("iRODS download is enabled")
    #irod environment file
    #make sure there is a field "irods_user_name" in this json file
    irods_env_file = "/home/coh22/.irods/irods_environment.json"

else:
    # Perform actions if irods is False
    print("Note: iRODS download is not enabled. Run with --irods if iRODS download is necessary")

##Modify information in this section ###################################
#number of simultaneous samples
simultaneousProcess = 3

#list of paths of input paired-end data files,
inputFileList1 = dataLocation['Path1'].tolist()
inputFileList2 = dataLocation['Path2'].tolist()

#list of sample Names
inputSampleList = dataLocation['Sample'].tolist()

#this directory is used for storing local file  and run pipeline
tmpDir = f"{rootDir}/tmp"

#result directory keep all final results
resultDir = f"{rootDir}/assembly_results"

# Prep dependencies
cmd = f"tar xzvf data/assembly_dependencies.tar.gz"
process = subprocess.call(cmd, shell=True) 
cmd = f"git clone https://bitbucket.org/tasseladmin/tassel-5-standalone.git"
process = subprocess.call(cmd, shell=True)                                  

#######################################################################

def main():
    global session

    if (not os.path.exists(tmpDir)):
        os.mkdir(tmpDir)
    if (not os.path.exists(resultDir)):
        os.mkdir(resultDir)

    logging.basicConfig(filename=f"{tmpDir}/run.log",level=logging.DEBUG)

    if args.irods:
        #start a iRods session
        session=iRODSSession(irods_env_file=irods_env_file)

    #number of input files
    inputFileCount = len(inputSampleList)
    sampleToFileList = []

    for i in range(0, inputFileCount):
        filePath1 = inputFileList1[i]
        filePath2 = inputFileList2[i]
        sampleName = inputSampleList[i]
        sampleToFileList.append((filePath1, filePath2, sampleName))

    pool = multiprocessing.Pool(processes= simultaneousProcess)
    pool.starmap(myFunction, sampleToFileList)
    pool.close()

    # delete temporary directories
    cmd = f"rm -fr tassel-5-standalone/"
    subprocess.call(cmd, shell=True)
    cmd = f"rm -fr assembly_dependencies"
    subprocess.call(cmd, shell=True)

def myFunction (filePath1, filePath2, sampleName):
    logging.info(f"Process {sampleName}: {filePath1}")
    try:
        #create a workDir for each sample, the directory name is the sample name
        workDir = f"{tmpDir}/{sampleName}"
        if (not os.path.exists(workDir)):
            os.mkdir(workDir)
        # Split multi-file paths into individual paths
        forwardList = filePath1.split(",")  # Assuming filePath1 and filePath2 are strings containing file paths separated by commas
        revList = filePath2.split(",")
        nFiles = len(forwardList)  # Assuming you want the number of files in forwardList
        if len(forwardList) == len(revList):
            logging.info(f"{sampleName}: Using {nFiles} paired-end file(s) for assembly")
        else:
            logging.info(f"{sampleName}: Error: number of forward and reverse read files do not match")
        
        # Download to workdir
        if args.irods:
            print(f"Getting iRODS data for {workDir}...")
        for file_index in range(nFiles):  # Fixed loop syntax
            forwardFileName = f"input1_{file_index+1}.fastq.gz"  # Corrected string formatting
            revFileName = f"input2_{file_index+1}.fastq.gz"
            forwardFilePath = forwardList[file_index]
            revFilePath = revList[file_index]

            session.data_objects.get(forwardFilePath, f"{workDir}/{forwardFileName}")
            session.data_objects.get(revFilePath, f"{workDir}/{revFileName}")

        #step 1 - assembly
        print("Running assembly step 1...")
        # set fractions of threads and memory to use
        threadFraction = math.floor(int(args.threads) / simultaneousProcess)
        #memFraction = round(1 / simultaneousProcess, 1)
        memFraction = 0.5
        # Find all forward read files and concatenate their real paths together for megahit
        forwardFiles = glob.glob(f"{workDir}/input1*.fastq.gz")
        forwardPaths = ",".join([os.path.realpath(file) for file in forwardFiles])
        print(forwardPaths)
        # Find all reverse read files and concatenate their real paths
        reverseFiles = glob.glob(f"{workDir}/input2*.fastq.gz")
        reversePaths = ",".join([os.path.realpath(file) for file in reverseFiles])
        print(reversePaths)
        if args.interleaved:
            #cmd = f"megahit --12 {workDir}/input1.fastq.gz -m {memFraction} -t {threadFraction} --k-min 31 -o {workDir}/megahit"
            cmd = f"megahit --12 {forwardPaths} -m {memFraction} -t {threadFraction} --k-min 31 -o {workDir}/megahit"

        else:
            #cmd = f"megahit -1 {workDir}/input1.fastq.gz -2 {workDir}/input2.fastq.gz -m {memFraction} -t {threadFraction} --k-min 31 -o {workDir}/megahit"
            cmd = f"megahit -1 {forwardPaths} -2 {reversePaths} -m {memFraction} -t {threadFraction} --k-min 31 -o {workDir}/megahit"
        returned_value = subprocess.call(cmd, shell=True)
        if (returned_value==0):
            logging.info(f"{sampleName}: Step 1 finished correctly.")
        else:
            logging.info(f"{sampleName}: Step 1 failed.")

        #step 2 assembly
        print("Running assembly step 2...")
        cmd = f"cp {workDir}/megahit/final.contigs.fa {workDir}/megahit/{sampleName}.final.contigs.fa"
        returned_value = subprocess.call(cmd, shell=True)
        if (returned_value==0):
            logging.info(f"{sampleName}: Step 2 finished correctly.")
        else:
            logging.info(f"{sampleName}: Step 2 failed.")

        #step 3 assembly
        print("Running assembly step 3...")
        cmd = f"sh {rootDir}/src/01_shortReadAssembly/tabasco_wrapper.sh {workDir}/megahit/{sampleName}.final.contigs.fa assembly_dependencies/transcript_seqs 0.6"
        returned_value = subprocess.call(cmd, shell=True)
        if (returned_value==0):
            logging.info(f"{sampleName}: Step 3 finished correctly.")
        else:
            logging.info(f"{sampleName}: Step 3 failed.")


        #step 4 copy result file delete intermediate files
        print("Running step 4...")
        cmd = f"mkdir -p {resultDir}/{sampleName}"
        returned_value = subprocess.call(cmd, shell=True)
        cmd = f"cp -r {workDir}/megahit/* {resultDir}/{sampleName}"
        returned_value = subprocess.call(cmd, shell=True)
        if (returned_value==0):
            logging.info(f"{sampleName}: copy results finished correctly.")
        else:
            logging.info(f"{sampleName}:copy results failed.")

        #delete temporary data files
        #cmd = f"rm -fr {workDir}"

    except Exception as e:
        logging.info(f"Error: {sampleName}")
        traceback.print_exc()
        print()
        raise e


if __name__=="__main__":
        main()
