#!/bin/bash

# Script to submit trident background data files to the grid
# M Wallbank, July 2018

## Submit using
## jobsub_submit --group dune --role=Analysis -N 100 --OS=SL6 --expected-lifetime=8h --memory=6000MB file:///pnfs/dune/persistent/path/to/file...
## Number of jobs should correspond to number of input files
## Right now, we have 100

## Get the background file number from the process
## Ensure it is zero-padded and two digits long
export BGFILENUM=${PROCESS}
printf -v BGFILENUMFORMAT "%02d" ${BGFILENUM}
export OUTFILE=TridentMuMuOut_BG${BGFILENUMFORMAT}.root

echo Making background file $OUTFILE

DATADIR="/pnfs/dune/persistent/users/jmalbos/Trident/data/sim/mumu/"
OUTDIR="/pnfs/dune/persistent/users/wallbank/trident/data/"

## Setup everything
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup root v6_12_04e -q e15:prof
setup ifdhc

## Copy the relevant code and data file to the worker node
ifdh cp ${DATADIR}/mumu.bkg.${BGFILENUMFORMAT}.g4.root mumu.bkg.${BGFILENUMFORMAT}.g4.root

## Run the root module
root -l -b -q "compile_trident_mumu.C(-1, 0, ${BGFILENUM}, \".\", \"${OUTFILE}\")"

## Copy the file back
ifdh cp ${OUTFILE} ${OUTDIR}/${OUTFILE}
