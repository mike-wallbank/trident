#!/bin/bash

# Script to submit trident background data files to the grid
# M Wallbank, July 2018

## Get the background file number from the process
## Ensure it is zero-padded and two digits long
printf -v BGFILE "%02d\n" ${PROCESS}
export OUTFILE=TridentMuMuOut_BG${BGFILE}.root

echo Making background file $OUTFILE

#USERDIR="/dune/app/users/wallbank/trident/"
USERDIR="/pnfs/dune/persistent/users/wallbank/trident/source/"
DATADIR="/pnfs/dune/persistent/users/jmalbos/Trident/data/sim/mumu/"
OUTDIR="/pnfs/dune/persistent/users/wallbank/trident/data/"

## Setup everything
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

echo Set up everything

## Copy the relevant code and data file to the worker node
ifdh cp ${USERDIR}/* .
echo Successfully copied source code
ifdh cp ${DATADIR}/mumu.bkg.${BGFILE}.g4.root .
echo Successfully copied input file

## Run the root module
root "compile_trident_mumu.C(-1, 0, ${BGFILE}, \".\", ${OUTFILE})"
echo Finished running code

## Copy the file back
ifdh cp ${OUTFILE} ${OUTDIR}
echo Successfully copied file back