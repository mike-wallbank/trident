#!/bin/bash

# Script to run the trident analyzer over just signal events
# M Wallbank, August 2018

DATADIR="/pnfs/dune/persistent/users/jmalbos/Trident/data/sim/mumu/"
OUTDIR="/pnfs/dune/persistent/users/wallbank/trident/data/"
OUTFILE="TridentMuMuOut_SIG.root"

## Setup everything
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup root v6_12_04e -q e15:prof

## Run the root module
root -l -b -q "compile_trident_mumu.C(-1, 1, -2, \"${DATADIR}\", \"${OUTFILE}\")"

## Copy the file back
mv ${OUTFILE} ${OUTDIR}
