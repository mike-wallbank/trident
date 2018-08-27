## Script to submit trident jobs to the grid

## First, copy everything which is required into the source directory

## Remove tarball if it exists
if [ -f source.tar ]
then
    rm -f source.tar
fi

## Make directory if it doesn't exist
if [ ! -d source ]
then
    mkdir source
else
    echo "Using existing source directory"
fi

## Copy everything into the source directory
cp Particle* source
cp compile_trident_mumu.C source
cp trident_mumu* source
cp -r dunend source &>/dev/null

## Make sure we have a valid proxy
kx509
voms-proxy-init -noregen -rfc -voms dune:/dune/Role=Analysis

## Submit the job
jobsub_submit --group dune --role=Analysis -N 100 --OS=SL6 --expected-lifetime=8h --memory=6000MB --tar_file_name=tardir:///dune/app/users/wallbank/trident/source/ file:///dune/app/users/wallbank/trident/grid_trident.sh