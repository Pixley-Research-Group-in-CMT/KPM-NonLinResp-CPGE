#!/bin/bash

# extract Julia tar.gz file
tar -xzf julia-1.7.1-linux-x86_64.tar.gz
tar -xzf julia-packages.tar.gz

# add Julia binary to PATH
export PATH=$_CONDOR_SCRATCH_DIR/julia-1.7.1/bin:$PATH
# add Julia packages to DEPOT variable
export JULIA_DEPOT_PATH=$_CONDOR_SCRATCH_DIR/julia-packages

# run Julia script
julia --project=julia-packages collectHdis.jl $1

# Transfer the output file to stash storage
#output_file="HdisL100gam08W2035dis$1.jld2"
#stash_destination="stash:///osgconnect/public/angkunwu/CPGEdata/$output_file"
#stashcp $output_file $stash_destination
