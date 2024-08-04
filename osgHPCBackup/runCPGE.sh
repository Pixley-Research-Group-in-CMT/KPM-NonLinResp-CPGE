#!/bin/bash

# extract Julia tar.gz file
tar -xzf julia-1.8.5-linux-x86_64.tar.gz
tar -xzf julia-packages.tar.gz

# add Julia binary to PATH
export PATH=$_CONDOR_SCRATCH_DIR/julia-1.8.5/bin:$PATH
# add Julia packages to DEPOT variable
export JULIA_DEPOT_PATH=$_CONDOR_SCRATCH_DIR/julia-packages

# run Julia script
julia --project=julia-packages test.jl $1

