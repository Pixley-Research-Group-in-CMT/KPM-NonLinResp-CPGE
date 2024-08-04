using Statistics
using LinearAlgebra, Arpack, SparseArrays
using KPM
using JLD2
using CUDA

@show CUDA.has_cuda()
@show KPM.whichcore()
# check number of threads
println("Number of threads: $(Threads.nthreads())")

disind = parse(Int64,ARGS[1])
pt = "/scratch/aw666/CPGEdata/"  # amarel
fn = pt*"HdisL100gam08W0PotDis"*"$disind"*".jld2"
H_norm = load(fn, "H_norm")
D = load(fn, "D")

NC = 2^15
NR = 10

mu = KPM.kpm_1d(H_norm, NC, NR)

fn = pt*"DOSmuL100gam08D6NC215NR10Sam"*"$disind"*".jld2"
jldsave(fn; NC=NC, D = D, mu=mu)
