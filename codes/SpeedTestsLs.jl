using Statistics
using LinearAlgebra, Arpack, SparseArrays
using KPM
using JLD2
using CUDA

include("model.jl")
@show CUDA.has_cuda()
@show KPM.whichcore()
# check number of threads
println("Number of threads: $(Threads.nthreads())")

function BuildHdis(L;W=0.0,twists=[0.0,0.0,0.0], D=6.0)
	Ls = [L,L,L]
	M = 2.0 # 2.0
	t = 1.0
	γ = 0.8  # 0.8
	boundary = exp.(1im.*twists)
	@time Ham = TwoBandModel(Ls, t, M, γ;boundary=boundary)
	Jx, Jy, Jz = velocityOperator(Ls,Ham)
	H_norm = Ham./D
	return H_norm,Jx,Jy,Jz
end

function obtainGamma(H_norm,Jx,Jy,Jz;NC=2^6,D=6.0,NR=1)
	NH = H_norm.n
	mu_3d_xyz = zeros(ComplexF64, NC, NC, NC)
	psi_in_l = exp.(2pi * 1im * rand(H_norm.n, NR));
	KPM.normalize_by_col(psi_in_l, NR)
	psi_in_r = psi_in_l
	KPM.kpm_3d!(H_norm, Jx, Jy, Jz, NC, NR, NH, mu_3d_xyz, psi_in_l, psi_in_r; verbose=1, arr_size=16)
	return mu_3d_xyz
end

L = 20
@time H_norm, Jx,Jy,Jz = BuildHdis(L);
runtime = @elapsed obtainGamma(H_norm,Jx,Jy,Jz);

runsam = parse(Int64,ARGS[1])

pt = "/scratch/aw666/CPGEdata/"  # amarel
fn = pt*"Gamma$(L)NC64CPUrunsam$(runsam)"*".jld2"
jldsave(fn; runtime=runtime,L=L)



