using KPM
using SuperLattice,SuperLattice.Util
using JLD2,HDF5

include("model.jl")

L = 18
Ls = [L,L,L]
M = 2 # 2.0 1.5
t = 1.0 
γ = 0.8  # 0.8 0
W = 0.5 

NC=32
NR=10

Emax = (M+2)+sqrt(1+γ^2)

boundary = [1,1,1]
H0 = TwoBandModel(Ls, t, M, γ;boundary=boundary)

# load disorders
f = h5open("CPGEdata/E2sL18M15W05gamma0Sam1050.h5","r")
E2s = read(f,"E2s")
disorders = read(f,"disorders")
close(f)

disind = argmin(E2s)[2]
Ham = H0 + spdiagm(disorders[:,k])#disind])
Jx, Jy, Jz = velocityOperator(Ls,Ham)
#D, _ = SuperLattice.Util.normalizeH(Ham;eps=0.01,setA=0.0)
D = 5.5245/(1-0.005) # choose the largest largest eigenvalues over disorder
H_norm = Ham./D

NH = H_norm.n
mu_3d_xyz = zeros(ComplexF64, NC, NC, NC)
psi_in_l = exp.(2pi * 1im * rand(H_norm.n, NR));
KPM.normalize_by_col(psi_in_l, NR)
psi_in_r = psi_in_l
KPM.kpm_3d!(H_norm, Jx, Jy, Jz, NC, NR, NH, mu_3d_xyz, psi_in_l, psi_in_r; verbose=1, arr_size=16);


mu_3d_yxz = zeros(ComplexF64, NC, NC, NC)
psi_in_l = exp.(2pi * 1im * rand(H_norm.n, NR));
KPM.normalize_by_col(psi_in_l, NR)
psi_in_r = psi_in_l
KPM.kpm_3d!(H_norm, Jy, Jx, Jz, NC, NR, NH, mu_3d_yxz, psi_in_l, psi_in_r; verbose=1, arr_size=16);


jldsave("GammamnpL18NC32W0NR10.jld2"; Ls=Ls, NC=NC, D = D, Gammaxyz=mu_3d_xyz, Gammayxz=mu_3d_yxz)

