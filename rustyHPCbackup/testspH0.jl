using SuperLattice
using SuperLattice.Util.Physics, SuperLattice.Util
using Statistics
using LinearAlgebra, Arpack, SparseArrays
using KrylovKit
using JLD2
include("model.jl")

#=
# model parameters
#=
L = 12
Ls = [L,L,L]
M = 2.0 # 2.0
t = 1.0
gamma = 0.0  # 0.8
=#
t=1.0
gamma=0.0
M=2
L=12 #30
sizes=[L,L,L]
W=0.0


OBC=false
tbc_theta = 0.5*ones(3) #zeros(3) #rand(3)
println("theta=$(tbc_theta)")

println("define model")
uc = UnitCell(;d=3)
addAtom(uc, :a, 2, [0,0,0]) # 2 band model


# potential term - disorder
Vr = function(x, y, v)
    v .= M * σz #+ chemmu * σ0 chemical potential enters in the final calculation
end


addPot(uc, :a, Vr)
# hopping term (only external)
addHop(uc, :a, :a, [1, 0, 0], -0.5im * t * σx - 0.5 * t * σz ; hc=true)
addHop(uc, :a, :a, [0, 1, 0], -0.5im * t * σy - 0.5 * t * σz; hc=true)
#addHop(uc, :a, :a, [0, 1, 0], -0.5im * t * σx - 0.5 * t * σz; hc=true)
addHop(uc, :a, :a, [0, 0, 1], 0.5im * gamma * σ0 - 0.5 * t * σz; hc=true)

# set up lattice
ltc = Lattice(d=3,pv=[[1.,0,0],[0.,1,0],[0.,0,1]],sizes=sizes,OBC=OBC);
uc_symb = addUC(ltc, uc);
# add magnetic field
# addMagneticField(ltc; B=[0,0,0.1])

populateUC(ltc, uc_symb);
Hsp_gen = get_operator_gen(ltc);

@time All_ops = SuperLattice.to_sparse_matrix(Hsp_gen; Ops=[
                                                       :H,
                                                       :Jx,
                                                       :Jy,
                                                       :Jz
                                                      ],
                                         #da=[1.0,0.0],db=[0,1.0,0.0],
                                         tbc_theta=tbc_theta,
                                         #save_ascii="ops-$(sample_n).dat"
                                        );
Hsp = All_ops[1]
#Jx = All_ops[2]
Jx = All_ops[2];Jy = All_ops[3];Jz = All_ops[4]
Hsp = (Hsp + Hsp')/2
#@time a, H_norm = SuperLattice.Util.normalizeH(Hsp)# if need a fixed bandwidth, use s
=#

fn = "/mnt/home/awu10/ceph/CPGEdata/spH0L80M2gam0.jld2"
H0 = load(fn,"H0")
Ls = load(fn, "Ls")
L = Ls[1]
@show nnz(H0)
@time H2 = H0*H0
@show nnz(H2)
@time es, vs = eigs(H0;nev=1,tol=0.001,maxiter=300)
@show es
@time es, vs = eigsolve(H0,1,:LM)
@show es
#@time es, vs = eigs(H0;nev=1,which=:SM,tol=0.001,maxiter=300)

@time es, vs = eigsolve(H2,1,:SR)
@show es
@time es, vs = eigs(H2;nev=1,which=:SR)#,tol=0.001,maxiter=300)
@show es

#@time es, vs = eigs(H2;nev=1,which=:SM)
#@show es
tks = [0.5,0.5,0.5]./L
Eminana = dispersion(tks[1],tks[2],0.25-tks[3],1,2,0)
E2min = Eminana[1]^2

#=
@time Es, _ = eigen(Matrix(Hsp))
@time E2s, _ = eigen(Matrix(Hsp*Hsp))

@show Esmin = minimum(abs.(Es))
@show E2smin = minimum(E2s)

tks = [0.5,0.5,0.5]./L
Eminana = dispersion(tks[1],tks[2],0.25-tks[3],1,2,0)
display(Eminana)

@time es, vs = eigsolve(Hsp*Hsp,1,:SR)
@show es
@time es, vs = eigs(Hsp*Hsp;nev=1,which=:SM)
@show es
=#
