using SuperLattice
using SuperLattice.Util.Physics, SuperLattice.Util
using Statistics
using LinearAlgebra, Arpack, SparseArrays
using KrylovKit
using JLD2

include("model.jl")

L = 12
Ls = [L,L,L]
M = 2.0 # 2.0
t = 1.0
γ = 0.0  # 0.8
twists = pi*ones(3)#[0.0,0,0] #[pi,pi,pi]
boundary = exp.(1im.*twists)
@time H0 = TwoBandModel(Ls, t, M, γ;boundary=boundary)

Eana = zeros(2*L^3)
tks = twists #[0.0,0.0,0.0] # twists
tks = tks./L/2/pi
for x in 1:L, y in 1:L, z in 1:L
        kx = 1.0*(x-1)/L
        ky = 1.0*(y-1)/L
        kz = 1.0*(z-1)/L
        pos = (x-1)*L^2+(y-1)*L+z
        E1, E2 = dispersion(kx+tks[1],ky+tks[2],kz+tks[3],t,M,γ)
        Eana[2*pos-1] = E1
        Eana[2*pos] = E2
end
Eana = sort(Eana)

@time Emy, _ = eigen(Matrix(H0))
@show maximum(Eana-Emy)



OBC=false
tbc_theta = [0.5,0.5,0.5] #pi*ones(3) #zeros(3) #rand(3)
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
addHop(uc, :a, :a, [0, 0, 1], 0.5im * γ * σ0 - 0.5 * t * σz; hc=true)

# set up lattice
ltc = Lattice(d=3,pv=[[1.,0,0],[0.,1,0],[0.,0,1]],sizes=Ls,OBC=OBC);
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

@time Eyx, _ = eigen(Matrix(Hsp))
@show maximum(Eana-Eyx)

es, vs = eigsolve(Hsp*Hsp,1,:SR)
@show es
es, vs = eigs(Hsp;nev=1,which=:SM)
@show es
es, vs = eigs(Hsp*Hsp;nev=1,which=:SR)
@show es

Emin = minimum(abs.(Eana))
Eminana = dispersion(tks[1],tks[2],0.25-tks[3],t,M,γ)

