using SuperLattice
using SuperLattice.Util.Physics, SuperLattice.Util
using Statistics
using LinearAlgebra, Arpack, SparseArrays
using JLD2

@show Threads.nthreads()
#include("model.jl")
# Should only use 1 thread!!!!!!!
L = 100
Ls = [L,L,L]
M = 2.0 # 2.0
t = 1.0
gamma = 0.1  # 0.8
#twists = [pi,pi,pi]
#boundary = exp.(1im.*twists)
#@time H0 = TwoBandModel(Ls, t, M, γ;boundary=boundary)

OBC=false
tbc_theta = zeros(3) #0.5*ones(3) zeros(3) rand(3) # twist angles in [0,1] \pi implicit
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
addHop(uc, :a, :a, [0, 0, 1], 0.5im * gamma * σ0 - 0.5 * t * σz; hc=true)

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

save("spH0L100M2gam01PBC.jld2","H0",Hsp,"Jx",Jx,"Jy",Jy,"Jz",Jz,"Ls",Ls,"M",M,"gamma",gamma,"boundary",tbc_theta)



