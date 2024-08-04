using SuperLattice
using SuperLattice.Util.Physics, SuperLattice.Util
using Statistics
using LinearAlgebra, Arpack, SparseArrays
using KPM
using JLD2
using CUDA

@show CUDA.has_cuda()
@show KPM.whichcore()
# check number of threads
println("Number of threads: $(Threads.nthreads())")

αx, αy, αz, β, id = gammaMatrices(3; rep=:Dirac)

# model parameters
t=1.0
gamma=0.8
M=2
L=30 #30
sizes=[L,L,L]
W=0.0

######## random function
rf = SuperLattice.Util.rand_gen(prod(sizes); rf=randn, offset=true) # random gaussian disorder that offset to 0
########

OBC=false
tbc_theta = zeros(3) #rand(3)
println("theta=$(tbc_theta)")

# KPM
NC= 32#512
Ntilde = NC*6
NR= 1 #25 # number of random vectors

N_rep = 10

# KPM running params
arr_size = NC ÷ 2 # may want to adjust this according to mem availability

println("define model")
uc = UnitCell(;d=3)
addAtom(uc, :a, 2, [0,0,0]) # 2 band model

# potential term - disorder
Vr = function(x, y, v)
    v .= W * rf() * σ0 + M * σz #+ chemmu * σ0 chemical potential enters in the final calculation
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


All_ops = SuperLattice.to_sparse_matrix(Hsp_gen; Ops=[
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
a, H_norm = SuperLattice.Util.normalizeH(Hsp)# if need a fixed bandwidth, use setA=a_fixed
#size(All_ops)



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





#jldsave("GammamnpL30NC512W0NR125.jld2"; L=L, NC=NC, D = a, Gammaxyz=mu_3d_xyz, Gammayxz=mu_3d_yxz)




