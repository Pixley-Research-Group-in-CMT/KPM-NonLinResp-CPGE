using HDF5
using Dates
using LinearAlgebra
@show BLAS.get_num_threads()
BLAS.set_num_threads(1)
@show BLAS.get_num_threads()

include("model.jl")
include("EDcpge.jl")
@show Threads.nthreads()

Ls = [5,5,100]
M = 2.0
t = 1.0
γ = 0.8
Emax = (M+2)+sqrt(1+γ^2)

σ = 0.01
η = 0.01
Ω = 10^(-2)

H0 = TwoBandModel(Ls,t,M,γ)
H = H0./Emax
hx, hy, hz = velocityOperator(Ls,H0)
@time Es, hxp, hyp, hzp = EDvelocity(Matrix(H), hx, hy, hz)

w1 = 0.8/Emax
w2 = -w1
w1 += Ω

h = 0.001
x_all = collect(-0.5:h:0.5)
y_xyz = complex(x_all)
npoints = size(x_all)[1]
#npoints = 25
Threads.@threads for k in 1:npoints
	#@show starttime = now()
	x = x_all[k]
	@time y_xyz[k] = EDcpge(x, w1, w2, Es, hxp, hyp, hzp;σ=σ,η=η)
	#@show endtime = now()
	@show k
end
y_yxz = complex(x_all)
Threads.@threads for k in 1:npoints
        #@show starttime = now()
        x = x_all[k]
        @time y_yxz[k] = EDcpge(x, w2, w1, Es, hyp, hxp, hzp;σ=σ, η=η)
        #@show endtime = now()
        @show k
end
y_all = y_xyz .+ y_yxz

y_int = cumsum(y_all)*h
volume = Ls[1]*Ls[2]*Ls[3]
cpge = imag.(y_int)./volume

f = h5open("EDcpgeLz100PBCOn2etan2.h5","w")
write(f, "xs", x_all)
write(f, "y_all",y_all)
write(f,"cpge",cpge)
write(f,"w1", w1)
write(f,"w2", w2)
write(f,"Omega",Ω)
close(f)

