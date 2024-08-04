using Random
using HDF5

include("model.jl")
include("EDcpge.jl")

function TBCcpge1(L, M, t, γ, samples;h=0.002, Ω=0.0,σ=0.01,η=0.01)
	xs = collect(-1:h:1)
	Nx = size(xs)[1]
	res = zeros(Nx,samples)
	Emax = (M+2)+sqrt(1+γ^2)
	w1 = 0.8/Emax
        w2 = -w1
        w1 += Ω
	Threads.@threads for k in 1:samples
		boundary = exp.(1im*2*pi*rand(3))
		H = TwoBandModel(L,t,M,γ;boundary=boundary)
		H = H./Emax
		#Es, Vs = eigen(Matrix(H))
		hx, hy, hz = velocityOperatorNew(H)
		Es, hxp, hyp, hzp = EDvelocity(Matrix(H), hx, hy, hz)
		temp = zeros(ComplexF64,Nx)
		for (i, x) in enumerate(xs)
        		temp[i] = EDcpge(x, w1, w2, Es, hxp, hyp, hzp;σ=σ,η=η)+EDcpge(x, w2, w1, Es, hyp, hxp, hzp;σ=σ, η=η)
		end
		y_int = cumsum(temp)*h
		res[:,k] = imag(y_int)
		println(k)
	end
	return xs, res
end

L = 8
M = 2.0
t = 1.0
γ = 0.8
samples = 100
xs, cpges = TBCcpge1(L, M, t, γ, samples)

Emax = (M+2)+sqrt(1+γ^2)
w1 = 0.8/Emax
w2 = -w1
w1 += 0.0
f = h5open("EDcpgeL8TBC100O0sigman2.h5","w")
write(f, "xs", xs)
write(f,"cpges",cpges)
write(f,"w1", w1)
write(f,"w2", w2)
close(f)
