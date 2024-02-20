include("model.jl")
include("EDcpge.jl")

L = 300
M = 2.0 # 2.0
t = 1.0
γ = 0.8 # 0.8

energyp,energyn = Float64[],Float64[]
for x in 1:L, z in 1:L, y in 1:L
	kx = 2*pi*(x-1)/L
        ky = 2*pi*(y-1)/L
        kz = 2*pi*(z-1)/L
	Ep,En = dispersion(kx,ky,kz, t,M,γ)
	if abs(Ep-En-0.8) < 0.01
		push!(energyp,Ep)
		push!(energyn,En)
        end
end
sort!(energyp)
sort!(energyn)
