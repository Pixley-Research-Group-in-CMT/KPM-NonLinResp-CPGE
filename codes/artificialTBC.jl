include("model.jl")

L = 18
samples = 100
t = 1.0
M = 2.0
γ = 0.4

#energies = zeros(2*L^3,samples)
Esx = zeros(2*L^3,samples)
Esz = zeros(2*L^3,samples)
thetas = collect(1:samples)/samples
for x in 1:L, y in 1:L, z in 1:L
	kx = 1.0*(x-1)/L
	ky = 1.0*(y-1)/L
	kz = 1.0*(z-1)/L
	pos = (x-1)*L^2+(y-1)*L+z
	for phi in 1:samples
		shift = (phi-1)/samples/L
		E1, E2 = dispersion(kx+shift,ky,kz,t,M,γ)
		Esx[2*pos-1, phi] = E1
		Esx[2*pos, phi] = E2
		E1, E2 = dispersion(kx,ky,kz+shift,t,M,γ)
                Esz[2*pos-1, phi] = E1
                Esz[2*pos, phi] = E2
	end
end

for k=1:samples
	Esx[:,k] = sort(Esx[:,k])
	Esz[:,k] = sort(Esz[:,k])
end

