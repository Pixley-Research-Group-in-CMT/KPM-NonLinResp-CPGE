using LinearAlgebra

"""
Compute the Green's function with KPM moments with Lorentz kernel
"""
function GreensFun(w::Float64, μ::Vector{Number};λ::Float64=0.5)
	NC = size(μ)[1]
	res = μ[1]
	for k = 2:NC
		res += μ[k]*sinh(λ(1-(k-1)/NC))/sinh(λ)
	end
	return res
end


function generatekstate(Ls::Vector{Number},ks::Vector{Number};spinor::Number=1)
	N = 2*Ls[1]*Ls[2]*Ls[3]
	V = N/2
	res = zeros(ComplexF64,N)
	if spinor > 2
		println("Wrong spinor!")
		return
	end
	for x in 1:Ls[1], y in 1:Ls[2], z in 1:Ls[3]
		pos = (x-1)*Ls[2]*Ls[3] + (y-1)*Ls[3] + z
		res[2*pos-spinor+1] = exp(-1im*ks.*[x,y,z])/sqrt(V)
	end
	return res
end




