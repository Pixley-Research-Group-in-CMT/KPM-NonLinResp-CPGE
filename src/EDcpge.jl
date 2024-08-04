"""
	Use Exact diagonalization for cpge
"""
function EDvelocity(Ham, hx, hy, hz)
	Es, Vs = eigen(Ham)
        hxp = Vs'*hx*Vs
        hyp = Vs'*hy*Vs
        hzp = Vs'*hz*Vs
	return Es, hxp, hyp, hzp
end

function DeltaH(ϵ, Es, η; broaden="Lorentzian")
	N = size(Es)[1]
        res = spzeros(N,N)
	if broaden == "Lorentzian"
		for k in 1:N
			res[k,k] = (η/pi)/((ϵ-Es[k])^2+η^2)
        	end
	elseif broaden == "Gaussian"
        	for k in 1:N
                	res[k,k] = exp(-((ϵ-Es[k])/η)^2)/(sqrt(2*pi)*η)
        	end
	else
		println("wrong type of broadening!")
		return
	end
        return res
end

function Green(ϵ, Es, η, type)
        N = size(Es)[1]
        if type == "A"
                flag = -1
        elseif type == "R"
                flag = 1
        else
                return println("wrong Green type!")
        end
        res = spzeros(ComplexF64, N, N)
        for k = 1:N
                res[k,k] = 1/(ϵ-Es[k]+flag*1im*η)
        end
        return res
end

function EDcpge(ϵ, ω1, ω2, Es, hxp, hyp, hzp; σ=0.01, η=0.01)
	#ω1, ω2 = ω, Ω-ω
	N = size(Es)[1]
	Ω = ω1 + ω2
	Mdelta = DeltaH(ϵ, Es, σ)
	G11, G12 = Green(ϵ+Ω,Es, η,"R"), Green(ϵ+ω2,Es, η,"R")
	G21, G22 = Green(ϵ+ω1,Es, η,"R"), Green(ϵ-ω2,Es, η,"A")
	G31, G32 = Green(ϵ-ω1,Es, η,"A"), Green(ϵ-Ω,Es, η,"A")
	resmat = spzeros(ComplexF64,N,N)
	resmat += hxp*G11*hyp*G12*hzp*Mdelta
	resmat += hxp*G21*hyp*Mdelta*hzp*G22
	resmat += hxp*Mdelta*hyp*G31*hzp*G32
	return tr(resmat)./(ω1*ω2)
end
