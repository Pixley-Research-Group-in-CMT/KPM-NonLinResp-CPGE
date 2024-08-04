using LinearAlgebra
using SparseArrays
using Random, Distributions

"""
	Generate the two band Weyl semimetal
	-t hopping parameter
	-M stagger onsite potential
	-γ z-dir hopping
	implicity parameters:
	-boundary [1,1,1] for PBC, [0,0,0] for OBC, and TBC with angles
	-disorder directly input onsite disorders so as to keep specific ones
"""
function TwoBandModel(Ls::Vector{Int64}, t::Number, M::Number, γ::Number;
                boundary::Vector{<:Number} = [1,1,1])
        Lx, Ly, Lz = Ls[1], Ls[2], Ls[3]
        Ix,Jy = Int64[], Int64[]
        Vals = ComplexF64[]
	t,γ = -t, -γ # to match Yixing's construction
        σ0 = Matrix(1.0I,2,2)
        σx, σy, σz = zeros(2,2), zeros(ComplexF64,2,2), zeros(2,2)
        σx[1,2], σx[2,1] = 1.0, 1.0
        σy[1,2], σy[2,1] = -1im, 1im
        σz[1,1], σz[2,2] = 1.0, -1.0

        function cellpos(x, y, z)
                return (x-1)*Ly*Lz + (y-1)*Lz + z
        end
        function addToList!(pos1,pos2,value,Ix,Jy,Vals)
                push!(Ix, pos1)
                push!(Jy, pos2)
                push!(Vals, value)
                return
        end
	function addByUnitCells!(cell1,cell2,mat,Ix,Jy,Vals)
                for k1 in 1:2, k2 in 1:2
                        if mat[k1,k2] != 0
                                addToList!(2*cell1+k1-2,2*cell2+k2-2,mat[k1,k2],Ix,Jy,Vals)
                        end
                end
                return
        end
	for x in 1:Lx, y in 1:Ly,z in 1:Lz
                pos = cellpos(x,y,z)
                addByUnitCells!(pos,pos,-M*σz,Ix,Jy,Vals)
                # x direction
                temp = 0.5*t
                if x + 1 > Lx
                        pos2 = cellpos(1,y,z)
                        temp *= boundary[1]
                else
                        pos2 = cellpos(x+1, y, z)
                end
                addByUnitCells!(pos2,pos,temp*(σz+1im*σx),Ix,Jy,Vals)
                addByUnitCells!(pos,pos2,(temp*(σz+1im*σx))',Ix,Jy,Vals)
                # y direction
                temp = 0.5*t
                if y + 1 > Ly
                        pos2 = cellpos(x,1,z)
                        temp *= boundary[2]
                else
                        pos2 = cellpos(x,y+1,z)
                end
                addByUnitCells!(pos2,pos,temp*(σz+1im*σy),Ix,Jy,Vals)
                addByUnitCells!(pos,pos2,(temp*(σz+1im*σy))',Ix,Jy,Vals)
		# z direction
                temp = 0.5
                if z + 1 > Lz
                        pos2 = cellpos(x,y,1)
                        temp *= boundary[3]
                else
                        pos2 = cellpos(x,y,z+1)
                end
                addByUnitCells!(pos2,pos,temp*(t*σz-1im*γ*σ0),Ix,Jy,Vals)
                addByUnitCells!(pos,pos2,(temp*(t*σz-1im*γ*σ0))',Ix,Jy,Vals)
        end
        return sparse(Ix,Jy,Vals)
end

function TwoBandModelSlow(Ls::Vector{Int64}, t::Number, M::Number, γ::Number;
		boundary::Vector{<:Number} = [1,1,1])
	#@show boundary
	Lx, Ly, Lz = Ls[1], Ls[2], Ls[3]
	N = 2*Lx*Ly*Lz
	Ham = spzeros(ComplexF64, N, N)
	
	σ0 = Matrix(1.0I,2,2)
	σx, σy, σz = zeros(2,2), zeros(ComplexF64,2,2), zeros(2,2)
	σx[1,2], σx[2,1] = 1.0, 1.0
	σy[1,2], σy[2,1] = -1im, 1im
	σz[1,1], σz[2,2] = 1.0, -1.0
	
	function cellpos(x, y, z)
		return (x-1)*Ly*Lz + (y-1)*Lz + z
	end
	for x in 1:Lx, y in 1:Ly, z in 1:Lz
		pos = cellpos(x,y,z)
		Ham[(2*pos-1):(2*pos),(2*pos-1):(2*pos)] += -M*σz
		# x direction
		temp = 0.5*t
		if x + 1 > Lx
			pos2 = cellpos(1,y,z)
			temp *= boundary[1]
		else
			pos2 = cellpos(x+1, y, z)
		end
		Ham[(2*pos2-1):(2*pos2),(2*pos-1):(2*pos)] += temp*(σz+1im*σx)
		Ham[(2*pos-1):(2*pos),(2*pos2-1):(2*pos2)] += (temp*(σz+1im*σx))'
		# y direction
		temp = 0.5*t
		if y + 1 > Ly
			pos2 = cellpos(x,1,z)
			temp *= boundary[2]
		else
			pos2 = cellpos(x,y+1,z)
		end
		Ham[(2*pos2-1):(2*pos2),(2*pos-1):(2*pos)] += temp*(σz+1im*σy)
                Ham[(2*pos-1):(2*pos),(2*pos2-1):(2*pos2)] += (temp*(σz+1im*σy))'
                # z direction
		temp = 0.5
		if z + 1 > Lz
			pos2 = cellpos(x,y,1)
			temp *= boundary[3]
		else
			pos2 = cellpos(x,y,z+1)
		end
		Ham[(2*pos2-1):(2*pos2),(2*pos-1):(2*pos)] += temp*(t*σz-1im*γ*σ0)
                Ham[(2*pos-1):(2*pos),(2*pos2-1):(2*pos2)] += (temp*(t*σz-1im*γ*σ0))'
	end
	return Ham
end

function dispersion(kx::Number, ky::Number, kz::Number, t::Number, M::Number, γ::Number)
	kx *= 2*pi
	ky *= 2*pi
	kz *= 2*pi
	res = sqrt((M-t*(cos(kx)+cos(ky)+cos(kz)))^2+t^2*(sin(kx)^2+sin(ky)^2))
	res2 = -res
	res += γ*sin(kz)
	res2 += γ*sin(kz)
	return res, res2
end

function spectrum(Ls::Vector{Int64}, t::Number, M::Number, γ::Number)
	Lx, Ly, Lz = Ls[1], Ls[2], Ls[3]
	N = 2*Lx*Ly*Lz
	res = zeros(N)
	for x in 1:Lx, y in 1:Ly, z in 1:Lz
		kx = 2*pi*(x-1)/Lx
		ky = 2*pi*(y-1)/Ly
		kz = 2*pi*(z-1)/Lz
		pos = (x-1)*Ly*Lz + (y-1)*Lz + z
		res[2*pos-1] = sqrt((M-t*(cos(kx)+cos(ky)+cos(kz)))^2+t^2*(sin(kx)^2+sin(ky)^2))
		res[2*pos] = -res[2*pos-1]
		res[2*pos-1] += γ*sin(kz)
		res[2*pos] += γ*sin(kz)
	end
	return sort(res)
end

function velocityOperator(Ls::Vector{Int64}, Ham)
        NZdictx, NZdicty, NZdictz = Dict(), Dict(), Dict()
        Lx, Ly, Lz = Ls[1], Ls[2], Ls[3]
        Ix, Jx, Iy, Jy, Iz, Jz = Int64[], Int64[], Int64[], Int64[], Int64[], Int64[]
        Valx, Valy, Valz = ComplexF64[],ComplexF64[],ComplexF64[]
        function cellpos(x, y, z)
                return (x-1)*Ly*Lz + (y-1)*Lz + z
        end
        function addToList!(pos1,pos2,value,I,J,Vals)
                push!(I, pos1)
                push!(J, pos2)
                push!(Vals, value)
                return
        end
        function addByUnitCells!(cell1,cell2,mat,NZdict)
                for k1 in 1:2, k2 in 1:2
                        if mat[k1,k2] != 0
                                x = 2*cell1+k1-2
                                y = 2*cell2+k2-2
                                if haskey(NZdict,(x,y))
                                        NZdict[(x,y)] += mat[k1,k2]
                                else
					NZdict[(x,y)] = mat[k1,k2]
                                end
                        end
                end
                return
        end
	for x in 1:Lx, y in 1:Ly,z in 1:Lz
                pos = cellpos(x,y,z)
                if x + 1 > Lx
                        pos2 = cellpos(1,y,z)
                else
                        pos2 = cellpos(x+1, y, z)
                end
                mat = Ham[(2*pos2-1):(2*pos2),(2*pos-1):(2*pos)]
                addByUnitCells!(pos2,pos,mat,NZdictx)
                addByUnitCells!(pos,pos2,-mat',NZdictx)
                if y + 1 > Ly
                        pos2 = cellpos(x,1,z)
                else
                        pos2 = cellpos(x, y+1, z)
                end
                mat = Ham[(2*pos2-1):(2*pos2),(2*pos-1):(2*pos)]
                addByUnitCells!(pos2,pos,mat,NZdicty)
                addByUnitCells!(pos,pos2,-mat',NZdicty)
		if z + 1 > Lz
                        pos2 = cellpos(x,y,1)
                else
                        pos2 = cellpos(x, y, z+1)
                end
                mat = Ham[(2*pos2-1):(2*pos2),(2*pos-1):(2*pos)]
                addByUnitCells!(pos2,pos,mat,NZdictz)
                addByUnitCells!(pos,pos2,-(mat'),NZdictz)
        end
        for (x, y) in keys(NZdictx)
                addToList!(x,y,NZdictx[(x,y)],Ix,Jx,Valx)
        end
        hx = -1im*sparse(Ix,Jx,Valx)
	for (x, y) in keys(NZdicty)
                addToList!(x,y,NZdicty[(x,y)],Iy,Jy,Valy)
        end
        hy = -1im*sparse(Iy,Jy,Valy)
        for (x, y) in keys(NZdictz)
                addToList!(x,y,NZdictz[(x,y)],Iz,Jz,Valz)
        end
        hz = -1im*sparse(Iz,Jz,Valz)
        return hx, hy, hz
end

#=
function velocityOperator(Ham)
	N = size(Ham)[1]
	L = Int(cbrt(N/2))
	function cellpos(x, y, z)
                return (x-1)*L^2 + (y-1)*L + z
        end
	rx, ry, rz = spzeros(ComplexF64, N, N), spzeros(ComplexF64, N, N), spzeros(ComplexF64, N, N)
        for x in 1:L, y in 1:L, z in 1:L
		pos = cellpos(x,y,z)
		rx[2*pos-1,2*pos-1] = x
		rx[2*pos, 2*pos] = x
		ry[2*pos-1,2*pos-1] = y
		ry[2*pos,2*pos] = y
		rz[2*pos-1,2*pos-1] = z
		rz[2*pos, 2*pos] = z
	end
	hx = -1im*(rx*Ham-Ham*rx)
	hy = -1im*(ry*Ham-Ham*ry)
	hz = -1im*(rz*Ham-Ham*rz)
	return hx, hy, hz
end

function velocityOperatorSlow(Ls::Vector{Int64}, Ham)
	Lx, Ly, Lz = Ls[1], Ls[2], Ls[3]
	N = size(Ham)[1]
        function cellpos(x, y, z)
                return (x-1)*Ly*Lz + (y-1)*Lz + z
        end
	hx, hy, hz = spzeros(ComplexF64, N, N), spzeros(ComplexF64, N, N), spzeros(ComplexF64, N, N)
	for x in 1:Lx, y in 1:Ly, z in 1:Lz
		pos = cellpos(x,y,z)
		if x + 1 > Lx
                        pos2 = cellpos(1,y,z)
                else
                        pos2 = cellpos(x+1, y, z)
                end
		hx[(2*pos2-1):(2*pos2),(2*pos-1):(2*pos)] += Ham[(2*pos2-1):(2*pos2),(2*pos-1):(2*pos)]
		hx[(2*pos-1):(2*pos),(2*pos2-1):(2*pos2)] -= Ham[(2*pos-1):(2*pos),(2*pos2-1):(2*pos2)]
		if y + 1 > Ly
                        pos2 = cellpos(x,1,z)
                else
                        pos2 = cellpos(x, y+1, z)
                end
                hy[(2*pos2-1):(2*pos2),(2*pos-1):(2*pos)] += Ham[(2*pos2-1):(2*pos2),(2*pos-1):(2*pos)]
                hy[(2*pos-1):(2*pos),(2*pos2-1):(2*pos2)] -= Ham[(2*pos-1):(2*pos),(2*pos2-1):(2*pos2)]
		if z + 1 > Lz
                        pos2 = cellpos(x,y,1)
                else
                        pos2 = cellpos(x, y, z+1)
                end
                hz[(2*pos2-1):(2*pos2),(2*pos-1):(2*pos)] += Ham[(2*pos2-1):(2*pos2),(2*pos-1):(2*pos)]
                hz[(2*pos-1):(2*pos),(2*pos2-1):(2*pos2)] -= Ham[(2*pos-1):(2*pos),(2*pos2-1):(2*pos2)]
	end
	return -1im*hx,-1im*hy,-1im*hz
end
=#
function pathInMomentum(cor1, cor2; points = 100)
	res = zeros(3, points)
	for k in 1:3
		res[k,:] = collect(range(cor1[k],cor2[k],points))
	end
	return res
end

function generateDisorder(Ls::Vector{Int64},W::Number)
	Lx, Ly, Lz = Ls[1], Ls[2], Ls[3]
        N = 2*Lx*Ly*Lz
	disorder = W*rand(Normal(0.0, 1.0),N)
	avg = mean(disorder)
	return (disorder .- avg)
end

function generatePotDisorder(Ls::Vector{Int64},W::Number)
        Lx, Ly, Lz = Ls[1], Ls[2], Ls[3]
        N = Lx*Ly*Lz
	L = 2*N
        disUnit = W*rand(Normal(0.0, 1.0),N)
        avg = mean(disUnit)
	disUnit = disUnit .- avg
	res = zeros(L)
	for k = 1:N
		res[2*k-1] = disUnit[k]
		res[2*k] = disUnit[k]
	end
        return res
end


