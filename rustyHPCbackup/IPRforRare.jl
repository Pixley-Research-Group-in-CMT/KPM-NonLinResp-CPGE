using Random
using HDF5

include("model.jl")

function IPR(Vs)
        n = size(Vs)[1]
        res = zeros(n,1)
        for k in 1:n
                res[k] = sum((abs.(Vs[:,k])).^4)
        end
        return res
end

L = 18
Ls = [L,L,L]
M = 2.0 # 2.0
t = 1.0
γ = 0.4 # 0.8
@show Emax = (M+2)+sqrt(1+γ^2)

f = h5open("/mnt/home/awu10/ceph/CPGEdata/E2sL18M15W05gamma0Sam1050.h5","r")
E2s = read(f,"E2s")
disorders = read(f,"disorders")
close(f)

disind = argmin(E2s)[2]
boundary = [1,1,1]
H1 = TwoBandModel(Ls, t, M, γ;boundary=boundary)
Ham = Matrix(H1 + spdiagm(disorders[:,disind]))
@time EsPBC, VsPBC = eigen(Ham)
@time iprPBC = IPR(VsPBC)

boundary = [-1,-1,-1]
H1 = TwoBandModel(Ls, t, M, γ;boundary=boundary)
Ham = Matrix(H1 + spdiagm(disorders[:,disind]))
@time EsTBC, VsTBC = eigen(Ham)
@time iprTBC = IPR(VsTBC)

f = h5open("EigenStatesRareM2g04TPBC.h5","w")
write(f,"EsPBC",EsPBC)
#write(f, "VsPBC",VsPBC)
write(f, "EsTBC",EsTBC)
#write(f,"VsTBC",VsTBC)
write(f,"iprPBC",iprPBC)
write(f,"iprTBC", iprTBC)
close(f)

