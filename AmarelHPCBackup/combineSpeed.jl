using JLD2, Statistics
#=
let
	pt = "/scratch/aw666/CPGEdata/"
	suffix = ".jld2"

	L = 100
	runtimes = zeros(100)
	for runsam in 1:100
		prefix = "Gamma$(L)NC27GPUrunsam$(runsam).jld2"
		fn = pt*prefix
		temp = load(fn, "runtime")
		runtimes[runsam] = temp
		println(temp," ",runsam)
	end
	pt = "/home/aw666/data/CPGEdata/"
	fn = pt*"Gamma$(L)NC27GPUrunsams100.jld2"
	jldsave(fn; L=L, runtimes=runtimes)
end
=#

let 
	pt = "/home/aw666/data/CPGEdata/"
	Ls = [20,40,60,80,100]
	runtimesLs = zeros(100,5)
	for Lind in 1:5
		L = Ls[Lind]
		prefix = "Gamma$(L)NC27GPUrunsams100.jld2"
		fn = pt*prefix
		temp = load(fn, "runtimes")
		runtimesLs[:,Lind] = temp
		println(sum(temp)/100," ",Lind)
	end
	fn = pt*"GammaLsNC27GPUrunsams100.jld2"
	jldsave(fn; Ls=Ls, runtimes=runtimesLs)
end

