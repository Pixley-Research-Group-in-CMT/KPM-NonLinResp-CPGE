using Statistics
using JLD2
using KPM

pt = "/scratch/aw666/CPGEdata/"
NC = 2^15

function simplesum(NC)
	musum = zeros(NC)
	for k = 1:100
		fn = pt*"DOSmuL100gam08D6NC215NR10Sam"*"$k"*".jld2"
		mu = load(fn,"mu")
		musum += mu
		println(k)
	end
	musum = musum/100
	return musum
end

mu = simplesum(NC)
D = 6.0

pt = "/home/aw666/data/CPGEdata/"
fn = pt*"DOSmuL100gam08D6NC215NR10TBC100.jld2"
jldsave(fn;NC=NC, D = D, mu=mu)

