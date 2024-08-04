using JLD2, Statistics

let
pt = "/scratch/aw666/CPGEdata/"
suffix = ".jld2"

L = 18#100
Ls = [L,L,L]
NC = 512
D = 6
W2 = "025"

mu_3d_xyz = zeros(ComplexF64, NC, NC, NC)
mu_3d_yxz = zeros(ComplexF64, NC, NC, NC)
#=
for disind in 1:26
	prefix = "GammamnpL100gam08D6NC512W01NR1Sam"
	fn = pt*prefix*"$disind"*suffix
	temp = load(fn, "Gammaxyz")
	mu_3d_xyz += temp
	temp = load(fn, "Gammayxz")
	mu_3d_yxz += temp
	println(disind)
end
=#
Nsam = 100
for disind in 1:Nsam
	GC.@preserve mu_3d_xyz mu_3d_yxz begin
		#prefix = "GammaxyzL100gam08D6NC512W2"*W2*"NR1Sam"
		prefix = "GammaxyzL$(L)gam03D6NC512W05NR100Sam"
		fn = pt*prefix*"$disind"*suffix
		temp = load(fn, "Gammaxyz")
		mu_3d_xyz += temp
		temp = nothing  # Free the temporary variable
		#prefix = "GammayxzL100gam08D6NC512W2"*W2*"NR1Sam"
		prefix = "GammayxzL$(L)gam03D6NC512W05NR100Sam"
		fn = pt*prefix*"$disind"*suffix
		temp = load(fn, "Gammayxz")
		mu_3d_yxz += temp
		temp = nothing  # Free the temporary variable
	end
	GC.gc()
	println(disind)
end
mu_3d_xyz = mu_3d_xyz ./ Nsam
mu_3d_yxz = mu_3d_yxz ./ Nsam

pt = "./CPGEdata/"  # path to save data
#pt = "/home/aw666/data/CPGEdata/"
#fn = pt*"GammamnpL100gam08D6NC512PotW2"*W2*"NR1Sams50.jld2"
fn = pt*"GammamnpL$(L)gam03D6NC512PotW05NR100TBCs$(Nsam).jld2"
jldsave(fn; Ls=Ls, NC=NC, D = D, Gammaxyz=mu_3d_xyz, Gammayxz=mu_3d_yxz)

end
