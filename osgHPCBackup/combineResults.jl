using JLD2, Statistics

let
pt = "/public/angkunwu/CPGEdata/"
suffix = ".jld2"

L = 100
Ls = [L,L,L]
NC = 512
D = 6

mu_3d_xyz = zeros(ComplexF64, NC, NC, NC)
mu_3d_yxz = zeros(ComplexF64, NC, NC, NC)


stash_destination="stash:///osgconnect/public/angkunwu/CPGEdata/"
pt = "/home/angkunwu/"
for disind in 0:49
	prefix = "GammaxyzL100gam08D6NC512W2002NR1Sam"
	fn = prefix*"$disind"*suffix
	path1 = stash_destination*fn
	path2 = pt*fn
	run(`stashcp $path1 $path2`)
	temp = load(fn, "Gammaxyz")
	mu_3d_xyz += temp
	rm("/home/angkunwu/"*fn)
	println(disind)
	prefix = "GammayxzL100gam08D6NC512W2002NR1Sam"
	fn = prefix*"$disind"*suffix
	path1 = stash_destination*fn
        path2 = pt*fn
        run(`stashcp $path1 $path2`)
	temp = load(fn, "Gammayxz")
	mu_3d_yxz += temp
	rm("/home/angkunwu/"*fn)
	println(disind)
end
mu_3d_xyz = mu_3d_xyz ./ 50
mu_3d_yxz = mu_3d_yxz ./ 50

fn = "GammamnpL100gam08D6NC512W2002NR1Sams50.jld2"
jldsave(fn; Ls=Ls, NC=NC, D = D, Gammaxyz=mu_3d_xyz, Gammayxz=mu_3d_yxz)

end
