using JLD2, Statistics
#=
let
	pt = "/scratch/aw666/CPGEdata/"
	
	Nsam,Nvs = 100,5
	L = 100
	disorders = zeros(2*L^3,Nsam)
	E2s = zeros(Nvs, Nsam)
	IPRs = zeros(Nvs, Nsam)
	for k = 1:Nsam
		prefix = "E2s5L100M2W05gamma0Sam1Ind$(k)MyTBC.jld2"
		fn = pt*prefix
		disorders[:,k] = load(fn,"disorders")[:,1]
		E2s[:,k] = load(fn,"E2s")[:,1]
		IPRs[:,k] = load(fn,"IPRs")[:,1]
		println(E2s[1,k]," ", k)
	end
	pt = "/home/aw666/data/CPGEdata/"
	fn = pt*"E2s5L100M2W05gamma0Sam$(Nsam)TBCnoDis.jld2"
	#jldsave(fn;disorders=disorders,E2s=E2s,IPRs=IPRs)
	jldsave(fn;E2s=E2s,IPRs=IPRs)
end
=#
let
	pt = "/scratch/aw666/CPGEdata/"
	Nsam,Nvs = 100,5
	E2sx = zeros(Nvs,Nsam)
	E2sy = zeros(Nvs,Nsam)
	E2sz = zeros(Nvs,Nsam)
	IPRsx,IPRsy,IPRsz = zeros(Nvs,Nsam),zeros(Nvs,Nsam),zeros(Nvs,Nsam)
	disind = zeros(1)
	for k = 1:Nsam
		fn = pt*"E2sIPRrareL100M2W05gamma01Twist$(k).jld2"
		disind[1] = load(fn,"disind")
		E2sxyz = load(fn,"E2sxyz")
		E2sx[:,k],E2sy[:,k],E2sz[:,k]=E2sxyz[:,1],E2sxyz[:,2],E2sxyz[:,3]
		IPRsxyz = load(fn,"IPRsxyz")
		IPRsx[:,k] = IPRsxyz[:,1]
		IPRsy[:,k] = IPRsxyz[:,2]
		IPRsz[:,k] = IPRsxyz[:,3]
		println(E2sx[1,k]," ",k)
	end
	pt = "/home/aw666/data/CPGEdata/"
	#fn = pt*"E2sIPRrareL100M2W05gamma0Twists100.jld2"
	fn = pt*"E2sIPRrareL100M2W05gamma01Twists100.jld2"
	jldsave(fn;disind=disind[1],E2sx=E2sx,E2sy=E2sy,E2sz=E2sz,
		IPRsx=IPRsx,IPRsy=IPRsy,IPRsz=IPRsz)
end


