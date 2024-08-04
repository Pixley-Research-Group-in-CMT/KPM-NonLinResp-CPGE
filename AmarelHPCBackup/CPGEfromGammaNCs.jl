using JLD2
using Distributed
#using ProgressMeter
using KPM
using MAT # write .mat file for matlab
using Dates
using CUDA

@show CUDA.has_cuda()
@show KPM.whichcore()
@show Threads.nthreads() # check threads

function plot_cpge(NC, w1, w2, Gammaxyz, Gammayxz; kernel=KPM.JacksonKernel, y_mult=1.0)
    #lbl = "$NC, $(round(w1+w2; digits=3)), $(round(abs(w1-w2)/2; digits=3))"
    h = 0.001;
    x_all = collect(-0.5:h:0.5)
    y_all = complex(x_all)
    mu_xyz_partial = KPM.maybe_to_device(Gammaxyz[1:NC, 1:NC, 1:NC])
    mu_yxz_partial = KPM.maybe_to_device(Gammayxz[1:NC, 1:NC, 1:NC])

    #updatetime = 1
    #@showprogress updatetime "Computing..." 
    for (i, x) in enumerate(x_all)
        y_all[i] = KPM.d_cpge(mu_xyz_partial, NC, w1, w2, x; δ=1e-9, kernel=kernel) + KPM.d_cpge(mu_yxz_partial, NC, w2, w1, x; δ=1e-9, kernel=kernel)
    end
    y_all = y_all / w1 / w2
    y_integral = cumsum(y_all) * h;
    mu_xyz_partial = nothing
    mu_yxz_partial = nothing
    y_all = nothing
    GC.gc()
    return (x_all, y_integral)
end

#disind = parse(Int64,ARGS[1])
pt = "/cache/home/aw666/data/CPGEdata/"
#fn = "/Users/angkunwu/CPGE/GammamnpL100NC256W0NR50.jld2"
#fn = "/cache/home/aw666/data/CPGEdata/GammamnpL100gam08D6NC512W01NR1Sams50.jld2" #in amarel
#fn = "/mnt/home/awu10/CPGE/GammamnpL100NC256W0NR50.jld2" #in rusty
#fn = "/mnt/home/awu10/ceph/CPGEdata/GammamnpL100gam08D6NC512W05maxNR1Sams100.jld2" #in rusty
#prefix = "/mnt/home/awu10/ceph/CPGEdata/GammamnpL100gam04D6NC512NR1S1Ws"
#suffix = ".jld2"
W2 = "03"
#fn = pt*"GammamnpL100gam08D6NC512PotW2"*W2*"NR1Sams50.jld2"
fn = pt*"GammamnpL18gam03D6NC512PotW05NR100TBCs100.jld2"
#fn = prefix*"$disind"*suffix
Ls = load(fn,"Ls")
NC = load(fn,"NC")
D = load(fn,"D")
Gammaxyz = load(fn,"Gammaxyz")
Gammayxz = load(fn,"Gammayxz")
#Ws = load(fn,"Ws")

w = 0.3/D #0.8/D
W = 0.0
NCs = collect(32:10:512) #collect(64:20:1024)
CPGEs = zeros(Complex{Float64},1001,49)
Efs = zeros(1001,49)
#@show typeof(Efs)
#@show typeof(CPGEs)
starttime = now()

#Threads.@threads 
for k = 1:49
    nNC = NCs[k]
    w1 = w+W
    w2 = -w
    #Efstemp, CPGEstemp = plot_cpge(NC, w1, w2,Gammaxyz,Gammayxz)
    #@show typeof(Efstemp) 
    #@show typeof(CPGEstemp)
    Efs[:,k], CPGEs[:,k] = plot_cpge(nNC, w1, w2,Gammaxyz,Gammayxz)
    #sleep(k)
    endtime = now()
    @show k, endtime-starttime
    #Efs[:,k] = Efstemp
    #CPGEs[:,k] = CPGEstemp
    GC.gc()
end
Efs = Efs[:,1]

#prefix = "/mnt/home/awu10/ceph/CPGEdata/cpgesw04L100D6NCsWs"
#suffix = ".mat"
#file = matopen(pt*"cpgesw08L100PotW2"*W2*"D6NCsMean50.mat", "w")
file = matopen(pt*"cpgesw03L18PotW05D6NC512sNR100TBCs100.mat", "w")
#file = matopen(prefix*"$disind"*suffix, "w")
write(file, "Efs", Efs)
write(file, "D", D)
write(file, "Ls", Ls)
write(file, "NCs", NCs)
write(file, "w", w)
write(file, "W", W)
#write(file,"Ws",Ws)
#write(file,"disind",disind)
write(file, "CPGEs", CPGEs)
close(file)

