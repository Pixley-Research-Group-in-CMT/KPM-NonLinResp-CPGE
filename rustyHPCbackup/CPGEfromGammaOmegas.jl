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
    return (x_all, y_integral)
end


#fn = "/Users/angkunwu/CPGE/GammamnpL100NC256W0NR50.jld2"
#fn = "/cache/home/aw666/CPGE/GammamnpL100NC256W0NR50.jld2" #in amarel
#fn = "/mnt/home/awu10/CPGE/GammamnpL100NC256W0NR50.jld2" #in rusty
fn = "/mnt/home/awu10/ceph/CPGEdata/GammamnpL100NC210W0NRs5.jld2" #in rusty
Ls = load(fn,"Ls")
NC = load(fn,"NC")
D = load(fn,"D")
Gammaxyz = load(fn,"Gammaxyz")
Gammayxz = load(fn,"Gammayxz")


NC = 1024
w = 0.8/D
Ωs = collect(LinRange(-4.6,-2.3,20))
Ws = exp.(Ωs)./D
#Ws = collect(LinRange(0.03, 0.063, 20))./D
CPGEs = zeros(Complex{Float64},1001,20)
Efs = zeros(1001,20)
#@show typeof(Efs)
#@show typeof(CPGEs)
starttime = now()

#Threads.@threads 
for k = 19:20
    W = Ws[k]
    w1 = w+W
    w2 = -w
    #Efstemp, CPGEstemp = plot_cpge(NC, w1, w2,Gammaxyz,Gammayxz)
    #@show typeof(Efstemp) 
    #@show typeof(CPGEstemp)
    Efs[:,k], CPGEs[:,k] = plot_cpge(NC, w1, w2,Gammaxyz,Gammayxz)
    #sleep(k)
    endtime = now()
    @show k, endtime-starttime
    #Efs[:,k] = Efstemp
    #CPGEs[:,k] = CPGEstemp
end
Efs = Efs[:,1]


file = matopen("cpgesw08L100NC210Ws1920.mat", "w")
write(file, "Efs", Efs)
write(file, "D", D)
write(file, "Ls", Ls)
write(file, "NC", NC)
write(file, "w", w)
write(file, "Ws", Ws)
write(file, "CPGEs", CPGEs)
close(file)





