using JLD2
using Distributed
#using ProgressMeter
using KPM
using MAT # write .mat file for matlab
using Dates
using CUDA
using HDF5

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

pt = "/home/angkunwu/"
fn = pt*"GammamnpL100gam08D6NC512W203NR1Sams50.jld2"
Ls = load(fn,"Ls")
NC = load(fn,"NC")
D = load(fn,"D")
Gammaxyz = load(fn,"Gammaxyz")
Gammayxz = load(fn,"Gammayxz")
#Ws = load(fn,"Ws")

ws = collect(0.1:0.1:0.8)
W = 0.0
CPGEs = zeros(Complex{Float64},1001,8)
Efs = zeros(1001,8)
starttime = now()
for k = 1:8
    w1 = ws[k]/D+W
    w2 = -ws[k]/D
    Efs[:,k], CPGEs[:,k] = plot_cpge(NC, w1, w2,Gammaxyz,Gammayxz)
    endtime = now()
    @show k, endtime-starttime
end
Efs = Efs[:,1]

file = matopen(pt*"cpgeswsL100W203D6NC512Mean50.mat", "w")
write(file, "Efs", Efs)
write(file, "D", D)
write(file, "Ls", Ls)
write(file, "NC", NC)
write(file, "ws", ws)
write(file, "W", W)
write(file, "CPGEs", CPGEs)
close(file)




