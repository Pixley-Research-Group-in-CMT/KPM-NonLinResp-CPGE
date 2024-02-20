using JLD2
using KPM

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

fn = "/Users/angkunwu/CPGE/GammamnpL100NC256W0NR50.jld2"
#fn = "/cache/home/aw666/CPGE/GammamnpL100NC256W0NR50.jld2" #in amarel
L = load(fn,"L")
NC = load(fn,"NC")
D = load(fn,"D")
Gammaxyz = load(fn,"Gammaxyz")
Gammayxz = load(fn,"Gammayxz")
@show D
NC = 128
Ω = 3*10^(-2)
@show w1 = 0.8/D+Ω
@show w2 = -0.8/D
Ef = -1.143/D
λ = 0*10^(-8)
@time Efs, cpge = plot_cpge(NC, w1, w2,Gammaxyz,Gammayxz)
#=
mu_xyz_partial = KPM.maybe_to_device(Gammaxyz[1:NC, 1:NC, 1:NC])
mu_yxz_partial = KPM.maybe_to_device(Gammayxz[1:NC, 1:NC, 1:NC])
@time v1 = KPM.d_cpge(mu_xyz_partial, NC, w1, w2, Ef;λ=λ)/(w1*w2)
@time v2 = KPM.d_cpge(mu_yxz_partial, NC, w2, w1, Ef;λ=λ)/(w1*w2)
=#
