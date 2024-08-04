using CUDA, Test

@show CUDA.has_cuda()

using KPM
KPM.whichcore() # this step will print a false notice saying GPU unavailable

#CUDA.devices()
@show CUDA.ones(2,2)

