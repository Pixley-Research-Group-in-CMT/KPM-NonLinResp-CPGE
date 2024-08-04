using CUDA, Test

@show CUDA.has_cuda()

using KPM
@show KPM.whichcore()

@show CUDA,devices()
@show CUDA.ones(2,2)



