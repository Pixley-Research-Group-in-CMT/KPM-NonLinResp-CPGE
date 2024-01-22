push!(LOAD_PATH,"../src/")
using CPGE

using Documenter

makedocs(
         sitename = "CPGE.jl",
         modules  = [CPGE],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="https://github.com/Pixley-Research-Group-in-CMT/KPM-NonLinResp-CPGE.git",
)

