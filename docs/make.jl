ENV["PLOTS_TEST"] = "true"
ENV["GKSwstype"] = "100"
using Documenter
using SplitClusterTest
using Literate
indir = joinpath(@__DIR__, "..", "examples")
outdir = joinpath(@__DIR__, "src", "examples")
files = ["demo-mirror.jl"]
for file in files
    Literate.markdown(joinpath(indir, file), outdir; credit = false)
end

makedocs(sitename="SplitClusterTest.jl",
        pages = [
            "Home" => "index.md",
            "Examples" => [
                "Mirror Statistics" => "examples/demo-mirror.md"
            ],
            "API" => "api.md"
        ]        
)

deploydocs(
    repo = "github.com/szcf-weiya/SplitClusterTest.jl"
)