ENV["PLOTS_TEST"] = "true"
ENV["GKSwstype"] = "100"
using Documenter
using SplitClusterTest
using Literate
indir = joinpath(@__DIR__, "..", "examples")
outdir = joinpath(@__DIR__, "src", "examples")
files = ["two-gaussians.jl", "two-poissons.jl", "cont-poissons.jl"]
for file in files
    Literate.markdown(joinpath(indir, file), outdir; credit = false)
end

makedocs(sitename="SplitClusterTest.jl",
        pages = [
            "Home" => "index.md",
            "Examples" => [
                "Two Gaussians" => "examples/two-gaussians.md",
                "Two Poissons" => "examples/two-poissons.md",
                "Continuous Poissons (Linear Pseduotime)" => "examples/cont-poissons.md"
            ],
            "API" => "api.md"
        ]        
)

deploydocs(
    repo = "github.com/szcf-weiya/SplitClusterTest.jl"
)