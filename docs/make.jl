using NetworkDynamicsInspector
using Documenter
using Literate
using OrdinaryDiffEq
using DiffEqCallbacks
using CairoMakie
using Graphs
using GraphMakie

DocMeta.setdocmeta!(NetworkDynamicsInspector, :DocTestSetup, :(using NetworkDynamicsInspector); recursive=true)

# generate examples
example_dir = joinpath(@__DIR__, "..", "examples")
outdir = joinpath(@__DIR__, "src", "generated")
isdir(outdir) && rm(outdir, recursive=true)
mkpath(outdir)

for example in filter(contains(r".jl$"), readdir(example_dir, join=true))
    Literate.markdown(example, outdir)
end

makedocs(;
    modules=[NetworkDynamicsInspector],
    authors="Hans Würfel <git@wuerfel.io> and contributors",
    repo="https://github.com/hexaeder/NetworkDynamicsInspector.jl/blob/{commit}{path}#{line}",
    sitename="NetworkDynamicsInspector.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://hexaeder.github.io/NetworkDynamicsInspector.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => ["Feature Walkthrough" => "generated/walkthrough.md"]
    ],
)

deploydocs(;
    repo="github.com/hexaeder/NetworkDynamicsInspector.jl",
    devbranch="main",
)
