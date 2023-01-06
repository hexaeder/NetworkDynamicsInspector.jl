using NetworkDynamicsInspector
using Documenter

DocMeta.setdocmeta!(NetworkDynamicsInspector, :DocTestSetup, :(using NetworkDynamicsInspector); recursive=true)

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
    ],
)

deploydocs(;
    repo="github.com/hexaeder/NetworkDynamicsInspector.jl",
    devbranch="main",
)
