#! julia

using Pkg
using REPL.TerminalMenus
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=dirname(@__DIR__))) # adds the package this script is called from
Pkg.instantiate()
Pkg.update()

port = isempty(ARGS) ? 8000 : parse(Int, ARGS[1])
@assert 8000 ≤ port ≤ 9000 "port has to be in range 8000..9000!"

using LiveServer
@async serve(dir=joinpath(@__DIR__, "build"))

menu = RadioMenu(["Run again!", "Quit!"])
while true
    revise()
    @info "Start building docs..."
    try
        include("make.jl")
    catch e
        @info "make.jl error" e
    end

    println("\nDocs are served at http://localhost:$port")

    if request("What now?", menu) != 1
        break
    end
end
