using NetworkDynamicsInspector
using Test

@testset "NetworkDynamicsInspector.jl" begin
    sol = include("testpowergrid.jl");

    pr = PRecord(p)

    sol.prob.f.f.graph
    gparguments(sol, pr, sol.prob.f.f.graph)
end
