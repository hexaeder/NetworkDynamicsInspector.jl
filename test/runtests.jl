using NetworkDynamicsInspector
using Test
using GraphMakie
using GLMakie
using NetworkDynamics

@testset "NetworkDynamicsInspector.jl" begin
    register_vstatelens!(r"_ω") do sol, p, idx, state
        u_r_lens = NetworkDynamics.vstatelens(sol, p, idx, :u_r)
        u_i_lens = NetworkDynamics.vstatelens(sol, p, idx, :u_i)
        ndcop = deepcopy(NetworkDynamics._get_nd(sol))
        (gd, t) -> begin
            u_r = u_r_lens(gd, t)
            u_i = u_i_lens(gd, t)
            dx = sol(t, Val{1})
            dgd = ndcop(dx, p(t), t, GetGD)
            u_dot_r = u_r_lens(dgd, t)
            u_dot_i = u_i_lens(dgd, t)
            return -(u_i*u_dot_r - u_r*u_dot_i)/(u_i^2 + u_r^2)
        end
    end

    sol = include("testpowergrid.jl");

    GLMakie.closeall()
    inspect_solution(sol)
end

@testset "Test FavSelect" begin
    fig = Figure()
    favsel1 = fig[1,1] = FavSelect(fig, [:a, :b], Observable([:x,:y,:z]))
    Label(fig[2,1], @lift(repr($(favsel1.selection))), tellwidth=false)
    favsel2 = fig[3,1] = FavSelect(fig, [:a, :b, :c], Observable([:x,:y,:z]); allowmulti=false)
    Label(fig[4,1], @lift(repr($(favsel2.selection))), tellwidth=false)

    fig = Figure()
    favsel1 = fig[1,1] = FavSelect(fig, [:a, :b], [:x,:y,:z])

    # symbol selector
    fig = Figure()
    selection = Observable(Symbol[:a,:b])
    tb = TBSelect(fig[1,1], selection; width=500)
    @test selection === tb.selection
    Label(fig[2,1], @lift(repr($(tb.selection))), tellwidth=false)

    selection[] = [:foo, :bar]

    # symbol selector for int
    fig = Figure()
    selection = Observable(Int[1,2,3])
    tb = TBSelect(fig[1,1], selection; width=500)
    @test selection === tb.selection
    Label(fig[2,1], @lift(repr($(tb.selection))), tellwidth=false)
end
