using NetworkDynamicsInspector
using Test
using GraphMakie
using GLMakie
using NetworkDynamics
using OrderedCollections

include("lenses_test.jl")

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

    fig = Figure()
    selection = Observable(OrderedSet([1,2,3]))
    tb = TBSelect(fig[1,1], selection; width=500)
    @test selection === tb.selection
    Label(fig[2,1], @lift(repr($(tb.selection))), tellwidth=false)

    fig = Figure()
    options = Observable([:x,:y,:z])
    favsel1 = fig[1,1] = FavSelect(fig, Symbol[:a, :b], options)

    options = Observable(Symbol[])
    favsel1 = fig[1,1] = FavSelect(fig, Symbol[:a, :b], options)
    options[] = [:x,:y,:z]
end

@testset "treesym style" begin
    d = Dict(:a => "fooo",
             2 => "bar",
             "string" => ('b', 'a', 'z'))
    using NetworkDynamicsInspector: treestyle_string
    treestyle_string(d)
end

@testset "PRecord" begin
    a = [1,2,3]
    b = [3,4,5]
    p = (a, b)

    pr = PRecord(p)
    p[1][1] = 4
    @test pr(0) == ([1,2,3], b)
    record!(pr, 1, p)
    @test pr(1.1) == p
    @test pr(0.9) == ([1,2,3], b)
end
