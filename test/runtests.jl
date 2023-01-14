using NetworkDynamicsInspector
using Test
using GraphMakie
using GLMakie
using NetworkDynamics
using OrderedCollections

@testset "NetworkDynamicsInspector.jl" begin
    empty_statestuff!()
    register_vstatelens!(r"^_ω$") do sol, p, idx, state
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
    @register_vstate :_ω => (:u_r, :u_i)

    register_vstatelens!(r"_i_([ri])$") do sol, p, idx, state
        m = match(r"_i_([ri])$", string(state))
        part = m[1] == "i" ? imag : real

        let part=part
            (gd, t) -> begin
                i = -total_current(get_dst_edges(gd, idx))
                part(i)
            end
        end
    end
    @register_vstate :_i_r :_i_i

    @register_vstatelens :_S => (:u_r + :u_i*im)*(:_i_r - :_i_i*im)
    @register_vstatelens :_P => real(:_S)
    @register_vstatelens :_Q => imag(:_S)

    register_estatelens!(r"^(srcv|dstv)_(.*)$") do sol, p, idx, state
        m = match(r"^(srcv|dstv)_(.*)$", string(state))
        vidx = if m[1] == "srcv"
            NetworkDynamics._get_nd(sol).f.graph_structure.s_e[idx]
        elseif m[1] == "dstv"
            NetworkDynamics._get_nd(sol).f.graph_structure.d_e[idx]
        end
        vlens = NetworkDynamics.vstatelens(sol, p, vidx, Symbol(m[2]))
        (gd, t) -> begin
            vlens(gd, t)
        end
    end
    @register_estate :srcv_u_i :srcv_u_r :dstv_u_i :dstv_u_r

    @register_estatelens r"^_(src|dst)_S$" => (s"\1v_u_r" + s"\1v_u_i"*im)*(s"\1_i_r" - s"\1_i_i"*im)
    @register_estatelens r"^_(src|dst)_P$" => real(s"_\1_S")
    @register_estatelens r"^_(src|dst)_Q$" => imag(s"_\1_S")

    @register_estatelens :_P => :_src_P

    @register_vstatelens r"^(.*)_mag$" => sqrt(s"\1_r"^2 + s"\1_i"^2)
    @register_vstatelens r"^(.*)_arg$" => atan(s"\1_i", s"\1_r")
    @register_estatelens r"^(.*)_mag$" => sqrt(s"\1_r"^2 + s"\1_i"^2)
    @register_estatelens r"^(.*)_arg$" => atan(s"\1_i", s"\1_r")

    @register_vstatefilter r"^(.*)_r" => (s"\1_arg", s"\1_mag")
    @register_vstatefilter r"^(.*)_i" => (s"\1_arg", s"\1_mag")
    @register_estatefilter r"^(.*)_r" => (s"\1_arg", s"\1_mag")
    @register_estatefilter r"^(.*)_i" => (s"\1_arg", s"\1_mag")

    NetworkDynamics.VSTATES
    NetworkDynamics.ESTATES

    sol = include("testpowergrid.jl");

    GLMakie.closeall()
    inspect_solution(sol)


    convert(Vector{Float64}, vec(a))
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

    fig = Figure()
    selection = Observable(OrderedSet([1,2,3]))
    tb = TBSelect(fig[1,1], selection; width=500)
    @test selection === tb.selection
    Label(fig[2,1], @lift(repr($(tb.selection))), tellwidth=false)
end

@testset "treesym style" begin
    d = Dict(:a => "fooo",
             2 => "bar",
             "string" => ('b', 'a', 'z'))
    using NetworkDynamicsInspector: treestyle_string
    treestyle_string(d)
end
