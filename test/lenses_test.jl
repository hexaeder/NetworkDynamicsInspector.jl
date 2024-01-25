using NetworkDynamics
using Graphs
using LinearAlgebra
using OrdinaryDiffEq
using Test
using NetworkDynamics: syms, observed_syms

@testset "Test lenses and state retrieval" begin
    # XXX: edge/vertex agnostic statelens is problematic, because complesx leneses require building of other lenses inside
    empty_lenses!()

    register_vstatelens!(r"^(.*)_(mag|arg)$") do sol, idx, state
        m = match(r"^(.*)_(mag|arg)$", string(state))
        if m[2] == "mag"
            transf = (r, i) -> sqrt(r^2 + i^2)
        elseif m[2] == "arg"
            transf = (r, i) -> atan(i, r)
        end
        rlens = NetworkDynamicsInspector.vstatelens(sol, idx, Symbol(m[1]*"_r"))
        ilens = NetworkDynamicsInspector.vstatelens(sol, idx, Symbol(m[1]*"_i"))
        let transf=transf
            (gd, p, t) -> begin
                r = rlens(gd, p, t)
                i = ilens(gd, p, t)
                transf(r,i)
            end
        end
    end

    register_vstatelens!(r"_ω") do sol, idx, state
        u_r_lens = NetworkDynamicsInspector.vstatelens(sol, idx, :u_r)
        u_i_lens = NetworkDynamicsInspector.vstatelens(sol, idx, :u_i)
        ndcop = deepcopy(NetworkDynamicsInspector._get_nd(sol))
        (gd, p, t) -> begin
            u_r = u_r_lens(gd, p, t)
            u_i = u_i_lens(gd, p, t)
            dx = sol(t, Val{1})
            dgd = ndcop(dx, p, t, GetGD)
            u_dot_r = u_r_lens(dgd, p, t)
            u_dot_i = u_i_lens(dgd, p, t)
            return -(u_i*u_dot_r - u_r*u_dot_i)/(u_i^2 + u_r^2)
        end
    end

    register_vstatelens!(r"_i_([ri])$") do sol, idx, state
        m = match(r"_i_([ri])$", string(state))
        part = m[1] == "i" ? imag : real

        let part=part
            (gd, p, t) -> begin
                i = -total_current(get_dst_edges(gd, idx))
                part(i)
            end
        end
    end

    @register_vstatelens r"^(.*)_mag_regex$" => sqrt(s"\1_r"^2 + s"\1_i"^2)
    @register_vstatelens :u_mag_simple => sqrt(:u_r^2 + :u_i^2)
    @register_vstatelens :_S => (:u_r + :u_i*im)*(:_i_r - :_i_i*im)
    @register_vstatelens :_P => real(:_S)
    @register_vstatelens :_Q => imag(:_S)

    include("testpowergrid.jl")
    nodep = [( 1, 1.1), # swing power and damping
             ( 1, 1.1),
             (-2, NaN)] # PQ load
    edgep = -5im        # line admittance
    p = (nodep, edgep)
    tspan = (0, 10)
    u0 = zeros(length(nd.syms))
    u0[2] = 0.1
    u0[idx_containing(nd, r"^u_r")] .= 1
    nd.syms .=> u0

    prob = ODEProblem(nd, u0, tspan, p)
    sol = solve(prob, Rodas4())
    # plot(sol)

    @test sol[:u_r_1] == vstatef(sol, p, 1, :u_r)(sol.t)
    @test sol[:u_i_1] == vstatef(sol, p, 1, :u_i)(sol.t)
    @test sol[:u_r_2] == vstatef(sol, p, 2, :u_r)(sol.t)
    @test sol[:u_i_2] == vstatef(sol, p, 2, :u_i)(sol.t)

    @test maximum(abs.(vstatef(sol,p,1,:ω)(sol.t) - vstatef(sol,p,1,:_ω)(sol.t))) < 1e-3
    @test maximum(abs.(vstatef(sol,p,2,:ω)(sol.t) - vstatef(sol,p,2,:_ω)(sol.t))) < 1e-3

    @test vstatef(sol, p, 1, :u_mag)(sol.t) ≈ sqrt.(sol[:u_r_1].^2 + sol[:u_i_1].^2)
    @test vstatef(sol, p, 2, :u_mag)(sol.t) ≈ sqrt.(sol[:u_r_2].^2 + sol[:u_i_2].^2)
    @test vstatef(sol, p, 3, :u_mag)(sol.t) ≈ sqrt.(sol[:u_r_3].^2 + sol[:u_i_3].^2)

    @test vstatef(sol, p, 1, :u_arg)(sol.t) ≈ atan.(sol[:u_i_1], sol[:u_r_1])
    @test vstatef(sol, p, 2, :u_arg)(sol.t) ≈ atan.(sol[:u_i_2], sol[:u_r_2])
    @test vstatef(sol, p, 3, :u_arg)(sol.t) ≈ atan.(sol[:u_i_3], sol[:u_r_3])

    @test vstatef(sol, p, 1, :u_mag_simple)(sol.t) == vstatef(sol, p, 1, :u_mag_regex)(sol.t) == vstatef(sol, p, 1, :u_mag)(sol.t)
    @test vstatef(sol, p, 2, :u_mag_simple)(sol.t) == vstatef(sol, p, 2, :u_mag_regex)(sol.t) == vstatef(sol, p, 2, :u_mag)(sol.t)
    @test vstatef(sol, p, 3, :u_mag_simple)(sol.t) == vstatef(sol, p, 3, :u_mag_regex)(sol.t) == vstatef(sol, p, 3, :u_mag)(sol.t)

    @test maximum(abs.(vstatef(sol, p, 1, :_P)(sol.t) + vstatef(sol, p, 2, :_P)(sol.t) + vstatef(sol, p, 3, :_P)(sol.t))) < 1e-10

    ## now for the edges
    register_estatelens!(r"^(srcv|dstv)_(.*)$") do sol, idx, state
        m = match(r"^(srcv|dstv)_(.*)$", string(state))
        vidx = if m[1] == "srcv"
            NetworkDynamicsInspector._get_nd(sol).f.graph_structure.s_e[idx]
        elseif m[1] == "dstv"
            NetworkDynamicsInspector._get_nd(sol).f.graph_structure.d_e[idx]
        end
        vlens = NetworkDynamicsInspector.vstatelens(sol, vidx, Symbol(m[2]))
        (gd, p, t) -> begin
            vlens(gd, p, t)
        end
    end

    estatef(sol, p, 1, :src_i_i)(sol.t)
    estatef(sol, p, 1, :src_i_r)(sol.t)
    @test estatef(sol, p, 1, :srcv_u_r)(sol.t) ≈ vstatef(sol, p, 1, :u_r)(sol.t)
    @test estatef(sol, p, 1, :dstv_u_r)(sol.t) ≈ vstatef(sol, p, 2, :u_r)(sol.t)
    @test estatef(sol, p, 3, :srcv_u_r)(sol.t) ≈ vstatef(sol, p, 2, :u_r)(sol.t)
    @test estatef(sol, p, 3, :dstv_u_r)(sol.t) ≈ vstatef(sol, p, 3, :u_r)(sol.t)

    @test estatef(sol, p, 3, :srcv_u_mag)(sol.t) ≈ vstatef(sol, p, 2, :u_mag)(sol.t)
    @test estatef(sol, p, 3, :dstv_u_mag)(sol.t) ≈ vstatef(sol, p, 3, :u_mag)(sol.t)

    @register_estatelens r"^_(src|dst)_S" => (s"\1v_u_r" + s"\1v_u_i"*im)*(s"\1_i_r" - s"\1_i_i"*im)
    @register_estatelens r"^_(src|dst)_P" => real(s"_\1_S")
    @register_estatelens r"^_(src|dst)_Q" => imag(s"_\1_S")

    estatef(sol, p, 1, :_src_P)(sol.t)
    estatef(sol, p, 1, :_src_Q)(sol.t)

    @register_estatelens :max_P => max(:_src_P, :_dst_P)

    @test map(eachrow(estatef(sol, p, 1, (:_src_P, :_dst_P, :max_P))(sol.t))) do (src, dst, m)
        m == max(src, dst)
    end |> all

    # test differen combinations of multiple arguments
    vstatef(sol, p, 1, :u_r)(1)
    vstatef(sol, p, 1:2, :u_r)(1)
    vstatef(sol, p, 1, (:u_r,:u_i))(1)
    vstatef(sol, p, 1:2, (:u_r,:u_i))(1)

    a = vstatef(sol, p, 1, :u_r)(0:0.1:1)
    b = vstatef(sol, p, 1:2, :u_i)(0:0.1:1)
    c = vstatef(sol, p, 1, (:u_r,:u_i))(0:0.1:1)
    d = vstatef(sol, p, 1:2, (:u_r,:u_i))(0:0.1:1)

    @test a == c[:, 1]
    @test a == d[:, 1, 1]
    @test b[:,2] == d[:,2,2]
    @test b[:,1] == c[:,2]
    @test b[:,1] == d[:,2,1]

    @test vstatef(sol, p, 1:3, :_ω)(sol.t) ≈ hcat(vstatef(sol, p, 1, :_ω)(sol.t), vstatef(sol, p, 2, :_ω)(sol.t), vstatef(sol, p, 3, :_ω)(sol.t))
    @test vstatef(sol, p, 1:3, :u_r)(sol.t) ≈ hcat(vstatef(sol, p, 1, :u_r)(sol.t), vstatef(sol, p, 2, :u_r)(sol.t), vstatef(sol, p, 3, :u_r)(sol.t))


    ## state registration
    empty!(NetworkDynamicsInspector.VSTATES)
    empty!(NetworkDynamicsInspector.VSTATEFILTER)
    listvstates(sol, 1)
    listvstates(sol, 3)
    @register_vstate :_ω => (:u_r, :u_i)
    @register_vstate s"\1_mag" => (r"(.+)_r", r"(.+)_i")
    @register_vstate s"\1_arg" => (r"(.+)_r", r"(.+)_i")
    listvstates(sol, 1)
    listvstates(sol, 3)
    @register_vstatefilter r"^(.+)_(r|i)$" => (s"\1_mag", s"\1_arg")
    listvstates(sol, 1)
    listvstates(sol, 3)

    empty!(NetworkDynamicsInspector.ESTATES)
    empty!(NetworkDynamicsInspector.ESTATEFILTER)
    listestates(sol, 1)
    @register_estate s"\1_mag" => (r"(.+)_r", r"(.+)_i")
    @register_estate s"\1_arg" => (r"(.+)_r", r"(.+)_i")
    @register_estate :_src_u_r
    @register_estate :_src_u_i
    @test listestates(sol, 1) == [:_src_u_arg,
                                  :_src_u_i,
                                  :_src_u_mag,
                                  :_src_u_r,
                                  :dst_i_arg,
                                  :dst_i_mag,
                                  :src_i_arg,
                                  :src_i_mag,
                                  :dst_i_r,
                                  :dst_i_i,
                                  :src_i_r,
                                  :src_i_i]
    @register_estatefilter r"^(.+)_(r|i)$" => (s"\1_mag", s"\1_arg")
    @test listestates(sol, 1) == [:_src_u_arg,
                                  :_src_u_mag,
                                  :dst_i_arg,
                                  :dst_i_mag,
                                  :src_i_arg,
                                  :src_i_mag]

    empty!(NetworkDynamicsInspector.VSTATES)
    empty!(NetworkDynamicsInspector.VSTATEFILTER)
    @test listvstates(sol, 1:3) == [:u_r, :u_i]
    @register_vstate :a :b
end

@testset "match conditions" begin
    using NetworkDynamicsInspector: _match_cond, _additional_states
    state = s"u_mag"
    conds = [r"u_r", r"u_i"]
    states = [:u_r, :u_i]
    @test _match_cond(state, conds, states) == ["u_mag"]
    @test _additional_states(states, Dict(state=>conds)) == Set([:u_mag])

    state = s"\1_arg"
    conds = [r"^(.+)_r", r"^(.+)_i"]
    states = [:u_r, :u_i, :i_r, :i_i]
    @test _match_cond(state, conds, states) == ["u_arg", "i_arg"]
    @test _additional_states(states, Dict(state=>conds)) == Set([:u_arg, :i_arg])

    state = s"\2_\1_foo"
    conds = [r"^(.+?)_(.+?)_r", r"^(.+?)_(.+?)_i"]
    states = [:u_r, :u_i, :i_r, :i_i, :u_x_r, :u_x_i, :u_y_r]
    @test _match_cond(state, conds, states) == ["x_u_foo"]
    @test _additional_states(states, Dict(state=>conds)) == Set([:x_u_foo])

    state = s"_\1_S"
    conds = [r"(src|dst)v_u_i", r"(src|dst)v_u_i", r"(src|dst)_i_r", r"(src|dst)_i_i",]
    states = [:srcv_u_r, :srcv_u_i, :dstv_u_r, :dstv_u_i, :src_i_r, :src_i_i, :dst_i_r, :dst_i_i]
    @test _match_cond(state, conds, states) == ["_src_S", "_dst_S"]
    @test _additional_states(states, Dict(state=>conds)) == Set([:_src_S, :_dst_S])

    state = s"_\1_S"
    conds = [r"(src|dst)v_u_i", r"(src|dst)v_u_i", r"(src|dst)_i_r", r"(src|dst)_i_i",]
    states = [:srcv_u_r, :srcv_u_i, :dstv_u_r, :dstv_u_i, :src_i_i, :dst_i_r ]
    @test _match_cond(state, conds, states) == String[]
    @test _additional_states(states, Dict(state=>conds)) == Set(Symbol[])
end

@testset "filter states" begin
    using NetworkDynamicsInspector: _filter_states!
    empty!(NetworkDynamicsInspector.ESTATES)
    empty!(NetworkDynamicsInspector.VSTATES)
    empty!(NetworkDynamicsInspector.ESTATEFILTER)
    empty!(NetworkDynamicsInspector.VSTATEFILTER)

    states = [:srcv_u_r, :srcv_u_i, :dstv_u_r, :dstv_u_i, :src_i_i, :src_i_r, :dst_i_r ]
    @register_vstate s"\1_mag" => (r"(.+)_r", r"(.+)_i")
    @register_vstate :src_i_arg => (:src_i_r, :src_i_i)
    add = _additional_states(states, NetworkDynamicsInspector.VSTATES)
    states = vcat(states,add...)
    # @register_vstatefilter :srcv_u_r
    @register_vstatefilter r"^(.+)_(u|i)$" => (s"\1_mag", s"\1_arg")
    @test setdiff(states, _filter_states!(copy(states), NetworkDynamicsInspector.VSTATEFILTER)) == [:src_i_i]
    @register_vstatefilter r".+_r$"
    @register_vstatefilter r".+_i$"
    @test setdiff(states, _filter_states!(copy(states), NetworkDynamicsInspector.VSTATEFILTER)) == [:srcv_u_r,
                                                                                           :srcv_u_i,
                                                                                           :dstv_u_r,
                                                                                           :dstv_u_i,
                                                                                           :src_i_i,
                                                                                           :src_i_r,
                                                                                           :dst_i_r]

end

@testset "test allocations and missing data" begin
    include("testpowergrid.jl")
    nodep = [( 1, 1.1), # swing power and damping
             ( 1, 1.1),
             (-2, NaN)] # PQ load
    edgep = -5im        # line admittance
    p = (nodep, edgep)
    tspan = (0, 10)
    u0 = zeros(length(nd.syms))
    u0[2] = 0.1
    u0[idx_containing(nd, r"^u_r")] .= 1
    nd.syms .=> u0

    prob = ODEProblem(nd, u0, tspan, p)
    sol = solve(prob, Rodas4())

    # single state
    lens = vstatef(sol, p, 1, :u_r); lens(1); lens(sol.t);
    @test 6 == @allocations lens(1)
    @test 5+length(sol.t) == @allocations lens(sol.t)
    buf = lens(sol.t);
    @allocations lens(sol.t)

    lens = vstatef(sol, p, 1, :u_x; failmode=:none); lens(1); lens(sol.t);
    @test 4 == @allocations lens(1)
    @test 4+length(sol.t) == @allocations lens(sol.t)

    # multiple indices
    lens = vstatef(sol, p, 1:3, :u_r); lens(1); lens(sol.t);
    @test 5 == @allocations lens(1)
    @test 5+length(sol.t) == @allocations lens(sol.t)
    buf = lens(0)

    lens = vstatef(sol, p, 1:3, :u_x; failmode=:none); lens(1); lens(sol.t);
    @test 4 == @allocations lens(1)
    @test 4+length(sol.t) == @allocations lens(sol.t)

    # multiple variables
    lens = vstatef(sol, p, 1, [:u_r, :u_i]); lens(1); lens(sol.t);
    @test 5 == @allocations lens(1)
    @test 5+length(sol.t) == @allocations lens(sol.t)

    # multiple variables
    lens = vstatef(sol, p, 1:3, [:u_r, :u_i]); lens(1); lens(sol.t);
    @test 2 == @allocations lens(1)
    @test 2+length(sol.t) == @allocations lens(sol.t)

    # missings
    lens = vstatef(sol, p, 1:2, [:u_x, :u_i]; failmode=:none); lens(1); lens(sol.t);
    @test 6 == @allocations lens(1)
    @test 295 == @allocations lens(sol.t)

    # composite variables
    lens = vstatef(sol, p, 1, :u_mag); lens(1); lens(sol.t);
    @test 6 == @allocations lens(1)
    @test 5+length(sol.t) == @allocations lens(sol.t)
end


@testset "Test autoamtic creation of states for lenses" begin
    using NetworkDynamicsInspector: _statedef_from_lensdef
    using NetworkDynamicsInspector.MacroTools: striplines

    @test _statedef_from_lensdef(:(:A => (:B - :C)^2)) == :(:A => (:B, :C))

    ex1 = :(r"^(.*)_mag_regex$" => sqrt(s"\1_r"^2 + s"\1_i"^2))
    ex2 = :(s"\1_mag_regex" => (r"^(.*)_r$", r"^(.*)_i$"))
    @test striplines(_statedef_from_lensdef(ex1)) == striplines(ex2)

    ex1 = :(r"^_(src|dst)_S$" => (s"\1v_u_r" + s"\1v_u_i"*im)*(s"\1_i_r" - s"\1_i_i"*im))
    ex2 = :(s"_\1_S" => (r"^(src|dst)v_u_r$", r"^(src|dst)v_u_i$", r"^(src|dst)_i_r$", r"^(src|dst)_i_i$"))
    @test striplines(_statedef_from_lensdef(ex1)) == striplines(ex2)

    ex1 = :(r"^_(.*)_S$" => (s"\1_i_i" + :symbol))
    ex2 = :(s"_\1_S" => (:symbol, r"^(.*)_i_i$"))
    @test striplines(_statedef_from_lensdef(ex1)) == striplines(ex2)

    ex = :(r"^(.+?)_(.+?)_r$" => s"\2_\1_foo" + s"\1_\2_bar")
    @test_throws AssertionError _statedef_from_lensdef(ex)

    ex = :(r"^(.+?)_(.+?)_r" => s"\2_\1_foo" + s"\1_\2_\4bar")
    @test_throws AssertionError _statedef_from_lensdef(ex)
end
