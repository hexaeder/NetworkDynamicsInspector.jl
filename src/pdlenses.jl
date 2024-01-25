function register_pd_lenses!()
    empty_lenses!()

    # get the estimated frequency
    register_vstatelens!(r"^_ω$") do sol, idx, state
        u_r_lens = vstatelens(sol, idx, :u_r)
        u_i_lens = vstatelens(sol, idx, :u_i)
        ndcop = deepcopy(_get_nd(sol))
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
    @register_vstate :_ω => (:u_r, :u_i)

    # get the total current at the end of the lines
    register_vstatelens!(r"^_i_([ri])$") do sol, idx, state
        m = match(r"_i_([ri])$", string(state))
        part = m[1] == "i" ? 2 : 1

        let part=part
            (gd, p, t) -> begin
                i = _total_current(get_dst_edges(gd, idx))
                i[part]
            end
        end
    end
    @register_vstate :_i_r :_i_i

    # power related states
    @register_vstatelens :_S => (:u_r + :u_i*im)*(:_i_r - :_i_i*im)
    @register_vstatelens :_P => real(:_S)
    @register_vstatelens :_Q => imag(:_S)

    @register_vstatelens :_Smeas => (:u_meas_r + :u_meas_i*im)*(:i_meas_r - :i_meas_i*im)
    @register_vstatelens :_Pmeas => real(:_Smeas)
    @register_vstatelens :_Qmeas => imag(:_Smeas)


    # a, b, c component
    register_vstatelens!(r"^.*_[abc]$") do sol, idx, state
        m = match(r"^(.*)_(.)$", string(state))
        state = m[1]
        phase = Dict("a"=>1,"b"=>2,"c"=>3)[m[2]]
        r_lens = vstatelens(sol, idx, Symbol(state*"_r"))
        i_lens = vstatelens(sol, idx, Symbol(state*"_i"))
        ω=2π*50

        (gd, p, t) -> begin
            r = r_lens(gd, p, t)
            i = i_lens(gd, p, t)
            abc = Tdqinv(ω*t) * [r, i]
            return abc[phase]
        end
    end
    # @register_vstate s"\1_a" => (r"^(.*)_r$", r"^(.*)_i$")
    # @register_vstate s"\1_b" => (r"^(.*)_r$", r"^(.*)_i$")
    # @register_vstate s"\1_c" => (r"^(.*)_r$", r"^(.*)_i$")

    @register_vstatelens r"^(.*)_mag$" => sqrt(s"\1_r"^2 + s"\1_i"^2)
    @register_vstatelens r"^(.*)_arg$" => atan(s"\1_i", s"\1_r")
    @register_estatelens r"^(.*)_mag$" => sqrt(s"\1_r"^2 + s"\1_i"^2)
    @register_estatelens r"^(.*)_arg$" => atan(s"\1_i", s"\1_r")

    # register transformation by pll
    @register_just_vstatelens r"^(.*)_d$" => cos(-:δ_pll)*s"\1_r" - sin(-:δ_pll)*s"\1_i"
    @register_just_vstatelens r"^(.*)_q$" => sin(-:δ_pll)*s"\1_r" + cos(-:δ_pll)*s"\1_i"

    ## now for the edges
    register_estatelens!(r"^(srcv|dstv)_(.*)$") do sol, idx, state
        m = match(r"^(srcv|dstv)_(.*)$", string(state))
        vidx = if m[1] == "srcv"
            _get_nd(sol).f.graph_structure.s_e[idx]
        elseif m[1] == "dstv"
            _get_nd(sol).f.graph_structure.d_e[idx]
        end
        vlens = vstatelens(sol, vidx, Symbol(m[2]))
        (gd, p, t) -> begin
            vlens(gd, p, t)
        end
    end

    @register_estatelens r"^_(src|dst)_Slcl" => (s"\1v_u_r" + s"\1v_u_i"*im)*(s"\1_i_r" - s"\1_i_i"*im)
    @register_estatelens r"^_(src|dst)_P" => real(s"_\1_S")
    @register_estatelens r"^_(src|dst)_Q" => imag(s"_\1_S")

    @register_estatelens :_P => max(:_src_P, :_dst_P)
    # @register_estate :_P :_src_P :_dst_P
    @register_estatelens :_Q => max(:_src_Q, :_dst_Q)
    # @register_estate :_Q :_src_Q :_dst_Q

    # dont show the _r and _i if ther are _arg and _mag
    # @register_vstatefilter r"^(.*)_r$" => (s"\1_arg", s"\1_mag")
    # @register_vstatefilter r"^(.*)_i$" => (s"\1_arg", s"\1_mag")
    # @register_estatefilter r"^(.*)_r$" => (s"\1_arg", s"\1_mag")
    # @register_estatefilter r"^(.*)_i$" => (s"\1_arg", s"\1_mag")
end

function _total_current(edges)
    # ND convention: all flows entering a node are positive
    current = 0.0im
    @inbounds for e in edges
        current += e[1] + e[2]*im
    end
    current
end
