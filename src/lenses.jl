#=
Wishilist:
- get single value
- get timeseires
- both for edges and nodes
- get full vectors
- discover available symbols
- define "virtual" patterns

@register_virtual_state?
=#
# TODO: should lenses be a struct? EStateLens <: StateLens and so on?
using MacroTools

function get_vertexf(nd::AbstractDiffEqFunction, idx)
    ndobj = nd.f
    _group = findfirst(group -> idx ∈ group, ndobj.unique_v_indices)
    ndobj.unique_vertices![_group]
end
get_vertexf(sol, idx) = get_vertexf(_get_nd(sol), idx)

function get_edgef(nd::AbstractDiffEqFunction, idx)
    ndobj = nd.f
    _group = findfirst(group -> idx ∈ group, ndobj.unique_e_indices)
    ndobj.unique_edges![_group]
end
get_edgef(sol, idx) = get_edgef(_get_nd(sol), idx)

_get_nd(sol) = sol.prob.f

export vstatef, estatef

estatef(sol, p, idx, states; failmode=:error) = _statef(sol, p, idx, states, estatelens; failmode)

vstatef(sol, p, idx, states; failmode=:error) = _statef(sol, p, idx, states, vstatelens; failmode)

function _statef(sol::ODESolution, p, idxs, states, statelens; failmode=:error)
    pr = PRecord(p)
    nd = _get_nd(sol)

    drops = Int[]
    if idxs isa Integer
        push!(drops, 3)
    end
    if states isa Union{Symbol, String}
        push!(drops, 2)
    end
    drops = Tuple(drops)

    states = states isa Union{Symbol,String} ? (Symbol(states),) : Symbol.(states)

    Ni = length(idxs)
    Ns = states isa Union{Symbol,String} ? 1 : length(states)

    missinglens(gd,p,t) = missing
    hasmissings = false

    lenses = Vector{Any}(undef, length(states))
    for (si, state) in enumerate(states)
        idxslenses = Vector{Any}(undef, length(idxs))
        for (ii, idx) in enumerate(idxs)
            idxslenses[ii] = try
                statelens(sol, idx, state)
            catch e
                if failmode===:error || !(e isa ArgumentError)
                    rethrow(e)
                else
                    hasmissings = true
                    missinglens
                end
            end
        end
        lenses[si] = Tuple(idxslenses)
    end

    if failmode === :warn && hasmissings
        @warn "Some lenses for idx ∈ $idxs and sym ∈ $states could not be resolved! Timeseries will produce Missings."
    end

    # bufT = failmode!=:error ? Union{eltype(sol),Missing} : eltype(sol)
    bufT = hasmissings ? Union{eltype(sol),Missing} : eltype(sol)
    bufTwrap = Ref{bufT}()

    let drops=Tuple(drops), lenses=Tuple(lenses)
        # (ts, buf_in=Array{eltype(bufTwrap), 3}(undef, length(ts), Ns, Ni)) -> begin
        (ts) -> begin
            buf = Array{eltype(bufTwrap), 3}(undef, length(ts), Ns, Ni)
            last_t = typemin(eltype(ts))

            for (ti, t) in enumerate(ts)
                localp = t==last_t ? pr(t;direction=:right) : pr(t;direction=:left)

                gd = nd(sol(t), localp, t, GetGD)
                # potentially improve performance by unrolling
                # for i in 1:Ni
                #     for s in 1:Ns
                #         len = lenses[s][i]
                #         val::eltype(bufTwrap) = len(gd,localp,t)
                #         buf[ti, s, i] = val
                #         # buf[ti, s, i] = lenses[s, i](gd,localp,t)
                #     end
                # end
                for (s, idxslenses) in enumerate(lenses)
                    for (i, len) in enumerate(idxslenses)
                        val = len(gd,localp,t)
                        buf[ti, s, i] = val
                    end
                end

                last_t = t
            end

            buf_out = dropdims(buf; dims=drops)
            if buf_out isa Vector && ts isa Number
                return only(buf_out)
            else
                return buf_out
            end
        end
    end
end

function vstatelens(sol, idx, state::Symbol)
    # @info "Crreate lens for $i $state"
    nd = _get_nd(sol)
    cf = get_vertexf(nd,idx)

    if state ∈ syms(cf)
        stateidx::Int = findfirst(s->s==state, syms(cf))
        return (gd, p, t) -> begin
            get_vertex(gd, idx)[stateidx]
        end
    elseif state ∈ observed_syms(cf)
        obsidx::Int = findfirst(s->s==state, observed_syms(cf))
        return (gd,p,t) -> begin
            statevec = get_vertex(gd, idx)
            input = get_dst_edges(gd, idx)
            observed_f(cf)(statevec, input, NetworkDynamics.p_v_idx(p, idx), t)[obsidx::Int]
        end
    elseif (f = _get_from_regexdict(VSTATELENSES, state)) !== nothing
        return f(sol, idx, state)
    else
        throw(ArgumentError("Could not match rule to symbol $state"))
    end
end
function estatelens(sol, idx, state::Symbol)
    nd = _get_nd(sol)
    cf = get_edgef(nd, idx)

    if state ∈ syms(cf)
        stateidx = findfirst(s->s==state, syms(cf))
        return (gd, p, t) -> begin
            get_edge(gd, idx)[stateidx]
        end
    # elseif state ∈ observed_syms(cf)
    #     stateidx = findfirst(s->s==state, observed_syms(cf))
    #     return (gd::GraphData, estate) -> begin
    #         input = get_dst_edges(gd, idx)
    #         observed_f(cf)(statevec, input, NetworkDynamics.p_v_idx(p, idx), t)[stateidx]
    #     end
    elseif (f = _get_from_regexdict(ESTATELENSES, state)) !== nothing
        return f(sol, idx, state)
    else
        throw(ArgumentError("Could not match rule to symbol $state"))
    end
end

function _get_from_regexdict(dict, state)
    regexes = collect(keys(dict))
    ridxs = findall(occursin(string(state)), regexes)
    if isempty(ridxs)
        return nothing
    end
    ridx = first(ridxs)
    if length(ridxs) > 1
        @warn "Abiguity detected: $state matches multiple rules $(regexes[ridxs]). First match $(regexes[ridx]) chosen."
    end
    return dict[regexes[ridx]]
end

const VSTATELENSES = Dict{Regex, Function}()
const ESTATELENSES = Dict{Regex, Function}()
export register_statelens!, register_vstatelens!, register_estatelens!

register_vstatelens!(f, key) = _register_statelens!(VSTATELENSES, key, f)
register_estatelens!(f, key) = _register_statelens!(ESTATELENSES, key, f)

_register_statelens!(dict, sym::Symbol, f) = _register_statelens!(dict, string(sym), f)
_register_statelens!(dict, str::String, f) = _register_statelens!(dict, Regex("^"*str*"\$"), f)
function _register_statelens!(dict, rex::Regex, f)
    dict[rex] = f
    dict
end

using MacroTools: postwalk, @capture

export @register_vstatelens, @register_estatelens

macro register_just_vstatelens(ex)
    key, fun = _generate_statelens_function(ex, :vertex)
    quote
        register_vstatelens!($fun, $key)
    end
end

macro register_vstatelens(ex)
    key, fun = _generate_statelens_function(ex, :vertex)
    statedef = _statedef_from_lensdef(ex)
    quote
        @register_vstate $statedef
        register_vstatelens!($fun, $key)
    end
end
macro register_estatelens(ex)
    key, fun = _generate_statelens_function(ex, :edge)
    statedef = _statedef_from_lensdef(ex)
    quote
        @register_estate $statedef
        register_estatelens!($fun, $key)
    end
end

function _generate_statelens_function(ex, type)
    ex = deepcopy(ex)
    @capture(ex, A_ => B_)

    matchhead = if isexpr(A) && A.head === :macrocall && A.args[1] === Symbol("@r_str")
        :(m = match($(esc(A)), string(state)))
    else
        nothing
    end

    symbols = []
    postwalk(B) do x
        if x isa QuoteNode
            push!(symbols, x)
        elseif isexpr(x) && x.head === :macrocall && x.args[1] === Symbol("@s_str")
            isnothing(matchhead) && error("Can't use substitution strings for lhs $A (regex needed).")
            push!(symbols, x)
        end
        return x
    end
    unique!(symbols)

    lenssyms = map(s -> gensym(string(s)*"_lens"), symbols)
    lenses = map(symbols, lenssyms) do s, lenssym
        if isexpr(s) && s.head === :macrocall
            s = :(Symbol(replace(string(state), $A=>$s)))
        end
        if type===:vertex
            :($lenssym = vstatelens(sol, idx, $s))
        else
            :($lenssym = estatelens(sol, idx, $s))
        end
    end
    lensconstruct = Expr(:block, lenses...)

    valuesyms = map(s -> gensym(string(s)), symbols)
    values = map(lenssyms, valuesyms) do lens, val
        :($val = $lens(gd, p, t))
    end
    lenseval = Expr(:block, values...)

    body = postwalk(B) do x
        idxs = findfirst(isequal(x), symbols)
        if isnothing(idxs)
            return x
        else
            idx = only(idxs)
            return valuesyms[idxs]
        end
    end

    fun = quote
        (sol, idx, state) -> begin
            $matchhead
            $lensconstruct
            (gd, p, t) -> begin
                $lenseval
                $body
            end
        end
    end
    return A, fun
end

function _statedef_from_lensdef(ex)
    ex = deepcopy(ex)
    @capture(ex, A_ => B_)
    regexrule = isexpr(A) && A.head===:macrocall && A.args[1]===Symbol("@r_str")

    symbols = QuoteNode[]
    s_strs = Expr[]
    postwalk(B) do x
        if x isa QuoteNode && x.value isa Symbol
            push!(symbols, x)
        elseif isexpr(x) && x.head===:macrocall && x.args[1]===Symbol("@s_str")
            @assert regexrule "Lhs needs to be regex to use $x"
            push!(s_strs, x)
        end
        return x
    end

    if regexrule
        # collect all catchgroups
        group_idx = Set{Int}()
        for sstr in s_strs
            str = sstr.args[end] # the raw string
            matches = eachmatch(r"\\([1-9]+)", str)
            lasti = 0
            for m in matches
                i = parse(Int, m[1])
                @assert i == lasti+1 "Replace rules musst appear in order!"
                lasti = i
                push!(group_idx, i)
            end
        end
        group_idx = sort!(collect(group_idx))
        @assert group_idx == group_idx[begin]:group_idx[end]

        grx = r"\(.+?\)"
        groups = [m.match for m in eachmatch(grx, A.args[end])]
        @assert length(groups) == group_idx[end] "Found $(group_idx[end]) groups in substitution strings but $(length(groups)) in lhs."

        newstr = A.args[end]
        for i in group_idx
            newstr = replace(newstr, grx=>"\\$i"; count=1)
        end
        newstr = replace(newstr, r"^\^"=>s"") # no ^ at begin
        newstr = replace(newstr, r"\$$"=>s"") # no $ at end
        A.args[1] = Symbol("@s_str")
        A.args[end] = newstr

        for s in s_strs
            s.args[1] = Symbol("@r_str")
            newstr = s.args[end]
            for i in group_idx
                newstr = replace(newstr, "\\$i" => groups[i])
            end
            newstr = "^"*newstr*"\$"
            s.args[end] = newstr
        end
    end


    rhs = Expr(:tuple, symbols..., s_strs...)
    return :($A => $rhs)
end


const VSTATES = Dict{SubstitutionString, Union{Nothing, Vector{Regex}}}()
const ESTATES = Dict{SubstitutionString, Union{Nothing, Vector{Regex}}}()
const VSTATEFILTER = Dict{Regex, Union{Nothing, Vector{SubstitutionString}}}()
const ESTATEFILTER = Dict{Regex, Union{Nothing, Vector{SubstitutionString}}}()

register_vstate!(s, cond=nothing) = setindex!(VSTATES, _prepare_state_cond(cond, s, SubstitutionString, Regex)...)
register_estate!(s, cond=nothing) = setindex!(ESTATES, _prepare_state_cond(cond, s, SubstitutionString, Regex)...)
register_vstatefilter!(s, cond=nothing) = setindex!(VSTATEFILTER, _prepare_state_cond(cond, s, Regex, SubstitutionString)...)
register_estatefilter!(s, cond=nothing) = setindex!(ESTATEFILTER, _prepare_state_cond(cond, s, Regex, SubstitutionString)...)
function _prepare_state_cond(conditions, state, ST, CT)
    if state isa ST
        _state = state
    elseif state isa Symbol
        _state = ST(string(state))
    elseif state isa AbstractString
        _state = ST(state)
    else
        error("Can't handle state $state.")
    end

    if conditions isa Tuple
        conditions = collect(conditions)
    end

    if conditions == nothing
        _cond = nothing
    elseif conditions isa Vector{CT}
        _cond = conditions
    elseif conditions isa Union{Vector, Tuple}
        _cond = CT.(string.(conditions))
    elseif conditions isa CT
        _cond = [conditions]
    elseif conditions isa Union{Symbol, String}
        _cond = [CT(string(conditions))]
    else
        error("Can't handle state conditions $conditions.")
    end
    return _cond, _state
end


export @register_vstate, @register_estate
export @register_vstatefilter, @register_estatefilter
macro register_vstate(exes...)
    args = Expr[]
    for ex in exes
        _state, _cond = _extract_state_cond(ex)
        push!(args, :(register_vstate!($_state, $_cond)))
    end
    Expr(:block, args...)
end
macro register_estate(exes...)
    args = Expr[]
    for ex in exes
        _state, _cond = _extract_state_cond(ex)
        push!(args, :(register_estate!($_state, $_cond)))
    end
    Expr(:block, args...)
end
macro register_vstatefilter(ex)
    _state, _cond = _extract_state_cond(ex)
    :(register_vstatefilter!($_state, $_cond))
end
macro register_estatefilter(ex)
    _state, _cond = _extract_state_cond(ex)
    :(register_estatefilter!($_state, $_cond))
end
function _extract_state_cond(ex)
    if @capture(ex, s_ => cond_)
        return s, cond
    else
        return ex, :(nothing)
    end
end

function _additional_states(states::Vector{Symbol}, dict, matching=Set{Symbol}())
    matchlength = length(matching)
    for (s, cond) in dict
        if cond == nothing
            if Symbol(s) ∉ states
                push!(matching, Symbol(s))
            end
        else
            matches = _match_cond(s, cond, vcat(states, collect(matching)))
            if !isempty(matches)
                new = setdiff(Symbol.(matches), states)
                for n in new
                    push!(matching, n)
                end
            end
        end
    end
    if matchlength != length(matching)
        matching = _additional_states(states, dict, matching)
    end
    return matching
end

function _match_cond(subs, conds, states)
    states_str = string.(states)
    captures = map(conds) do rx
        matches = String[]
        for str in states_str
            if contains(str, rx)
                push!(matches, replace(str, rx=>subs))
            end
        end
        matches
    end
    ## for each condition the capture has to be identical
    reduce(intersect, captures)
end

function _filter_states!(states, dict)
    states_str = String.(states)
    del = Set{Int}()
    for (key, cond) in dict, (stateidx, state) in pairs(states_str)
        m = match(key, state)
        m == nothing && continue
        if isnothing(cond)
            push!(del, stateidx)
        else
            conditions = map(c->replace(state, key => c), cond)
            if conditions ⊆ states_str
                push!(del, stateidx)
            end
        end
    end
    deleteat!(states, sort!(collect(del)))
end

export listvstates, listestates
function listvstates(sol, i::Int)
    nd = _get_nd(sol)
    cf = get_vertexf(nd, i)
    states = vcat(syms(cf), observed_syms(cf))
    add = _additional_states(states, NetworkDynamicsInspector.VSTATES)
    fullstates = vcat(sort!(collect(add)), states)
    _filter_states!(fullstates, VSTATEFILTER)
end

function listestates(sol, i::Int)
    nd = _get_nd(sol)
    cf = get_edgef(nd, i)
    states = syms(cf)
    add = _additional_states(states, NetworkDynamicsInspector.ESTATES)
    fullstates = vcat(sort!(collect(add)), states)
    _filter_states!(fullstates, ESTATEFILTER)
end

listvstates(sol, idxs) = isempty(idxs) ? Symbol[] : mapreduce(i -> listvstates(sol, i), intersect, idxs)
listestates(sol, idxs) = isempty(idxs) ? Symbol[] : mapreduce(i -> listestates(sol, i), intersect, idxs)
listvstates(sol) = listvstates(sol, nv(_get_nd(sol).f.graph))
listestates(sol) = listestates(sol, ne(_get_nd(sol).f.graph))

export empty_lenses!
empty_lenses!() = (empty_statelenses!(), empty_statefilter!(), empty_states!())
empty_statelenses!() = (empty_vstatelenses!(); empty_estatelenses!())
empty_vstatelenses!() = empty!(VSTATELENSES)
empty_estatelenses!() = empty!(ESTATELENSES)
empty_statefilter!() = (empty_vstatefilter!(); empty_estatefilter!())
empty_vstatefilter!() = empty!(VSTATEFILTER)
empty_estatefilter!() = empty!(ESTATEFILTER)
empty_states!() = (empty_vstates!(); empty_estates!())
empty_vstates!() = empty!(VSTATES)
empty_estates!() = empty!(ESTATES)


export printvstates, printestates
printvstates(args...) = _reduce_states(string.(listvstates(args...)))
printestates(args...) = _reduce_states(string.(listestates(args...)))

function _reduce_states(allstates)
    endings = ["_arg", "_mag", "_i", "_r"]
    allstates = map(allstates) do state
        m = match(r"(.*?)(_arg|_mag|_i|_r)$", state)
        if m !== nothing
            m[1]*"_*"
        else
            state
        end
    end
    unique!(allstates)
    sort!(allstates)
end
