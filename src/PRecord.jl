struct PRecord{PT}
    t::Vector{Float64}
    p::Vector{PT}
end
Base.eltype(::PRecord{PT}) where {PT} = PT

PRecord(prob::ODEProblem) = PRecord([Float64(prob.tspan[1])], [prob.p])
PRecord(pr::PRecord) = pr
PRecord(p) = PRecord([-Inf], [p])

function record!(pr::PRecord, integrator)
    push!(pr.t, integrator.t)
    push!(pr.p, deepcopy(integrator.p))
end

function (pr::PRecord)(t::Number; direction=:right)
    fun = direction==:right ? tr->tr>t : tr->tr≥t
    idx = findfirst(fun, pr.t)
    if isnothing(idx)
        pr.p[end]
    elseif idx == 1
        pr.p[begin]
    else
        pr.p[idx - 1]
    end
end

function (pr::PRecord{PT})(ts) where {PT}
    p = Vector{PT}(undef, length(ts))
    lastt = -Inf
    for (i, t) in enumerate(ts)
        if t == lastt
            p[i] = pr(t; direction=:right)
        else
            p[i] = pr(t; direction=:left)
        end
        lastt = t
    end
    p
end
