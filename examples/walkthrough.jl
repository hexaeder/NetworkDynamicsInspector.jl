#=
# Feature walkthrough

In this example, we'll go through the main features of the
NetworkDynamicsInspector package.

This package is all about inspecting the dynamic states of a network. To do so,
we start by defining dynamic problem. The problem is stated in normal
`PowerDynamics.jl` convention. This means, the edge states visible to the nodes are
the complex currents `(i_r, i_i)` flowing into the node, the node states visible to the
edges are the complex currents `(u_r, u_i)` established at the node.
=#
using NetworkDynamics
using NetworkDynamicsInspector
using Graphs
using GLMakie
using GraphMakie
using OrdinaryDiffEq

g = SimpleGraph(5)
add_edge!(g, 1, 2); add_edge!(g, 1, 4); add_edge!(g, 2, 3); add_edge!(g, 2, 4);
add_edge!(g, 2, 5); add_edge!(g, 3, 5); add_edge!(g, 4, 5);

ntypes = [:gen, :load, :gen, :load, :load]

node_color = map(t->t==:gen ? Makie.wong_colors()[5] : Makie.wong_colors()[2], ntypes)
fig, ax, p = graphplot(g; ilabels=repr.(1:5), node_color)
hidedecorations!(ax); hidespines!(ax)
fig #hide

#=
Here, the blue nodes represent generators while the orange nodes represent loads.

Next we can define our node models. The swing node is typical represenation of
the swing quation and has two parameters, mechanical turbine power `P` and
damping `D`. It has three states, the internal angle `δ` and both components of
the complex voltage.
=#

# du = rand(3)
# u = rand(3)
# edges = [[1.0,2.2], [1.2,2.2]]
# p = (1,2,3)

# @code_warntype nd.f.unique_vertices![1].f(du, u, edges, p, 0.0)
# @allocated nd.f.unique_vertices![1].f(du, u, edges, p, 0.0)
# @btime $(nd.f.unique_vertices![1].f)($du, $u, $edges, $p, 0.0)

function swing_vertex!(dv, v, edges, (P, M, D), t)
    u_r, u_i, ω = v
    i = total_current(edges)
    u = u_r + u_i * im
    δ = angle(u)

    Pel = real(u*conj(i))
    dω = M * (P + Pel - D*ω)

    du_r = -sin(δ)*ω
    du_i =  cos(δ)*ω

    dv .= du_r, du_i, dω
    nothing
end
swing = ODEVertex(; f=swing_vertex!, dim=3, mass_matrix=1.0, sym=[:u_r, :u_i, :ω], psym=[:P, :M, :D])

#=
The load vertices is implemented as a constraint. It just forces `u_r` and `u_i` such that `P=u⋅i*`:
=#
function pq_vertex!(dv, v, edges, (P,), t)
    current = total_current(edges)
    voltage = v[1] + v[2] * im
    residual = P + voltage * conj(current)
    dv[1] = real(residual)
    dv[2] = imag(residual)
    nothing
end
load = ODEVertex(; f=pq_vertex!, dim=2, mass_matrix=0.0, sym=[:u_r, :u_i], psym=[:P])

#=
Both vertices need a helper function to find the total current inflow at the node:
=#

function total_current(edges)
    ## ND convention: all flows entering a node are positive
    current = 0.0im
    @inbounds for e in edges
        current += e[1] + e[2]*im
    end
    current
end

#=
For the edge model, we take the RMS represenation (static edge) of an edge with
complex admittance (RL line), where the admittance `Y` is the onlye parameter.
=#
function complex_admittance_edge!(e, v_s, v_d, (Y,), t)
    src_voltage = v_s[1] + v_s[2] * im
    dst_voltage = v_d[1] + v_d[2] * im
    ## If current is flowing away from the source, it is negative at the source.
    complex_current = Y * (dst_voltage - src_voltage)
    e[1] = real(complex_current)
    e[2] = imag(complex_current)
    nothing
end
edge = StaticEdge(; f=complex_admittance_edge!, dim=2, sym=[:i_r, :i_i], psym=[:Y], coupling=:antisymmetric)

#=
Now we can finally create the nd object, the parameters and simulate:
=#

vertex_list = map(t->t==:gen ? swing : load, ntypes)
nd = network_dynamics(vertex_list, edge, g)

p_vert = [( 2.0, 1.0, 0.5), # inj and damping gen 1
          (-1.0, NaN, NaN),  # draw load 2
          ( 1.0, 1.0, 0.5), # inj and damping gen 3
          (-1.0, NaN, NaN),  # draw load 2
          (-1.0, NaN, NaN)]  # draw load 2
## fillig up with `NaN` here helps us to achive concret parameter array eltype which is important for performance
@assert isconcretetype(eltype(p_vert))

## the edges are homogenious
p_edge = 5im

## network apraemters are the combination of vertex and edge parameters
p = (p_vert, p_edge)

# For the initial condition, we first initialize everything to 0 and than rais all the real parts of the voltage to 0
u0 = zeros(length(nd.syms)) #+ 0.1*randn(length(u0))
u0[idx_containing(nd, r"^u_r")] .= 1

prob = ODEProblem(nd, u0, (0, 100), p)
sol = solve(prob, Rodas5P())


register_pd_lenses!()
inspect_solution(sol)

# nd.syms sol[:ω_1]


# NetworkDynamicsInspector.register_pd_lenses!()
# inspect_solution(sol)

# ax.linestyle
# ax.palette
# propertynames(ax).
# ax.scene.palette

# fig, ax, p = lines(rand(10))
# ax.scene.theme.palette.linestyle
# ax.scene.theme.linestyle


# ax.theme
