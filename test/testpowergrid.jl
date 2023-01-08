using NetworkDynamics
using Graphs
using Random
using OrdinaryDiffEq

Ngen  = 10
Nload = 10

g = watts_strogatz(Ngen + Nload, 4, 0.05; seed=1)
# graphplot(g; layout=Stress())


function complex_admittance_edge!(e, v_s, v_d, p, t)
    source_voltage = v_s[1] + v_s[2] * im
    destination_voltage = v_d[1] + v_d[2] * im
    # If current is flowing away from the source, it is negative at the source.
    complex_current = p[1] * (destination_voltage - source_voltage)
    e[1] = real(complex_current)
    e[2] = imag(complex_current)
    nothing
end

function total_current(edges)
    # Keeping with the convention of negative sign for outging current
    current = 0.0im
    for e in edges
        current -= e[1] + e[2] * im
    end
    current
end

function swing_vertex!(dv, v, edges, (P, D), t)
    i = total_current(edges)
    u = v[1] + v[2] * im
    dv[3] = P - D * v[3] + real(u * conj(i))
    dvolt = 1.0im * v[3] * u - (abs(u) - 1) * u
    dv[1] = real(dvolt)
    dv[2] = imag(dvolt)
    nothing
end

function pq_vertex!(dv, v, edges, p, t)
    current = total_current(edges)
    voltage = v[1] + v[2] * im
    residual = p[1] + voltage * conj(current)
    dv[1] = real(residual)
    dv[2] = imag(residual)
    nothing
end

pq    = ODEVertex(; f=pq_vertex!, dim=2, mass_matrix=0.0, sym=[:u_r, :u_i])
swing = ODEVertex(; f=swing_vertex!, dim=3, mass_matrix=1.0, sym=[:u_r, :u_i, :ω])
edge = StaticEdge(; f=complex_admittance_edge!, dim=2, sym=[:i_r, :i_i], coupling=:antisymmetric)

# vertex_list = swing;#[swing, swing, swing]
vertex_list = vcat(fill(swing, Ngen), fill(pq,Nload))

nd = network_dynamics(vertex_list, edge, g)

Pgen  = 1.0
Pload = - Ngen*Pgen / Nload

nodep = vcat(fill((Pgen, 1.1), Ngen),
             fill((Pload, NaN), Nload))

edgep = -5im        # line admittance
p = (nodep, edgep)
tspan = (0, 10)

Random.seed!(1)
u0 = 0.1*randn(length(nd.syms))
u0[idx_containing(nd, r"^u_r")] .+= 1
# nd.syms .=> u0

prob = ODEProblem(nd, u0, tspan, p)
sol = solve(prob, Rodas4())

# Plots.plot(sol)
