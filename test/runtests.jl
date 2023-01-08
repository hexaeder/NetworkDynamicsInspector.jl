using NetworkDynamicsInspector
using Test

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

    u0 = zeros(20)
    u0 .= @views vstatef(sol, p, 1:20, :u_r)(0)[1, :]

    Plots.plot(sol.t, vstatef(sol, p, 1, :_ω)(sol.t))
    Plots.plot!(sol.t, vstatef(sol, p, 3, :_ω)(sol.t))

    Plots.plot(sol.t, vstatef(sol, p, 1:20, :_ω)(sol.t))

    Plots.plot(sol.t, vstatef(sol, p, 1:20, "_ω")(sol.t))

    vstatef(sol, p, 1:20, "_ω")(sol.t);

    values = vstatef(sol, p, 1:20, "_ω")(sol.t)

    values[1,:]

    filter(!isnan, r)
    size(values, 2)
    eachindex
    axes(values,2)


end
