module NetworkDynamicsInspector

using SciMLBase: AbstractDiffEqFunction, ODESolution, ODEProblem
using NetworkDynamics: syms, observed_syms
using MacroTools
using Printf

export PRecord, record!, register_pd_lenses!

include("utils.jl")
include("PRecord.jl")
include("lenses.jl")
include("FavSelect.jl")
include("interactive.jl")
include("pdlenses.jl")


end
