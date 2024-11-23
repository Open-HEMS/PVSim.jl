module PVSim

using ModelingToolkit, ModelingToolkitStandardLibrary
using ModelingToolkit: t_nounits as t
using OrdinaryDiffEq
using DataInterpolations

include("PVModule/PVModule.jl")

end # module