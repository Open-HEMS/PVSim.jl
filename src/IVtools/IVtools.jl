"""
Utilities and models for constructing IV curves.
"""
module IVtools

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Thermal
using ModelingToolkitStandardLibrary.Blocks

using DataInterpolations
using OrdinaryDiffEq

export SingleDiode, VariableResistor
include("single_diode.jl")

export iv_curve
include("iv_curve.jl")

end # module