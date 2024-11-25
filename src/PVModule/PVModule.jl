"""
Library of DC cell / module models.
"""
module PVModule

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Thermal
using ModelingToolkitStandardLibrary.Blocks

using DataInterpolations

export SingleDiode, VariableResistor
include("single_diode.jl")

end