using PVSim
using Documenter

DocMeta.setdocmeta!(PVSim, :DocTestSetup, :(using PVSim); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
    file for file in readdir(joinpath(@__DIR__, "src")) if
    file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [PVSim],
    authors = "Stefan de Lange",
    repo = "https://github.com/Open-HEMS/PVSim.jl/blob/{commit}{path}#{line}",
    sitename = "PVSim.jl",
    format = Documenter.HTML(; canonical = "https://Open-HEMS.github.io/PVSim.jl"),
    # pages = ["index.md"; numbered_pages],
    pages = ["index.md"],
)

deploydocs(; repo = "github.com/Open-HEMS/PVSim.jl")
