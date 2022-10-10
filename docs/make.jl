using Sensemakr
using Documenter

DocMeta.setdocmeta!(Sensemakr, :DocTestSetup, :(using Sensemakr); recursive=true)

makedocs(;
    modules=[Sensemakr],
    authors="Rodrigo Grijalba, Alexander Quispe",
    repo="https://github.com/d2cml-ai/Sensemakr.jl/blob/{commit}{path}#{line}",
    sitename="Sensemakr.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://d2cml-ai.github.io/Sensemakr.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/d2cml-ai/Sensemakr.jl",
    devbranch="master",
)
