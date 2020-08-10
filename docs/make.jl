using Ising2D
using Documenter

makedocs(;
    modules=[Ising2D],
    authors="genkuroki <genkuroki@gmail.com> and contributors",
    repo="https://github.com/genkuroki/Ising2D.jl/blob/{commit}{path}#L{line}",
    sitename="Ising2D.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://genkuroki.github.io/Ising2D.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/genkuroki/Ising2D.jl",
)
