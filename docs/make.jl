using Documenter, MinFEM

Home = "Home" => "index.md"

Introduction = "The Finite Element Method" => "fem.md"

GettingStarted = "gettingstarted.md"

Visualization = "Visualization" => "paraview.md"

Examples = "Examples" => [
    "A Poisson Problem" => "examples/poisson.md",
    "A Boundary Source Problem" => "examples/boundary_source.md",
    "A Vector-Valued Problem" => "examples/elasticity.md",
    "A Semi-Linear Problem" => "examples/semilinear.md",
    "A Time-Dependent Problem" => "examples/parabolic.md"
    ]

Library = "Library" => [
    "Public" => "lib/public.md",
    "Internals" => "lib/internals.md"
    ]

License = "License" => "license.md"

PAGES = [
    Home,
    Introduction,
    GettingStarted,
    Visualization,
    Examples,
    Library,
    License
    ]

FORMAT = Documenter.HTML(
            prettyurls = true,
            assets = ["assets/favicon.ico"]
        )

makedocs(
    modules = [MinFEM],
    sitename = "MinFEM.jl",
    authors = "Martin Siebenborn, Henrik Wyschka",
    format = FORMAT,
    pages = PAGES
)

deploydocs(
    repo = "github.com/MinFEM/MinFEM.jl.git"
)
