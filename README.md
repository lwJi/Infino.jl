# Infino.jl

Infrastructure to test Subcycling

[![Build Status](https://github.com/lwJi/Infino.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lwJi/Infino.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Installation

1. Clone the repo: `git clone git@github.com:lwJi/Infino.jl.git`

2. Add **Infino** to your Pkg environment: `] dev \path\to\the\repo`

3. Import **Infino**: `using Infino`

## Run an example

```bash
julia run/Subcycling.jl test/integration/scalarwave_3levels.toml
julia tool/generate_gifs.jl test/integration/scalarwave_3levels/
```

## Refs

* [rkprol branch](https://bitbucket.org/cactuscode/cactusnumerical/src/11f6c32fcacc9c5e1f7fce0c49b94159e27957b2/?at=ianhinder%2Frkprol)

* [wiki](https://docs.einsteintoolkit.org/et-docs/Prolongation)
