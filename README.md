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

## Generate julia function from Symbolic expression

```bash
cd symbolic
julia generate_functions.jl
```

## Refs

* [rkprol branch](https://bitbucket.org/cactuscode/cactusnumerical/src/11f6c32fcacc9c5e1f7fce0c49b94159e27957b2/?at=ianhinder%2Frkprol)

* [wiki](https://docs.einsteintoolkit.org/et-docs/Prolongation)

* [pseudocode](https://docs.google.com/document/d/1M65w_8keIf6Ypas43YXf1CA-tbuhOssFPpBPlGS1tD0/edit)

## RK4 Steps

* Initialise

* March the first substep for all levels (coarse to fine)

```bash
 ------

        ------
               ------
 ------ ------ ------
```

* March the other substeps to the same time slice

```julia
for s = 2:2^(lmax-1)  # from second to final substeps of the finest level
    for l = 2:lmax  # march all levels execpt the coarest (coarset to fine)
        if l == lmax || (levs[l].time == levs[l+1].time)
            if l < lmax
                Restriction  # from l+1 to l
            Prolongation  # from l-1 to l
            rk4(levs[l])
```

* Restriction all levels
