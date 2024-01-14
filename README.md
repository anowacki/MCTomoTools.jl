# MCTomoTools

MCTomoTools is a Julia package to deal with the input to and output
from the [MCTomo](https://github.com/xin2zhang/MCTomo) Bayesian
tomography software.


## Installation and requirements

MCTomoTools.jl requires Julia v1.8+.

Install it using Julia's inbuilt package manager like so:

```julia
julia> import Pkg

julia> Pkg.add(url="https://github.com/anowacki/MCTomoTools.jl")
```


## MCTomo file structure

MCTomoTools.jl requires that you have the input and output from a previous
MCTomo run.  This means that within the directory where the program was run,
you have the following files and folders for `n` chains:

```
.
├── MCTomo.inp
├── bsources.dat  <-- File name set in MCTomo.inp
└── Results
    ├── ini
    │   ├── InitialSample_1.dat … InitialSample_n.dat
    │   └── InitialSigma_1.dat … InitialSigma_n.dat
    ├── samples_1.dat … samples_n.dat
    │ # Files below are not needed by MCTomoTools.jl
    ├── likelihood_1.dat … likelihood_n.dat
    ├── noise_body_1.dat … noise_body_n.dat
    └── temperatures_1.dat … temperatures_n.dat
```

Other files will be present but are not needed by this package.


## Sampling from a previous MCTomo run

MCTomoTools allows you to use the chains of proposed steps generated during a
previous MCTomo run to sample the posterior probability density function of
the velocity field and source locations.

### Loading run parameters
To do this, first load the run parameters from the `MCTomo.inp` file:

```julia
julia> using MCTomoTools

julia> settings = read_settings("MCTomo.inp");
```

Here we assume that the commands are being run in the same directory as
the MCTomo run, but this is not required.
You may pass a relative or absolute path to `read_settings`; the path
is stored in `settings`, making it convenient to read other files, so long
as you have not changed the default file structure of an MCTomo run.

### Reading a chain
Sampling involves reading an initial model, then perturbing it by each
accepted proposal in the chain.  To read a set of proposals, do this:

```julia
julia> chain_index = 2;

julia> steps = read_raw_samples(settings, chain_index)
1979000-element Vector{MCTomoTools.RawSample}:
 MCTomoTools.RawSample(7, 0, 0, 0x000000000000001c, 532.9068706355182, 2.155027912402892, -752.3633358306661, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
 ⋮
```

`chain_index` is an integer giving the chain number.  By default, all proposals
are read in, but this can be limited with the `n` keyword argument to only the
first `n` proposals.

### Reading an initial model
To sample from a chain, one must start with the initial model.  This is read in
like so:

```julia
julia> model = Model(settings, chain_index)
Model{Float64}:
 nodes: (63 nodes)
 body_noise0: [0.0, 0.0]
 body_noise1: [0.4985924041826726, 0.4662786807087899]
 surfacewave_noise0: [0.0010113237269714762]
 surfacewave_noise1: [0.012812444780211069]
 locations: (639 sources)
```

Optionally, you can construct a `Chain` from the initial model and the update
steps:

```julia
julia> chain = Chain(model, steps)
```

### Sampling
Now, for each chain used, you can sample from the chain.  There are currently
three functions to do this:

1. `sample_grid`.  This can be used to reconstruct the velocity field on an
   arbitrary grid of points.
2. `sample_source`.  This reconstructs the positions (in space and time) of
   an individual event throughout the chain.
3. `sample_chain`.  This runs a function of your choice on each sample of the
   chain.

(As well as the copying versions, there are also in-place versions of the above
with trailing exclamation marks `!`.)

#### Sampling velocities with `sample_grid`
To use `sample_grid`, first create a grid of points covering the area you want,
and with the sampling density you require.  If you are unsure of this, you
can easily create a new, empty grid by passing `settings` to `Grid` like so:

```julia
julia> grid = Grid(settings)
101×101×51 Grid{Float64}:
 x: 460.0:0.13:473.0
 y: 365.0:0.15:380.0
 z: 0.0:0.1:5.0
 data: [values in range (NaN, NaN)]
```

If you know which area you want to sample, however, you can call the other
`Grid` constructors—see the docstring for help.

Armed with a `Grid`, you then pass this in alongside the chain of steps and
the initial model to generate a set of velocity grids.

```julia
julia> vp_grids = sample_grid(model, grid, steps, :vp, thin=10_000);
```

As well as `:vp` for P-wave velocity, you can also pass `:vs` for S-wave velocity
(both in km/s) or `:density` (in g/cm³), though note density is simply scaled from
_V_<sub>P</sub> and is not independently sampled.

#### Sampling source locations with `sample_source`
`sample_source` takes the index of the source of interest and returns a set of
locations:

```julia
julia> source_index = 1;

julia> source1_locations = sample_source(model, steps, source_index; thin=1000)
1158-element Vector{NamedTuple{(:x, :y, :z, :t), NTuple{4, Float64}}}:
 (x = 465.861, y = 370.003, z = 0.9, t = 0.0)
 ⋮
```

Note that each element of the returned vector is a named tuple of:
- `x`: Easting in km
- `y`: Northing in km
- `z`: Depth in km
- `t`: Difference from original origin time in s

You can also return a `Vector{Vector{NamedTuple}}` of all events using
`sample_sources`.

#### General sampling with `sample_chain`
`sample_chain` allows you to perform any action at each sample.  At each proposed
step, the function passed as the first argument to `sample_chain` is called with
a single argument, `sample`.  This is a named tuple
`(; isample::Int, model::Model)`, where `isample` is the number of the proposal,
and `model` is a `MCTomoTools.Model` giving the state of the model after proposal
number `isample`.

For example, to calculate the average number of cells across the chain, you
could do:

```julia
julia> mean_nnodes = let sum_nnodes = 0
           nsamples = sample_chain!(model, steps; thin=1) do sample
               sum_nnodes += length(sample.model.nodes)
           end
           sum_nnodes/nsamples
       end
```


## Relationship to MCTomo

MCTomoTools.jl is meant to make post-processing MCTomo runs easier for
Julia users.  It does not perform any inference itself.  Although we are
users of MCTomo, we are not its primary authors, and any questions about
MCTomo itself should be directed via the
[MCTomo GitHub repo](https://github.com/xin2zhang/MCTomo).

## Licensing

MCTomoTools.jl is distributed under the MIT licence, and is able to be so
because it is not a derivative work of MCTomo, nor does it link to it.

The primary author is Andy Nowacki (@anowacki).

## Contributions

Bug reports and improvements suggestions can be made at the package's GitHub
home by
[creating an issue](https://github.com/anowacki/MCTomoTools.jl/issues/new/choose),
and code contributions can be made by
[forking the repo](https://github.com/anowacki/MCTomoTools.jl/fork) and opening a
[pull request](https://github.com/anowacki/MCTomoTools.jl/compare) against your
branch with the new feature or fix.
