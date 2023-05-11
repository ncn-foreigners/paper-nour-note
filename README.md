# Readme

Repository for the paper "A note on Nour (1982) dual-system estimator". See paper [here](paper/note_on_nour.pdf)

## Setup

1. Install Julia ([download here](https://julialang.org/downloads/))
2. Clone repository with `git clone https://github.com/ncn-foreigners/paper-nour-note.git`
3. Run `cd paper-nour-note`
4. Run `julia` in this folder (you may check the working directory using `pwd()`)
5. Run the following code

```julia
using Pkg
Pkg.activate(".")

using IJulia
notebook(dir=".")
```

## Reproducing

### Environment

```julia
Julia Version 1.9.0
Commit 8e630552924 (2023-05-07 11:25 UTC)
```

### Packages

```julia
  [a93c6f00] DataFrames v1.5.0
⌃ [31c24e10] Distributions v0.25.87
  [7073ff75] IJulia v1.24.0
  [b964fa9f] LaTeXStrings v1.3.0
⌅ [23fbe1c1] Latexify v0.15.18
  [442fdcdd] Measures v0.3.2
⌃ [91a5bcdd] Plots v1.38.9
  [9a3f8284] Random
  [10745b16] Statistics v1.9.0
```

