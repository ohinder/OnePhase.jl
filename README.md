# A one-phase interior point method for nonconvex optimization

This package is the implementation of a one-phase interior point method that finds KKT points of optimization problems of the form:
<!--
$$
\min f(x)  \quad \text{s.t.} \quad a(x) \le 0.
$$
-->

![min f(x) s.t. a(x) < 0](misc/problem-statement.gif)

where the functions ![f : R^n -> R](misc/f.gif) and ![a : R^n -> R^m](misc/a.gif) are twice differentiable. The one-phase algorithm also handles bound constraints and nonlinear equalities.

*Currently, the package is in development.* Although you are welcome to try it out. Please let me know if there are any bugs etc. Note that the code is generally significantly slower than Ipopt in terms of raw runtime, particularly on small problems (the iteration count is competitive). However, we recommend trying our one-phase IPM if: Ipopt is failing to solve, the problem is very large or might be infeasible.

## How to install

```julia
Pkg.clone("https://github.com/ohinder/advanced_timer.jl")
Pkg.clone("https://github.com/ohinder/OnePhase.git")
```

## How to use with JuMP

Here is a simple example where a [JuMP](http://www.juliaopt.org/JuMP.jl/0.18/JuMP) model is passed to the one-phase solver

```julia
using OnePhase, JuMP

m = Model(OnePhase())
@variable(m, x, start=-3)
@objective(m, Min, x)
@NLconstraint(m, x^2 >= 1.0)
@NLconstraint(m, x >= -1.0)

status = solve(m)
```

To interpret the solver output see [here](solver_output.md)

## Example using CUTEst

Install [CUTEst](http://juliasmoothoptimizers.github.io/CUTEst.jl/latest/) then run
```julia
using OnePhase, CUTEst
nlp = CUTEstModel("CHAIN")
iter = one_phase_solve(nlp);
```

## Feedback?

If you have found some bug or think there is someway I can improve the code feel free to contact me! My webpage is https://stanford.edu/~ohinder/.
