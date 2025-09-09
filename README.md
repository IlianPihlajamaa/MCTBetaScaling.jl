# MCTBetaScaling.jl
[![Build status (Github Actions)](https://github.com/IlianPihlajamaa/MCTBetaScaling.jl/workflows/CI/badge.svg)](https://github.com/IlianPihlajamaa/MCTBetaScaling.jl/actions)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://IlianPihlajamaa.github.io/MCTBetaScaling.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://IlianPihlajamaa.github.io/MCTBetaScaling.jl/dev)

## Beta-scaling equation

The beta-scaling model is implemented to make it easier to find critical exponents of MCT. The equation is

$$\sigma - \delta t + \lambda (g(t))^2 = \partial_t∫g(t-\tau)g(\tau)d\tau.$$

Here, $\sigma$ is the distance from the critical point, $\lambda$ is the relevant eigenvalue of the stability matrix. $g(t)$ describes the deviation of the order parameter from the plateau. $\delta$ is an optional hopping parameter, defaulting to 0 if not specified. Each of the parameters have to be floating point numbers.  

### Installation

The package is not registered. This means you need to add it from this github repository directly.

```julia
import Pkg
Pkg.activate(".") # activate current folder as virtual environment
Pkg.add("ModeCouplingTheory") # dependency
Pkg.add("Plots") # plotting
Pkg.add(url="https://github.com/IlianPihlajamaa/MCTBetaScaling.jl")
```


### Example
In order to solve the beta-scaling equation, we have to specify the parameters defining the equation and a time-scale `t0` that shifts the results. 
```julia
using ModeCouplingTheory, Plots
λ = 0.7; ϵ = -0.1; t0 = 0.001
equation = BetaScalingEquation(λ, ϵ, t0)
sol = solve(equation, TimeDoublingSolver(t_max=10^4.))
plot(log10.(sol.t), log10.(abs.(sol.F)), ylabel="log_{10}(|g(t)|)", xlabel="log_{10}(t)", label="g(t)")
```

![image](docs/src/images/beta.png)


In the figure, the slopes of the straight lines are given by the parameters $-a$ and $b$, which describe the relaxation towards and away from the plateau value of the correlator. These exponents are automatically computed, and are stored in `equation.coeff.a` and `equation.coeff.b`.

## Stochastic Beta-relaxation

The code also implements a numerical solver for the SBR equation:

$$\sigma(x) +\alpha \nabla^2g(x,t)- \delta t + \lambda (g(x,t))^2 = \partial_t∫g(x,t-\tau)g(x,\tau)d\tau.$$

in 1, 2, and 3 dimensions with periodic boundary conditions. See the [documentation](https://ilianpihlajamaa.github.io/MCTBetaScaling.jl/dev/#Stochastic-Beta-Relaxation-(SBR)).



### References
Götze, J Phys Condens Matter 2, 8485 (1990)

## To do:
Add code to identify the eigenvalue `\lambda`. (update: code is written, but no tests but docs yet.)
