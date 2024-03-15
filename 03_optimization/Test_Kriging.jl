
#include(joinpath("/home/fynn/Repositories/BladeWakeOPT/03_optimization/Test_Kriging.jl"))

using Surrogates
using Plots
#default()

# JULIA EXAMPLE 1D optimization
# https://www.sfu.ca/~ssurjano/forretal08.html
# Forrester et al. (2008) Function
#f(x) = (6 * x - 2)^2 * sin(12 * x - 4)
function f(x)
    value = (6 * x - 2)^2 * sin(12 * x - 4)
    return value
end

n_samples = 4
lower_bound = 0.0
upper_bound = 1.0

xs = lower_bound:0.001:upper_bound

x = sample(n_samples, lower_bound, upper_bound, SobolSample())
y = f.(x)

scatter(x, y, label="Sampled points", xlims=(lower_bound, upper_bound), ylims=(-7, 17))
plot!(xs, f.(xs), label="True function", legend=:top)

kriging_surrogate = Kriging(x, y, lower_bound, upper_bound);

plot(x, y, seriestype=:scatter, label="Sampled points", xlims=(lower_bound, upper_bound), ylims=(-7, 17), legend=:top)
plot!(xs, f.(xs), label="True function", legend=:top)
plot!(xs, kriging_surrogate.(xs), label="Surrogate function", ribbon=p->std_error_at_point(kriging_surrogate, p), legend=:top)


@show surrogate_optimize(f, SRBF(), lower_bound, upper_bound, kriging_surrogate, SobolSample())

scatter(x, y, label="Sampled points", ylims=(-7, 7), legend=:top)
plot!(xs, f.(xs), label="True function", legend=:top)
plot!(xs, kriging_surrogate.(xs), label="Surrogate function", ribbon=p->std_error_at_point(kriging_surrogate, p), legend=:top)



# JULIA EXAMPLE 2D optimization
function branin(x)
    x1=x[1]
    x2=x[2]
    a=1;
    b=5.1/(4*pi^2);
    c=5/pi;
    r=6;
    s=10;
    t=1/(8pi);
    a*(x2-b*x1+c*x1-r)^2+s*(1-t)*cos(x1)+s
  end
  
  n_samples = 10
  lower_bound = [-5.0, 0.0]
  upper_bound = [10.0, 15.0]
  
  xys = sample(n_samples, lower_bound, upper_bound, GoldenSample())
  zs = branin.(xys);
  println(xys)
  println(zs)
  
  kriging_surrogate = Kriging(xys, zs, lower_bound, upper_bound, p=[2.0, 2.0], theta=[0.03, 0.003])
  
  surrogate_optimize(branin, SRBF(), lower_bound, upper_bound, kriging_surrogate, SobolSample(); maxiters = 100, num_new_samples = 10)