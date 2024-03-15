#=##############################################################################
# DESCRIPTION
    Functions used for optimization.

# AUTHORSHIP
  * Author          : Fynn Gerhardy
  * Email           : fygerh@gmail.com
  * Created         : Mar 2024
  * Last updated    : Mar 2024
  * License         : -
=###############################################################################

# Function to create sampling points via Latin Hypercube Sampling
function get_solution_points(method, Nsamples, Ndimensions, dim_ranges, generations)

    if method == "random"
        lhs = scaleLHC(randomLHC(Nsamples, Ndimensions), dim_ranges)
    
    elseif method == "optimized"
        optimized_plan, fitness_history = LHCoptim(Nsamples, Ndimensions, generations) 
        lhs = scaleLHC(optimized_plan, dim_ranges)
        #println(fitness_history)

    end
    
    return lhs
end

"
Calculate the Bernstein polynomial B_i,n(t).
"
function bernstein_poly(i, n, t)
    return binomial(n, i) * (t ^ i) * ((1 - t) ^ (n - i))
end


"
Calculate the value of a Bézier curve for a specific parameter t.
"
function bezier_curve(control_points, t)
    
    n = length(control_points) - 1
    curve_point = zeros(size(control_points[1]))
    for (i, point) in enumerate(control_points)
        curve_point .+= point * bernstein_poly(i, n, t)
    end
    return curve_point
end