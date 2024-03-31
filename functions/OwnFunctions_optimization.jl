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

"
Function to create sampling points via Latin Hypercube Sampling.
"
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
Calculate the value of a Bézier curve for a specific parameter t.
"
function bezier_point(t::Float64, control_points::Array{Array{Float64,1},1})
    n = length(control_points) - 1
    result = [0.0, 0.0]
    for i in 0:n
        result .+= binomial_coefficient(n, i) * (1 - t)^(n - i) * t^i * control_points[i + 1]
    end
    return result
end

function binomial_coefficient(n, k)
    return factorial(n) / (factorial(k) * factorial(n - k))
end


"
Return a scaled twist distribution at x_values.
"
function get_scaled_twist(x_values, scaling_curve_x, scaling_curve_y, twist_ref, ref_LB, ref_UB)
    # interpolate scalig factors from scaling curve at twist locations
    itp = LinearInterpolation(scaling_curve_x, scaling_curve_y, extrapolation_bc=Line())
    scaling_factors = itp(x_values)

    scaled_twist = []
    
    for i in 1:length(scaling_factors)
        SF = scaling_factors[i]
        if SF < 0
            push!(scaled_twist, (twist_ref[i] + (twist_ref[i] - ref_LB[i]) * SF))
            #scaled_twist[i] = twist_ref[i] + (twist_ref[i] - ref_LB[i]) * SF
        elseif SF >= 0
            push!(scaled_twist, (twist_ref[i] + (ref_UB[i] - twist_ref[i]) * SF))
            #scaled_twist[i] = twist_ref[i] + (ref_UB[i] - twist_ref[i]) * SF
        end
    end
    
    return scaled_twist
end