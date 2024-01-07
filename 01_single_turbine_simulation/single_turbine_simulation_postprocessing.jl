#=##############################################################################
# DESCRIPTION
    Postprocessing single turbine simulations
    Plots created for:
        => Blade loads
        => 
=###############################################################################


function single_turbine_simulation_postprocessing(save_path::String, save_path_post::String, run_name::String, R::Float64, AOA::Float64; 
                                                  # ----- POSTPROCESSING EXECUTION ----------------
                                                  plot_bladeloads::Bool=true,     # postprocessing the bladeloads?
                                                  postprocess_fdom::Bool=true,    # postprocessing the fluiddomain?
                                                  # ----- SETTINGS FOR POSTPROCESSING -------------
                                                  rev_to_average_idx=1,           # Revolution to wich the postprocessing should be applied on
                                                  nrevs_to_average=1,             # number of Revolutions to average for postprocessing the bladeloads
                                                  num_elements::Int64=1,          # number of blade elements per blade
                                                  debug::Bool=false,              # postprocess dimensionless coefficients too? => NOTE: debug statement must be set to true for uns.run_simulation. Otherwise the simulation files will not contain the coefficient data.
                                                  suppress_plots::Bool=true       # suppresses the plots to show up on the display
                                                  )

  println("\nPostprocessing...\n")
  
  # postprocess the blade loadings if desired
  if plot_bladeloads
    # postprocess the statistics first by applying function uns.postprocess_statistics(...)
    OwnFunctions.postprocess_statistics(save_path, run_name;
                                        rev_to_average_idx=rev_to_average_idx, 
                                        nrevs_to_average=nrevs_to_average)
    
    #statistics_save_path = joinpath(save_path, run_name)

    # create a folder to save the plots in
    if isdir(save_path_post)
      rm(save_path_post, recursive=true)
      mkdir(save_path_post)
    else
      mkdir(save_path_post)
    end

    # postprocess statistics of the "_vlm" statistics file regarding blade 1
    OwnFunctions.plot_blade_loading(save_path, save_path_post, run_name, R; 
                                    file_marker="_vlm",
                                    rev_to_average_idx=rev_to_average_idx, 
                                    nrevs_to_average=nrevs_to_average,
                                    num_elements=num_elements,
                                    suppress_plots=suppress_plots)
    # postprocess statistics of the "_loft" statistics file regarding blade 1
    OwnFunctions.plot_blade_loading(save_path, save_path_post, run_name, R; 
                                    file_marker="_loft",
                                    rev_to_average_idx=rev_to_average_idx, 
                                    nrevs_to_average=nrevs_to_average,
                                    num_elements=num_elements,
                                    debug=debug,
                                    suppress_plots=suppress_plots)

    println("\nBlade load plots saved under $(save_path_post)\n")
  end

  if postprocess_fdom

    # read simulation data
    #simdata = CSV.read(joinpath(save_path, run_name*"_convergence.csv"), DataFrame)
    # Calculate nsteps_per_rev
    #nsteps_per_rev[run_name] = ceil(Int, 360 / (simdata[2, 1] - simdata[1, 1]))

    timesteps_to_evaluate = [140]
    OwnFunctions.postprocess_fluiddomain(save_path, run_name, R, AOA, timesteps_to_evaluate;
                                         # ----- GRID OPTIONS ----------------  
                                         x_resolution    = 50,                   # discretization of grid in x-direction
                                         y_resolution    = 50,                   # discretization of grid in y-direction
                                         z_resolution    = 50,                   # discretization of grid in z-direction
                                         x_bound_min     = -0.1,                    # minimum bounds in x-direction (bound_factor*2*R in meters)
                                         y_bound_min     = -0.6,                  # minimum bounds in y-direction (bound_factor*2*R in meters)
                                         z_bound_min     = -0.6,                  # minimum bounds in z-direction (bound_factor*2*R in meters)
                                         x_bound_max     = 1,                     # maximum bounds in x-direction (bound_factor*2*R in meters)
                                         y_bound_max     = 0.6,                   # maximum bounds in y-direction (bound_factor*2*R in meters)
                                         z_bound_max     = 0.6,                   # maximum bounds in z-direction (bound_factor*2*R in meters)
                                         )
  end

end

