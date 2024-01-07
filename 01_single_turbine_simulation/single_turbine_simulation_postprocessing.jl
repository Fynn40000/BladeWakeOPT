#=##############################################################################
# DESCRIPTION
    Postprocessing single turbine simulations
    Plots created for:
        => Blade loads
        => 
=###############################################################################


function single_turbine_simulation_postprocessing(save_path::String, save_path_post::String, run_name::String, R::Float64; 
                                                  # ----- POSTPROCESSING EXECUTION -------------
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
    #OwnFunctions.postprocess_fluiddomain()
  end

end

