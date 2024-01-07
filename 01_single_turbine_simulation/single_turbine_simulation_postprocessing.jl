#=##############################################################################
# DESCRIPTION
    Postprocessing single turbine simulations
    Plots created for:
        => Blade loads
        => 
=###############################################################################


function single_turbine_simulation_postprocessing(save_path::String, save_path_post::String, run_name::String, R::Float64; 
                                                  plot_bladeloads::Bool=true, 
                                                  postprocess_fdom::Bool=true, 
                                                  rev_to_average_idx=1, 
                                                  nrevs_to_average=1,
                                                  num_elements::Int64=1,
                                                  debug::Bool=false)

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
                                    num_elements=num_elements)
    # postprocess statistics of the "_loft" statistics file regarding blade 1
    OwnFunctions.plot_blade_loading(save_path, save_path_post, run_name, R; 
                                    file_marker="_loft",
                                    rev_to_average_idx=rev_to_average_idx, 
                                    nrevs_to_average=nrevs_to_average,
                                    num_elements=num_elements,
                                    debug=debug)
  end

  if postprocess_fdom
    #OwnFunctions.postprocess_fluiddomain()
  end

end

