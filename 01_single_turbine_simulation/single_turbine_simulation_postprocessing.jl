#=##############################################################################
# DESCRIPTION
    Postprocessing single turbine simulations
    Plots created for:
        => Blade loads
        => 
=###############################################################################


function single_turbine_simulation_postprocessing(save_path::String, save_path_post::String, run_name::String)

  println("\nPostprocessing...\n")
  
  
  OwnFunctions.plot_blade_loading(save_path, save_path_post, run_name)

  
  
end

