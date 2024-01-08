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


  file_suffixes = []
  if postprocess_fdom

    # read simulation data
    #simdata = CSV.read(joinpath(save_path, run_name*"_convergence.csv"), DataFrame)
    # Calculate nsteps_per_rev
    #nsteps_per_rev[run_name] = ceil(Int, 360 / (simdata[2, 1] - simdata[1, 1]))

    timesteps_to_evaluate = [140]



    # --------------- postprocess x-y-plane -------------------

    z_loc = 0   # z coordinate location of plane = z_loc*2*R in meters
    file_suffix = "_x_y_atz$(z_loc)D"
    push!(file_suffixes, file_suffix)
    vol_thickness = 0 # plane thickness in meters => set this value to 0 to get a 2D solution (plane)

    # grid resolution
    x_res = 100
    y_res = 100
    z_res = 1
    # grid minimum boundaries (bound_factor*2*R in meters)
    x_b_min = -0.2
    y_b_min = -0.6
    z_b_min = (-(vol_thickness/2)/(2*R))+z_loc # choose volume that has a thickness of vol_thickness meter in the z-dimension (-0.5m under z_loc and +0.5m over z_loc)
    # grid maximum boundaries (bound_factor*2*R in meters)
    x_b_max = 10
    y_b_max = 0.6
    z_b_max = ((vol_thickness/2)/(2*R))+z_loc # choose volume that has a thickness of vol_thickness meter in the z-dimension (-0.5m under z_loc and +0.5m over z_loc)

    OwnFunctions.postprocess_fluiddomain(save_path, run_name, file_suffix, R, AOA, timesteps_to_evaluate;
                                         # ----- GRID OPTIONS ----------------  
                                         x_resolution    = x_res,                   # discretization of grid in x-direction
                                         y_resolution    = y_res,                   # discretization of grid in y-direction
                                         z_resolution    = z_res,                   # discretization of grid in z-direction
                                         x_bound_min     = x_b_min,                 # minimum bounds in x-direction (bound_factor*2*R in meters)
                                         y_bound_min     = y_b_min,                 # minimum bounds in y-direction (bound_factor*2*R in meters)
                                         z_bound_min     = z_b_min,                 # minimum bounds in z-direction (bound_factor*2*R in meters)
                                         x_bound_max     = x_b_max,                 # maximum bounds in x-direction (bound_factor*2*R in meters)
                                         y_bound_max     = y_b_max,                 # maximum bounds in y-direction (bound_factor*2*R in meters)
                                         z_bound_max     = z_b_max,                 # maximum bounds in z-direction (bound_factor*2*R in meters)
                                         prompt          = false,
                                         verbose         = false,
                                         debug           = false,#debug                    # saves fdom grid as a file
                                         first_fdom      = false                    # prevents question if fluiddomain postprocessing folder should be removed or not
                                         )


                                         
    # --------------- postprocess y-z-plane -------------------

    for x in 1:2:10 # x represents the stations a y-z-plane will be added
      x_loc = x   # x coordinate location of plane = x_loc*2*R in meters
      file_suffix = "_y_z_atx$(x_loc)D"
      push!(file_suffixes, file_suffix)
      vol_thickness = 0 # plane thickness in meters => set this value to 0 to get a 2D solution (plane)

      # grid resolution
      x_res = 1
      y_res = 100
      z_res = 100
      # grid minimum boundaries (bound_factor*2*R in meters)
      x_b_min = (-(vol_thickness/2)/(2*R))+x_loc # choose volume that has a thickness of vol_thickness meter in the x-dimension (-0.5m before x_loc and +0.5m behind x_loc)
      y_b_min = -0.6
      z_b_min = -0.6
      # grid maximum boundaries (bound_factor*2*R in meters)
      x_b_max = ((vol_thickness/2)/(2*R))+x_loc # choose volume that has a thickness of vol_thickness meter in the x-dimension (-0.5m before x_loc and +0.5m behind x_loc)
      y_b_max = 0.6
      z_b_max = 0.6

      OwnFunctions.postprocess_fluiddomain(save_path, run_name, file_suffix, R, AOA, timesteps_to_evaluate;
                                          # ----- GRID OPTIONS ----------------  
                                          x_resolution    = x_res,                   # discretization of grid in x-direction
                                          y_resolution    = y_res,                   # discretization of grid in y-direction
                                          z_resolution    = z_res,                   # discretization of grid in z-direction
                                          x_bound_min     = x_b_min,                 # minimum bounds in x-direction (bound_factor*2*R in meters)
                                          y_bound_min     = y_b_min,                 # minimum bounds in y-direction (bound_factor*2*R in meters)
                                          z_bound_min     = z_b_min,                 # minimum bounds in z-direction (bound_factor*2*R in meters)
                                          x_bound_max     = x_b_max,                 # maximum bounds in x-direction (bound_factor*2*R in meters)
                                          y_bound_max     = y_b_max,                 # maximum bounds in y-direction (bound_factor*2*R in meters)
                                          z_bound_max     = z_b_max,                 # maximum bounds in z-direction (bound_factor*2*R in meters)
                                          prompt          = false,
                                          verbose         = false,
                                          debug           = false,#debug                    # saves fdom grid as a file
                                          first_fdom      = false                    # prevents question if fluiddomain postprocessing folder should be removed or not
                                          )
    end


    # --------------- save the suffixes of the _fdom files that will be opened in paraview -------------------
    if length(timesteps_to_evaluate) == 1
      for i in 1:length(file_suffixes)
        file_suffixes[i] = file_suffixes[i]*"_fdom.$(timesteps_to_evaluate[1]).xmf;"
      end
    elseif length(timesteps_to_evaluate) > 1
      for i in 1:length(file_suffixes)
        file_suffixes[i] = file_suffixes[i]*"_fdom...xmf;"
      end
    end

  end

  return file_suffixes
end

