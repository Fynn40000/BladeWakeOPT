#=##############################################################################
# DESCRIPTION
    Postprocessing single turbine simulations
    => calculates the bladeloads (averaged over time and number of blades)
    => calculates the velocity field of the fluiddomain as:
        => x-y plane at hub height
        => y-z planes at several axial stations x behind the turbine

    => provides all necessary functions to postprocess the data of a single turbine simulation
=###############################################################################


function single_turbine_simulation_postprocessing(simulation_path::String, save_path_post::String, run_name::String, R::Float64, AOA::Float64; 
                                                  # ----- POSTPROCESSING EXECUTION ----------------
                                                  plot_bladeloads::Bool=true,     # postprocessing the bladeloads?
                                                  postprocess_fdom::Bool=true,    # postprocessing the fluiddomain?
                                                  # ----- SETTINGS FOR POSTPROCESSING -------------
                                                  Vinf = (X, t)->zeros(3),        # Freestream Velocity
                                                  magVinfx = 1,                   # Freestream Velocity in x direction
                                                  sim_time = 1,                   # Overall (real) simulation time in seconds
                                                  rev_to_average_idx=1,           # Revolution to wich the postprocessing should be applied on
                                                  nrevs_to_average=1,             # number of Revolutions to average for postprocessing the bladeloads
                                                  num_elements::Int64=1,          # number of blade elements per blade
                                                  tsteps = [1],                   # time steps to be postprocessed
                                                  debug::Bool=false,              # postprocess dimensionless coefficients too? => NOTE: debug statement must be set to true for uns.run_simulation. Otherwise the simulation files will not contain the coefficient data.
                                                  suppress_plots::Bool=true,      # suppresses the plots to show up on the display
                                                  gridsize_x_y::Float64=0.5,      # grid size of x-y fluid domain plane in meters
                                                  gridsize_y_z::Float64=0.5,      # grid size of y-z fluid domain plane in meters
                                                  cylindrical_grid = false,       # if true, the y-z plane will be calculated as a cylindrical grid and the wake velocity profiles will be saved within a .csv file
                                                                                  # this grid will be set automatically with the turbine diameter as its diameter
                                                  verbose = true
                                                  )

  
  # postprocess the blade loadings if desired
  if plot_bladeloads
    # postprocess the statistics first by applying function uns.postprocess_statistics(...)
    OwnFunctions.postprocess_statistics(simulation_path, run_name;
                                        rev_to_average_idx=rev_to_average_idx, 
                                        nrevs_to_average=nrevs_to_average)
    
    #statistics_save_path = joinpath(simulation_path, run_name)

    # create a folder to save the plots in
    if isdir(save_path_post)
      rm(save_path_post, recursive=true)
      mkdir(save_path_post)
    else
      mkdir(save_path_post)
    end

    # postprocess statistics of the "_vlm" statistics file regarding blade 1
    OwnFunctions.postprocess_blade_loading(simulation_path, save_path_post, run_name, R; 
                                    file_marker="_vlm",
                                    rev_to_average_idx=rev_to_average_idx, 
                                    nrevs_to_average=nrevs_to_average,
                                    num_elements=num_elements,
                                    suppress_plots=suppress_plots)
    # postprocess statistics of the "_loft" statistics file regarding blade 1
    OwnFunctions.postprocess_blade_loading(simulation_path, save_path_post, run_name, R; 
                                    file_marker="_loft",
                                    rev_to_average_idx=rev_to_average_idx, 
                                    nrevs_to_average=nrevs_to_average,
                                    num_elements=num_elements,
                                    debug=debug,
                                    suppress_plots=suppress_plots)

    #println("\nBlade load plots saved under $(save_path_post)\n")
  end


  file_suffixes = []
  if postprocess_fdom

    # read simulation data
    #simdata = CSV.read(joinpath(simulation_path, run_name*"_convergence.csv"), DataFrame)
    # Calculate nsteps_per_rev
    #nsteps_per_rev[run_name] = ceil(Int, 360 / (simdata[2, 1] - simdata[1, 1]))
    
    # calculate number of planes to postprocess 
    approx_x_range = magVinfx * sim_time
    x_max_calc = approx_x_range / (2*R)       # calculated maximum x boundary based on wind velocity and simulation time

    # set time steps to be postprocessed via the postprocess_fluiddomain function
    timesteps_to_evaluate = tsteps



    # --------------- postprocess x-y-plane -------------------

    z_loc = 0   # z coordinate location of plane = z_loc*2*R in meters
    file_suffix = "_x_y_atz$(z_loc)D"
    push!(file_suffixes, file_suffix)
    vol_thickness = 0                         # plane thickness in meters => set this value to 0 to get a 2D solution (plane)
    gridsize_m = gridsize_x_y                 # grid size in meters => 0.5 equals 0.5m

    # grid resolution
    x_res = R*(1/gridsize_m)
    y_res = R*(1/gridsize_m)
    z_res = 1
    # grid minimum boundaries (bound_factor*2*R in meters)
    x_b_min = -0.2
    y_b_min = -0.7
    z_b_min = (-(vol_thickness/2)/(2*R))+z_loc # choose volume that has a thickness of vol_thickness meter in the z-dimension (-0.5m under z_loc and +0.5m over z_loc)
    # grid maximum boundaries (bound_factor*2*R in meters)
    x_b_max = x_max_calc#15
    y_b_max = 0.7
    z_b_max = ((vol_thickness/2)/(2*R))+z_loc # choose volume that has a thickness of vol_thickness meter in the z-dimension (-0.5m under z_loc and +0.5m over z_loc)

    OwnFunctions.postprocess_fluiddomain(simulation_path, run_name, file_suffix, R, AOA, timesteps_to_evaluate;
                                         # ----- FREESTREAM VELOCITY ---------
                                         Vinf            = Vinf,
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
                                         verbose         = verbose,
                                         debug           = false,#debug                    # saves fdom grid as a file
                                         first_fdom      = false                    # prevents question if fluiddomain postprocessing folder should be removed or not
                                         )



    # --------------- postprocess y-z-plane -------------------

    x_max_calc = floor(x_max_calc)
    for x in 1:1:x_max_calc # x represents the x-coordinate stations a y-z-plane will be calculated for
      x_loc = x   # x coordinate location of plane = x_loc*2*R in meters
      file_suffix = "_y_z_atx$(x_loc)D"
      push!(file_suffixes, file_suffix)
      vol_thickness = 0                         # plane thickness in meters => set this value to 0 to get a 2D solution (plane)
      gridsize_m = gridsize_y_z                 # grid size in meters => 0.5 equals 0.5m

      
      if cylindrical_grid
        # NOTE: The plane is created in the x_y plane first, then rotated into the y-z plane
              # => therefore the x_loc is used to define the z boundaries (location in the z plane) and then rotated to -90 deg around the y axis
        
        # grid resolution
        x_res = R*(1/gridsize_m)
        y_res = 49                                                # Angular Step = 360/y_res in (°) !!!For some reason an error occurs when y_res >= 50!!! Nice would be: y_res = 2*pi*R*(1/gridsize_m) to have the minimum gridsize at the outer cylindrical plane section
        z_res = R*(1/gridsize_m)                                  # THIS WILL DO NOTHING IF vol_thickness == 0 (meaning z_b_min==z_b_max)
        # grid minimum boundaries (bound_factor*2*R in meters)
        x_b_min = 0.0
        #y_b_min = 0.0
        z_b_min = (-(vol_thickness/2)/(2*R))+x_loc
        # grid maximum boundaries (bound_factor*2*R in meters)
        x_b_max = 1.0
        #y_b_max = 2*pi
        z_b_max = ((vol_thickness/2)/(2*R))+x_loc

      else
        # grid resolution
        x_res = 1
        y_res = R*(1/gridsize_m)
        z_res = R*(1/gridsize_m)
        # grid minimum boundaries (bound_factor*2*R in meters)
        x_b_min = (-(vol_thickness/2)/(2*R))+x_loc # choose volume that has a thickness of vol_thickness meter in the x-dimension (-0.5m before x_loc and +0.5m behind x_loc)
        y_b_min = -0.7
        z_b_min = -0.7
        # grid maximum boundaries (bound_factor*2*R in meters)
        x_b_max = ((vol_thickness/2)/(2*R))+x_loc # choose volume that has a thickness of vol_thickness meter in the x-dimension (-0.5m before x_loc and +0.5m behind x_loc)
        y_b_max = 0.7
        z_b_max = 0.7

      end



      OwnFunctions.postprocess_fluiddomain(simulation_path, run_name, file_suffix, R, AOA, timesteps_to_evaluate;
                                          # ----- FREESTREAM VELOCITY ---------
                                          Vinf            = Vinf,
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
                                          cylindrical_grid=cylindrical_grid,         # a cylindrical grid will be calculated
                                          prompt          = false,
                                          verbose         = verbose,
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




"
Applies function uns.postprocess_statistics
"
function postprocess_statistics(save_path::String, run_name::String;
                                rev_to_average_idx=1, nrevs_to_average=1)

  rotor_axis = [-1.0, 0.0, 0.0]       # Rotor centerline axis

  nsteps_per_rev = Dict()

  # read simulation data
  simdata = CSV.read(joinpath(save_path, run_name*"_convergence.csv"), DataFrame)

  # Calculate nsteps_per_rev
  nsteps_per_rev[run_name] = ceil(Int, 360 / (simdata[2, 1] - simdata[1, 1]))

  # Generate statistics (mean and deviation of loading)
  r_path = save_path
  s_path = joinpath(r_path, run_name*"-statistics")

  nums      = range(rev_to_average_idx-nrevs_to_average, rev_to_average_idx; length = nsteps_per_rev[run_name]+1)# Process outputs between revolutions 3 and 4!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nums      = ceil.(Int, nums * nsteps_per_rev[run_name])

  uns.postprocess_statistics(r_path, s_path, nums;
                             cyl_axial_dir = rotor_axis,
                             prompt = false,
                             verbose=false)

end


"
Postprocesses and plots the blade loads of a single turbines simulation.
"
function postprocess_blade_loading(save_path::String, save_path_post::String, run_name::String, R::Float64; 
                            file_marker::String="_vlm", rev_to_average_idx=1, 
                            nrevs_to_average=1, 
                            num_elements::Int64=1, 
                            debug::Bool=false,
                            suppress_plots::Bool=true # suppresses the plots to show up on the display
                            )
  
  # Plot styles
  stl = "-" # style
  clr = "dodgerblue" # color
  alpha = 1.0 # alpha
  lbl = run_name # label

  rotor_axis = [-1.0, 0.0, 0.0]       # Rotor centerline axis

  # Read and plot FLOWUnsteady mean blade loading and deviation
  r_path = joinpath(save_path, run_name*"-statistics")

  if file_marker == "_vlm"
    (rs, Gamma, Np, Tp) = uns.postprocess_bladeloading_turbine(r_path;
                                                                O           = zeros(3),
                                                                rotor_axis  = rotor_axis,
                                                                filename    = run_name*"_Rotor_Blade1"*file_marker*"-statistics.vtk",
                                                                fieldsuff   = "-mean",
                                                                num_elements= num_elements
                                                                )
  elseif file_marker == "_loft"
    (roR1, Gamma1, Np1, Tp1, L1, D1,
     cn1, ct1, cl1, cd1,
     _, _,#a, a_tangential,
     aoa1, twist1, flowangle1, 
     Vx1, Vy1, w1,
     F1, loc_solidity1) =        uns.postprocess_bladeloading_turbine(r_path;
                                                              O           = zeros(3),
                                                              rotor_axis  = rotor_axis,
                                                              filename    = run_name*"_Rotor_Blade1"*file_marker*"-statistics.vtk",
                                                              fieldsuff   = "-mean",
                                                              num_elements= num_elements,
                                                              debug=debug
                                                              )

     (roR2, Gamma2, Np2, Tp2, L2, D2,
     cn2, ct2, cl2, cd2,
     _, _,#a, a_tangential,
     aoa2, twist2, flowangle2, 
     Vx2, Vy2, w2,
     F2, loc_solidity2) =        uns.postprocess_bladeloading_turbine(r_path;
                                                              O           = zeros(3),
                                                              rotor_axis  = rotor_axis,
                                                              filename    = run_name*"_Rotor_Blade2"*file_marker*"-statistics.vtk",
                                                              fieldsuff   = "-mean",
                                                              num_elements= num_elements,
                                                              debug=debug
                                                              )

     (roR3, Gamma3, Np3, Tp3, L3, D3,
     cn3, ct3, cl3, cd3,
     _, _,#a, a_tangential,
     aoa3, twist3, flowangle3, 
     Vx3, Vy3, w3,
     F3, loc_solidity3) =        uns.postprocess_bladeloading_turbine(r_path;
                                                              O           = zeros(3),
                                                              rotor_axis  = rotor_axis,
                                                              filename    = run_name*"_Rotor_Blade3"*file_marker*"-statistics.vtk",
                                                              fieldsuff   = "-mean",
                                                              num_elements= num_elements,
                                                              debug=debug
                                                              )

    roR = ((roR1 .+ roR2 .+ roR3) ./ 3)
    rs = roR .* R
    Gamma = ((Gamma1 .+ Gamma2 .+ Gamma3) ./ 3)
    Np = ((Np1 .+ Np2 .+ Np3) ./ 3)
    Tp = ((Tp1 .+ Tp2 .+ Tp3) ./ 3)
    L = ((L1 .+ L2 .+ L3) ./ 3)
    D = ((D1 .+ D2 .+ D3) ./ 3)
    cn_oneblade = ((cn1 .+ cn2 .+ cn3) ./ 3)
    ct_oneblade = ((ct1 .+ ct2 .+ ct3) ./ 3)
    cn = (cn1 .+ cn2 .+ cn3)
    ct = (ct1 .+ ct2 .+ ct3)
    cl = ((cl1 .+ cl2 .+ cl3) ./ 3)
    cd = ((cd1 .+ cd2 .+ cd3) ./ 3)
    aoa = ((aoa1 .+ aoa2 .+ aoa3) ./ 3)
    twist = ((twist1 .+ twist2 .+ twist3) ./ 3)
    flowangle = ((flowangle1 .+ flowangle2 .+ flowangle3) ./ 3)
    Vx = ((Vx1 .+ Vx2 .+ Vx3) ./ 3)
    Vy = ((Vy1 .+ Vy2 .+ Vy3) ./ 3)
    w = ((w1 .+ w2 .+ w3) ./ 3)
    F = ((F1 .+ F2 .+ F3) ./ 3)
    loc_solidity = ((loc_solidity1 .+ loc_solidity2 .+ loc_solidity3) ./ 3)

    a = 1 ./ (1 .+ ((4 .* sin.(flowangle) .* sin.(flowangle)) ./ (loc_solidity .* cl .* cos.(flowangle))) )
    a_tangential = 1 ./ (((4 .* cos.(flowangle)) ./ (loc_solidity .* cl)) .- 1)
  end

  #--------------------------------------------
  # Start plotting - Gamma
  xlabel=L"Radial position $r/R$" 
  ylabel=L"Circulation $\Gamma$ (m$^2$/s)"
  linelabel=L"$\Gamma$ distribution"
  filename="Gamma_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx)"*file_marker
  plot_var(rs/R, Gamma, xlabel, ylabel, linelabel,save_path_post, filename;
           xlims=[0,1],
           dx=0.1,
           suppress_plots=suppress_plots
           )

  #--------------------------------------------
  # Start plotting - Np
  xlabel=L"Radial position $r/R$" 
  ylabel=L"Thrust force ($\mathrm{N/m}$)"
  linelabel=L"$F_{T}$ distribution"
  filename="Np_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx)"*file_marker
  plot_var(rs/R, Np, xlabel, ylabel, linelabel,save_path_post, filename;
           xlims=[0,1],
           dx=0.1,
           suppress_plots=suppress_plots
           )

  #--------------------------------------------
  # Start plotting - Tp
  xlabel=L"Radial position $r/R$" 
  ylabel=L"Tangential force ($\mathrm{N/m}$)"
  linelabel=L"$F_{t}$ distribution"
  filename="Tp_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx)"*file_marker
  plot_var(rs/R, Tp, xlabel, ylabel, linelabel,save_path_post, filename;
           xlims=[0,1],
           dx=0.1,
           suppress_plots=suppress_plots
           )

  if debug

    #--------------------------------------------
    # Start plotting - L
    xlabel=L"Radial position $r/R$" 
    ylabel=L"Lift $F_{L}$ ($\mathrm{N/m}$)"
    linelabel=L"$F_{L}$ distribution"
    filename="L_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx)"*file_marker
    plot_var(rs/R, L, xlabel, ylabel, linelabel,save_path_post, filename;
            xlims=[0,1],
            dx=0.1,
            suppress_plots=suppress_plots
            )

    #--------------------------------------------
    # Start plotting - D
    xlabel=L"Radial position $r/R$" 
    ylabel=L"Drag $F_{D}$ ($\mathrm{N/m}$)"
    linelabel=L"$F_{D}$ distribution"
    filename="D_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx)"*file_marker
    plot_var(rs/R, D, xlabel, ylabel, linelabel,save_path_post, filename;
            xlims=[0,1],
            dx=0.1,
            suppress_plots=suppress_plots
            )

    #--------------------------------------------
    # Start plotting - cn_oneblade
    xlabel=L"Radial position $r/R$" 
    ylabel=L"Thrust coefficient $C_{T}$ (-)"
    linelabel=L"$C_{T}$ distribution"
    filename="cn_oneblade_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx)"*file_marker
    plot_var(rs/R, cn_oneblade, xlabel, ylabel, linelabel,save_path_post, filename;
            xlims=[0,1],
            dx=0.1,
            suppress_plots=suppress_plots
            )

    #--------------------------------------------
    # Start plotting - ct_oneblade
    xlabel=L"Radial position $r/R$" 
    ylabel=L"Tangential force coefficient $C_{t}$ (-)"
    linelabel=L"$C_{t}$ distribution"
    filename="ct_oneblade_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx)"*file_marker
    plot_var(rs/R, ct_oneblade, xlabel, ylabel, linelabel,save_path_post, filename;
            xlims=[0,1],
            dx=0.1,
            suppress_plots=suppress_plots
            )

    #--------------------------------------------
    # Start plotting - cn
    xlabel=L"Radial position $r/R$" 
    ylabel=L"Thrust coefficient $C_{T}$ (-)"
    linelabel=L"$C_{T}$ distribution"
    filename="cn_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx)"*file_marker
    plot_var(rs/R, cn, xlabel, ylabel, linelabel,save_path_post, filename;
            xlims=[0,1],
            dx=0.1,
            suppress_plots=suppress_plots
            )

    #--------------------------------------------
    # Start plotting - ct
    xlabel=L"Radial position $r/R$" 
    ylabel=L"Tangential force coefficient $C_{t}$ (-)"
    linelabel=L"$C_{t}$ distribution"
    filename="ct_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx)"*file_marker
    plot_var(rs/R, ct, xlabel, ylabel, linelabel,save_path_post, filename;
            xlims=[0,1],
            dx=0.1,
            suppress_plots=suppress_plots
            )

    #--------------------------------------------
    # Start plotting - cl
    xlabel=L"Radial position $r/R$" 
    ylabel=L"Lift coefficient $C_{L}$ (-)"
    linelabel=L"$C_{L}$ distribution"
    filename="cl_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx)"*file_marker
    plot_var(rs/R, cl, xlabel, ylabel, linelabel,save_path_post, filename;
            xlims=[0,1],
            dx=0.1,
            suppress_plots=suppress_plots
            )

    #--------------------------------------------
    # Start plotting - cd
    xlabel=L"Radial position $r/R$" 
    ylabel=L"Drag coefficient $C_{D}$ (-)"
    linelabel=L"$C_{D}$ distribution"
    filename="cd_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx)"*file_marker
    plot_var(rs/R, cd, xlabel, ylabel, linelabel,save_path_post, filename;
            xlims=[0,1],
            dx=0.1,
            suppress_plots=suppress_plots
            )

    #--------------------------------------------
    # Start plotting - a
    xlabel=L"Radial position $r/R$" 
    ylabel=L"Axial induction factor $a$ (-)"
    linelabel=L"$a$ distribution"
    filename="a_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx)"*file_marker
    plot_var(rs/R, a, xlabel, ylabel, linelabel,save_path_post, filename;
            xlims=[0,1],
            dx=0.1,
            suppress_plots=suppress_plots
            )

    #--------------------------------------------
    # Start plotting - a_tangential
    xlabel=L"Radial position $r/R$" 
    ylabel=L"Tangential induction factor $a'$ (-)"
    linelabel=L"$a'$ distribution"
    filename="a_tangential_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx)"*file_marker
    plot_var(rs/R, a_tangential, xlabel, ylabel, linelabel,save_path_post, filename;
            xlims=[0,1],
            dx=0.1,
            suppress_plots=suppress_plots
            )
    
    #--------------------------------------------
    # Start plotting - aoa
    xlabel=L"Radial position $r/R$" 
    ylabel=L"angle of attack $\alpha$ (°)"
    linelabel=L"$\alpha$ distribution"
    filename="aoa_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx)"*file_marker
    plot_var(rs/R, aoa, xlabel, ylabel, linelabel,save_path_post, filename;
            xlims=[0,1],
            dx=0.1,
            suppress_plots=suppress_plots
            )

    #--------------------------------------------
    # Start plotting - twist
    xlabel=L"Radial position $r/R$" 
    ylabel=L"twist angle $\Theta$ (°)"
    linelabel=L"$\Theta$ distribution"
    filename="twist_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx)"*file_marker
    plot_var(rs/R, twist, xlabel, ylabel, linelabel,save_path_post, filename;
            xlims=[0,1],
            dx=0.1,
            suppress_plots=suppress_plots
            )

    #--------------------------------------------
    # Start plotting - flowangle
    xlabel=L"Radial position $r/R$" 
    ylabel=L"flow angle $\beta$ (°)"
    linelabel=L"$\beta$ distribution"
    filename="flowangle_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx)"*file_marker
    plot_var(rs/R, flowangle, xlabel, ylabel, linelabel,save_path_post, filename;
            xlims=[0,1],
            dx=0.1,
            suppress_plots=suppress_plots
            )

    #--------------------------------------------
    # Start plotting - Vx
    xlabel=L"Radial position $r/R$" 
    ylabel=L"local flow velocity x-direction $V_{x}$ ($\mathrm{m/s}$)"
    linelabel=L"$V_{x}$ distribution"
    filename="Vx_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx)"*file_marker
    plot_var(rs/R, Vx, xlabel, ylabel, linelabel,save_path_post, filename;
            xlims=[0,1],
            dx=0.1,
            suppress_plots=suppress_plots
            )

    #--------------------------------------------
    # Start plotting - Vy
    xlabel=L"Radial position $r/R$" 
    ylabel=L"local flow velocity y-direction $V_{y}$ ($\mathrm{m/s}$)"
    linelabel=L"$V_{y}$ distribution"
    filename="Vy_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx)"*file_marker
    plot_var(rs/R, Vy, xlabel, ylabel, linelabel,save_path_post, filename;
            xlims=[0,1],
            dx=0.1,
            suppress_plots=suppress_plots
            )

    #--------------------------------------------
    # Start plotting - w
    xlabel=L"Radial position $r/R$" 
    ylabel=L"local relative flow velocity $w$ ($\mathrm{m/s}$)"
    linelabel=L"$w$ distribution"
    filename="w_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx)"*file_marker
    plot_var(rs/R, w, xlabel, ylabel, linelabel,save_path_post, filename;
            xlims=[0,1],
            dx=0.1,
            suppress_plots=suppress_plots
            )

    #--------------------------------------------
    # Start plotting - F
    xlabel=L"Radial position $r/R$" 
    ylabel=L"Prandtl loss factor (-)"
    linelabel=L"$F$ distribution"
    filename="F_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx)"*file_marker
    plot_var(rs/R, F, xlabel, ylabel, linelabel,save_path_post, filename;
            xlims=[0,1],
            dx=0.1,
            suppress_plots=suppress_plots
            )

    #--------------------------------------------
    # Start plotting - loc_solidity
    xlabel=L"Radial position $r/R$" 
    ylabel=L"local solidity (-)"
    linelabel=L"$\sigma'$ distribution"
    filename="loc_solidity_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx)"*file_marker
    plot_var(rs/R, loc_solidity, xlabel, ylabel, linelabel,save_path_post, filename;
            xlims=[0,1],
            dx=0.1,
            suppress_plots=suppress_plots
            )

  end

end


"
Creates and saves a plot of a given x and y data row.
"
function plot_var(x_data, y_data, xlabel, ylabel, linelabel,
                  save_path_post::String, filename::String;
                  stl::String="-",              # linestyle
                  alpha::Float64=1.0,           # opacity (1.0 -> solid)
                  color::String="dodgerblue",   # line color
                  linewidth::Float64=2.0,       # line thickness
                  xlims =[0, 0],                # x axis limits => if not set by user, the axis is set automatically
                  dx=0.1,                       # x axis step => if not set by user, the axis is set automatically
                  ylims =[0,0],                 # y axis limits => if not set by user, the axis is set automatically
                  dy=5,                         # y axis step => if not set by user, the axis is set automatically
                  legend_loc="center left",
                  legend_bbox_to_anchor=(1, 0.5),
                  legend_frameon=false,
                  legend_fontsize=10,
                  plotYzero::Bool=true,         # plot dashed black line at y=0
                  filetype::String=".pdf",      # plot file type the plot is saved as
                  x_name::String="x_data",      # name of the x data
                  y_name::String="y_data",      # name of the y data
                  suppress_plots::Bool=true,    # suppresses the plots to show up on the display
                  save_csv::Bool=true           # if true, a .csv file will be created containing the x-axis and y-axis values
                  )

  fig = plt.figure(figsize=[7*1.5, 5] * 2/3)
  ax = fig.gca()

  # Plot y=0
  if plotYzero
    ax.plot([0, 1], zeros(2), ":k", linewidth=1)
  end

  # Plot var
  ax.plot(x_data, y_data, stl; alpha=alpha, label=linelabel, color=color, linewidth=linewidth)
  
  # Set axes
  if xlims != [0,0]
    ax.set_xlim(xlims)
    ax.set_xticks(xlims[1]:dx:xlims[2])
  end
  if ylims != [0,0]
    ax.set_ylim(ylims)
    ax.set_yticks(ylims[1]:dy:ylims[2])
  end

  # Set axes labels
  ax.set_xlabel(xlabel)
  ax.set_ylabel(ylabel)

  # Set legend
  ax.legend(loc=legend_loc, bbox_to_anchor=legend_bbox_to_anchor, frameon=legend_frameon, fontsize=legend_fontsize)

  ax.spines["right"].set_visible(false)
  ax.spines["top"].set_visible(false)

  fig.tight_layout()

  # save the plot
  fig.savefig(joinpath(save_path_post, filename*filetype))
  if suppress_plots
    close(fig)
  end

  if save_csv
    _create_csv(x_data, y_data, save_path_post, filename;
                x_name=x_name,
                y_name=y_name
                )
  end

end


"
Postprocesses the fluiddomain of a single turbines simulation.
"

function postprocess_fluiddomain(# ---- ESSENTIAL ARGUMENTS ---------
                                simulation_path::String,                      # path to folder where simulation is stored in
                                run_name::String,                       # name of simulation
                                file_suffix::String,                    # postprocessed fluid domain file suffix
                                R::Float64,                             # Rotor tip radius
                                AOA::Float64,                           # Angle of attack
                                nums::Vector{Int64};                    # Time steps to process
                                # ----- OPTIONAL ARGUMENTS ----------
                                save_path       = nothing,              # folder to store data in (if nothing, data will be stored under simulation_path)
                                # ----- FREESTREAM VELOCITY ---------
                                Vinf            = (X, t)->zeros(3),
                                # ----- GRID OPTIONS ----------------  
                                x_resolution    = 50,                   # discretization of grid in x-direction
                                y_resolution    = 50,                   # discretization of grid in y-direction
                                z_resolution    = 50,                   # discretization of grid in z-direction
                                x_bound_min     = -1,                   # minimum bounds in x-direction (bound_factor*2*R in meters)
                                y_bound_min     = -0.5,                 # minimum bounds in y-direction (bound_factor*2*R in meters)
                                z_bound_min     = -0.5,                 # minimum bounds in z-direction (bound_factor*2*R in meters)
                                x_bound_max     = 1,                    # maximum bounds in x-direction (bound_factor*2*R in meters)
                                y_bound_max     = 0.5,                  # maximum bounds in y-direction (bound_factor*2*R in meters)
                                z_bound_max     = 0.5,                  # maximum bounds in z-direction (bound_factor*2*R in meters)
                                cylindrical_grid= false,                # is the grid a cylindrical grid?
                                # ----- VPM OPTIONS ----------------    
                                maxsigma        = nothing,              # Particles larger than this get shrunk to this size (this helps speed up computation)
                                maxmagGamma     = Inf,                  # Any vortex strengths larger than this get clipped to this value
                                include_staticparticles = true,         # Whether to include the static particles embedded in the solid surfaces
                                # ----- OUTPUT OPTIONS ---------------- 
                                prompt          = true,                 # Whether to prompt the user
                                verbose         = true,                 # Enable verbose
                                v_lvl           = 0,                    # Verbose indentation level
                                debug           = false,                # saves fdom grid as a file
                                # ----- OTHER ------------------------- 
                                first_fdom      = true                  # if true, creates folder to store the postprocessed fluiddomains in
                                )


  pfield_prefix   = run_name*"_pfield"               # Prefix of particle field files to read
  staticpfield_prefix = run_name*"_staticpfield"     # Prefix of static particle field files to read

  r_path = simulation_path                           # path to read data from (path to folder all the simulation data is stored in)
  if save_path !== nothing                           # path to folder all postprocessed fluiddomain files will be stored in
    s_path = save_path
  else
    s_path = simulation_path                         
  end

  output_prefix   = run_name                # Prefix of output files

  # Grid
  L               = R                                               # (m) reference length = rotor tip radius
  

  if cylindrical_grid
    spacetransform  = uns.gt.cylindrical3D
    O               = zeros(3)                                        # translation and re-orientation to the given origin and orientation  
    Oaxis           = uns.gt.rotation_matrix2(0, -90, AOA)            # Orientation of grid (flip the grid aroud and locate it behind the turbine in x direction)

    Pmin            = [L*x_bound_min, 0.0, L*z_bound_min]
    Pmax            = [L*x_bound_max, 2*pi, L*z_bound_max]
    dradial_annulus, dpolar_angle, dplane = L/x_resolution, y_resolution, L/z_resolution
    NDIVS           = ceil.(Int, [(Pmax[1] - Pmin[1])/dradial_annulus, dpolar_angle, (Pmax[3] - Pmin[3])/dplane])

  else
    spacetransform  = nothing
    O               = zeros(3)                                        # translation and re-orientation to the given origin and orientation  
    Oaxis           = uns.gt.rotation_matrix2(0, 0, AOA)              # Orientation of grid

    Pmin            = 2*L*[x_bound_min, y_bound_min, z_bound_min]     # (m) minimum bounds (bound_factor * rotor diameter)
    Pmax            = 2*L*[x_bound_max, y_bound_max, z_bound_max]     # (m) maximum bounds (bound_factor * rotor diameter)
    dx, dy, dz      = L/x_resolution, L/y_resolution, L/z_resolution  # (m) cell size in each direction
    NDIVS           = ceil.(Int, (Pmax .- Pmin)./[dx, dy, dz])        # Number of cells in each direction

  end
  
  nnodes          = prod(NDIVS .+ 1)                                # Total number of nodes

  
  

  # VPM settings
  maxparticles    = Int(1.0e6 + nnodes)                             # Maximum number of particles
  fmm             = uns.vpm.FMM(; p=4, ncrit=50, theta=0.4, phi=0.3)# FMM parameters (decrease phi to reduce FMM noise)
  scale_sigma     = 1.00                                            # Shrink smoothing radii by this factor
  f_sigma         = 0.5                                             # Smoothing of node particles as sigma = f_sigma*meansigma

  if maxsigma === nothing
    maxsigma        = L/10                                          # Particles larger than this get shrunk to this size (this helps speed up computation)
  end

  other_file_prefs = include_staticparticles ? [staticpfield_prefix] : []
  other_read_paths = [r_path for i in 1:length(other_file_prefs)]

  if verbose
    println("\t"^(v_lvl)*"Fluid domain grid")
    println("\t"^(v_lvl)*"NDIVS =\t$(NDIVS)")
    println("\t"^(v_lvl)*"Number of nodes =\t$(nnodes)")
  end

  # --------------- PROCESSING SETUP ---------------------------------------------
  if verbose
    println("\t"^(v_lvl)*"Getting ready to process $(r_path)")
    println("\t"^(v_lvl)*"Results will be saved under $(s_path)")
  end

  # Create save path
  if first_fdom
    if s_path != r_path
      uns.gt.create_path(s_path, prompt)
    end
  end

  # Copy this driver file
  #cp(@__FILE__, joinpath(s_path, splitdir(@__FILE__)[2]); force=true)

  # Generate function to process the field clipping particle sizes
  preprocessing_pfield = uns.generate_preprocessing_fluiddomain_pfield(maxsigma, maxmagGamma;
                                                                       verbose=verbose, v_lvl=v_lvl+1)
                                                                    
  # --------------- PROCESS SIMULATION -------------------------------------------

  nthreads        = 1                         # Total number of threads
  nthread         = 1                         # Number of this thread
  dnum = floor(Int, length(nums)/nthreads)    # Number of time steps per thread
  threaded_nums = [view(nums, dnum*i+1:(i<nthreads-1 ? dnum*(i+1) : length(nums))) for i in 0:nthreads-1]

  for these_nums in threaded_nums[nthread:nthread]
    Xdummy = zeros(3)
    U = t->Vinf(Xdummy, t)
    #uns.computefluiddomain(Pmin, Pmax, NDIVS,
    computefluiddomain(Pmin, Pmax, NDIVS,
                           maxparticles,
                           these_nums, r_path, pfield_prefix;
                           spacetransform=spacetransform,
                           O=O,
                           Oaxis=Oaxis,
                           fmm=fmm,
                           pfield_optargs=[:Uinf => U],
                           f_sigma=f_sigma,
                           save_path=s_path,
                           file_pref=output_prefix, grid_names=[file_suffix*"_fdom"],
                           other_file_prefs=other_file_prefs,
                           other_read_paths=other_read_paths,
                           userfunction_pfield=preprocessing_pfield,
                           add_J=true, add_Uinf=true,
                           verbose=verbose, v_lvl=v_lvl, debug=debug)

  end

end



"
Creates and saves a csv file of a given x and y data row.
"
function _create_csv(x_data, y_data, save_path_post::String, filename::String;
                     x_name::String="x_data",      # name of the x data
                     y_name::String="y_data",      # name of the y data
                     )

  fname = joinpath(save_path_post, filename*".csv")             # file path to store data in

  #f = open(fname, "w")

  data_matrix = hcat(x_data, y_data)                            # matrix with all data
  #data_table = Tables.table(data_matrix)                        # matrix converted into table
  df = DataFrame(data_matrix, :auto)                            # matrix converted into dataframe

  CSV.write(fname, df, header=[x_name, y_name])                 # Create and save csv file

end
