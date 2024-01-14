#=##############################################################################
# DESCRIPTION
    Postprocessing single turbine simulations
    Plots created for:
        => Blade loads
        => 
=###############################################################################

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
Plots the blade loads of a single turbines simulation.
"
function plot_blade_loading(save_path::String, save_path_post::String, run_name::String, R::Float64; 
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
    (roR, Gamma, Np, Tp, L, D,
     cn, ct, cl, cd,
     a, a_tangential, 
     aoa, twist, flowangle, 
     Vx, Vy, w) =        uns.postprocess_bladeloading_turbine(r_path;
                                                              O           = zeros(3),
                                                              rotor_axis  = rotor_axis,
                                                              filename    = run_name*"_Rotor_Blade1"*file_marker*"-statistics.vtk",
                                                              fieldsuff   = "-mean",
                                                              num_elements= num_elements,
                                                              debug=debug
                                                              )
    rs = roR .* R

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
  ylabel=L"Loading ($\mathrm{N/m}$)"
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
  ylabel=L"Loading ($\mathrm{N/m}$)"
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
    ylabel=L"Lift $F_{L}$ ($\mathrm{N/m}$))"
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
    ylabel=L"Drag $F_{D}$ ($\mathrm{N/m}$))"
    linelabel=L"$F_{D}$ distribution"
    filename="D_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx)"*file_marker
    plot_var(rs/R, D, xlabel, ylabel, linelabel,save_path_post, filename;
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
    ylabel=L"Thangential force coefficient $C_{t}$ (-)"
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
    ylabel=L"Induction factor $a$ (-)"
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
