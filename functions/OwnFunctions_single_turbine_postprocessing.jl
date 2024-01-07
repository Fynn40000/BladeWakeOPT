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
    (roR, Gamma, Np, Tp, 
     cn, ct, cl, cd,
     a, a_tangential) = uns.postprocess_bladeloading_turbine(r_path;
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

  end

end

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
                  plotYzero::Bool=true,          # plot dashed black line at y=0
                  filetype::String=".pdf",       # plot file type the plot is saved as
                  suppress_plots::Bool=true      # suppresses the plots to show up on the display
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

end
