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
  s_path = r_path*"-statistics"

  nums      = range(rev_to_average_idx-nrevs_to_average, rev_to_average_idx; length = nsteps_per_rev[run_name]+1)# Process outputs between revolutions 3 and 4!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nums      = ceil.(Int, nums * nsteps_per_rev[run_name])

  uns.postprocess_statistics(r_path, s_path, nums;
                             cyl_axial_dir = rotor_axis,
                             prompt = false)

end


"
Plots the blade loads of a single turbines simulation.
"

function plot_blade_loading(save_path::String, save_path_post::String, run_name::String, R::Float64; 
                            file_marker::String="_vlm", rev_to_average_idx=1, nrevs_to_average=1, num_elements::Int64=1, debug::Bool=false)
  
  # Plot styles
  stl = "-" # style
  clr = "dodgerblue" # color
  alpha = 1.0 # alpha
  lbl = run_name # label

  rotor_axis = [-1.0, 0.0, 0.0]       # Rotor centerline axis

  # Read and plot FLOWUnsteady mean blade loading and deviation
  r_path = joinpath(save_path, "..", run_name*"-statistics")

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
    #println("CALCULATED DATA OF loft FILE")
    #println("length of rs")
    #println(length(rs))
    #println("length of cn")
    #println(length(cn))
    #println(" ")
    #println("rs")
    #println(rs)
    #println(" ")
    #println("Gamma")
    #println(Gamma)
    #println("Np")
    #println(Np)
    #println("Tp")
    #println(Tp)
    #println("cn")
    #println(cn)
    #println("ct")
    #println(ct)
    #println("cl")
    #println(cl)
    #println("cd")
    #println(cd)
    #println("a")
    #println(a)
    #println("a_tangential")
    #println(a_tangential)
  end

  #--------------------------------------------
  # Start plotting - Np
  fig = plt.figure(figsize=[7*1.5, 5] * 2/3)
  ax = fig.gca()
  
  # Plot y=0
  ax.plot([0, 1], zeros(2), ":k", linewidth=1)

  # Plot Np
  ax.plot(rs/R, Np, stl; alpha=alpha, label=lbl*" - Np", color=clr, linewidth=2.0)#"-", "dodgerblue", 1.0, "rVPM - high fidelity"
  
  # Format plot
  xlims, dx = [[0, 1], 0.2]
  ylims, dy = [[-100, 100], 5]

  ax.set_xlim(xlims)
  ax.set_xticks(xlims[1]:dx:xlims[2])
  ax.set_xlabel(L"Radial position $r/R$")

  #ax.set_ylim(ylims)
  #ax.set_yticks(0:dy:ylims[2])
  ax.set_ylabel(L"Loading ($\mathrm{N/m}$)")

  ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), frameon=false, fontsize=10)

  ax.spines["right"].set_visible(false)
  ax.spines["top"].set_visible(false)

  fig.tight_layout()

  # Save plot
  filename = "Np_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx).pdf"
  save_plot(save_path_post, filename; fig=fig)


  #--------------------------------------------
  # Start plotting - Tp
  fig2 = plt.figure(figsize=[7*1.5, 5] * 2/3)
  ax2 = fig2.gca()
  
  # Plot y=0
  ax2.plot([0, 1], zeros(2), ":k", linewidth=1)

  # Plot Tp
  ax2.plot(rs/R, Tp, stl; alpha=alpha, label=lbl*" - Tp", color=clr, linewidth=2.0)

  # Format plot
  xlims, dx = [[0, 1], 0.2]
  ylims, dy = [[-100, 100], 5]

  ax2.set_xlim(xlims)
  ax2.set_xticks(xlims[1]:dx:xlims[2])
  ax2.set_xlabel(L"Radial position $r/R$")

  #ax2.set_ylim(ylims)
  #ax2.set_yticks(0:dy:ylims[2])
  ax2.set_ylabel(L"Loading ($\mathrm{N/m}$)")

  ax2.legend(loc="center left", bbox_to_anchor=(1, 0.5), frameon=false, fontsize=10)

  ax2.spines["right"].set_visible(false)
  ax2.spines["top"].set_visible(false)

  fig2.tight_layout()
  typeof(fig2)

  # Save plot
  filename = "Tp_over_rtoR-$(nrevs_to_average)RevsAveragedToRevNo$(rev_to_average_idx).pdf"
  save_plot(save_path_post, filename; fig=fig2)

end



function save_plot(save_path_post::String, filename::String; fig)

  if isdir(save_path_post)
    rm(save_path_post, recursive=true)
    mkdir(save_path_post)
  else
    mkdir(save_path_post)
  end

  fig.savefig(joinpath(save_path_post, filename))

end











  "
  Plots the blade loads of a single turbines simulation.
  "

function plot_KOEFFIZIENT(save_path_post::String, run_name::String)

  # ------------ INITIALIZE PLOT -------------------------
  uns.formatpyplot()

  fig = plt.figure(figsize=[7*1.5, 5*0.75] * 2/3)
  ax = fig.gca()

  xlims, dx = [[0, 10], 2]
  ylims, dy = [[0.04, 0.13], 0.02]

  # ------------ PLOT SIMULATION DATA -----------------------------------


  # ----------- FORMAT -------------------------------------------------------
  ax.set_xlim(xlims)
  ax.set_xticks(xlims[1]:dx:xlims[2])
  ax.set_ylim(ylims)
  ax.set_yticks(ylims[1]:dy:ylims[2])

  ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), frameon=false, fontsize=10)

  ax.set_xlabel("Number of revolutions")
  ax.set_ylabel(L"Thrust coefficient $C_T$")
  ax.set_title(L"$C_T = \frac{T}{\rho n^2 d^4}$", color="gray")

  ax.spines["right"].set_visible(false)
  ax.spines["top"].set_visible(false)


  fig.tight_layout()

  # Save plot
  fig.savefig("dji9443-CTcomparison.png", dpi=300, transparent=true)

end
