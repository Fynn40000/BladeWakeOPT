#=##############################################################################
# DESCRIPTION
    Postprocessing single turbine simulations
    Plots created for:
        => Blade loads
        => 
=###############################################################################


"
Plots the blade loads of a single turbines simulation.
"

function plot_blade_loading(save_path::String, save_path_post::String, run_name::String)
  
  # Plot styles
  stl = "-" # style
  clr = "dodgerblue" # color
  alpha = 1.0 # alpha
  lbl = run_name # label

  rotor_axis = [-1.0, 0.0, 0.0]       # Rotor centerline axis
  R          = 63.0                   # (m) rotor radius
  rev_to_average_idx = 4              # This is the last revolution that is used for averaging the values
  nrevs_to_average = 1                # Number of revolutions to average


  nsteps_per_rev = Dict()
  loading_mean = Dict()
  Loading_std = Dict()

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

  # Start plotting
  fig = plt.figure(figsize=[7*1.5, 5] * 2/3)
  ax = fig.gca()
  
  # Plot y=0
  ax.plot([0, 1], zeros(2), ":k", linewidth=1)

  # Read and plot FLOWUnsteady mean blade loading and deviation
  r_path = joinpath(save_path, "..", run_name*"-statistics")

  (rs, Gamma,
        Np, Tp) = uns.postprocess_bladeloading(r_path;
                                                O           = zeros(3),
                                                rotor_axis  = rotor_axis,
                                                filename    = run_name*"_Rotor_Blade1_vlm-statistics.vtk",
                                                fieldsuff   = "-mean"
                                                )

  ax.plot(rs/R, Np, stl; alpha=alpha, label=lbl, color=clr, linewidth=2.0)#"-", "dodgerblue", 1.0, "rVPM - high fidelity"


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

  filename = "testplot.pdf"
  save_plot(save_path_post, filename; fig=fig)

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
