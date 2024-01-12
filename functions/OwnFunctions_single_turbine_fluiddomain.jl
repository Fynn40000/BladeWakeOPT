#=##############################################################################
# DESCRIPTION
    Postprocessing single turbine simulations fluid domain
=###############################################################################

"
Postprocesses the fluiddomain of a single turbines simulation.
"

function postprocess_fluiddomain(# ---- ESSENTIAL ARGUMENTS ---------
                                save_path::String,                      # path to folder where simulation is stored in
                                run_name::String,                       # name of simulation
                                file_suffix::String,                    # postprocessed fluid domain file suffix
                                R::Float64,                             # Rotor tip radius
                                AOA::Float64,                           # Angle of attack
                                nums::Vector{Int64};                    # Time steps to process
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

  r_path = save_path                        # path to read data from (path to folder all the simulation data is stored in)
  s_path = save_path#joinpath(r_path, run_name*"-fdom")# path to folder all postprocessed fluiddomain files will be stored in

  output_prefix   = run_name                # Prefix of output files

  # Grid
  L               = R                                               # (m) reference length = rotor tip radius
  dx, dy, dz      = L/x_resolution, L/y_resolution, L/z_resolution  # (m) cell size in each direction
  Pmin            = 2*L*[x_bound_min, y_bound_min, z_bound_min]     # (m) minimum bounds (bound_factor * rotor diameter)
  Pmax            = 2*L*[x_bound_max, y_bound_max, z_bound_max]     # (m) maximum bounds (bound_factor * rotor diameter)
  NDIVS           = ceil.(Int, (Pmax .- Pmin)./[dx, dy, dz])        # Number of cells in each direction
  nnodes          = prod(NDIVS .+ 1)                                # Total number of nodes

  Oaxis           = uns.gt.rotation_matrix2(0, 0, AOA)                  # Orientation of grid
  

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
    uns.computefluiddomain(Pmin, Pmax, NDIVS,
                           maxparticles,
                           these_nums, r_path, pfield_prefix;
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