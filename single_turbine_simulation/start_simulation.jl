#=##############################################################################
# DESCRIPTION
    Simulation of a single wind turbine rotor.

# AUTHORSHIP
  * Author          : Fynn Gerhardy
  * Email           : fygerh@gmail.com
  * Created         : Dec 2023
  * Last updated    : Dec 2023
  * License         : -
=###############################################################################


import FLOWUnsteady as uns
import FLOWVLM as vlm
import FLOWVPM as vpm

run_name        = "single_turbine_simulation"       # Name of this simulation

start_simulation_path = splitdir(@__FILE__)[1]
save_path       = joinpath(start_simulation_path, "data_out", run_name) # Where to save this simulation #splitdir(@__FILE__)[1] * "/data_out" * "/" * run_name                 
paraview        = true                      # Whether to visualize with Paraview

# ----------------- GEOMETRY PARAMETERS ----------------------------------------

# Rotor geometry
rotor_file      = "apc10x7.csv"                                             # Rotor geometry
data_path       = joinpath(start_simulation_path, "..", "database")         # Path to rotor database
pitch           = 0.0                                                       # (deg) collective pitch of blades
CW              = false                                                     # Clock-wise rotation
xfoil           = true                                                      # Whether to run XFOIL
ncrit           = 9                                                         # Turbulence criterion for XFOIL

# NOTE: If `xfoil=true`, XFOIL will be run to generate the airfoil polars used
#       by blade elements before starting the simulation. XFOIL is run
#       on the airfoil contours found in `rotor_file` at the corresponding
#       local Reynolds and Mach numbers along the blade.
#       Alternatively, the user can provide pre-computer airfoil polars using
#       `xfoil=false` and pointing to polar files through `rotor_file`.

# Discretization
n               = 20                        # Number of blade elements per blade
r               = 1/5                       # Geometric expansion of elements

# NOTE: Here a geometric expansion of 1/5 means that the spacing between the
#       tip elements is 1/5 of the spacing between the hub elements. Refine the
#       discretization towards the blade tip like this in order to better
#       resolve the tip vortex.

# Read radius of this rotor and number of blades
R, B            = uns.read_rotor(rotor_file; data_path=data_path)[[1,3]]

# ----------------- SIMULATION PARAMETERS --------------------------------------

# Operating conditions
RPM             = 9200                      # RPM
J               = 0.4                       # Advance ratio Vinf/(nD)
AOA             = 0                         # (deg) Angle of attack (incidence angle)

rho             = 1.225                     # (kg/m^3) air density
mu              = 1.81e-5                   # (kg/ms) air dynamic viscosity
speedofsound    = 342.35                    # (m/s) speed of sound

magVinf         = J*RPM/60*(2*R)
Vinf(X, t)      = magVinf*[cosd(AOA), sind(AOA), 0] # (m/s) freestream velocity vector

ReD             = 2*pi*RPM/60*R * rho/mu * 2*R      # Diameter-based Reynolds number
Matip           = 2*pi*RPM/60*R / speedofsound      # Tip Mach number

println("""
    RPM:    $(RPM)
    Vinf:   $(Vinf(zeros(3), 0)) m/s
    Matip:  $(round(Matip, digits=3))
    ReD:    $(round(ReD, digits=0))
""")

# ----------------- SOLVER PARAMETERS ------------------------------------------

# Aerodynamic solver
VehicleType     = uns.UVLMVehicle           # Unsteady solver
# VehicleType   = uns.QVLMVehicle           # Quasi-steady solver
const_solution  = VehicleType==uns.QVLMVehicle  # Whether to assume that the
                                                # solution is constant or not
# Time parameters
nrevs           = 4                         # Number of revolutions in simulation
nsteps_per_rev  = 36                        # Time steps per revolution
nsteps          = const_solution ? 2 : nrevs*nsteps_per_rev # Number of time steps
ttot            = nsteps/nsteps_per_rev / (RPM/60)       # (s) total simulation time

# VPM particle shedding
p_per_step      = 2                         # Sheds per time step
shed_starting   = true                      # Whether to shed starting vortex
shed_unsteady   = true                      # Whether to shed vorticity from unsteady loading
max_particles   = ((2*n+1)*B)*nsteps*p_per_step + 1 # Maximum number of particles

# Regularization
sigma_rotor_surf= R/40                      # Rotor-on-VPM smoothing radius
lambda_vpm      = 2.125                     # VPM core overlap
                                            # VPM smoothing radius
sigma_vpm_overwrite = lambda_vpm * 2*pi*R/(nsteps_per_rev*p_per_step)

# Rotor solver
vlm_rlx         = 0.7                       # VLM relaxation <-- this also applied to rotors
hubtiploss_correction = vlm.hubtiploss_nocorrection # Hub and tip loss correction

# VPM solver
vpm_viscous     = vpm.Inviscid()            # VPM viscous diffusion scheme

# NOTE: In most practical situations, open rotors operate at a Reynolds number
#       high enough that viscous diffusion in the wake is negligible.
#       Hence, it does not make much of a difference whether we run the
#       simulation with viscous diffusion enabled or not.


if VehicleType == uns.QVLMVehicle
    # NOTE: If the quasi-steady solver is used, this mutes warnings regarding
    #       potential colinear vortex filaments. This is needed since the
    #       quasi-steady solver will probe induced velocities at the lifting
    #       line of the blade
    uns.vlm.VLMSolver._mute_warning(true)
end