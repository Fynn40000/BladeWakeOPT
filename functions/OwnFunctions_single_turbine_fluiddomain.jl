#=##############################################################################
# DESCRIPTION
    Postprocessing single turbine simulations fluid domain
=###############################################################################

"
Postprocesses the fluiddomain of a single turbines simulation.
"

function postprocess_fluiddomain(# INPUT OPTIONS
                                 simulation_name = "rotorhover-example",               # Simulation to read
                                 read_path       = "/home/fynn/Repositories/BladeWakeOPT/"*simulation_name, # Where to read simulation from
                                
                                 pfield_prefix   = "singlerotor_pfield",      # Prefix of particle field files to read
                                 staticpfield_prefix = "singlerotor_staticpfield", # Prefix of static particle field files to read
                                 
                                 nums            = [719],              # Time steps to process)

    # OUTPUT OPTIONS
save_path       = joinpath(read_path, "..", simulation_name*"-fdom"),  # Where to save fluid domain
output_prefix   = "singlerotor",             # Prefix of output files
prompt          = true,                      # Whether to prompt the user
verbose         = true,                      # Enable verbose
v_lvl           = 0                         # Verbose indentation level
  
)

end