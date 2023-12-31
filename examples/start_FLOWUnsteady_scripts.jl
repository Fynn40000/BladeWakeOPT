#import FLOWUnsteady as uns
#include(joinpath(uns.examples_path, "wing", "wing.jl"))
#include(joinpath(uns.examples_path, "tetheredwing.jl"))
#include(joinpath(uns.examples_path, "rotorhover", "rotorhover.jl"))
#include(joinpath(uns.examples_path, "rotorhover", "rotorhover_fluiddomain.jl"))
#include(joinpath(uns.examples_path, "propeller", "propeller.jl"))

# VPM Leapfrog example
import FLOWVPM as vpm
include(joinpath(vpm.examples_path, "vortexrings", "run_leapfrog.jl"))
