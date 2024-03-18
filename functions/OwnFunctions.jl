


module OwnFunctions


# ------------ GENERIC MODULES -------------------------------------------------
import PyPlot as plt
import CSV
import DataFrames
import DataFrames: DataFrame
import Printf: @printf
import PyPlot: @L_str

using LatinHypercubeSampling
using Interpolations

# ------------ FLOW CODES ------------------------------------------------------
# NOTE: Unregistered packages available at https://github.com/byuflowlab
import FLOWUnsteady as uns
import FLOWVLM as vlm
import FLOWVPM as vpm

import AirfoilPrep
ap = AirfoilPrep

import GeometricTools
gt = GeometricTools

import CCBlade
ccb = CCBlade

# ------------ HEADERS ---------------------------------------------------------
for header_name in ["start_single_turbine_simulation", "single_turbine_simulation_postprocessing", "optimization", 
                     "utilities_fluiddomain", "files_to_zip", "UnsteadyTools", "FLOWUnsteady_monitors", "FLOWUnsteady_rotor", "FLOWVLM_rotor"]#, "FLOWUnsteady_monitors", "FLOWUnsteady_postprocessing"
#"single_turbine_postprocessing"
include("OwnFunctions_"*header_name*".jl")

end



end # END OF MODULE