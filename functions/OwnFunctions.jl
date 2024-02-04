


module OwnFunctions


# ------------ GENERIC MODULES -------------------------------------------------
import PyPlot as plt
import CSV
import DataFrames: DataFrame
import Printf: @printf
import PyPlot: @L_str

# ------------ FLOW CODES ------------------------------------------------------
# NOTE: Unregistered packages available at https://github.com/byuflowlab
import FLOWUnsteady as uns
import FLOWVLM as vlm
import FLOWVPM as vpm


# ------------ HEADERS ---------------------------------------------------------
for header_name in ["single_turbine_simulation_postprocessing", "single_turbine_fluiddomain", 
                    "start_single_turbine_simulation", "utilities_fluiddomain"]#, "FLOWUnsteady_monitors", "FLOWUnsteady_postprocessing", "FLOWUnsteady_rotor", "FLOWVLM_rotor"
#"single_turbine_postprocessing"
include("OwnFunctions_"*header_name*".jl")

end



end # END OF MODULE