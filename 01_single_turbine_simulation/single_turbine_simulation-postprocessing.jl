#=##############################################################################
# DESCRIPTION
    Postprocessing single turbine simulations
    Plots created for:
        => Blade loads
        => 
=###############################################################################


import FLOWUnsteady as uns
import PyPlot as plt
import CSV
import DataFrames: DataFrame
import Printf: @printf
import PyPlot: @L_str

# declare own functions folder to acces all own functions
start_simulation_path = splitdir(@__FILE__)[1]
include(start_simulation_path, "..", "functions")

uns.formatpyplot()

println("\nPostprocessing...\n")