
#include(joinpath("/home/fynn/Repositories/BladeWakeOPT/05_files_to_zip/files_to_zip.jl"))

using CSV
using Base
using Base.Iterators: filter
this_file_path = splitdir(@__FILE__)[1]
include(joinpath(this_file_path, "..", "functions", "OwnFunctions.jl"))        # Read file with module definition first
using .OwnFunctions


# ----------------- ENTER ALL SIMULATION FOLDER PATHS TO BE COPIED --------------------------------------------------------  
#     # BLADE ELEMENT STUDY
# study_folder = joinpath(this_file_path, "..", "02_parameter_study", "Convergence", "ParameterStudy-Blade_Elements")
# BEs = [170]
# sim_folders = []
# for BE in BEs
#     push!(sim_folders, joinpath(study_folder, "NREL5MW_BladeElement_$(BE)"))
# end
# last_steps_number = 100 # COPY THIS MANY LAST TIME STEPS

#     # TIME STEPS STUDY
# study_folder = joinpath(this_file_path, "..", "02_parameter_study", "Convergence", "ParameterStudy-Time_Steps")
# tsteps = [80, 90, 103, 120, 144, 180]
# sim_folders = []
# for tstep in tsteps
#     push!(sim_folders, joinpath(study_folder, "NREL5MW_TimeSteps_$(tstep)"))
# end
# last_steps_number = 200 # COPY THIS MANY LAST TIME STEPS

    # TWIST STUDY
study_folder = joinpath(this_file_path, "..", "02_parameter_study", "ParameterSpace", "ParameterStudy-Twist_Angle")
sim_folders = [joinpath(study_folder, "NREL5MW_TwistAngle_Knauer"),
                joinpath(study_folder, "NREL5MW_TwistAngle_Kelley"),
                joinpath(study_folder, "NREL5MW_TwistAngle_LLR"),
                joinpath(study_folder, "NREL5MW_TwistAngle_HLR"),
                joinpath(study_folder, "NREL5MW_TwistAngle_CBR1_2BC"),
                joinpath(study_folder, "NREL5MW_TwistAngle_CBR2_2BC"),
                joinpath(study_folder, "NREL5MW_TwistAngle_CBR1_4BC"),
                joinpath(study_folder, "NREL5MW_TwistAngle_CBR2_4BC")]
last_steps_number = 130 # COPY THIS MANY LAST TIME STEPS

#     # OPTIMIZATION STUDY
# study_folder = joinpath(this_file_path, "..", "03_optimization", "solutions-Knauer_opt")
# sim_folders = [joinpath(study_folder, "NREL5MW_solution_1"),
#                 joinpath(study_folder, "NREL5MW_solution_2"),
#                 joinpath(study_folder, "NREL5MW_solution_3"),
#                 joinpath(study_folder, "NREL5MW_solution_4"),
#                 joinpath(study_folder, "NREL5MW_solution_5"),
#                 joinpath(study_folder, "NREL5MW_solution_6"),
#                 joinpath(study_folder, "NREL5MW_solution_7"),
#                 joinpath(study_folder, "NREL5MW_solution_8"),
#                 joinpath(study_folder, "NREL5MW_solution_9"),
#                 joinpath(study_folder, "NREL5MW_solution_10")]
# last_steps_number = 130 # COPY THIS MANY LAST TIME STEPS



# ----------------- CREATE "_to_zip" FOLDERS --------------------------------------------------------  
for source_folder in sim_folders
    destination_folder = joinpath(this_file_path, basename(source_folder) * "_to_zip")
    OwnFunctions.files_to_zip(source_folder, destination_folder, last_steps_number)     # Copy files into according "_to_zip" folder
end

