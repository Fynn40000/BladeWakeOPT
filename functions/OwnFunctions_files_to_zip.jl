#=##############################################################################
# DESCRIPTION
    Functions to copy files used for postprocessing in folder that needs to be zipped.

# AUTHORSHIP
  * Author          : Fynn Gerhardy
  * Email           : fygerh@gmail.com
  * Created         : Mar 2024
  * Last updated    : Mar 2024
  * License         : -
=###############################################################################




function extract_number_from_filename(filename::String)
    match = Base.match(r".(\d+).xmf", filename)
    if match !== nothing
        return parse(Int, match.captures[1])
    else
        return 0
    end
end

function find_highest_number_in_filenames(folder::String)
    max_number = 0
    for file in readdir(folder)
        number = extract_number_from_filename(file)
        max_number = max(max_number, number)
    end
    return max_number
end


function files_to_zip(source_folder::String, destination_folder::String, last_steps_number::Int64)

    # create destination folder if it does not exist
    if !isdir(destination_folder)
        mkdir(destination_folder)
    end

    # copy paste postprocessing folder
    postprocessing_folder = joinpath(source_folder, "postprocessing")
    if isdir(postprocessing_folder)
        Base.Filesystem.cptree(postprocessing_folder, joinpath(destination_folder, "postprocessing"))
        println("Der Ordner 'postprocessing' und dessen Inhalt wurden nach $destination_folder kopiert.")
    end

    

    # Iteration through all files in source folder
    highest_number = find_highest_number_in_filenames(source_folder)
    for file in readdir(source_folder)
        source_file_path = joinpath(source_folder, file)

        # Überprüfen, ob die Datei eine CSV-Datei ist
        if (isfile(source_file_path) && endswith(file, ".csv")) || 
            (isfile(source_file_path) && endswith(file, ".txt"))
            # Konstruieren des Ziel-Pfades
            destination_file_path = joinpath(destination_folder, file)

            # Kopieren der CSV-Datei in den Zielordner
            cp(source_file_path, destination_file_path)

        # Überprüfen, ob die Datei dem Muster entspricht und die Zahl im angegebenen Bereich liegt
        elseif isfile(source_file_path) && occursin(r"_pfield\.(\d+)\.", file)
            number = parse(Int, match(r"_pfield\.(\d+)\.", file).captures[1])
            if highest_number- last_steps_number <= number <= highest_number
                # Konstruieren des Ziel-Pfades
                destination_file_path = joinpath(destination_folder, file)

                # Kopieren der Datei in den Zielordner
                cp(source_file_path, destination_file_path)
            end
        elseif isfile(source_file_path) && occursin(r"_staticpfield\.(\d+)\.", file)
            number = parse(Int, match(r"_staticpfield\.(\d+)\.", file).captures[1])
            if highest_number- last_steps_number -1 <= number <= highest_number
                # Konstruieren des Ziel-Pfades
                destination_file_path = joinpath(destination_folder, file)

                # Kopieren der Datei in den Zielordner
                cp(source_file_path, destination_file_path)
            end
        end

    end
end