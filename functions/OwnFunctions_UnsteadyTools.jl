#=##############################################################################
# NOTE:
    This code is copied from the byuflowlab package "UnsteadyTools"
=###############################################################################




"""
chopWake!(pfield::vpm.ParticleField, cut::String, location)

Chops the wake according to the directions specified.

If `cut == "plane"`, removes particles at all locations satisfying:

* x_particle[px, py, pz] s.t. a*px + b*py + c*pz > 1, where
* location : [a, b, c] describes the equation for the plane ax + by + cz = 1

E.g., to remove all particles above the plane x = 2.0, run:

```julia
chopWake!(pfield, "plane", [0.5, 0.0, 0.0])
```

If `cut == "cylinder"`, removes particles at all locations that don't lie inside a cylinder defined by:

* location : [radius_rotor1::Number, radius_rotor2::Number...] describes the radii of the cylinders created by any of the rotors (must be length Nrotors)

E.g., to remove all particles outside a cylinder with radius 0.25R_rotor aligned with the system's single rotor

```julia
chopWake!(pfield, "cylinder", [0.25])
```

If `cut == "timeshed"`, removes all particles shed at times within the inclusive bounds:

* location : [lowerboundtimestep, upperboundtimestep]

E.g., to remove all particles shed from steps 5 to 50, run:

```julia
chopWake!(pfield, "timeshed", [5,50])
```
"""
function chopWake!(pfield::vpm.ParticleField, cut::String, location; rotors=[], verbose=true, v_lvl=0, reverselogic=false)
    #println("\t"^v_lvl * "Chopping wake...")
    if cut == "plane"
        #println("\t"^v_lvl,"Begin plane chop...")
        # println("\t"^v_lvl * "Sherlock! commencing `plane` cut...")
        # println("\t"^v_lvl * "Sherlock! N particles = ",pfield.np)
        a = location[1]
        b = location[2]
        c = location[3]
        if a < 0 || b < 0 || c < 0
            throw("negative parameters a, b, c not supported")
        end
        # println("\t\tpfield._p_field[pfield.np,:] = ",pfield._p_field[pfield.np,:])
        ineq = reverselogic ? "<" : ">"
        #println("\t"^v_lvl * "Chopping all px " * ineq * "$(round(1/a, digits=2))")
        numdeleted = 0
        for  pi = pfield.np:-1:1
            x = vpm.get_X(pfield, pi)
            chopaway = (a*x[1] + b*x[2] + c*x[3]) > 1.0
            reverselogic ? chopaway = !chopaway : nothing
            if chopaway
                #println("\t"^(v_lvl+1),"Deleting particle n=",pi)
                #vpm.delparticle(pfield, pi)
                vpm.remove_particle(pfield, pi)
                numdeleted += 1
            end
        end
        #println("\t"^v_lvl,"particles deleted: ",numdeleted)
    elseif cut == "cylinder"
        println("\t"^v_lvl,"Begin cylinder chop...")
        numdeleted = 0
        if length(rotors) != length(location); throw("Location vector for cylinder chop inconsistent with number of rotors."); end
        for pi = pfield.np:-1:1
            incylinders = Bool[]
            x = vpm.get_X(pfield, pi)            
            for rotori = 1:length(rotors)
                point = rotors[rotori]._wingsystem.O
                δx = x .- point
                orientation = rotors[rotori]._wingsystem.Oaxis[:,1]
                radius = location[rotori]
                incylinder = norm(δx .- dot(δx, orientation)) < radius * rotors[rotori].rotorR
                reverselogic ? chopaway = !chopaway : nothing
                push!(incylinders,incylinder)
            end
            if !(true in incylinders)
                println("\t"^(v_lvl+1), "Deleting particle n=",pi)
                #vpm.delparticle(pfield, pi)
                vpm.remove_particle(pfield, pi)
                numdeleted += 1
            end
        end
    elseif cut == "timeshed"
        println("Sherlock! \n\tpfield._timeshed = $(pfield._timeshed)\n")
        chopthese = find(x -> x >= location[1] && x <= location[2], pfield._timeshed)
        for pi = pfield.np:-1:1
            if pi in chopthese
                println("\t"^(v_lvl+1), "Deleting particle n=",pi)
                #vpm.delparticle(pfield, pi)
                vpm.remove_particle(pfield, pi)
            end
        end
        println("\tpfield._timeshed = $(pfield._timeshed)\n")
    elseif cut == "maxgamma"
        println("Sherlock! removing high circulation particles")
        chopthese = find(x -> x >= location, [maximum(pfield._p_field[i,4:6]) for i in 1:size(pfield._p_field)[1]])
        for pi = pfield.np:-1:1
            if pi in chopthese
                println("\t"^(v_lvl+1), "Deleting particle n=",pi)
                #vpm.delparticle(pfield, pi)
                vpm.remove_particle(pfield, pi)
            end
        end
    else
        throw("Desired `cut` not found.")
    end
    #println("\t"^v_lvl * "Finished chopping wake.")
    #println("")
end

"""
test the logic of the `reverselogic` variable in chopWake!()
"""
function testChopWake(x,location,reverselogic)
    a = location[1]
    b = location[2]
    c = location[3]
    if (a*x[1] + b*x[2] + c*x[3]) * (-1)^reverselogic > (-1)^reverselogic
        return true
    else
        return false
    end
end