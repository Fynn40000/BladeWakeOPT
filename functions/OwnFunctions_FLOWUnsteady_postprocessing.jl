#=##############################################################################
# NOTE:
    This code has to be copied and pasted into the file "src/FLOWUnsteady_postprocessing.jl" 
    in the FLOWUnsteady Repository under the Julia packages
=###############################################################################



"""
    postprocess_bladeloading_turbine(read_path;
                                O           = zeros(3),     # Rotor center
                                rotor_axis  = [-1, 0, 0],   # Rotor centerline axis
                                Ftot_axis   = nothing,      # Use a different centerline axis for forces if given
                                filename    = "singlerotor_Rotor_Blade1_vlm-statistics.vtk", # File name
                                fieldsuff   = "-mean"       # Suffix of fields "Gamma" and "Ftot", if any
                                num_elements::Int64= 1,     # Number of blade elements used for discretization (only used when "_loft.vtk" files are evaluated)
                                debug::Bool=false
                                )

Read a blade VTK file `filename` under directory `read_path` and returns
the circulation and force components of the load distribution along the blade.

Return: `rs, Gamma, Np, Tp, Rp, Zhat, Rhat, That, Ftot` for "_vlm" files 
and
Return: `roR, Gamma, Np, Tp, L, D, cn, ct, cl, cd, a, a_tangential, aoa, twist, flowangle, Vx, Vy, w, F, loc_solidity, Zhat` for "_loft" files
"""
function postprocess_bladeloading_turbine(read_path;
                                    O           = zeros(3),     # Rotor center
                                    rotor_axis  = [-1, 0, 0],   # Rotor centerline axis
                                    Ftot_axis   = nothing,      # Use a different centerline axis for forces if given
                                    filename    = "singlerotor_Rotor_Blade1_vlm-statistics.vtk", # File name
                                    fieldsuff   = "-mean",      # Suffix of fields "Gamma" and "Ftot", if any
                                    num_elements::Int64= 1,     # Number of blade elements used for discretization (only used when "_loft.vtk" files are evaluated)
                                    debug::Bool=false
                                    )

    @assert isapprox(norm(rotor_axis), 1.0) "Non-unitary rotor axis $(rotor_axis)"

    _Ftot_axis = Ftot_axis==nothing ? rotor_axis : Ftot_axis

    points, cells, cell_types, data = gt.read_vtk(filename; path=read_path)

    if occursin("_vlm", filename)
        maxind = Int(size(cells, 1)/2)

        Xs = [mean([points[:, pi+1] for pi in cell]) for cell in cells[1:maxind]]
        Rs = [(X-O) - rotor_axis*dot(X-O, rotor_axis) for X in Xs]

        # Blade centerline
        Rhat = mean(Rs)
        Rhat /= norm(Rhat)

        rs = [dot(R, Rhat) for R in Rs]

        # Tangent vector
        That = cross(_Ftot_axis, Rhat)

        # Grabs only the rectangular cells
        Gamma = data["CELL_DATA"]["Gamma"*fieldsuff][1:maxind]
        Ftot = data["CELL_DATA"]["Ftot"*fieldsuff][:, 1:maxind]
      
        nr = length(rs)
        Np = zeros(nr)              # Normal component
        Rp = zeros(nr)              # Radial component
        Tp = zeros(nr)              # Tangential component

        for i in 1:nr
          F = Ftot[:, i]
          Np[i] = dot(F, _Ftot_axis)
          Rp[i] = dot(F, Rhat)
          Tp[i] = dot(F, That)
        end

        return rs, Gamma, Np, Tp, Rp, _Ftot_axis, Rhat, That, Ftot

    elseif occursin("_loft", filename)
        #maxind = Int(floor(size(cells, 1)/2))

        roR = data["POINT_DATA"]["roR"*fieldsuff][1:end]
        Gamma = data["POINT_DATA"]["Gamma"*fieldsuff][1:end]
        Np = data["POINT_DATA"]["Np"*fieldsuff][1:end] # Normal (Thrust) force
        Tp = data["POINT_DATA"]["Tp"*fieldsuff][1:end] # Tangential force
        if debug
            L = data["POINT_DATA"]["Lift"*fieldsuff][1:end] # Lift
            D = data["POINT_DATA"]["Drag"*fieldsuff][1:end] # Drag

            cn = data["POINT_DATA"]["cn"*fieldsuff][1:end] # Thrust coefficient
            ct = data["POINT_DATA"]["ct"*fieldsuff][1:end] # Tangential coefficient NOT Thrust coefficient
            cl = data["POINT_DATA"]["cl"*fieldsuff][1:end] # Lift coefficient
            cd = data["POINT_DATA"]["cd"*fieldsuff][1:end] # Drag coefficient
            
            aoa = data["POINT_DATA"]["ThetaEffDeg"*fieldsuff][1:end] # effective angle of attack
            twist = data["POINT_DATA"]["twistdeg"*fieldsuff][1:end] # twist angle of the 
            flowangle = data["POINT_DATA"]["flowangledeg"*fieldsuff][1:end] # flow angle

            Vx = data["POINT_DATA"]["Vx_out"*fieldsuff][1:end] # local flow velocity at controlpoint (x direction)
            Vy = data["POINT_DATA"]["Vy_out"*fieldsuff][1:end] # local flow velocity at controlpoint (y direction)

            F = data["POINT_DATA"]["F_out"*fieldsuff][1:end]   # Prandtl loss factor
            loc_solidity = data["POINT_DATA"]["loc_solidity"*fieldsuff][1:end]   # local solidity
        else
            L = nothing
            D = nothing
            cn = nothing
            ct = nothing
            cl = nothing
            cd = nothing
            aoa = nothing
            twist = nothing
            flowangle = nothing
            Vx = nothing
            Vy = nothing
            F = nothing
            loc_solidity = nothing
            w = nothing
            a = nothing
            a_tangential = nothing
        end


        # Clear the arrays up...
        roR = _clear_arr(roR, num_elements)
        Gamma = _clear_arr(Gamma, num_elements)
        Np = _clear_arr(Np, num_elements)
        Tp = _clear_arr(Tp, num_elements)

        if debug
            L = _clear_arr(L, num_elements)
            D = _clear_arr(D, num_elements)

            cn = _clear_arr(cn, num_elements)
            ct = _clear_arr(ct, num_elements)
            cl = _clear_arr(cl, num_elements)
            cd = _clear_arr(cd, num_elements)

            aoa = _clear_arr(aoa, num_elements)
            twist = _clear_arr(twist, num_elements)
            flowangle = _clear_arr(flowangle, num_elements)
            Vx = _clear_arr(Vx, num_elements)
            Vy = _clear_arr(Vy, num_elements)
            F = _clear_arr(F, num_elements)
            loc_solidity = _clear_arr(loc_solidity, num_elements)

            # Calculate relative flow velocity based on given flow triangle
            w = sqrt.(Vx .* Vx .+ Vy .* Vy)

            # Calculate induction factor distribution based on thrust coeff distribution
            a = (1 .- sqrt.(1 .- cn)) ./ 2
            a_tangential = (1 .- sqrt.(1 .- ct)) ./ 2
        end

        return roR, Gamma, Np, Tp, L, D, cn, ct, cl, cd, a, a_tangential, aoa, twist, flowangle, Vx, Vy, w, F, loc_solidity, _Ftot_axis
    end

end

function _clear_arr(arr::Vector{Float64}, num_elements::Int64=1)
    step = Int(length(arr)/(num_elements+2))
    arr = arr[1:step:end]
    popfirst!(arr)
    pop!(arr)
    return arr
end