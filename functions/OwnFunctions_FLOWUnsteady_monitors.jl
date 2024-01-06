#=##############################################################################
# NOTE:
    This code has to be copied and pasted into the file "src/FLOWUnsteady_monitors.jl" 
    in the FLOWUnsteady Repository under the Julia packages
=###############################################################################

"""
    generate_monitor_turbines(rotors::Array{vlm.Rotor}, J_ref::Real,
                                rho_ref::Real, RPM_ref::Real, nsteps_sim::Int, magVinfx::Real, turbine_flag::Bool;
                                save_path=nothing)

Generate a turbine monitor plotting the aerodynamic performance and blade loading
at every time step.

The aerodynamic performance consists of thrust coefficient
\$C_T = \\frac{T}{q*A}\$, torque coefficient
\$C_Q = \\frac{Q}{q*R*A}\$, and power coefficient
(((((((((\$C_P = \\frac{P}{q*A*u_\\infty}\$)))))))))

with
\$q = 0.5*rho*u_\\infty^2\$, free stream kinetic energy
\$A = pi*R^2\$, rotor area
\$u_\\infty = magVinfx\$

* `J_ref` and `rho_ref` are the reference advance ratio and air density used for calculating propulsive efficiency and coefficients. The advance ratio used here is defined as \$J=\\frac{u_\\infty}{n d}\$ with \$n = \\frac{\\mathrm{RPM}}{60}\$.
* `RPM_ref` is the reference RPM used to estimate the age of the wake.
* `nsteps_sim` is the number of time steps by the end of the simulation (used for generating the color gradient).
* Use `save_path` to indicate a directory where to save the plots. If so, it will also generate a CSV file with \$C_T\$, \$C_Q\$, and \$C_P\$.

"""
function generate_monitor_turbines( rotors::Array{vlm.Rotor, 1},
                                    J_ref::Real, rho_ref::Real, RPM_ref::Real,
                                    nsteps_sim::Int,
                                    magVinfx::Real,                 # Freestream velocity
                                    turbine_flag::Bool;
                                    t_scale=1.0,                    # Time scaling factor
                                    t_lbl="Simulation time (s)",    # Time-axis label
                                    # OUTPUT OPTIONS
                                    out_figs=[],
                                    out_figaxs=[],
                                    save_path=nothing,
                                    run_name="rotor",
                                    figname="monitor_rotor",
                                    disp_conv=true,
                                    conv_suff="_convergence.csv",
                                    save_init_plots=true,
                                    figsize_factor=5/6,
                                    nsteps_plot=1,
                                    nsteps_savefig=10,
                                    colors="rgbcmy"^100,
                                    stls="o^*.px"^100, )

    fcalls = 0                  # Number of function calls

    # Name of convergence file
    if save_path!=nothing
        fname = joinpath(save_path, run_name*conv_suff)
    end

    # Call figure
    if disp_conv
        formatpyplot()
        fig = plt.figure(figname, figsize=[7*3, 5*2]*figsize_factor)
        axs = fig.subplots(2, 3)
        axs = [axs[6], axs[2], axs[4], axs[1], axs[3], axs[5]]

        push!(out_figs, fig)
        push!(out_figaxs, axs)
    end
    
    # Function for run_vpm! to call on each iteration
    function extra_runtime_function(sim::Simulation{V, M, R},
                                    PFIELD::vpm.ParticleField,
                                    T::Real, DT::Real; optargs...
                                   ) where{V<:AbstractVLMVehicle, M, R}

        # rotors = vcat(sim.vehicle.rotor_systems...)
        angle = T*360*RPM_ref/60
        t_scaled = T*t_scale
        
        if fcalls==0
            # Format subplots
            if disp_conv
                ax = axs[1]
                ax.set_title("Circulation distribution", color="gray")
                ax.set_xlabel("Element index")
                ax.set_ylabel(L"Circulation $\Gamma$ (m$^2$/s)")

                ax = axs[2]
                ax.set_title("Normal force distribution", color="gray")
                ax.set_xlabel("Element index")
                ax.set_ylabel(L"Normal load $N_p$ (N/m)")

                ax = axs[3]
                ax.set_title("Tangential force distribution", color="gray")
                ax.set_xlabel("Element index")
                ax.set_ylabel(L"Tangential load $T_p$ (N/m)")

                ax = axs[4]
                ax.set_title(L"$C_T = \frac{T}{0.5 \rho A u_\infty^2}$", color="gray")
                ax.set_xlabel(t_lbl)
                ax.set_ylabel(L"Thrust Coefficient $C_T$")

                ax = axs[5]
                ax.set_title(L"$C_Q = \frac{Q}{0.5 \rho A u_\infty^2 R}$", color="gray")
                ax.set_xlabel(t_lbl)
                ax.set_ylabel(L"Torque Coefficient $C_Q$")

                ax = axs[6]
                ax.set_title(L"$C_P = \frac{P}{0.5 \rho A u_\infty^3}$", color="gray")
                ax.set_xlabel(t_lbl)
                ax.set_ylabel(L"Power Coefficient $C_P$")


                for ax in axs
                    ax.spines["right"].set_visible(false)
                    ax.spines["top"].set_visible(false)
                    # ax.grid(true, color="0.8", linestyle="--")
                end

                fig.tight_layout()
            end

            # Convergence file header
            if save_path!=nothing
                f = open(fname, "w")
                print(f, "ref age (deg),T,DT")
                for (i, rotor) in enumerate(rotors)
                    print(f, ",RPM_$i,CT_$i,CQ_$i,CP_$i")
                end
                print(f, "\n")
                close(f)
            end

            # Save initialization plots (including polars)
            if save_init_plots && save_path!=nothing
                for fi in plt.get_fignums()
                    this_fig = plt.figure(fi)
                    this_fig.savefig(joinpath(save_path, run_name*"_initplot$(fi).png"),
                                                            transparent=false, dpi=300)
                end
            end
        end

        # Write rotor position and time on convergence file
        if save_path!=nothing
            f = open(fname, "a")
            print(f, angle, ",", T, ",", DT)
        end


        # Plot circulation and loads distributions
        if  PFIELD.nt%nsteps_plot==0 && disp_conv

            cratio = PFIELD.nt/nsteps_sim
            cratio = cratio > 1 ? 1 : cratio
            clr = fcalls==0 && false ? (0,0,0) : (1-cratio, 0, cratio)
            stl = fcalls==0 && false ? "o" : "-"
            alpha = fcalls==0 && false ? 1 : 0.5

            # Circulation distribution
            this_sol = []
            for rotor in rotors
                this_sol = vcat(this_sol, [vlm.get_blade(rotor, j).sol["Gamma"] for j in 1:rotor.B]...)
            end
            axs[1].plot(1:size(this_sol,1), this_sol, stl, alpha=alpha, color=clr)

            # Np distribution
            this_sol = []
            for rotor in rotors
                this_sol = vcat(this_sol, rotor.sol["Np"]["field_data"]...)
            end
            axs[2].plot(1:size(this_sol,1), this_sol, stl, alpha=alpha, color=clr)

            # Tp distribution
            this_sol = []
            for rotor in rotors
                this_sol = vcat(this_sol, rotor.sol["Tp"]["field_data"]...)
            end
            axs[3].plot(1:size(this_sol,1), this_sol, stl, alpha=alpha, color=clr)
        end

        # Plot performance parameters
        for (i, rotor) in enumerate(rotors)

            #CT, CQ = vlm.calc_thrust_torque_coeffs(rotor, rho_ref, magVinfx, turbine_flag)
            thrust, torque = vlm.calc_thrust_torque(rotor)
            power = torque * (2*pi*rotor.RPM)/60
            q = 0.5*rho_ref*magVinfx^2
            A = pi*rotor.rotorR^2
            CT = thrust/(q*A)
            CQ = torque/(q*rotor.rotorR*A)
            CP = power/(q*A*magVinfx)

            roR = rotor.sol["roR"]["field_data"][1]
            #println("roR")
            #println(roR)
            #cn = rotor.sol["cn"]["field_data"][1]
            #println("cn")
            #println(cn)
            #ct = rotor.sol["ct"]["field_data"][1]
            #println("ct")
            #println(ct)
            #println(" ")

            if PFIELD.nt%nsteps_plot==0 && disp_conv
                axs[4].plot([t_scaled], [CT], "$(stls[i])", alpha=alpha, color=clr)
                axs[5].plot([t_scaled], [CQ], "$(stls[i])", alpha=alpha, color=clr)
                axs[6].plot([t_scaled], [CP], "$(stls[i])", alpha=alpha, color=clr)
            end

            if save_path!=nothing
                print(f, ",", rotor.RPM, ",", CT, ",", CQ, ",", CP)
            end
        end

        if disp_conv
            # Save figure
            if PFIELD.nt%nsteps_savefig==0 && fcalls!=0 && save_path!=nothing
                fig.savefig(joinpath(save_path, run_name*"_convergence.png"),
                                                    transparent=false, dpi=300)
            end
        end

        # Close convergence file
        if save_path!=nothing
            print(f, "\n")
            close(f)
        end

        fcalls += 1

        return false
    end

    return extra_runtime_function
end