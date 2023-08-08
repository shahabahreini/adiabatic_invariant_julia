module Simulation2D

# Import necessary modules and functions
include("SharedFunctions.jl")
using .SharedFunctions
using DifferentialEquations
using ArbNumerics
using Distributed
using DataFrames
using CSV

export CylindricalProblem!, SolvingtheProblem, ParticleMotion

"""
    calculateOmega(z_, rho_, type_)

	Calculate the omega values for the given z and rho values based on the type.
	- If type is "rho", omega is calculated based on rho and z values.
	- Otherwise, omega is calculated based only on z values.
"""
function calculateOmega(z_, rho_, type_)
    omega = similar(z_)

    if type_ == "rho"
        omega .= @. -1 / 2 * rho_^2 / z_^6
    else
        omega .= @. 1 / (4 * z_^4)
    end

    return omega
end

"""
    exportData(diffrentialSol, endTime)

	Export the solution data to a CSV file.
	- Renames columns for clarity.
	- Calculates omega values for rho and z.
	- Constructs a filename based on various parameters.
	- Writes the DataFrame to a CSV file.
"""
function exportData(diffrentialSol, endTime)
    df = DataFrame(diffrentialSol)
    rename!(df, :"value1" => "drho", :"value2" => "dz", :"value3" => "rho", :"value4" => "z")

    rho = diffrentialSol[3, :]
    z = diffrentialSol[4, :]
    omega_rho = calculateOmega(z, rho, "rho")
    omega_z = calculateOmega(z, rho, "z")

    df[!, :"omega_rho"] = omega_rho
    df[!, :"omega_z"] = omega_z

    round_endTime = round(endTime[2], digits=1)
    filename = "2D_export-eps$epsilon-epsphi$eps_phi-kappa$kappa-deltas$delta_star-beta$(round(rad2deg(beta_0)))-alpha$(round(rad2deg(alpha_0)))-theta$(round(rad2deg(theta_0)))-time$round_endTime.csv"

    CSV.write(filename, df)
    return nothing
end

"""
    CylindricalProblem!(ddu, du, u, p, t)

	Calculate the differential equations for the cylindrical problem in 2D.
	- Uses global variables for initial conditions.
	- Computes the differential equations based on the current state.
"""
function CylindricalProblem!(ddu, du, u, p, t)
    global z_0, rho_0, dphi0
    rho, z = u

    # ------------------------------ Exact Equations ----------------------------- #
    # Calculate the l0, frac0, and frac1 terms for efficiency
    #= l0 = epsilon * sin(alpha_0) * sin(beta_0) / sin(theta_0) * rho_0 - z_0 / sqrt(z_0^2 + rho_0^2)
    frac0 = (l0 + z / sqrt(z^2 + rho^2))
    frac1 = (rho^2 + z^2)^(3 / 2)

    # Update the differential equations
    ddu[1] = (1 / rho^3) * frac0^2 + z / (rho * frac1) * frac0
    ddu[2] = -(1 / frac1) * frac0 =#

    # --------------------------- Approximate Equation --------------------------- #
    eps_p = epsilon * sin(alpha_0) * sin(beta_0) / sin(theta_0)
    ddu[1] = 1 / (4 * rho^3) * (rho_0^4 * (1 + 2 * eps_p)^2 - (rho / z)^4)
    ddu[2] = -1 / (2 * z^3) * (rho_0^2 * (1 + eps_p) - (rho / z)^2)

    return nothing
end

"""
    SolvingtheProblem(CylindricalProblem, du0, u0, tspan)

	Solve the differential equations for the given initial conditions.
	- Uses the CylindricalProblem! function to compute the differential equations.
	- Extracts the solved position values.
	- Exports the solution data to a CSV file.
"""
function SolvingtheProblem(CylindricalProblem, du0, u0, tspan)
    problem = SecondOrderODEProblem(CylindricalProblem, du0, u0, tspan)
    sol = solve(problem, Feagin14(), reltol=1e-35, abstol=1e-40)

    # Extract the solved position values
    rho = sol[3, :]
    z = sol[4, :]

    # Export the solution data
    exportData(sol, tspan)

    return rho, z
end

"""
    ParticleMotion()

	Simulate the motion of a particle based on the initial conditions.
	- Loads the initial conditions.
	- Sets up the initial state and time span.
	- Solves the differential equations using the SolvingtheProblem function.
"""
function ParticleMotion()
    global alpha_0, theta_0, beta_0, phi_0, epsilon, eps_phi, delta_star, kappa, delta_star, rho_0, z_0, time_end, drho0, dz0, dphi0
    alpha_0, theta_0, beta_0, phi_0, epsilon, eps_phi, delta_star, kappa, delta_star, rho_0, z_0, time_end, drho0, dz0, dphi0 = load_initial_conditions()

    du0 = [drho0; dz0]
    u0 = [rho_0; z_0]
    tspan = (0.0, time_end)
    rho, z = SolvingtheProblem(CylindricalProblem!, du0, u0, tspan)

    return nothing
end

end