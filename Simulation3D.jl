module Simulation3D

# Import necessary modules and functions
include("SharedFunctions.jl")
using .SharedFunctions
using DifferentialEquations
using ArbNumerics
using Distributed
using DataFrames
using CSV

# Declare global variables without initializing them
global alpha_0, theta_0, beta_0, phi_0, epsilon, eps_phi, delta_star, kappa, rho_0, z_0, time_end, drho0, dz0, dphi0

export CylindricalProblem!, SolvingtheProblem, ParticleMotion, Converting

# Initialization function to load initial conditions
function __init__()
	global alpha_0, theta_0, beta_0, phi_0, epsilon, eps_phi, delta_star, kappa, rho_0, z_0, time_end, drho0, dz0, dphi0
	alpha_0, theta_0, beta_0, phi_0, epsilon, eps_phi, delta_star, kappa, rho_0, z_0, time_end, drho0, dz0, dphi0 = load_initial_conditions()
end

"""
	CylindricalProblem!(ddu, du, u, p, t)

    Calculate the differential equations for the cylindrical problem in 3D.
    - Uses global variables for constants.
    - Computes the differential equations based on the current state.
"""
function CylindricalProblem!(ddu, du, u, p, t)
    global kappa, delta_star, eps_phi
	rho, phi, z = u
	drho, dphi, dz = du

	inv_rho2_z2 = 1 / (rho^2 + z^2)
	log_term = log(inv_rho2_z2)

	dPhi_dR = rho * (z^2 * (1 - 2 * kappa) - 2 * rho^2 * kappa + 2 * kappa * (rho^2 - z^2 * log_term) * delta_star^2) * inv_rho2_z2^2
	dPhi_dz = -z * (2 * z^2 * kappa + rho^2 * (1 + 2 * kappa) - 2 * rho^2 * kappa * (1 + log_term) * delta_star^2) * inv_rho2_z2^2

	Brho = rho * inv_rho2_z2^(3 / 2)
	Bz = z * inv_rho2_z2^(3 / 2)

	ddu[1] = rho * dphi^2 + rho * dphi * Bz - eps_phi * dPhi_dR
	ddu[2] = (dz * Brho - drho * Bz - 2 * drho * dphi) / rho
	ddu[3] = -rho * dphi * Brho - eps_phi * dPhi_dz
end

"""
	SolvingtheProblem(CylindricalProblem!, du0, u0, tspan)

    Solve the differential equations for the given initial conditions.
    - Uses the CylindricalProblem! function to compute the differential equations.
    - Extracts the solved position values.
    - Exports the solution data to a CSV file.
"""
function SolvingtheProblem(CylindricalProblem!, du0, u0, tspan)
	problem = SecondOrderODEProblem(CylindricalProblem!, du0, u0, tspan)
	sol = solve(problem, Feagin14(), reltol = 1e-35, abstol = 1e-40)
	exportData(sol, tspan)
	return sol[4, :], sol[5, :], sol[6, :]
end

"""
	Converting(rho_, phi_, z_)

    Convert cylindrical coordinates to cartesian coordinates.
    - rho, phi, and z are the cylindrical coordinates.
    - Returns x, y, and z as cartesian coordinates.
"""
function Converting(rho_, phi_, z_)
	x = rho_ .* cos.(phi_)
	y = rho_ .* sin.(phi_)
	return x, y, z_
end

"""
	exportData(diffrentialSol, endTime)

    Export the solution data to a CSV file.
    - Renames columns for clarity.
    - Calculates magnetic and electric fields for rho and z.
    - Constructs a filename based on various parameters.
    - Writes the DataFrame to a CSV file.
"""
function exportData(diffrentialSol, endTime)
	df = DataFrame(diffrentialSol)
	rho = diffrentialSol[4, :]
	z = diffrentialSol[6, :]
	MagRho, MagZ = calculateMagnetigField(z, rho)
	ERho, EZ = calculateElectricField(z, rho)
	df[!, :"Magnetic_rho"] = MagRho
	df[!, :"Magnetic_z"] = MagZ
	df[!, :"Electric_rho"] = ERho
	df[!, :"Electric_z"] = EZ

	filename = "3D_export-eps$(epsilon)-epsphi$(eps_phi)-kappa$(kappa)-deltas$(delta_star)-beta$(round(rad2deg(beta_0)))-alpha$(round(rad2deg(alpha_0)))-theta$(round(rad2deg(theta_0)))-time$(round(endTime[2], digits=1)).csv"
	CSV.write(filename, df)
end

"""
	calculateMagnetigField(z_, rho_)

    Calculate the magnetic field for the given z and rho values.
    - Returns magnetic field components for rho and z.
"""
function calculateMagnetigField(z_, rho_)
	inv_rho2_z2 = 1 ./ (rho_ .^ 2 + z_ .^ 2)
	Brho_ = rho_ .* inv_rho2_z2 .^ (3 / 2)
	Bz_ = z_ .* inv_rho2_z2 .^ (3 / 2)
	return Brho_, Bz_
end

"""
	calculateElectricField(z_, rho_)

    Calculate the electric field for the given z and rho values.
    - Uses global variables for constants.
    - Returns electric field components for rho and z.
"""
function calculateElectricField(z, rho)
	global kappa, delta_star

	@assert length(z) == length(rho) "Vectors z and rho must have the same length."

	inv_rho2_z2 = 1.0 ./ (rho .^ 2 + z .^ 2)
	log_term = log.(inv_rho2_z2)

	dPhi_dR = rho .* (z .^ 2 .* (1 .- 2 * kappa) - 2 * rho .^ 2 * kappa + 2 * kappa * (rho .^ 2 - z .^ 2 .* log_term) * delta_star^2) .* inv_rho2_z2 .^ 2
	dPhi_dz = -z .* (2 * z .^ 2 * kappa + rho .^ 2 * (1 .+ 2 * kappa) - 2 * kappa * rho .^ 2 .* (1 .+ log_term) .* delta_star^2) .* inv_rho2_z2 .^ 2

	return dPhi_dR, dPhi_dz
end

"""
	ParticleMotion()

    Simulate the motion of a particle in a magnetic field based on the initial conditions.
    - Loads the initial conditions.
    - Sets up the initial state and time span.
    - Solves the differential equations using the SolvingtheProblem function.
"""
function ParticleMotion()
    global alpha_0, theta_0, beta_0, phi_0, epsilon, eps_phi, delta_star, kappa, delta_star, rho_0, z_0, time_end, drho0, dz0, dphi0
    alpha_0, theta_0, beta_0, phi_0, epsilon, eps_phi, delta_star, kappa, delta_star, rho_0, z_0, time_end, drho0, dz0, dphi0 = load_initial_conditions()

	du0 = [drho0; dphi0; dz0]
	u0 = [rho_0; phi_0; z_0]
	tspan = (0.0, time_end)
	rho, phi, z = SolvingtheProblem(CylindricalProblem!, du0, u0, tspan)
end

end
