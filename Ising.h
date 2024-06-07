#include <random>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#pragma once
class Ising
{
private:
	int N; //Lattice Dimensions
	int num_steps; // Number of Steps to iterate for Monte Carlo Simulation
	double T; // Temperature of Ferromagnet
	double J; //Coupling Constant
	double Kb; //Boltzmann Constant
	long double current_energy; // variables for holding the current energy and magnetization of the lattice 
	long double current_mag;
	std::vector <std::vector <int>> lattice; // Lattice Vector of Vectors array
	std::vector <double> magnetization_per_spin_average; // Average Magnetization per Spin Vector
	std::vector <double> energy_per_spin; // Vector to store Energy per Site 
	std::vector <double> magnetization_per_spin_vector; // Vector to store Mag per spin at each iteration
	std::vector <double> _E2_; // Vectors for Energy per Spin variance
	std::vector <double> _E_2;
	std::vector <double> _M2_; // Vectors for Magnetization per Spin variance
	std::vector <double> _M_2;
	std::vector <double> specific_heat_capacity_vector; // Vector to store Cv at each iteration
	std::vector <double> magnetic_susceptibility; // Vector to store MS
public:
	Ising(int dimension, double temp, double coupling_const, double boltzmann)
	{
		N = dimension;
		T = temp; 
		J = coupling_const;
		Kb = boltzmann;
		lattice.resize(N);
	} // Contstructor
	void initialise_lattice(); // Randomly initialize lattice spins
	void display_lattice(); // Show the lattice, helper method intially just to ensure above method was correctly implemented
	int calculate_hamiltonian(); // Calculate the Energy of the Lattice
	void flip_site(); // Implement the Metropolis Step
	double calc_average_magnetization_per_spin(int iter, double sample_rate); // get <M> for the lattice
	void lattice_to_text(int iter); // Print the lattice to text for visualization of magnetization domains
	void mag_per_spin_to_vec(); // Write <M> for each step of the metropolis algorithm to a text file for analysis in Python
	void energy_mag_vals(); // At each Temperature step of the Algorithm, the Energy and Magnetization must be calculated per iteration step  
	double calculate_specific_heat_capacity(); // Specific Heat Capacity Method to get <Cv> at each Temp value
	double calculate_magnetic_suscepibility(); // Magnetic Susceptibility Method to get <Xm> at each Temp value
	double get_average_energy(); // Methods to get <E> and <M> after the system has reached equilibrium
	double get_average_magnetization();
};

