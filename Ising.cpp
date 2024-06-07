#include "Ising.h"

void Ising::initialise_lattice()
{
	std::vector <std::vector <int>> lattice(N);
	// We use the Mersenne Twister to randomly assign the lattice points
	std::random_device rd;
	std::mt19937 gen(rd());
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			std::uniform_int_distribution<> dis2(0, 1);
			int plus_minus = (dis2(gen) == 0) ? -1 : 1;
			this->lattice[i].push_back(plus_minus);
		}
	}
}
void Ising::display_lattice()
{
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			std::cout << lattice[i][j] << "\t";
		}
		std::cout << "\n";
	}
}
int Ising::calculate_hamiltonian()
{
	int energy = 0;
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			int up = ((i - 1) + N) % N ;
			int down = (i + 1) % N;
			int left = ((j - 1) + N) % N;
			int right = (j+1) % N;
			energy += -J * lattice[i][j] * (lattice[up][j] + lattice[down][j] + lattice[i][left] + lattice[i][right]);
		}
	}
	energy = energy / 2;
	energy_per_spin.push_back(energy);
	return energy;
}
void Ising::flip_site()
{
	// We use the Mersenne Twister to choose a random lattice point
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis1(0, N-1);
	int random_x = dis1(gen);
	int random_y = dis1(gen); 
	// get current energy before we flip the spin at the site
	int old_energy = calculate_hamiltonian();
	// flip the spin and check if the energy has reduced
	lattice[random_x][random_y] = -1 * lattice[random_x][random_y];
	int new_energy = calculate_hamiltonian();
	if (new_energy - old_energy < 0)
	{
		// If it has, then leave the site flipped
		return;
	}
	else
	{
		// If it hasn't, check if a random number alpha < our Metropolis condition
		std::random_device rd1;
		std::mt19937 gen1(rd1());
		std::uniform_int_distribution<> dis100(0, 10000);
		double check = dis100(gen1);
		check = check / 10000;
		if (check < exp(-(new_energy - old_energy)/T))
		{
			return;
		}
		else
		{
			lattice[random_x][random_y] = -1 * lattice[random_x][random_y];
			return;
		}
	}

}
double Ising::calc_average_magnetization_per_spin(int iter, double sample_rate)
{
	double magnetization_sum = 0;
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			magnetization_sum += lattice[i][j];
		}
	}
	int point_record = 1 / sample_rate;
	magnetization_sum = magnetization_sum / pow(N, 2);
	if (iter % point_record == 0)
	{
		magnetization_per_spin_average.push_back(magnetization_sum);
	}
	magnetization_per_spin_vector.push_back(magnetization_sum);
	return magnetization_sum;
}
void Ising::mag_per_spin_to_vec()
{
	std::string filename = "Magnetization_per_spin_T_" + std::to_string(T) + ".txt";
	std::ofstream outFile(filename);
	outFile << "X \t" << "Value \n";
	for (int i = 0; i < magnetization_per_spin_average.size(); ++i)
	{
		outFile << i << "\t" << magnetization_per_spin_average[i] << "\n";
	}
	outFile.close();
}
void Ising::lattice_to_text(int iter)
{
	std::string filename = "Lattice_at_" + std::to_string(iter) + "_iters_for_T_" + std::to_string(T) + ".txt";
	std::ofstream outFile(filename);
	outFile << "T : \t" << T << "\t Iter : \t" << iter << "\n";
	for (int i = 0; i < lattice.size(); ++i)
	{
		for (int j = 0; j < lattice.size(); ++j)
		{
			outFile << lattice[i][j] << "\t";
		}
		outFile << "\n";
	}
	outFile.close();
}
void Ising::energy_mag_vals()
{
	double val_e_1 = 0;
	double val_e_2 = 0;
	double val_m_1 = 0;
	double val_m_2 = 0;
	for (int i = 0; i < N; ++i)
	{
		for (int j =0; j < N; ++j)
		{
			//double energy = 0;
			int up = ((i - 1) + N) % N;
			int down = (i + 1) % N;
			int left = ((j - 1) + N) % N;
			int right = (j + 1) % N;
			val_e_1 += -J * lattice[i][j] * (lattice[up][j] + lattice[down][j] + lattice[i][left] + lattice[i][right]);
			val_e_2 += -J * lattice[i][j] * (lattice[up][j] + lattice[down][j] + lattice[i][left] + lattice[i][right]);
			val_m_1 += lattice[i][j];
			val_m_2 += lattice[i][j];
		}
	}
	val_e_1 = val_e_1 / pow(N, 2);
	val_e_2 = val_e_2 / pow(N, 2);
	val_m_1 = val_m_1 / pow(N, 2);
	val_m_2 = val_m_2 / pow(N, 2);

	val_e_1 = val_e_1 / (2);
	val_e_2 = val_e_2 / (2);

	val_e_2 = pow(val_e_2, 2);
	val_m_2 = pow(val_m_2, 2);

	_E2_.push_back(val_e_2);
	_E_2.push_back(val_e_1);
	_M2_.push_back(val_m_2);
	_M_2.push_back(val_m_1);
}
double Ising::calculate_specific_heat_capacity()
{
	double prefix = 1 / (Kb * T * T);
	double sum_term_1 = 0;
	double sum_term_2 = 0;
	for (int i = floor((_E2_.size() * 0.90)); i < _E2_.size(); ++i)
	{
		sum_term_1 += _E2_[i];
		sum_term_2 += _E_2[i];
	}
	double denominator = int(_E2_.size() - floor((_E2_.size() * 0.90)));
	sum_term_1 = sum_term_1 / denominator;
	sum_term_2 = sum_term_2 / denominator;
	double Cv = prefix * (sum_term_1 - pow(sum_term_2,2));
	specific_heat_capacity_vector.push_back(Cv);
	return Cv;
}
double Ising::calculate_magnetic_suscepibility()
{
	double prefix = 1 / (Kb * T * N * N);
	double sum_term_1 = 0;
	double sum_term_2 = 0;
	for (int i = floor((_M2_.size()*0.90)); i < _M2_.size(); ++i)
	{
		sum_term_1 += _M2_[i];
		sum_term_2 += _M_2[i];
	}
	double denominator = int(_E2_.size() - floor((_E2_.size() * 0.90)));
	sum_term_1 = sum_term_1 / denominator;
	sum_term_2 = sum_term_2 / denominator;
	double chi = prefix * (sum_term_1 - pow(sum_term_2, 2));
	magnetic_susceptibility.push_back(chi);
	return chi;
}
double Ising::get_average_energy()
{
	double energy_sum = 0;
	for (int i = floor(_E_2.size()*0.90); i < _E_2.size(); ++i )
	{
		energy_sum += _E_2[i];
	}
	double energy_denominator = int(_E_2.size() - floor(_E_2.size() * 0.90));
	energy_sum = energy_sum / energy_denominator;
	return energy_sum;
}
double Ising::get_average_magnetization()
{
	double magnetization_sum = 0;
	for (int i = floor(_M_2.size() * 0.90); i < _M_2.size(); ++i)
	{
		magnetization_sum += _M_2[i];
	}
	double magnetization_denominator = int(_M_2.size() - floor(_M_2.size() * 0.90));
	magnetization_sum = magnetization_sum / magnetization_denominator;
	return magnetization_sum;
}