#include "Ising.h"

using namespace std; 

int main()
{
	/*{
	Ising obj(200, 0.1, 1, 1);
	obj.initialise_lattice();
	obj.lattice_to_text(0);
	for (int i = 0; i < 500001; ++i)
	{
		obj.flip_site();
		obj.calc_average_magnetization_per_spin(i, 0.01);
		if (i % 100 == 0)
		{
			cout << i << "\t" << obj.calc_average_magnetization_per_spin(i, 0.01) << "\n";
		}
		if (i % 25000 == 0)
		{
			obj.lattice_to_text(i);
		}
	}
	obj.mag_per_spin_to_vec();
	}*/
	{
		double temp = 0.1;
		for (int iteration = 0; iteration < 200; ++iteration)
		// We repeat each Temp value simulation 200 times, in an attempt to reduce statistical noise
		{


			//Write our Observables to text, just to see how they're looking, with iteration value in name for uniqueness
			std::string filename = "Lattice_Observables_at_T_" + std::to_string(temp) + "_iter_" + std::to_string(iteration) + ".txt";

			// Open a file to write the results
			std::ofstream outFile(filename);
			outFile << "Iteration : \t" << iteration << "\n";
			outFile << "Temperature \t" << "Energy \t" << "Magnetization \t" << "Cv \t" << "X\n";
			for (int i = 0; i < 40; ++i)
			// Went to T = 4K which is passed the Curie Temp of 2.27K
			{
				// For each Temp, initialise object with random lattice config. Then proceed with Metropolis Algorithm
				Ising obj(15, temp, 1, 1);
				obj.initialise_lattice();
				temp = 0.1 + i * 0.1;
				cout << "Iteration :  \t" << iteration << "\tTemperature : \t" << temp << "\n";
				for (int i = 0; i < 70000; ++i)
				// This value would be more accurate larger, but computational/time restraints impose limits
				{
					obj.flip_site();
					obj.calc_average_magnetization_per_spin(i, 0.01);
					obj.energy_mag_vals();
				}
				outFile << temp << "\t" << obj.get_average_energy() << "\t" << obj.get_average_magnetization() << "\t" << obj.calculate_specific_heat_capacity() << "\t" << obj.calculate_magnetic_suscepibility() << "\n";
			}
			outFile.close();
		}
	}
	return 0;
}