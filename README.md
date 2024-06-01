# Metropolis-Algorithm-for-the-Ising-Model
This repository showcases a detailed simulation of the Ising model using the Metropolis algorithm, implemented with an object-oriented programming (OOP) approach in C++. The simulation investigates ferromagnetism and physical observables as functions of temperature, providing insights into phase transitions and critical phenomena.

## Project Overview
### Object-Oriented Implementation in C++:
The simulation begins by instantiating an Ising object with a lattice attribute. This approach allows for a modular and scalable design, making the simulation easy to extend and maintain.
### Magnetization Analysis
We plotted the average magnetization over 2 million iterations at several temperatures, ranging from 0.5K past the Curie Temperature to 2.5K. The results clearly show the phase transition occurring at the Curie Temperature.
### Lattice Snapshots
Snapshots of the lattice were taken at various iteration values to visually depict the formation of magnetic domains within the lattice. These snapshots provide a clear representation of the microstates of the system during the simulation.
### Extended Temperature Range Simulation 
The temperature was increased from 0.5K to 4K, with the Metropolis algorithm applied 75,000 times per iteration for a 15x15 lattice. Physical observables such as energy per spin, magnetization per spin, specific heat capacity, and magnetic susceptibility were calculated.
### Repetition for Concreteness 
The simulation from 0.5K to 4K was repeated 200 times to ensure statistical significance and robustness of the results. This extensive repetition helps in reducing statistical noise and improving the reliability of the observed phenomena.
### Data Analysis and Visualization
The results were analyzed and graphed using Python and Jupyter Notebook. Error bars were calculated using the Median Absolute Deviation (MAD) values from the StatsModels library, providing a robust measure of variability.
### Oberving Phase Transition 
The phase transition was evident in the results, as shown by the critical behaviors in the specific heat capacity and magnetic susceptibility.

## Conclusion 
This simulation provides a comprehensive look at the behavior of the Ising model using the Metropolis algorithm. By leveraging an OOP approach in C++ and thorough data analysis in Python, we have created a robust framework for studying phase transitions and ferromagnetism. The results offer valuable insights and a clear visualization of the critical phenomena within the model.
