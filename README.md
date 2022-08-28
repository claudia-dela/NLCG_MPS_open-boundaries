# Variational NLCG on matrix product states with open boundary conditions.
	
The code in this repository is part of my PhD work. The code consists in a variational nonlinear conjugate gradient method, constrained to the variety of matrix product states with open boundary conditions. The method is applied to solve the expectation value minimization problem of the AKLT Hamiltonian. The vector which realizes the minimum value of the functional, which coincides with the lowest eigenvalue of the Hamiltonian, provide an approximation of the ground state. Moreover, a variation of the algorithm that modify the line search routine is implemented. The details are presented in the manuscript of my Phd thesis. 
	
The main file runs three algorithm: the NLCG on matrix product states with open boundary conditions with (1) standard restart, (2) reduced restart 
and (3) a variation of the NLCG (the details are presented in the manuscript of my Phd thesis). 
The code generates data files containing the functional values reached by the algorithms, other relevant informations (runtime, mean time of the line search, number of line searches performed) and it collects the respective vectors that realize the functional values. 
