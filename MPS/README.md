# Variational NLCG on MPS with open boundary conditions.
The main file runs three algorithm: the NLCG on general MPS with open boundary conditions with (1) standard restart, (2) reduced restart 
and (3) a variation of the NLCG (the details are presented in the manuscript of my Phd thesis). 
The code generates data files containing the functional values reached by the algorithms, other relevant informations (runtime, mean time of the line search, number of line searches performed) and it collects the respective vectors that realize the functional values. 

The code makes use of the general tensor network contraction routine ncon() available at https://arxiv.org/abs/1402.0939
