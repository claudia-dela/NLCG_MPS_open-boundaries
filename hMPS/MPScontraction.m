% Input:
% N-  number of sites
% W - MPO Hamiltonian
% Alist - hMPS tensors
% Blist - structural tensor (reshape of the identity matrix)
% Output: 
% MPS contraction: vector in C^{n^d} 

function [V] = MPScontraction(N,Alist)
    Nsiti=N-2;
	Ai1=Alist{1};	
	AiN=Alist{3};
	    
	for i=1:Nsiti 
        Ai1=ContractAi(Ai1,Alist,i);
    end
    list = {Ai1,AiN};
    ind = {[-1,1],[1]};
    con = [1];
    finO = [-1]; 
    dx = ncon(list,ind,con,finO);

    V=dx;
 	
end
