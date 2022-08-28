% Input:
% N - number of site
% A - MPS tensors
% W - MPO Hamiltonian
% Output:
% expectation value of the Hamiltonian on the image of A (via the MPS map)
function [ev] = ExpVal(W,A)
	M=W{2};
	sites = size(A);
	N= sites(2);
	j=1;
	if N==2
		E = ncon({HLeft(W,A,j),HRight(W,A,j+1)},{[1,2,3],[1,2,3]},[1,2,3],[]);
	    	Norm = ncon({NL(A,j),NR(A,j+1)},{[1,2],[1,2]},[1,2],[]);
	    	ev= E/Norm;
	else
	    	E = SandWichAMA(HLeft(W,A,j),HRight(W,A,j+2),M,A,j+1);
	    	Norm = SandWichAA(NL(A,j),NR(A,j+2),A,j+1);
	    	ev= E/Norm;
	end
end
