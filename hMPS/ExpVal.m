% Input:
% N - number of site
% A - hMPS tensors
% W - MPO Hamiltonian
% Output:
% expectation value of the Hamiltonian on the image of A (via the hMPS map)
function [ev] = ExpVal(N,W,A)
	M=W{2};
	j=1;
	if N==2
		E = ncon({HLeft(N,W,A,j),HRight(N,W,A,j+1)},{[1,2,3],[1,2,3]},[1,2,3],[]);
    		Norm = ncon({NL(N,A,j),NR(N,A,j+1)},{[1,2],[1,2]},[1,2],[]);
    		ev= E/Norm;
	else
    		E = SandWichAMA(HLeft(N,W,A,1),HRight(N,W,A,3),M,A);
    		Norm = SandWichAA(NL(N,A,1),NR(N,A,3),A);
    		ev= E/Norm;
	end
end
