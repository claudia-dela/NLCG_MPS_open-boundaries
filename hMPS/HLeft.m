% Input:
% N - number of sites 
% W - MPO representation of the Hamiltonian
% Alist - hMPS tensors
% site - integer refferring to the site (in the chain)
% Output:
% contraction of A-M-A* for site from 1 to site 

function [hlj] = HLeft(N,W,Alist,site)
	dims = size(Alist{2});
	d= dims(1);
	D= dims(2);

	ML = W{1};
	M = W{2};
	MR = W{3};
	
	Ai=Alist{2};
	orderM =size(M);
	mpo = orderM(1);

	A1=Alist{1};
	Ai1 = A1;%reshape(A1(:,:),[D,1]);
	CAi1=conj(Ai1);
    	list = {Ai1,ML,CAi1};
    	ind = {[-1],[-2],[-3]};
    	con = [];
    	finO = [-1,-2,-3];
    	SP = reshape(ncon(list,ind,con,finO),[D,mpo,D]);
	if site== 1
    		hlj = SP;
    	else
        	for i=2:(site) %compute Hleft
            		SP = Pleft(N,SP,Alist,M,MR,i); 
        	end
	end
	hlj=SP;
end
