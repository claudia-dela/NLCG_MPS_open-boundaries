% Input:
% N - number of sites 
% W - MPO representation of the Hamiltonian
% Alist - hMPS tensors
% site - integer refferring to the site (in the chain)
% Output:
% contraction of A-M-A* for site from site to N 

function [hrj] = HRight(N,W,Alist,site)
	dims = size(Alist{2});
	d= dims(1);
	D= dims(2);

	ML = W{1};
	M = W{2};
	MR = W{3};

	Ai=Alist{2};
	orderM = size(M);
	mpo = orderM(1);

	AN=Alist{3};
	AiN = AN;    
	CAiN=conj(AiN);
    	list = {AiN,MR,CAiN};
    	ind = {[-1],[-2],[-3]};
    	con = [];
    	finO = [-1,-2,-3];
    	SP = reshape(ncon(list,ind,con,finO),[D,mpo,D]);
    
	if site== N
    		hrj = SP;
    	else
        	for i=1:(N-site) %compute Hleft
            		SP = Pright(N,SP,Alist,M,ML,N-i); 
        	end
	end
	hrj=SP;
end
