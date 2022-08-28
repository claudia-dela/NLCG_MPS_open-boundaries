% Input: 
% W - MPO representation of the Hamiltonian
% Alist - MPS tensors
% site - integer refferring to the site (in the chain)
% Output:
% contraction of A(site)-M-A*(site) for site from site to N 

function [hrj] = HRight(W,Alist,site)
	sites = size(Alist);
	N= sites(2);
	dims = size(Alist{2});
	d= dims(1);
	D= dims(2);

	ML=W{1};
	M=W{2};
	MR=W{3};
	ordM=size(M);
	mpo=ordM(1);
	
	Ai=Alist{site};
	AN=Alist{N};
	AiN = reshape(AN(:,:),[d,D]);
	CAiN=conj(AiN);
	list = {AiN,M,CAiN,MR};
	ind = {[1,-1],[-2,2,1,3],[3,-3],[2]};
	con = [2,1,3];
	finO = [-1,-2,-3];
	SP = reshape(ncon(list,ind,con,finO),[D,mpo,D]);
	    
	if site== N
		hrj = SP;
	else
		for i=1:(N-site) %compute Hleft
			SP = Pright(SP,Alist,M,ML,N-i); 
		end
	end
	hrj=SP;
end
