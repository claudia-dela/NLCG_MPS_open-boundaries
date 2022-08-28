% Input: 
% W - MPO representation of the Hamiltonian
% Alist - MPS tensors
% site - integer refferring to the site (in the chain)
% Output:
% contraction of A(site)-M-A*(site) for site from 1 to site 

function [hlj] = HLeft(W,Alist,site)
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
	A1=Alist{1};
	Ai1 = reshape(A1(:,:),[d,D]);
	CAi1=conj(Ai1);
	list = {ML,Ai1,M,CAi1};
	ind = {[1],[2,-1],[1,-2,2,3],[3,-3]};
	con = [1,2,3];
	finO = [-1,-2,-3];
	SP = reshape(ncon(list,ind,con,finO),[D,mpo,D]);
	if site== 1
		hlj = SP;
	else
		for i=2:(site) 
			SP = Pleft(SP,Alist,M,MR,i); 
		end
	end
	hlj=SP;
end
