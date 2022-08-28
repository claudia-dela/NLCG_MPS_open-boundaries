% Input: 
% N - number of sites
% HR - right block
% Alist - MPS tensors
% M - MPO Hamiltonian (inner tensor)
% ML - MPO Hamiltonian (left tensor)
% site - integer refferring to the site (in the chain)
% Output:
% contraction between A(site)-M-A*(site) and the right block HR

function [right] = Pright(HR,Alist,M,ML,site)
	sites = size(Alist);
	N= sites(2);
	dims = size(Alist{2});
	d= dims(1);
	D= dims(2);

	Asite=Alist{site};

	if site==1 
		Ai = reshape(Asite(:,:),[d,D]);
	    	CAi=conj(Ai);
	    	list = {ML,Ai,M,CAi,HR};
	    	ind = {[1],[2,3],[1,4,2,5],[5,6],[3,4,6]};
	    	con = [1,2,3,4,5,6];
	    	finO = [];
	    	right = ncon(list,ind,con,finO);
	else
	    	Ai = reshape(Asite(:,:,:),[d,D,D]);
	    	CAi=conj(Ai);
	    	list = {Ai,M,CAi,HR};
	    	ind = {[1,-1,2],[-2,3,1,4],[4,-3,5],[2,3,5]};
	    	con = [2,1,3,4,5];
	    	finO = [-1,-2,-3];
	    	right = ncon(list,ind,con,finO);
	end
end
