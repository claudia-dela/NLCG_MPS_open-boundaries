% Input: 
% N - number of sites
% HL - left block
% Alist - MPS tensors
% M - MPO Hamiltonian (inner tensor)
% MR - MPO Hamiltonian (right tensor)
% site - integer refferring to the site (in the chain)
% Output:
% contraction between the left block HL and A(site)-M-A*(site)

function [left] = Pleft(HL,Alist,M,MR,site)
	sites = size(Alist);
	N= sites(2);
	dims = size(Alist{2});
	d= dims(1);
	D= dims(2);

	Asite=Alist{site};

	if site==N 
		Ai = reshape(Asite(:,:),[d,D]);
	    	CAi=conj(Ai);
	    	list = {HL,Ai,M,CAi,MR};
	    	ind = {[1,2,3],[4,1],[2,5,4,6],[6,3],[5]};
	    	con = [1,4,2,6,3,5];
	    	finO = [];
	    	left = ncon(list,ind,con,finO);
	else
	    	Ai = reshape(Asite(:,:,:),[d,D,D]);
	    	CAi=conj(Ai);
	    	list = {HL,Ai,M,CAi};
	    	ind = {[1,2,3],[4,1,-1],[2,-2,4,5],[5,3,-3]};
	    	con = [1,4,2,5,3];
	    	finO = [-1,-2,-3];
	    	left = ncon(list,ind,con,finO);
	end
end
