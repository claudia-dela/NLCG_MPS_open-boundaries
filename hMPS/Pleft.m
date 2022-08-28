% Input: 
% N - number of sites
% HL - left block
% Alist - hMPS tensors
% M - MPO Hamiltonian (inner tensor)
% MR - MPO Hamiltonian (right tensor)
% site - integer refferring to the site (in the chain)
% Output:
% contraction between the left block HL and A-M-A*

function [left] = Pleft(N,HL,Alist,M,MR,site)
	dims = size(Alist{2});
	d= dims(1);
	D= dims(2);

	Asite=Alist{2};

	if site==N % the site N is strange
	    	Ai=Alist{3};
	    	CAi=conj(Ai);
	    	list = {HL,Ai,MR,CAi};
	    	ind = {[1,2,3],[1],[2],[3]};
	    	con = [1,2,3];
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
