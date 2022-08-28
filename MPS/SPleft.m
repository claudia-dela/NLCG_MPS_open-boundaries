% Input: 
% N - number of sites
% Nl - right block
% Alist - MPS tensors
% site - integer refferring to the site (in the chain)
% Output:
% contraction between the left block Nl and A(site)-A(site)*

function [SP] = SPleft(Nl,Alist,site)
	sites = size(Alist);
	N= sites(2);
	dims = size(Alist{2});
	d= dims(1);
	D= dims(2);
	
	A=Alist{site};

	if site==N
		Ai = reshape(A(:,:),[d,D]);
	    	CAi=conj(Ai);
	    	list = {Nl,Ai,CAi};
	    	ind = {[1,2],[3,1],[3,2]};
	    	con = [1,3,2];
	    	finO = [];
	    	SP = ncon(list,ind,con,finO);
	else
	    
	Ai = reshape(A(:,:,:),[d,D,D]);
	CAi=conj(Ai);
	list = {Nl,Ai,CAi};
	ind = {[1,2],[3,1,-1],[3,2,-2]};
	con = [1,3,2];
	finO = [-1,-2];
	SP = reshape(ncon(list,ind,con,finO),[D,D]);
	end
end
