% Input: 
% N - number of sites
% Nl - right block
% Alist - hMPS tensors
% site - integer refferring to the site (in the chain)
% Output:
% contraction between the left block Nl and A-A*

function [SP] = SPleft(N, Nl,Alist,site)
	dims = size(Alist{2});
	d= dims(1);
	D= dims(2);
	A=Alist{2};

	if site==N
		Ai=Alist{3};
		Ai = reshape(Ai(:,:),[1,D]);
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
