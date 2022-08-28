% Input: 
% N - number of sites
% Nl - right block
% Alist - hMPS tensors
% site - integer refferring to the site (in the chain)
% Output:
% contraction between the left block Nl and A

function [SP] = SCleft(N, Nl,Alist,site)
	dims = size(Alist{2});
	d= dims(1);
	D= dims(2);
	A=Alist{2};

	if site==N
		Ai=Alist{3};
		Ai = reshape(Ai(:,:),[1,D]);
		list = {Nl,Ai};
		ind = {[-1,-2,1],[-3,1]};
		con = [1];
		finO = [-1,-2,-3];
		SP = ncon(list,ind,con,finO);
	else
		    
		Ai = reshape(A(:,:,:),[d,D,D]);
		list = {Nl,Ai};
        S= size(Nl);
        dime = S(1); %potenza di 3
		ind = {[-1,1],[-2,1,-3]};
		con = [1];
		finO = [-1,-2,-3];
		SP = reshape(ncon(list,ind,con,finO),[d*dime,D]);
	end
end
