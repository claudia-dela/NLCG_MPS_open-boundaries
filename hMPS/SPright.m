% Input: 
% N - number of sites
% Nr - right block
% Alist - hMPS tensors
% site - integer refferring to the site (in the chain)
% Output:
% contraction between A-A* and the right block Nr

function [SP] = SPright(N,Nr,Alist,site)
	dims = size(Alist{2});
	d= dims(1);
	D= dims(2);
	A=Alist{2};

	if site==1
    		Ai=Alist{1};
    		Ai = reshape(Ai(:,:),[1,D]);
    		CAi=conj(Ai);
    		list = {Ai,CAi,Nr};
    		ind = {[1,2],[1,3],[2,3]};
    		con = [1,3,2];
    		finO = [];
    		SP = (ncon(list,ind,con,finO));
	else
    		Ai = reshape(A(:,:,:),[d,D,D]);
		CAi=conj(Ai);
		list = {Ai,CAi,Nr};
		ind = {[1,-1,2],[1,-2,3],[2,3]};
		con = [1,3,2];
		finO = [-1,-2];
		SP = reshape(ncon(list,ind,con,finO),[D,D]);
	end

end
