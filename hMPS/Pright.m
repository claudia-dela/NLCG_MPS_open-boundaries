% Input: 
% N - number of sites
% HR - right block
% Alist - hMPS tensors
% M - MPO Hamiltonian (inner tensor)
% ML - MPO Hamiltonian (left tensor)
% site - integer refferring to the site (in the chain)
% Output:
% contraction between A-M-A* and the right block HR

function [right] = Pright(N,HR,Alist,M,ML,site)
	dims = size(Alist{2});
	d= dims(1);
	D= dims(2);

	Asite=Alist{2};

	if site==1 % the site N is strange
    		Ai= Alist{1};
    		CAi=conj(Ai);
    		list = {Ai,ML,CAi,HR};
    		ind = {[1],[2],[3],[1,2,3]};
    		con = [1,2,3];
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
