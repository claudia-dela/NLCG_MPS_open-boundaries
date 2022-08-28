% Input: 
% N - number of sites
% Alist - MPS tensors
% site - integer refferring to the site (in the chain)
% Output:
% contraction of the left side of the scalar product, untile site = site

function [nlj] = NL(Alist,site)
	sites = size(Alist);
	N= sites(2);
	dims = size(Alist{2});
	d= dims(1);
	D= dims(2);
	A1=Alist{1};
	A=Alist{site};

	Ai1 = reshape(A1(:,:),[d,D]);
	CAi1=conj(Ai1);
	list = {Ai1,CAi1};
	ind = {[1,-1],[1,-2]};
	con = [1];
	finO = [-1,-2];
	SP1 = reshape(ncon(list,ind,con,finO),[D,D]);
	
	if site== 1
		nlj = SP1;
	else
		for i=2:(site)
			SP1 = SPleft(SP1,Alist,i);   
		end
	end
	nlj = SP1;
end
