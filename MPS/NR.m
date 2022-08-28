% Input: 
% N - number of sites
% Alist - MPS tensors
% site - integer refferring to the site (in the chain)
% Output:
% contraction of the right side of the scalar product, untile site = site

function [nrj] = NR(Alist,site)
	sites = size(Alist);
	N= sites(2);
	dims = size(Alist{2});
	d= dims(1);
	D= dims(2);
	AN=Alist{N};
	A=Alist{site};

	AiN = reshape(AN(:,:),[d,D]);
	CAiN=conj(AiN);
	list = {AiN,CAiN};
	ind = {[1,-1],[1,-2]};
	con = [1];
	finO = [-1,-2];
	SPN = reshape(ncon(list,ind,con,finO),[D,D]);
	
	if site== N
		nrj = SPN;
	else
		for i=1:(N-site)
			SPN = SPright(SPN,Alist,N-i);   
		end
	end
	nrj = SPN;

end
