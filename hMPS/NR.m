% Input: 
% N - number of sites
% Alist - hMPS tensors
% site - integer refferring to the site (in the chain)
% Output:
% contraction of the right side of the scalar product, until site = site

function [nrj] = NR(N,Alist,site)
	dims = size(Alist{2});
	d= dims(1);
	D= dims(2);
	AN=Alist{3};
	A=Alist{2};

    	AiN = reshape(AN(:,:),[1,D]);
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
        		SPN = SPright(N,SPN,Alist,N-i);   
        	end
	end
	nrj = SPN;
end
