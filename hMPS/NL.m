% Input: 
% N - number of sites
% Alist - hMPS tensors
% site - integer refferring to the site (in the chain)
% Output:
% contraction of the left side of the scalar product, until site = site

function [nlj] = NL(N,Alist,site)
	dims = size(Alist{2});
	d= dims(1);
	D= dims(2);
	A1=Alist{1};
	A=Alist{2};

    	Ai1 = reshape(A1(:,:),[1,D]);
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
                	SP1 = SPleft(N,SP1,Alist,i);   
        	end
	end
	nlj = SP1;
end
