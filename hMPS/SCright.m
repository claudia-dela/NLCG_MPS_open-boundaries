% Input: 
% N - number of sites
% Nr - right block
% Alist - hMPS tensors
% site - integer refferring to the site (in the chain)
% Output:
% contraction between A and the right block Nr

function [SP] = SCright(N,Nr,Alist,site) %P->C
	dims = size(Alist{2});
	d= dims(1);
	D= dims(2);
	A=Alist{2};

	if site==1
    		Ai=Alist{1};
    		Ai = reshape(Ai(:,:),[1,D]);
    		list = {Ai,Nr};
            %size(Nr)
    		ind = {[-1,1],[-2,1,-3]};
    		con = [1];
    		finO = [-1,-2,-3];
    		SP = (ncon(list,ind,con,finO));
	else
    	Ai = reshape(A(:,:,:),[d,D,D]);
		list = {Ai,Nr};
        S= size(Nr);
        dime = S(1); %potenza di 3
		ind = {[-1,-2,1],[-3,1]};
		con = [1];
		finO = [-1,-3,-2];
		SP = reshape(ncon(list,ind,con,finO),[d*dime,D]);
	end

end
