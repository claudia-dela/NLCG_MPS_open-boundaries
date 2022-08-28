% Input:
% N-  number of sites
% Al - left block
% Alist - hMPS tensors
% Output: 
% contraction between Al-A 

function [dx] = ContractAi(Al,Alist,N)
    Ai = Alist{2};
    dims = size(Ai);
	d= dims(1);
	D= dims(2);
        
if N==1
    list = {Al,Ai};
    ind = {[1],[-1,1,-2]};
    con = [1];
    finO = [-1,-2]; 
    dx = reshape(ncon(list,ind,con,finO),[d,D]);
    
else
    S = size(Al);
    dimen = S(1); 
    list = {Al,Ai};
    ind = {[-1,1],[-2,1,-3]};
    con = [1];
    finO = [-1,-2,-3]; 
    dx = reshape(ncon(list,ind,con,finO),[d*dimen,D]);
end
end
