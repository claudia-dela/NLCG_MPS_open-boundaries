% Input:
% N-  number of sites
% Hl1 - left block
% Alist - hMPS tensors
% W - MPO Hamiltonian representation
% mpo - bond dimension MPO
% Output: 
% contraction between Hl1-(A-M) 

function [dx] = MAi(Hl1,W,Alist,N)
%costruisco tensore i
    Ai = Alist{2};
    dims = size(Ai);
	d= dims(1);
	D= dims(2);

	M=W{2};
	orderM =size(M);
	mpo = orderM(1);

    
    list = {Ai,M};
    ind = {[1,-1,-2],[-3,-4,1,-5]};
    con = [1];
    finO = [-1,-2,-3,-4,-5]; %ultimo quello sotto libero
    AiTensor = reshape(ncon(list,ind,con,finO),[D,D,mpo,mpo,d]);
    
if N==1
    %contrai diretto
    list = {Hl1,AiTensor};
    ind = {[1,2],[1,-1,2,-2,-3]};
    con = [1,2];
    finO = [-1,-2,-3]; %ultimo quello sotto libero
    dx = reshape(ncon(list,ind,con,finO),[D,mpo,d]);
    
else
    S = size(Hl1);
    dimen = S(3); %potenza di 3
    list = {Hl1,AiTensor};
    ind = {[1,2,-1],[1,-2,2,-3,-4]};
    con = [1,2];
    finO = [-2,-3,-1,-4]; %ultimo quello sotto libero
    dx = reshape(ncon(list,ind,con,finO),[D,mpo,d*dimen]);
end
end
