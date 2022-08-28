% Input:
% N-  number of sites
% W - MPO Hamiltonian representation
% Alist - hMPS tensors
% Output: 
% MPO- MPS contraction

function [V] = MPOMPS(N,W,Alist)
    Nsiti=N-2; 
	dims = size(Alist{2});
	d= dims(1);
	D= dims(2);

	ML = W{1};
	M = W{2};
	MR = W{3};
	orderM =size(M);
	mpo = orderM(1);

	Al1=Alist{1};	
	AlN=Alist{3};
	
    Ai1 = Al1;
    list = {Ai1,ML};
    ind = {[-1],[-2]};
    con = [];
    finO = [-1,-2];
    Hl1 = reshape(ncon(list,ind,con,finO),[D,mpo]); 

	
	AiN = AlN;    
    list = {AiN,MR};
    ind = {[-1],[-2]};
    con = [];
    finO = [-1,-2];
    HrN = reshape(ncon(list,ind,con,finO),[D,mpo]);
    
	for i=1:Nsiti 
        Hl1=MAi(Hl1,W,Alist,i);
    end
    S=size(Hl1);
    dimen = S(3); % power of local dimension
    list = {Hl1,HrN};
    ind = {[1,2,-1],[1,2]};
    con = [1,2];
    finO = [-1]; 
    dx = ncon(list,ind,con,finO);

    V=dx;
	
end
