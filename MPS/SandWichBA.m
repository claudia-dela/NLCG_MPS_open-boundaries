% Input: 
% HL, HR - left and right blocks
% M - MPO Hamiltonian
% Alist - MPS tensors
% Blist - structural tensor                        	    
% i - site                              	 
% Output: 
% "sandwich" = contraction  HL-(B(i)-A(i)*)-HR

function [sandw] = SandWichBA(HL,HR,Alist,Blist,i)
	sites = size(Blist);
	N= sites(2);
	dims = size(Blist{2});
	d= dims(1);
	D= dims(2);
	s= dims(4);

	A=Alist{i};
	B=Blist{i};

	Ai = reshape(A(:,:,:),[d,D,D]);
	CAi=conj(Ai);
	Bi = reshape(B(:,:,:,:),[d,D,D,s]);
	
	list = {HL,Bi,CAi,HR};
	index = {[1,2],[3,1,4,-1],[3,2,5],[4,5]};
	con = [1,3,2,4,5];
	finO = [-1];
	sandw = ncon(list,index,con,finO);
	
end
