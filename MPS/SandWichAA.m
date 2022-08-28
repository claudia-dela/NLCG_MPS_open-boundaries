% Input: 
% HL, HR - left and right blocks
% M - MPO Hamiltonian
% Alist - MPS tensors
% Blist - structural tensor
% i - site                               	                                  	 
% Output: 
% "sandwich" = contraction  HL-(A(i)-A(i)*)-HR

function [sandw] = SandWichAA(HL,HR,Alist,i)
	sites = size(Alist);
	N= sites(2);
	dims = size(Alist{2});
	d= dims(1);
	D= dims(2);

	A=Alist{i};
	Ai = reshape(A(:,:,:),[d,D,D]);
	CAi=conj(Ai);

	list = {HL,Ai,CAi,HR};
	index = {[1,2],[3,1,4],[3,2,5],[4,5]};
	con = [1,3,2,4,5];
	finO = [];
	sandw = ncon(list,index,con,finO);
end
