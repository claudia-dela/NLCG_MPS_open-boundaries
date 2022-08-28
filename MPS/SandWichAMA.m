% Input: 
% HL, HR - left and right blocks
% M - MPO Hamiltonian
% Alist - MPS tensors
% i - site
% Output: 
% "sandwich" = contraction  HL-(A(i)-M-A(i)*)-HR

function [sandw] = SandWichAMA(HL,HR,M,Alist,i)
	sites = size(Alist);
	N= sites(2);
	dims = size(Alist{2});
	d= dims(1);
	D= dims(2);

	A=Alist{i};
	Ai = reshape(A(:,:,:),[d,D,D]);
	CAi=conj(Ai);
	list = {HL,Ai,M,CAi,HR};
	index = {[1,2,3],[4,1,5],[2,6,4,7],[7,3,8],[5,6,8]};
	con = [1,4,2,7,3,5,6,8];
	finO = [];
	sandw = ncon(list,index,con,finO);
end
