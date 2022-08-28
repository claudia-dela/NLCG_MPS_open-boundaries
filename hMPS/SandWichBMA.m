% Input: 
% HL, HR - left and right blocks
% M - MPO Hamiltonian (inner tensor)
% Alist - hMPS tensors
% Blist - structural tensor
% Output: 
% "sandwich" = contraction  HL(B-M-A*)M-HR

function [sandw] = SandWichBMA(HL,HR,M,Alist,Blist)
	dims = size(Blist{2});
	d= dims(1);
	D= dims(2);
	s= dims(4);

	A=Alist{2};
	B=Blist{2};
	Ai = reshape(A(:,:,:),[d,D,D]);
	CAi=conj(Ai);
	Bi = reshape(B(:,:,:,:),[d,D,D,s]);
	list = {HL,Bi,M,CAi,HR};
	index = {[1,2,3],[4,1,5,-1],[2,6,4,7],[7,3,8],[5,6,8]};
	con = [1,4,2,7,3,5,6,8];
	finO = [-1];
	sandw = ncon(list,index,con,finO);
end
