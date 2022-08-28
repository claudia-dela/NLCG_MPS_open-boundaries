% Input: 
% HL, HR - left and right blocks
% M - MPO Hamiltonian (inner tensor)
% Alist - hMPS tensors
% Output: 
% "sandwich" = contraction  HL-(A-M-A*)-HR

function [sandw] = SandWichAMA(HL,HR,M,Alist)
	dims = size(Alist{2});
	d= dims(1);
	D= dims(2);

	A=Alist{2};
	Ai = reshape(A(:,:,:),[d,D,D]);
	CAi=conj(Ai);
	list = {HL,Ai,M,CAi,HR};
	index = {[1,2,3],[4,1,5],[2,6,4,7],[7,3,8],[5,6,8]};
	con = [1,4,2,7,3,5,6,8];
	finO = [];
	sandw = ncon(list,index,con,finO);
end
