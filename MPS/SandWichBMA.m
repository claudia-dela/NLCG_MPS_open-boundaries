% Input: 
% HL, HR - left and right blocks
% M - MPO Hamiltonian
% Alist - MPS tensors
% Blist - structural tensor
% i - site
% Output: 
% "sandwich" = contraction  HL-(B(i)-M-A(i)*)-HR

function [sandw] = SandWichBMA(HL,HR,M,Alist,Blist,i)
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
	list = {HL,Bi,M,CAi,HR};
	index = {[1,2,3],[4,1,5,-1],[2,6,4,7],[7,3,8],[5,6,8]};
	con = [1,4,2,7,3,5,6,8];
	finO = [-1];
	sandw = ncon(list,index,con,finO);
end
