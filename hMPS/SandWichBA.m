% Input: 
% HL, HR - left and right blocks
% Alist - hMPS tensors
% Blist - structural tensor
% Output: 
% "sandwich" = contraction  HL-(B-A*)-HR

function [sandw] = SandWichBA(HL,HR,Alist,Blist)
	dims = size(Blist{2});
	d= dims(1);
	D= dims(2);

	if length(dims)==4
		s= dims(4);
	else
    		s=1;
	end
	A=Alist{2};
	B=Blist{2};

	Ai = reshape(A(:,:,:),[d,D,D]);
	CAi=conj(Ai);
	Bi = reshape(B(:,:,:,:),[d,D,D,s]);

	list = {HL,Bi,CAi,HR};
	index = {[1,2],[3,1,4,-1],[3,2,5],[4,5]};
	con = [1,3,2,4,5];
	finO = [-1];
	sandw = ncon(list,index,con,finO);

end
