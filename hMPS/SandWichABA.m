% Input: 
% HL, HR - left and right blocks
% Alist - hMPS tensors
% Blist - structural tensor
% Output: 
% "sandwich" = contraction  HL-(B-A*)-HR

function [sandw] = SandWichABA(HL,HR,Alist,Blist) %era BA
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
	Bi = reshape(B(:,:,:,:),[d,D,D,s]);

	list = {HL,Bi,HR};
    S1= size(HL);
    dime1 = S1(1); %potenza di 3
    S2= size(HR);
    dime2 = S2(1); %potenza di 3
	
    index = {[-1,1],[-2,1,2,-3],[-4,2]};
	con = [1,2];
	finO = [-1,-2,-4,-3];
	sandw =reshape(ncon(list,index,con,finO),[d*dime1*dime2,s]);

end
