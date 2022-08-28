% Input:
% N-  number of sites
% W - MPO Hamiltonian
% Alist - hMPS tensors
% Blist - structural tensor
% Output: 
% gradient of the expectation value formula in Figure 
% (5.5) for the standard algoritm
% (6.4) for the variation

function [gradi,Hdx,Ndx] = MPSgrad(N,W,Alist,Blist)
	dims = size(Blist{2});
	d= dims(1);
	D= dims(2);
	s= dims(4);

	ML = W{1};
	M = W{2};
	MR = W{3};
	
    if ~isreal(Alist{2})
        %complex case
        Bconj={conj(Blist{1}),conj(Blist{2}),conj(Blist{3})};
        Blist=Bconj;
        % complex case: partial derivative on conjugated variables
        Aconj={conj(Alist{1}),conj(Alist{2}),conj(Alist{3})};
        Alist=Aconj;
    end
    
	if N==2
		E = ncon({HLeft(N,W,Alist,1),HRight(N,W,Alist,2)},{[1,2,3],[1,2,3]},[1,2,3],[]);
	    	Norm = ncon({NL(N,Alist,1),NR(N,Alist,2)},{[1,2],[1,2]},[1,2],[]);
	else
		E = SandWichAMA(HLeft(N,W,Alist,1),HRight(N,W,Alist,3),M,Alist);
	    	Norm = SandWichAA(NL(N,Alist,1),NR(N,Alist,3),Alist);
	end

	gradi = zeros(1,s);

	Al1=Alist{1};
	Bl1=Blist{1};
	A1 = Al1;
	CA1=conj(A1);
	B1 = reshape(Bl1(:,:,:,:),[D,s]);
	
	AlN=Alist{3};
	BlN=Blist{3};
	AN = AlN;
	CAN=conj(AN);
	BN = reshape(BlN(:,:,:,:),[D,s]);
	
	Hl = HLeft(N,W,Alist,1); %first site
	Nl = NL(N,Alist,1);

	Hr1 = HRight(N,W,Alist,N);
	Nr1=NR(N,Alist,N);

	Hdx={Hr1};
	Ndx={Nr1}; %cells of tensors to store

	for i=2:(N-1) 
		Hr1=Pright(N,Hr1,Alist,M,MR,N-i+1);
	    	Nr1=SPright(N,Nr1,Alist,N-i+1);
	    	Hdx{end+1}= Hr1;
	    	Ndx{end+1}= Nr1;
	end

	for i=1:N
        if i == 1
			di = ncon({B1,CA1,Ndx{N-1}},{[1,-1],[2],[1,2]},[1,2],[-1]).';
			ci = ncon({B1,ML,CA1,Hdx{N-1}},{[1,-1],[2],[3],[1,2,3]},[1,2,3],[-1]).';
        elseif i==N
			di = ncon({BN,CAN,Nl},{[1,-1],[2],[1,2]},[1,2],[-1]).';
			ci = ncon({BN,MR,CAN,Hl},{[1,-1],[2],[3],[1,2,3]},[1,2,3],[-1]).';
	    else    
			di = SandWichBA(Nl,Ndx{N-i},Alist,Blist).'; 
			ci = SandWichBMA(Hl,Hdx{N-i},M,Alist,Blist).';  
			Hl=Pleft(N,Hl,Alist,M,MR,i);
			Nl=SPleft(N,Nl,Alist,i);
        end 
	    gradi = gradi + (ci * Norm - E * di)/Norm^2;
    end
    if ~isreal(gradi)
        gradi=2*gradi;
    end
end
