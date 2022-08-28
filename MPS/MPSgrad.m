% Input:
% W - MPO representation of the Hamiltonian 
% N - number of sites
% Alist - MPS tensors
% Blist - structural tensor
% Output: 
% gradient of the expectation value formula in Figure 
% (5.5) for the standard algoritm
% (6.4) for the variation

function [gradi,Hdx,Ndx] = MPSgrad(W,Alist,Blist) 
	sites = size(Blist);
	N= sites(2);
	dims = size(Blist{2});
	d= dims(1);
	D= dims(2);
	s= dims(4);
	
	ML=W{1};
	M=W{2};
	MR=W{3};
	ordM=size(M);
	mpo=ordM(1);
	
    if ~isreal(Alist{2}) %complex case
        %complex case
        Bconj={};
        Aconj={};
        for i=1:N
            Bconj{i}= conj(Blist{i});
            Aconj{i}= conj(Alist{i});
        end
        Blist=Bconj;
        % complex case: partial derivative on conjugated variables
        Alist=Aconj;
    end
    
	if N==2
		E = ncon({HLeft(W,Alist,1),HRight(W,Alist,2)},{[1,2,3],[1,2,3]},[1,2,3],[]);
	    	Norm = ncon({NL(Alist,1),NR(Alist,2)},{[1,2],[1,2]},[1,2],[]);
	else
		E = SandWichAMA(HLeft(W,Alist,1),HRight(W,Alist,3),M,Alist,2);
	  	Norm = SandWichAA(NL(Alist,1),NR(Alist,3),Alist,2);
	end

	gradi = zeros(1,s);
	
	Al1=Alist{1};
	Bl1=Blist{1};
	A1 = reshape(Al1(:,:),[d,D]);
	CA1=conj(A1);
	B1 = reshape(Bl1(:,:,:,:),[d,D,s]);
	
	AlN=Alist{N};
	BlN=Blist{N};
	AN = reshape(AlN(:,:),[d,D]);
	CAN=conj(AN);
	BN = reshape(BlN(:,:,:,:),[d,D,s]);
	   
	Hl = HLeft(W,Alist,1); %first site
	Nl = NL(Alist,1);

	Hr1 = HRight(W,Alist,N);
	Nr1=NR(Alist,N);

	Hdx={Hr1};
	Ndx={Nr1}; %%cells of tensors to store

	for i=2:(N-1) 
		Hr1=Pright(Hr1,Alist,M,MR,N-i+1);
	    	Nr1=SPright(Nr1,Alist,N-i+1);
	    	Hdx{end+1}= Hr1;
	    	Ndx{end+1}= Nr1;
	end

	for i=1:N
		if i == 1
			di = ncon({B1,CA1,Ndx{N-1}},{[1,2,-1],[1,3],[2,3]},[1,2,3],[-1]).'; 
			ci = ncon({ML,B1,M,CA1,Hdx{N-1}},{[1],[2,3,-1],[1,4,2,5],[5,6],[3,4,6]},[1,2,5,3,4,6],[-1]).';
		elseif i==N
			di = ncon({BN,CAN,Nl},{[1,2,-1],[1,3],[2,3]},[1,2,3],[-1]).';
			ci = ncon({MR,BN,M,CAN,Hl},{[1],[2,3,-1],[4,1,2,5],[5,6],[3,4,6]},[1,2,5,3,4,6],[-1]).';
	    	else    
			di = SandWichBA(Nl,Ndx{N-i},Alist,Blist,i).'; 
			ci = SandWichBMA(Hl,Hdx{N-i},M,Alist,Blist,i).';  
			Hl=Pleft(Hl,Alist,M,MR,i);
			Nl=SPleft(Nl,Alist,i);
	    	end    
	    	gradi = gradi + (ci * Norm - E * di)/Norm^2;
    end
    if ~isreal(gradi) % if not real
        gradi=2*gradi;
    end
end
