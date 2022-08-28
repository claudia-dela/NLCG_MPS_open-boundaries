% Input:
% N-  number of sites
% Alist - hMPS tensors
% Blist - structural tensor
% Output: 
% differential of the hMPS map 

function [gradi,Hdx,Ndx] = MPSdiff(N,Alist,Blist)
	dims = size(Blist{2});
	d= dims(1);
	D= dims(2);
	s= dims(4);

    if ~isreal(Alist{2})
        %complex case
        Bconj={conj(Blist{1}),conj(Blist{2}),conj(Blist{3})};
        Blist=Bconj;
        % complex case: partial derivative on conjugated variables
        Aconj={conj(Alist{1}),conj(Alist{2}),conj(Alist{3})};
        Alist=Aconj;
    end
    
	gradi = zeros(s,d^(N-2)); 

	Al1=Alist{1};
	Bl1=Blist{1};
	A1 = Al1;
	B1 = reshape(Bl1(:,:,:,:),[D,s]);
	
	AlN=Alist{3};
	BlN=Blist{3};
	AN = AlN;
	BN = reshape(BlN(:,:,:,:),[D,s]);
	
	Hl = A1.'; %first site
	Hr1 = AN.'; %last site

	Hdx={Hr1};

	for i=2:(N-1) 
        Hr1=SCright(N,Hr1,Alist,N-i+1);
	    	Hdx{end+1}= Hr1;
	end

	for i=1:N
        if i == 1
			di = ncon({B1,Hdx{N-1}},{[1,-1],[-2,1]},[1],[-2,-1]).'; 
        elseif i==N
			di = ncon({BN,Hl},{[1,-1],[-2,1]},[1],[-2,-1]).'; 
	    else    
			di = SandWichABA(Hl,Hdx{N-i},Alist,Blist).'; 
			Hl=SCleft(N,Hl,Alist,i);
        end 
        
	    gradi = gradi + di ; 
        
    end

end
