% Input:
% N - number of sites
% W - MPO Hamiltonian
% B - structural tensor
% v0 - starting point (vector), dimension campatible with the MPS domain
% nmax - bound for the line search
% toll - tolerance of the gradient 
% nstepLS -  bound of the number of iterations of the zoom function
% boundruns - iteration bound

% Output:
% rhok - expectation value of the functional
% cond - convergence condition: 0=yes, 1=no
% count - count of the iteration of NLCG
% TLS - list of time of the linsearches
function [rhok,cond,count,TLS,timeKO,ground_state] = NLCG_variation(N , W, B, v0, nmax, toll,nstepLS,boundruns)
	dims = size(B{2});
	d= dims(1);
    	D= dims(2);
    	s= dims(4);
    
    	timeKO=0;
    	TLS=[];
    	cond=1;
    	vk = v0;
    	A = BuildA(B,vk); %point  (A(1),...,A(N)) = Q in D_{MPS(D,d,N)}
    
    	rhok = ExpVal(W,A);
    
    	[MM,time] = giveOrthMM(A); %vector space orthogonal to the tangent directions of the gauge orbit
        %MM = [ v1 ... vs ] in C^(D_{MPS} x dim(MPS)) 
    	
        timeKO=timeKO+time;
    	Borth = TensorBnew(N,d,D,MM); %tensor of the affine linear space Q+< v1,...,v_s > 
    	gk = MPSgrad(W,A,Borth);
    
    	iter = 0;
    	count=0;
    	while ((norm(gk) > toll) && (length(TLS) < boundruns)) %&& (count < boundruns))
    		iter = iter + 1;
        	if iter== nmax+1
            		count=count+1;
            		iter=1;
        	end
        	if iter == 1 
            		pk = - gk;
            	else
                       beta = ((conj(gk)*gk.')/(conj(gkprec)*gkprec.'));
            		pk = -gk + beta * pk;
            	end
		%%%
		%LINE-SEARCH
		of = rhok;  
		slope0 = conj(gk)* pk.';   
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		lstempo=tic;
		Ak = LineSearchOrthOne(N,W,Borth,A,pk,slope0,of,nstepLS);
		TLS=[TLS,toc(lstempo)];
		
		gkprec = gk; %old gradient
		
	   	AB = BuildA(Borth,pk);
		AA={};
		for i=1:N
			AA{i}= A{i}+ Ak *AB{i};
		end
		A=AA;
		
		[MM,time] = giveOrthMM(A); %vector space orthogonal to the tangent direstions of the fiber
		timeKO=timeKO+time;
		Borth = TensorBnew(N,d,D,MM); %tensor of the affine linear space Q+< v1,...,v_s > 
		
		gk = MPSgrad(W,A,Borth); %update gradient!
		rhok = ExpVal(W,A);
				
    	end 
    	if (norm(gk) <= toll)
        	cond=0;
        end
        ground_state= Buildv(A);
end

