% Input:
% N - number of sites
% W - MPO Hamiltonian
% B - structural tensor
% v0 - starting point (vector), dimension campatible with the hMPS domain
% nmax - bound for the line search
% toll - tolerance of the gradient 
% nstepLS - bound of the number of iterations of the zoom function

% boundruns - iteration bound

% Output:
% rhok - expectation value of the functional
% cond - convergence condition: 0=yes, 1=no
% count - count of the iteration of NLCG
% TLS - list of time of the linsearches
function [rhok,cond,count,TLS, ground_state] = NLCG_hom_variation(N , W, B, v0, nmax, toll,nstepLS,remark_bound)
	dims = size(B{2});
	d= dims(1);
    	D= dims(2);
    	s= dims(4);
    	
    	timels=0;
    	TLS=[];
    	cond=1;
    
    	vk = v0;
    	A = BuildA(B,vk); %point  (v_L,A,v_R) = Q in C^{dim(D_{hMPS})}
    	
    	rhok = ExpVal(N,W,A);
    
    	MM = giveOrthMM(A,N); % < v1,...,v_s > = vector space orthogonal to the tangent directions of the gauge orbit
                        
    	Borth = TensorBnew(d,D,MM); %tensor of the affine linear space Q+< v1,...,v_s > 
    	gk = MPSgrad(N,W,A,Borth);
 
 	   iter = 0;
 	   count=0;
    	while ((norm(gk) > toll) && length(TLS)< remark_bound ) %&& (count < boundruns)
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
        	timels=tic;
        	Ak = LineSearchOrthOne(N,W,B,Borth,A,pk,slope0,of,nstepLS);
        	TLS=[TLS,toc(timels)];
        	
	        gkprec = gk; %previous gradient
           
        	AB = BuildA(Borth,pk);
	        AA={};
        	for i=1:3
        		AA{i}= A{i}+ Ak *AB{i};
        	end
        	A=AA;
        
        	MM = giveOrthMM(A,N); %< v1,...,v_s > = vector space orthogonal to the tangent directions of gauge orbit
        	Borth = TensorBnew(d,D,MM); %tensor of the affine linear space Q+< v1,...,v_s > 
        	gk = MPSgrad(N,W,A,Borth); %update gradient!
        	rhok = ExpVal(N,W,A);
        
    	end
         
    	if (norm(gk) <= toll)
        	cond=0;
        end
        ground_state = Buildv(A);
end

