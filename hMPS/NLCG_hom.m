% Input:
% N - number of sites
% W - MPO Hamiltonian
% B  - structural tensor
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

function [rhok,cond,count,TLS,ground_state] = NLCG_hom(N,W,B,a0, nmax, toll,nstepLS,remark_bound)
	tstart=0;
	TLS=[];
    	cond=1;
    	
    	ak = a0;
    	A = BuildA(B,ak); %tensor=point  
    	rhok = ExpVal(N,W,A);
    	gk = MPSgrad(N,W,A,B);
    	
    	iter = 0;
    	count=0;
    	
    	while ((norm(gk) > toll) && length(TLS)< remark_bound) %&& (count < boundruns)
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
        
        	%LINE-SEARCH
        	of =rhok;
        	slope0 = conj(gk)*pk.';
        	%time line search
        	tstart=tic;
        	Ak = LineSearch(N,W,B,ak,pk,slope0,of,nstepLS);
        	TLS=[TLS,toc(tstart)];
	
        	gkprec = gk;
        	ak = ak + Ak * pk;
        	
        	A = BuildA(B,ak);
        	gk = MPSgrad(N,W,A,B);
        	rhok = ExpVal(N,W,A);
	
    	end
	
	if (norm(gk) <= toll)
	        cond=0;
    end
    ground_state = ak;
end

