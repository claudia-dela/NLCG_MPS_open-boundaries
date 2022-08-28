% Strong Wolfe line search
% Input:
% N -          number of sites
% B, BORTH -   structural tensors
% W  -         MPO Hamiltonian
% A  -         current point
% dir -        search direction
% slope0 -     gradient of functional at A in direction dir
% of -         functional value in A
% nstep - bound of the number of iterations of the zoom function

% FUNCTIONS
% ExpVal() -   functional
% MPSgrad() -  gradient

% Output: 
% alpha - lenght step from A along BORT(dir), Formula (6.21),
%         which minimize the funtional value

function alpha = LineSearchOrthOne(N,W,B,BORTH,A,dir,slope0,of,nstep)
	% parameter for sufficient decrease condition
	c1 = 0.0001;
	% parameter for curvature condition
	c2=0.1;

	alphaMax = 100;
	alpha    = 1;

	alpha_0  = 0;
	alpha_1  = alpha;

	of_x = of;
	of_0 = of;
	iter = 0;
	while 1
		AAA = BuildA(BORTH, alpha_1 *dir);
	    	XC={};
	    	for i=1:3
			XC{i}= A{i}+ AAA{i};
	    	end
	    	of = ExpVal(N,W,XC);
	    	slopec = conj(MPSgrad(N,W,XC,BORTH))*dir.';
       
    		if  (of > of_0 + slope0*c1*alpha_1) || ((of >= of_x ) && (iter > 0))
			alpha = Zoom(N,W,B,BORTH,A,dir,slope0, alpha_0, alpha_1,of_0,of_x ,c1,c2,nstep);
        		break;
    		end
		if(abs(slopec) <= -c2*slope0)
			alpha = alpha_1;
			break;
		end
		if (slopec >= 0)
			alpha = Zoom(N,W,B,BORTH,A,dir,slope0,alpha_1 , alpha_0,of_0,  of,c1,c2,nstep);
			break;
		end
		
		alpha_0 = alpha_1;
		alpha_1 = min(alphaMax, alpha_1*3);
		of_x = of;
		iter = iter + 1;
	end
end

function alpha = Zoom(N,W,B,BO,A,dir,slope0,alphaLo,alphaHi,of_0,ofLo,c1,c2,nstep)

	for i=1:nstep 
    		alpha = (alphaLo+alphaHi)/2;
    		AAA = BuildA(BO, alpha *dir);
    		XC={};
    		for i=1:3
        		XC{i}= A{i}+ AAA{i};
    		end
    		of = ExpVal(N,W,XC);

    		if of > of_0 + c1*alpha*slope0 || of >= ofLo
        		alphaHi = alpha;       
    		else        
        		slopec = conj(MPSgrad(N,W,XC,BO))*dir.';
             		if abs(slopec) <= -c2*slope0
            			return;
        		end
        
        		if slopec*(alphaHi-alphaLo) >= 0  
            			alphaHi = alphaLo;
            			alphaLo = alpha;
            			ofLo    = of;
        		end
    		end
    
	end
end
