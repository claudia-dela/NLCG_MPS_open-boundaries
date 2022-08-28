%Input:
% B - tensor
% P - vector
% with compatible dimension (see Output)
% Output:
% matrix reshape of B and B-P product 
function [A] = BuildA(B,P)
	dims = size(B{2});
	d= dims(1);
	D= dims(2);
	s= dims(4);
	ND=3;
	A={};
	for j=1:ND
	    	C=B{j};
	    	if (j==1 || j==ND)
			d1=1;
	    	else
			d1=d;
	    	end
	    	for i=1:d1
			if (j==1)
			    	D1=1;
			    	D2=D;
			elseif (j==ND)
		    		D1=D;
		    		D2=1;
			else
		    		D1=D;
		    		D2=D;
			end
			for a=1:D1
		    		for b=1:D2
		        	       T(i,a,b) = (reshape(C(i,a,b,:),[1,s]))*P.';
		    		end
			end
	    	end
	    	if (j==1 || j==ND)
			T=reshape(T,[D,1]);
	    	end
	    	A{j}=T;
	    	clear T;
	end
end
