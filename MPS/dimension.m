% Input:
% N - number of sites
% d1 - physical dimension
% D1 - bond dimension
% Output:
% dimension of the variety MPS(D1,d1,N)

function [s]=dimension(N,d1,D1)
	b1=min(d1,D1);
	b2=min(d1,D1^2);
	if b1<d1 
		bb1 = b1*(d1-b1);
	else 
	    	bb1=0;
	end
	if b2<d1 
	    	bb2 = b2*(d1-b2);
	else 
	    	bb2=0;
	end
	s = min((N-2) * (D1^2 * b2 -1)+ 2*(D1*b1-1)+ 2*bb1+(N-2)*bb2-(N-1)*(D1^2-1)+1,d1^N);
end
