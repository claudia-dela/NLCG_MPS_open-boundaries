%Input: 
% N,d,D - integers
% Mat - matrix nxs, with n=dim(D_{MPS})=2*d*D +(N-2)*d*D^2
% Output:
% Collection of tensors: {{d×1×D×s}, (N-2) tensors of order {d×D×D×s},{d×D×1×s}}
% built from Formula (6.19) with PA = Mat 
function [listB] = TensorBnew(N,d,D,Mat)
	n = 2*d*D +(N-2)*d*D^2;
	k=1;
	j=1;
	i=1;
	a=1;
	b=1;
	l=1;
	listB={};
	while j<= N 
		i=1;
	    	while i<= d 
			a=1;
			if (j==1)
		    		D1=1;
		    		D2=D;
			elseif (j==N)
		    		D1=D;
		    		D2=1;
			else
		    		D1=D;
		    		D2=D;
			end
			while a<= D1 
		    		b=1;
		    		while b <= D2 
		        		B(i,a,b,:)= Mat(k,:);
		        		k=k+1;
		        		b = b+1;
		    		end
		    		a = a+1;
			end
			i = i+1;
	    	end
	    	j=j+1;
	    	listB{end+1}=B;
	    	clear B;
	end

end
