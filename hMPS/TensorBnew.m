%Input: 
% d,D - integers
% Mat - matrix nxs, n=dim(D_{hMPS})
% Output:
% Collection of tensors: {{1×1×D×s},{d×D×D×s},{1×D×1×s}}
% built from Formula (6.19) with PA = Mat 

function [listB] = TensorBnew(d,D,Mat)
	n = 2*D + d*D^2;
	ND=3;
	k=1;
	j=1;
	i=1;
	a=1;
	b=1;
	l=1;
	listB={};
	while j<= ND 
		i=1;
	    	if (j==1 || j==ND)
			d1=1;
	    	else
			d1=d;
	    	end
	    	while i<= d1 
			a=1;
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
