%Input: 
% d,D,s - integers
% Output:
% Collection of tensors: {{d×1×D×s}, (N-2) tensors of order {d×D×D×s},{d×D×1×s}}
% if s=n=dim(D_{MPS}): Formula (6.19) with PA = Id_n
% if s<n=dim(D_{MPS}): Formula (6.19) with PA = matrix sxn = [Id_s|random(s x (n-s))] 

function [listB,Bmat] = TensorB(N,d,D,s) 
	n =2*d*D +(N-2)*d*D^2 ;
	listB={};
	k=1;
	j=1;
	i=1;
	a=1;
	b=1;
	l=1;
	Bmat = zeros(s,n);
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
		           		if k <= n-s
		               		vector = rand(1,s);
		              		else
		               		vector = zeros(1,s);
		               		vector(l)=1;
		               		l=l+1;
		           		end
		           		B(i,a,b,:)=vector;
		           		Bmat(:,k)=vector;
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

