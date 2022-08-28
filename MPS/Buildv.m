%Input:
% A - tensor
% Output:
% vector reshape of the tensor A
function [v] = Buildv(A)
	sites = size(A);
	N= sites(2);
	dims = size(A{2});
	d= dims(1);
	D= dims(2);

	v=[];
	for j=1:N
		B=A{j};
	    	for i=1:d
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
			for a=1:D1
		    		for b=1:D2
		        	        B=reshape(B,[d,D1,D2]);
		        		v=[v,B(i,a,b)];
		    		end
			end
	    	end
	end
end
