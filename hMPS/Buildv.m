%Input:
% A - tensor
% Output:
% vector reshape of the tensor A
function [v] = Buildv(A)
	sites = size(A);
	ND= sites(2);
	dims = size(A{2});
	d= dims(1);
	D= dims(2);

	v=[];
	for j=1:ND
    		B=A{j};
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
                	               B=reshape(B,[d1,D1,D2]);
                			v=[v,B(i,a,b)];
            			end
        		end
    		end
	end
end
