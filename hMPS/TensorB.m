%Input: 
% d,D,s integers
% Output:
% Collection of tensors: {{1×1×D×s},{d×D×D×s},{1×D×1×s}}
% if s=n=dim(D_{hMPS}): Formula (6.19) with PA = Id_n
% if s<n=dim(D_{hMPS}): Formula (6.19) with PA = matrix sxn = [Id_s|random(s x (n-s))] 

function [listB,Bmat] = TensorB(d,D,s) 
	n =2*D + d*D^2 ; %boundary vectors + 1 tensor 
	ND=3; 
	listB={};
	k=1;
	j=1;
	i=1;
	a=1;
	b=1;
	l=1;
	Bmat = zeros(s,n);
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

