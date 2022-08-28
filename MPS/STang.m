% Input: 
% Alist - MPS tensors 
% Output:
% Matrix whose columns span the tangent space to the gauge orbit, Formula (6.13)
function [tang] = STang(Alist)
	sites = size(Alist);
	N= sites(2);
	dims = size(Alist{2});
	d= dims(1);
	D= dims(2);
	mL=[];
	for i=1:d
    		mL = vertcat(mL,kron(Alist{1}(i,:),eye(D)));
	end
	z1=zeros(2*d,4*(N-2));
	mL=horzcat(mL,z1);

	Ain=[];
	for j=2:(N-1)
    		Aj=[];
    		for k=1:d
        		Aj=vertcat(Aj,...
        		horzcat(...
        		kron( -eye(D),(reshape((Alist{j}(k,:,:)),[D,D])).'),...
            		kron(  reshape((Alist{j}(k,:,:)),[D,D]),eye(D))...
            			)...
            		);
        
    		end
    		if j==2
        		zjsx=[];
    		else
        		zjsx=zeros(4*d,4*(j-2));
    		end
    		if j==N-1
        		zjdx=[];
    		else
        		zjdx=zeros(4*d,4*(N-2-j+1));
    		end
    		Aj=horzcat(zjsx,Aj,zjdx);
       	Ain=vertcat(Ain,Aj);
	end
	mR=[];
	for i=1:d
    		mR = vertcat(mR,kron(-eye(D),Alist{N}(i,:)));
	end   
	zN=zeros(2*d,4*(N-2));
	mR=horzcat(zN,mR);

	tang=vertcat(mL,Ain,mR);
end
