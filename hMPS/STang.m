% Input: 
% Alist - hMPS tensora 
% Output:
% Matrix whose columns span the tangent space to the gauge orbit, Formula (6.16)
function [tang] = STang(Alist)
	dims = size(Alist{2});
	d= dims(1);
	D= dims(2);

	vL=Alist{1};
	A=Alist{2};
	vR=Alist{3};

	mL = horzcat(vL,[vL(1);-vL(2)],[0;vL(1)],[vL(2);0]);

	Ain=[];
	for i=1:d
		m=[0,0,-A(i,2,1),A(i,1,2);...
		   0,-2*A(i,1,2),A(i,1,1)-A(i,2,2),0;...
		   0,2*A(i,2,1),0,A(i,2,2)-A(i,1,1);...
		   0,0,A(i,2,1),-A(i,1,2)];
	    	Ain=vertcat(Ain,m);
	end

	mR= horzcat(-vR,[-vR(1);vR(2)],[-vR(2);0],[0;-vR(1)]);
	tang=vertcat(mL,Ain,mR);
end
