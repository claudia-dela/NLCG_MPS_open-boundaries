%% Test tangent space to the gauge orbit and orthogonal complement
d=3;                % physical dimension
D=2;                % bond dimension
N=5;                % number of sites

n = 2*D + d*D^2 ;   % domain dim(D_{hMPS})
s= 2*D+d*D^2 -D^2;  % dim(D_{hMPS})-dim(G), G = gauge group

% Structural tensor: Formula (6.19) with PA = Id_n
% {{1×1×D×n}, one tensors of order {d×D×D×n},{1×D×1×n}}
Bfull = TensorB(d,D,n);

v0=randn(1,n);%+1i*randn(1,n);  % random starting vector 
A = BuildA(Bfull,v0); %point  (v_L,A,v_R) = Q in C^{dim(D_{hMPS})}

[Q,R]=qr(STang(A,N),0);
[MM] = giveOrthMM(A,N); %vector space orthogonal to the tangent direstions of the gauge orbit
        
disp(['Norm of Q^{\dag}MM: ',num2str(norm(Q'*MM))])
