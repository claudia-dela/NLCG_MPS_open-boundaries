%% Test tangent and orthogonal complement
d=3;                % physical dimension
D=2;                % bond dimension
N=3;

n = 2*D*d + (N-2)*d*D^2 ; % domain dim(D_{MPS})
s= dimension(N,d,D);      % dimension of the variety s=dim(MPS)

% Structural tensor: Formula (6.19) with PA = Id_n
% {{1×1×D×n},(N-2) tensors of order {d×D×D×n},{1×D×1×n}}
Bfull = TensorB(N,d,D,n);

v0=randn(1,n)+1i*randn(1,n);  % random starting vector 
A = BuildA(Bfull,v0); %point  (A(1),...,A(N)) = Q in C^{dim(D_{MPS})}

[Q,R]=qr(STang(A),0);
[MM,time] = giveOrthMM(A); %vector space orthogonal to the tangent direstions of the fiber
        
disp(['Norm of Q^{\dag}MM: ',num2str(norm(Q'*MM))])