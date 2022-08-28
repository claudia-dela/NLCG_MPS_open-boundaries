%remember 5sites= 2 boundaries+3inner
clear all
clc
for N=[5]
%Hamiltonian
W= MpoAKLT();

d=3;                % physical dimension fixed to 3 because AKLT Hamiltonian models a spin-1 chain
D=2;                % bond dimension
n = 2*D + d*D^2 ;   % domain dim(D_{hMPS})
Bfull = TensorB(d,D,n);


%aklt
AKLT=[];
AKLT(1,:,:)= [0,sqrt(2/3);0,0];                   
AKLT(2,:,:)= [-(1/sqrt(3)),0;0,1/(sqrt(3))];
AKLT(3,:,:)= [0,0;-(sqrt(2/3)),0];
aklt={[1,0].',AKLT,[1,0].'};

exact=ExpVal(N,W,aklt);

PRODUCT = MPOMPS(N,W,aklt);
LINE = MPScontraction(N,aklt);
eigLINE = exact*LINE;

disp(['norm(MA-aA) = 0 ? ',num2str(norm(PRODUCT-eigLINE))])
disp(['differential is zero? ',num2str(norm(MPSdiff(N,aklt,Bfull)))])
disp(['vector in kernel? ',num2str(norm(PRODUCT*LINE.'))])
 
end
