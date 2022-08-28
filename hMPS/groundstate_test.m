clear all
clc
filelist=["N3","N4","N5","N6","N7",...
    "N8","N9","N10","N11","N12","N13","N14","N15","N16",...
    "N17","N18","N19","N20","N21","N22",...
    "N23","N24","N25","N26","N27","N28",...
    "N29","N30",...
    "N31","N32","N33","N34","N35","N36","N37","N38","N39","N40"...
    ];
format long

W= MpoAKLT();
N=5;
for file=["testN3"] %,"N4","N5","N6","N7","N8"] %filelist
    disp(['d = ',num2str(N)])
Ni=importdata(strjoin([file,".dat"],'')) 
Nipoint=importdata(strjoin([file,"stationary_points.dat"],'')) 

point = str2num(Nipoint{2})
lineinfile = str2num(Ni{2})
eigval = real(lineinfile(3))

disp(['eigvalue = ',num2str(eigval)])


d=3;%dim fisica
%N=3;%numero siti
D=2;%bond dim

n = 2*D + d*D^2 ;   % domain dim(D_{hMPS})
c= 2*D+d*D^2 -D^2;  % c = dim(D_{hMPS})-dim(G), G = gauge group
s= c-1;             % dimension of the variety s = c-1 = dim(hMPS)

Bfull = TensorB(d,D,n);
A=BuildA(Bfull,point);
% Hamiltonian 


 %%test MPOMPS
 Alist=A;
 %'products'
 PRODUCT = MPOMPS(N,W,Alist).'
 %'lines'
 LINE = MPScontraction(N,Alist).';
 eigLINE = eigval*LINE;
 
  disp(['norm(MA-aA) = 0 ? ',num2str(norm(PRODUCT-eigLINE))])
  disp(['differential is zero? ',num2str(norm(MPSdiff(N,Alist,Bfull)))])
  disp(['vector in kernel? ',num2str(norm(PRODUCT*LINE.'))])
 
 N=N+1;
end