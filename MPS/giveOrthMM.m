%Input: 
% A - MPS tensors
% Output: 
% MM - kernel of the matrix TA, contructed via Formula (6.13)
% The columns of TA span the tangent space to the gauge orbit.
% MM is the orthogonal complement in the domain of the
% parametrization
% timeQR - time of computation of QR and null functions. 
function [MM,timeQR] = giveOrthMM(A)
	TimeStart1=tic;
	[Q,R]=qr(STang(A),0);
	kk=Q;
	MM=null(kk');
	timeQR=toc(TimeStart1);
end
