%Input: 
% A - hMPS tensors
% N - number of sites
% Output: 
% MM - kernel of the matrix TA, contructed via Formula (6.16)
% The columns of TA span the tangent space to the gauge orbit.
% MM is the orthogonal complement in the domain of the
% parametrization
function [MM] = giveOrthMM(A,N)
	[Q,R]=qr(STang(A,N),0);
	kk=Q;
	MM=null(kk');
end
