% MPO Hamiltonian
% Output:
% W = list of:
% ML - reshape([1;0;0;0;0;0;0;0;0;0;0],[chi,1,1]);
% M - inner tensor of order (d x d x chi x chi) of the MPO representation
% Formula (5.20)
% MR - reshape([0;0;0;0;0;0;0;0;0;0;1],[chi,1,1]);

function W = MpoAKLTnew() %[ML,M,MR] = MpoAKLTnew()
	d=3;
	S1=(1/sqrt(2))*[0 1 0;1 0 1; 0 1 0];
	S2=(1/sqrt(2))*[0 -1j 0;1j 0 -1j; 0 1j 0];
	S3 =[1 0 0;0 0 0; 0 0 -1];
	SI = eye(d);
	c = 1/3;

	M = zeros(11,11,d,d);

	M(1,1,:,:)= SI;
	M(11,11,:,:) = SI;

	M(1,2,:,:)= S1^2;
	M(1,3,:,:)= S2^2;
	M(1,4,:,:)= S3^2;
	M(1,5,:,:)= S1*S2;
	M(1,6,:,:)= S1*S3;
	M(1,7,:,:)= S2*S1;
	M(1,8,:,:)= S2*S3;
	M(1,9,:,:)= S3*S1;
	M(1,10,:,:)= S3*S2;

	a= -2/3;

	M(2,11,:,:)= c*S1^2;
	M(3,11,:,:)= c*S2^2;
	M(4,11,:,:)= c*S3^2;
	M(5,11,:,:)= a*S1*S2 + S2*S1;
	M(6,11,:,:)= a*S1*S3 + S3*S1;
	M(7,11,:,:)= a*S2*S1 + S1*S2;
	M(8,11,:,:)= a*S2*S3 + S3*S2;
	M(9,11,:,:)= a*S3*S1 + S1*S3;
	M(10,11,:,:)= a*S3*S2 + S2*S3;

	ML = reshape([1;0;0;0;0;0;0;0;0;0;0],[11,1,1]);
	MR = reshape([0;0;0;0;0;0;0;0;0;0;1],[11,1,1]);
	
	W = {ML,M,MR};
end

