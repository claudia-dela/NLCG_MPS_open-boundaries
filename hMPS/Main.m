% MAIN.
% Comparison between 3 NLCG methods: 
% Algorithm 5 with restart after n=dim(D_{hMPS})
% Algorithm 5 with restart after dim=dim(hMPS)
% Algorithm 6
clear all
close all 
clc
% Notation: 
d=3;                % physical dimension fixed to 3 because AKLT Hamiltonian models a spin-1 chain
D=2;                % bond dimension
n = 2*D + d*D^2 ;   % domain dim(D_{hMPS})
c= 2*D+d*D^2 -D^2;  % c = dim(D_{hMPS})-dim(G), G = gauge group
s= c-1;             % dimension of the variety s = c-1 = dim(hMPS)

%% display 
disp(['physical dimension fixed to ',num2str(d),' because AKLT Hamiltonian describes a spin-1 chain'])
disp(['bond dimension ',num2str(D)])
disp(['domain dim(D_{hMPS}) = ',num2str(n)])
disp(['dimension of the variety dim(hMPS) = ',num2str(s)])
%%%%%%%%%%%%%%%%%%%%% AKLT MPO Hamiltonian
W = MpoAKLT();

% Structural tensor: Formula (6.19) with PA = Id
%cells: {1×1×D×n},{d×D×D×n},{1×D×1×n}
Bfull = TensorB(d,D,n);

%AKLT point:
AKLT=[];
AKLT(1,:,:)= [0,sqrt(2/3);0,0];                   
AKLT(2,:,:)= [-(1/sqrt(3)),0;0,1/(sqrt(3))];
AKLT(3,:,:)= [0,0;-(sqrt(2/3)),0];
aklt={[1,0].',AKLT,[1,0].'};

ND=3; %Number of tensors in the domain. In the uniform case we have always ND=3, i.e. (vl,A,vr)

toll = 10^(-6); % gradient tolerance
bound=3*s;        % max number of NLCG iterations
launch=2;       % number of runs 
nmaxn=n;
nmaxdim=s;
nstep = 50;  
disp(['gradient tolerance = ', num2str(toll)])
disp(['iteration bound = ', num2str(bound)])
disp(['number of runs per site = ', num2str(launch)])
% number of sites
list=[5];

for N=list
	if N < 5
        	warning('almeno 5!');
        	return
    	end
    	disp(['Number of sites: ',num2str(N-2)])
    	% expectation value on the AKLT point
    	exact=ExpVal(N,W,aklt);

    	MAT=[];
    	GS=[];
    	ground_states=[];
    
    	for j=1:launch
		disp(['launch: ',num2str(j)])
		%v_aklt = Buildv(aklt);
		%v0=v_aklt
		v0=randn(1,n); %+1i*randn(1,n); % random starting vector  
		MAT=vertcat(MAT,[j,toll,bound,zeros(1,6-2)]); %toll,bound,
		GS= vertcat(GS,[j,toll,bound,zeros(1,n-1-2)]);
		%ground_states=vertcat(ground_states,[j,zeros(1,n-1)]);  %launch
        
		It={};
		Time={};
		Rho={};
		v0Save={};
		Cond={}; 
		LS={};
		LSmedia={};
		LSnum={};
		ground_states={};
        
		mat=[];
        	vec=[];

        	%% 
            	disp(['Algorithm 5 with restart after dim(D_{hMPS}) = ',num2str(n)])
	    	timefullnStart=tic;
		[rhon,condFulln,countn,LSfulln,gs_standard] = NLCG_hom(N,W, Bfull, v0, nmaxn , toll,nstep,bound);
	    	timefullnEnd= toc(timefullnStart);
	    	if condFulln == 1
			disp(['No convergence']);
		elseif condFulln == 0
			disp(['Convergence in ',num2str(timefullnEnd),' seconds']);
	    	end

            	%% 
            	disp(['Algorithm 5 with restart after dim(hMPS) = ',num2str(s)])
	    	timefullsStart=tic;
	    	[rhos,condFulls,counts,LSfulls,gs_restart] = NLCG_hom(N,W, Bfull, v0, nmaxdim , toll,nstep,bound);
	    	timefullsEnd=  toc(timefullsStart);
            	if condFulls == 1
                	disp(['No convergence']);
            	elseif condFulls == 0
                	disp(['Convergence in ',num2str(timefullsEnd),' seconds']);
            	end

            	%% 
            	disp('Algorithm 6: variation of NLCG')
            	timeKOdimStart=tic;
            	[rhoKOdim,condKOdim,countKOdim,LSKOdim,gs_variation] = NLCG_hom_variation(N,W, Bfull, v0, nmaxdim, toll,nstep,bound);
		timeKOdimEnd=  toc(timeKOdimStart);
            	if condKOdim == 1
                	disp(['No convergence']);
            	elseif condKOdim == 0
                	disp(['Convergence in ',num2str(timeKOdimEnd),' seconds']);
            	end
            
            	%% ground state 
            	ground_states=[gs_standard;gs_restart;gs_variation];
            
            	%% save data %%%%%%%%%%%%%%%%%%%%%%%
            	% total time line search
	    	lsfulln=sum(LSfulln,'all');
	    	lsfulls=sum(LSfulls,'all');
	    	lsKOdim=sum(LSKOdim,'all');
	    	LS=[lsfulln,lsfulls,lsKOdim];
            
            	% mean time line search 
            	LSmedia=[lsfulln/length(LSfulln),lsfulls/length(LSfulls),lsKOdim/length(LSKOdim)];
            
            	% mean number of line searches
            	LSnum=[length(LSfulln),length(LSfulls),length(LSKOdim)];
            
            	% time of algorithms
            	Time=[timefullnEnd,timefullsEnd,timeKOdimEnd];
            
            	% expectation value reached
            	Rho=[rhon,rhos,rhoKOdim];
            
            	% condition of convergence: 0 = yes, 1 = no
            	Cond=[condFulln,condFulls,condKOdim];
            
            	% number of iterations
            	It=[countn+1,counts+1,countKOdim+1];

        	%end
        
            	mat=vertcat(mat,[It.',Cond.',Rho.',Time.',LS.',LSmedia.',LSnum.']);
            	vec=vertcat(vec,[ground_states]);
        
        	MAT=vertcat(MAT,mat);
        	GS= vertcat(GS,vec);
    	end
    	writematrix(MAT,['testN',num2str(N-2),'.dat'])
    	writematrix(GS,['testN',num2str(N-2),'stationary_points','.dat'])

end


