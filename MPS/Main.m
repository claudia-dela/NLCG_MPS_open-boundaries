%MAIN 
% Comparison between 3 NLCG methods: 
% Algorithm 5 with restart after n=dim(D_{MPS})
% Algorithm 5 with restart after s=dim(MPS)
% Algorithm 6
clear all
close all 
clc

% Notation: 
d=3;                % physical dimension fixed to 3 because AKLT Hamiltonian models a spin-1 chain
D=2;                % bond dimension

%% display 
disp(['physical dimension fixed to ',num2str(d),' because AKLT Hamiltonian models a spin-1 chain'])
disp(['bond dimension ',num2str(D)])

%%%%%%%%%%%%%%%%%%%%% AKLT MPO Hamiltonian
W = MpoAKLT();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
launch=2; 	 % number of runs
nstep = 50;

% number of sites
list=[3];
for N=list
	n = 2*D*d + (N-2)*d*D^2 ; % domain dim(D_{MPS})
    	s= dimension(N,d,D);      % dimension of the variety s=dim(MPS)
        exp_toll=8;
        toll=10^(-exp_toll);
        bound=3*s; 	% max number of NLCG iterations
        
        disp(['Number of sites: ',num2str(N)])
        disp(['gradient tolerance = ', num2str(toll)])
        disp(['iteration bound = ', num2str(bound)])        
        disp(['number of runs per site = ', num2str(launch)])
        disp(['domain dim(D_{MPS}) = ',num2str(n)])
        disp(['dimension of the variety dim(MPS) = ',num2str(s)])

    	nmaxn=n;
    	nmaxdim=s;
    
    	% Structural tensor: Formula (6.19) with PA = Id_n
        % {{1×1×D×n},(N-2) tensors of order {d×D×D×n},{1×D×1×n}}
    	Bfull = TensorB(N,d,D,n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AKLT point:
        for i=1:N
            	Aklt(1,:,:)= [0,sqrt(2/3);0,0];                   
        	Aklt(2,:,:)= [-(1/sqrt(3)),0;0,1/(sqrt(3))];
        	Aklt(3,:,:)= [0,0;-(sqrt(2/3)),0];
        end

    	AKLT={};
    	for i=2:(N-1)
        	AKLT{i}=Aklt;
    	end
    	vlr=[1,0].';
    	vl=[1,0].';
    	vr=[0,1].';
    	AKLT{1}=ncon({vl,Aklt},{[1],[-1,1,-2]},[1],[-1,-2]);
    	AKLT{N}=ncon({Aklt,vr},{[-1,-2,1],[1]},[1],[-1,-2]);
        
        % expectation value on the AKLT point
    	exact=ExpVal(W,AKLT);
        
    	MAT=[];
        GS=[];
        ground_states=[];
        
    	for j=1:launch
        	disp(['launch: ',num2str(j)])
        	v0=randn(1,n);%+1i*randn(1,n);  % random starting vector
            	MAT=vertcat(MAT,[j,toll,bound,zeros(1,6-2)]); %toll,bound,
            	GS= vertcat(GS,[j,toll,bound,zeros(1,n-1-2)]);
        
        	LS={};
        	LSmedia={};
        	LSnum={};
        	It={};
        	Time={};
        	Rho={};
        	v0Save={};
        	Cond={}; 
        	Last={};
        
        	mat=[];
            	vec=[];
        	   
 			%% 
            	disp(['Algorithm 5 with restart after dim(D_{MPS}) = ',num2str(n)])
 		timefullnStart=tic;
            	[rhon,condFulln,countn,LSfulln,gs_standard] = NLCG(W, Bfull, v0, nmaxn , toll,nstep,bound);
            	timefullnEnd= toc(timefullnStart);
            	if condFulln == 1
                	disp(['No convergence']);
            	elseif condFulln == 0
                	disp(['Convergence in ',num2str(timefullnEnd),' seconds']);
            	end
			%% 
            	disp(['Algorithm 5 with restart after dim(MPS) = ',num2str(s)])
		timefullsStart=tic;
		[rhos,condFulls,counts,LSfulls,gs_restart] = NLCG(W, Bfull, v0, nmaxdim , toll,nstep,bound);
		timefullsEnd=  toc(timefullsStart);
            	if condFulls == 1
                	disp(['No convergence']);
            	elseif condFulls == 0
                	disp(['Convergence in ',num2str(timefullsEnd),' seconds']);
            	end
			%% 
            	disp('Algorithm 6: variation of NLCG')
           	timeKOdimStart=tic;
		[rhoKOdim,condKOdim,countKOdim,LSKOdim,timeKO,gs_variation] = NLCG_variation(N,W, Bfull, v0, nmaxdim, toll,nstep,bound);
		timeKOdimEnd=  toc(timeKOdimStart);
            	if condKOdim == 1
                	disp(['No convergence']);
            	elseif condKOdim == 0
                	disp(['Convergence in ',num2str(timeKOdimEnd),' seconds']);
            	end
            
            	%% ground state 
            	ground_states=[gs_standard;gs_restart;gs_variation];
          
		%% Data
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
 	        It=[countn,counts,countKOdim];
            
        	%end
        
            	mat=vertcat(mat,[It.',Cond.',Rho.',Time.',LS.',LSmedia.',LSnum.']);
            	vec=vertcat(vec,[ground_states]);

        	MAT=vertcat(MAT,mat);
            	GS= vertcat(GS,vec);

        end
        %% writematrix
    	writematrix(MAT,['testN',num2str(N),'.dat'])
        writematrix(GS,['testN',num2str(N),'stationary_points','.dat'])

end


