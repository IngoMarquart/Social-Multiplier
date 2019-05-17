%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONFIGURATION %%%%%%%%%%%%%%%%%%%%
% Graph firms
% Set to one to graph firms
% Note this saves all firms in a cell, set to 0 for large sims
graphIt = 1;
% toGraph denotes the desired graph
% "NW" graphs the attention network in the last period
% "AVGNW" graphs the attention network averaged over all periods
% "SM" graphs x-theta averages over time periods
toGraph = "NW";
saveIt = 0; % Save simulation results to a new folder
maxT = 1; % Time periods to run
avgOverT = 0; % Average results over all T - uncertainty sample
% Symmetric handles whether the simulation runs checking
% left+right skew and slackers<>climbers distributions
% symmetric=1 checks both left+right skew and opposite C/S ratios
% symmetric=0 checks only given C/S ratio and right skew
% symmetric=-1 checks only given C/S ratio and left skew
symmetric = 1;
% Concave Utility switches between focal groups and focal peers
% conUtil=1 gives concave benefit
% conUtil=0 gives linear benefit
% conUtil=-1 checks both cases
conUtil=0;
% Concavity parameter if needed
conParam=0.5;



%% Baseline parameters
paramsDefault.maxT = maxT; % Time periods for the company to run
paramsDefault.minEqmT = 2; % Minimum periods to achieve equilibrium convergence
paramsDefault.maxEqmT = 10; % Maximum period to achieve eqm convergence
paramsDefault.globalsearch = -1; % Unused
paramsDefault.thetaRepShockVar = 0; % Standard deviation of shock to theta representations
paramsDefault.rationality = 0; % First-stage rationality of agents - can anticipate A
paramsDefault.ceoAct = 0; % whether the CEO adapts over time
paramsDefault.pn = 0; % P parameter for G network. Set to 0 for full G!
paramsDefault.mn = 0; % M parameter for G network (Jackson&Rogers 2014 algorithm)

%% List of parameters to run
% % State Space
%nList=[50:10:80];
%mList=[1:30];
%eList=[0.10,  0.30,  0.50,  0.70, 0.90, 1.00, 5.00, 50,100,500];
%consList=[-1,-0.5,0,0.5,1];
% % Archetypes
% nList=[15:5:70];
% mList=[1:30];
% eList=[0.10,  0.30,  0.50,  0.70, 0.90, 1.00, 5.00, 50,100,500];
% consList=[-1,-0.5,0,0.5,1];
% % Single firm
nList = [40];
mList = [4];
eList = [1];
consList = [0];

%% Type settings
% Probabilities of climbers relative to slackers.
% Simulation will check symmetrically for slackers
% Archetypes
% PCscale=0.25:0.25:0.5; Wscale=[2/9,2/3];
% State Space
%PCscale=0.15:0.05:0.5; Wscale=[1/3,2/3];
% Single firm
 PCscale=0.5; Wscale=1/3;

%% Theta settings
% Archetypes
% thetascale=[2,5];
% State Space
%thetascale = [2:0.5:7];
% Single firm
 thetascale=[2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF CONFIG %%%%%%%%%%%%%%%%%%%%
