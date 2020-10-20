%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONFIGURATION %%%%%%%%%%%%%%%%%%%%
% Graph firms
% Set to one to graph firms
% Note this saves all firms in a cell, set to 0 for large sims
% Please disable all parfor loops in main.m
graphIt = 0;
% toGraph denotes the desired graph
% "NW" graphs the attention network in the last period
% "SM" graphs x-theta averages over time periods
toGraph = "NW";
% This saves the results into the folder ../Datasave/
saveIt = 1; % Save simulation results to a new folder

%% List of parameters to run
% How many random seeds? Each seed will be put through all parameter
% combinations.
mListTotal=[1:1000];
% How many "m" do we allow per block? Reduce if memory overflow
mIterations=4;
%% Per firm
% comments and uncomment desired
thetaVar=1;
thetaMean=0;
% 
%% Single firm
%  nList = [5]; % Number of workers
%  mList = [5]; % Random Seeds / Iterations
%  eList = [1000]; % Embedding levels
%  consList = [-1]; % Consolidation levels
%  PCscale=[0.5]; Wscale=1/4; % Ratio of Improvers/Enhancers and P(Assessors)
%  thetascale=[2]; % Parameters of Beta Distribution
% % The learning Rate delta
% learningRates=[0];
%% Probability of type switching
% probTypeSwitch=0.5;

%% Smaller State Space
mListTotal=[1:100];
mIterations=1;
identifier='PFT-3x3v02'; % Use "runtime" for a timestamp
nList=[10,20,30];
eList=[0.1,0.4,0.7,0.9];
consList=[0];
PCscale=[0.25,0.5]; Wscale=[1/3]; % "2/3"
thetascale=[2,5000000];
learningRates=[0,0.05,0.4,1];
% Dynamics
maxT = 200; % Time periods to run
% Probability of type switching
probTypeSwitch=0;
%% V Parameter
% SETTING THIS TO 1 REQUIRES GLOBAL OPTIMIZATION TOOLPACKAGE
% Concave Utility switches between focal groups and focal peers
% conUtil=1 gives concave benefit
% conUtil=0 gives linear benefit
% conUtil=-1 checks both cases
conUtil=1;
% Concavity parameter if needed
conParam=0.2;


%% Detailed Settings

%% Skew handling
% Symmetric handles whether the simulation runs checking
% left+right skew and slackers<>climbers distributions
% symmetric=1 checks both left+right skew and opposite C/S ratios
% symmetric=0 checks only given C/S ratio and right skew
% symmetric=-1 checks only given C/S ratio and left skew
symmetric = 1;






%% Baseline parameters
maxCellLength=10000; % Split simulation runs into blocks for memory
avgOverT = 0; % Average results over all T - uncertainty sample
paramsDefault.maxT = maxT; % Time periods for the company to run
paramsDefault.minEqmT = 2; % Minimum periods to achieve equilibrium convergence
paramsDefault.maxEqmT = 100; % Maximum period to achieve eqm convergence
paramsDefault.globalsearch = -1; % Unused
paramsDefault.thetaRepShockVar = 0; % Standard deviation of shock to theta representations
paramsDefault.rationality = 0; % First-stage rationality of agents - can anticipate A
paramsDefault.pn = 0; % P parameter for G network. Set to 0 for full G!
paramsDefault.mn = 0; % M parameter for G network (Jackson&Rogers 2014 algorithm)
paramsDefault.maxDegree = 5; % Maximum number of peers to monitor
paramsDefault.thetaVar=thetaVar;
paramsDefault.thetaMean=thetaMean;
paramsDefault.probTypeSwitch=probTypeSwitch;
paramsDefault.ceoAct="Off";
%paramsDefault.learningRate=learningRate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF CONFIG %%%%%%%%%%%%%%%%%%%%
