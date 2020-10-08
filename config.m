%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONFIGURATION %%%%%%%%%%%%%%%%%%%%
% Graph firms
% Set to one to graph firms
% Note this saves all firms in a cell, set to 0 for large sims
% Please disable all parfor loops in main.m
graphIt = 0;
% toGraph denotes the desired graph
% "NW" graphs the attention network in the last period
% "SM" graphs x-theta averages over time periods
toGraph = "SM";
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


%% Archetype State Space
identifier='PFT-3x3Figure3'; % Use "runtime" for a timestamp
%nList=[20:10:40];
nList=[60];
eList=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99,1];
%eList=[0.999,0.9999,0.99999,0.999999];
%eList=[0.10,  0.50, 1.00, 5.00, 50,100,1000,10000];
consList=[0];
PCscale=[0.25,0.5]; Wscale=[1/3]; % "2/3"
thetascale=[2,5000000];
learningRates=[0,0.05,0.4,0.7,1];
%% Dynamics
maxT = 200; % Time periods to run


% %% Figure 1 and 2
% identifier='PFT-3x3Figure1ad'; % Use "runtime" for a timestamp
% nList=[60];
% %eList=[0,0.01,0.001,0.1,0.5,0.99];
% eList=[0.01,0.05];
% consList=[0];
% PCscale=[0,0.25,0.35,0.45,0.5]; Wscale=[1/3]; % "2/3"
% thetascale=[2,3,4,5,8,10,250,5000,500000000];
% learningRates=[0.5];
% %% Dynamics
% maxT = 100; % Time periods to run


%% Detailed Settings

% Whether or not the CEO doubles or halves e at T/2
% Random, Double, Half, Zero, Off
ceoAct="Off";
% Times where the CEO intervenes
ceoActStartT=[1];

%% Skew handling
% Symmetric handles whether the simulation runs checking
% left+right skew and slackers<>climbers distributions
% symmetric=1 checks both left+right skew and opposite C/S ratios
% symmetric=0 checks only given C/S ratio and right skew
% symmetric=-1 checks only given C/S ratio and left skew
symmetric = 1;

%% V Parameter
% SETTING THIS TO 1 REQUIRES GLOBAL OPTIMIZATION TOOLPACKAGE
% Concave Utility switches between focal groups and focal peers
% conUtil=1 gives concave benefit
% conUtil=0 gives linear benefit
% conUtil=-1 checks both cases
conUtil=0;
% Concavity parameter if needed
conParam=0;

%% G Network
%% Choose method for underlying network G
% Normally: Use Full
% For robustness check: Use JR
% Options: Task, JR, Full
gMethod="Full";

%% Shuffle positions:
% Essentially determines if G is independent of theta by shuffling
% Options: Random, Mu, Theta, Shuffled(default)
ShufflePositions="Shuffled";


%% Baseline parameters
maxCellLength=10000; % Split simulation runs into blocks for memory
avgOverT = 0; % Average results over all T - uncertainty sample
paramsDefault.gMethod=gMethod;
paramsDefault.maxT = maxT; % Time periods for the company to run
paramsDefault.minEqmT = 2; % Minimum periods to achieve equilibrium convergence
paramsDefault.maxEqmT = 100; % Maximum period to achieve eqm convergence
paramsDefault.globalsearch = -1; % Unused
paramsDefault.thetaRepShockVar = 0; % Standard deviation of shock to theta representations
paramsDefault.rationality = 0; % First-stage rationality of agents - can anticipate A
paramsDefault.ceoAct = ceoAct; % whether the CEO adapts over time
paramsDefault.pn = 0; % P parameter for G network. Set to 0 for full G!
paramsDefault.mn = 0; % M parameter for G network (Jackson&Rogers 2014 algorithm)
paramsDefault.maxDegree = 40; % Maximum number of peers to monitor
paramsDefault.shufflePositions=ShufflePositions;
paramsDefault.thetaVar=thetaVar;
paramsDefault.thetaMean=thetaMean;

%paramsDefault.learningRate=learningRate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF CONFIG %%%%%%%%%%%%%%%%%%%%
