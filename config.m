%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONFIGURATION %%%%%%%%%%%%%%%%%%%%
% Graph firms
% Set to one to graph firms
% Note this saves all firms in a cell, set to 0 for large sims
graphIt = 0;
% toGraph denotes the desired graph
% "NW" graphs the attention network in the last period
% "AVGNW" graphs the attention network averaged over all periods
% "SM" graphs x-theta averages over time periods
toGraph = "NW";
saveIt = 1; % Save simulation results to a new folder

%% Dynamics
maxT = 150; % Time periods to run
avgOverT = 0; % Average results over all T - uncertainty sample
maxCellLength=100000; % Split simulation runs into blocks for memory
% Whether or not the CEO doubles or halves e at T/2
% Random, Double, Half, Zero, Off
ceoAct="Random";
ceoActStartT=[2,4,15,50];
learningRate=0.25;

%% Skew handling
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
conParam=0;

%% Choose method for underlying network G
% Options: Task, JR, Full
gMethod="Full";

%% Shuffle positions:
% Essentially determines if G is independent of theta by shuffling
% Options: Random, Mu, Theta, Shuffled(default)
ShufflePositions="Random";


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
% Archetypes Dynamics
nList=[40];
mList=[1];
eList=[0,0.5,1];
consList=[-1,0,1];
% % Archetypes v
% nList=[10,15,20];
% mList=[1:30];
% eList=[0.5,1,5,10000];
% consList=[-1,0,1];
% % % Single firm
% nList = [10];
% mList = [4];
% eList = [1];
% consList = [0];

%% Type settings
% Probabilities of climbers relative to slackers.
% Simulation will check symmetrically for slackers
% Archetypes
PCscale=0.25:0.25:0.5; Wscale=[2/9,2/3];
PCscale=0.25:0.25:0.5; Wscale=[1/3];

% State Space
%PCscale=0.15:0.05:0.5; Wscale=[1/3,2/3];
% % % % Single firm
% PCscale=[0.5]; Wscale=1/4;

%% Theta settings
% Archetypes
thetascale=[2,5];
% State Space
%thetascale = [2:0.5:7];
% % % % Single firm
% thetascale=[2];



%% Baseline parameters
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
paramsDefault.learningRate=learningRate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF CONFIG %%%%%%%%%%%%%%%%%%%%
