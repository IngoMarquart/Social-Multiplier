%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONFIGURATION %%%%%%%%%%%%%%%%%%%%
%% ROBUSTNESS CHECK: Loss Aversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This configuration file replicates the k-parameter robustness check
%% Saving and graphing
% This saves the results into the folder ../Datasave/
saveIt = 1; % Save simulation results to a new folder
% Graph firms - recommended for only a single firm
% ATTN: If set to 1, please disable all parfor loops in main.m
graphIt = 0;


%% V Parameter

% Concavity parameter "v"
% We check 0.2,0.5,2 (main results),1.5
conParam=1;
%% List of parameters to run
% How many random seeds? Each seed will be put through all parameter
% combinations.
mListTotal=[1:12];
% How many "m" do we allow per block? Reduce if memory overflow
mIterations=1;
% Identifier for the folders being set
identifier='PFT-3x3-LA'; % Use "runtime" for a timestamp
% Number of employees in the firm
nList=[10,20,30];
% Embedding levels to check
eList=[0.1,0.4,0.7,0.9];
% PCscale - Ratio of improvers to enhancers. Simulation also runs the
% equivalent ratio of enhancers to improvers.
% Wscale - Overall probability of assessors
PCscale=[0.25,0.5]; Wscale=[1/3]; % "2/3"
% Skew parameter: Second parameter of Beta distribution
% Simulation also checks the other direction
thetascale=[2,5000000];
% Alpha levels
learningRates=[0,0.05,0.4,1];
% Correlation of P and F (unused)
consList=[0];
% Dynamics
maxT = 200; % Time periods to run
% Probability of type switching
probTypeSwitch=0;
% Loss Aversion
kList=[0.01,0.1,0.3,0.5,0.7,0.9,1];


%% Further Settings
% Symmetry handling
% Symmetric handles whether the simulation runs checking
% left+right skew and slackers<>climbers distributions
% symmetric=1 checks both left+right skew and opposite C/S ratios
% symmetric=0 checks only given C/S ratio and right skew
% symmetric=-1 checks only given C/S ratio and left skew
symmetric = 1;
% Baseline parameters for theta: Variance and Mean
thetaVar=1;
thetaMean=0;
% toGraph denotes the desired graph
% "NW" graphs the attention network, averaged over all periods
% "SM" graphs x-theta averages over time periods
% "SM-stacked" visualizes x-theta for several firms
toGraph = "NW";
% The following explicitly allows for concave social benefit
% and is required for this robustness check
% SETTING THIS TO 1 REQUIRES GLOBAL OPTIMIZATION TOOLPACKAGE
% Concave Utility switches between focal groups and focal peers
% conUtil=1 gives concave benefit
% conUtil=0 gives linear benefit
% conUtil=-1 checks both cases
conUtil=0;

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
paramsDefault.ceoAct="Off"; % Embedding level set in first period, other options unused
%paramsDefault.learningRate=learningRate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF CONFIG %%%%%%%%%%%%%%%%%%%%
