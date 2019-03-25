%%
%   initFirm.m
%   For given parameter ranges, initiate a firm
%%
% @param: n - Number of actors
% @param: T - Dynamic time periods
% @param: gamma - Vector of type probabilities
% @param: thetaD - Parameters of Beta distribution
% @param: e - Embedding
% @param: m - Firm Seed
% @param: minEqmT - Minimum time periods to run (independent of convergence)
% @param: maxEqmT - Minimum time periods to run (independent of convergence)
% @param: globalsearch - Search method
% @param: pn,mn - Exogenous social matrix configuration
% @param: cons - Consolidation
% @return: firm - a firm
%%
function firm=initFirm(params)

% Set the random stream
s = RandStream('mcg16807','Seed',params.m);
RandStream.setGlobalStream(s);

% First period is initialization
params.maxT=params.maxT+1;

%% Generation of types and quality vectors
% Matrix to hold identity and thetas
TIVec = zeros(params.n,2);
pd = makedist('Beta',params.thetaD(1),params.thetaD(2));

% If exogenous scaling is desired use this and set fourth and third parameters
%TIVec(:,1) = (random(pd,n,1)-pd.mean).*thetaD(4)+thetaD(3);

% We scale automatically by fixing variance.
TIVec(:,1) = random(pd,params.n,1);
c=sqrt(1/var(TIVec(:,1)));
% This fixes variance to 1
TIVec(:,1) =TIVec(:,1).*c;
% De-mean theta
TIVec(:,1)=TIVec(:,1)-mean(TIVec(:,1));
% Ensure theta is positive
% TIVec(:,1)=TIVec(:,1)-min(TIVec(:,1))+1;

%% Generate identities
R = mnrnd(params.n,params.gamma);
TIVec(:,2) = [1.*ones(1,R(1)), 0.*ones(1,R(2)), -1.*ones(1,R(3))]';

%% Prepare consolidation
% We will sort a number of thetas by their size.
nrSort=round(abs(params.cons)*params.n);
% Randomly select thetas
[~,idx]=datasample(TIVec(:,1),nrSort,'Replace',false);
% Sort this part of the sample
if params.cons < 0 % Consolidation < 0: We sort ascending, since slackers at the bottom
    vec=sort(TIVec(idx,1));
    sidx=sort(idx);
    TIVec(sidx,1)=vec;
else % Consolidation > 0: We sort descending, since climbers at the top
    vec=sort(TIVec(idx,1),'descend');
    sidx=sort(idx);
    TIVec(sidx,1)=vec;
end

% Prepare vectors
theta=TIVec(:,1);
identity=TIVec(:,2);

%% Create helper functions
% Calculate range of theta
firm.thetaRange = abs(max(theta)-min(theta));
% Perceived Benefit functions
gemA=1;
gemL=1;
firm.psiWatcher=@(theta_i,theta_j) -gemL.*abs(theta_j-theta_i)+firm.thetaRange;
firm.psiClimber=@(theta_i,theta_j) -gemA.*(theta_i-theta_j);
firm.psiSlacker=@(theta_i,theta_j) -gemA.*(theta_j-theta_i);

%% Create G matrix
%% If pn is set, we create an underlying G Matrix
if params.pn==0
    firm.gMat=ones(params.n,params.n)-eye(params.n,params.n);
else
    firm.gMat=JacksonRogersNW(params.n,params.mn,params.pn, params.m);
end

%% Populate firm-level variables.
% aMat saves the attention choices per time period
firm.aMat=cell(1,params.maxT);
% In these matrices we will save vectors for each t in T
firm.xMat=repmat(theta,1,params.maxT);
% A matrix containing the thetas for each time period
firm.thetaMat=repmat(theta,1,params.maxT);
% A matrix containing the current representations of theta
firm.thetaRep=repmat(theta',params.n,1);
% Type matrix
firm.muMat=repmat(identity,1,params.maxT);
% Embedding matrix
firm.eMat=repmat(params.e,1,params.maxT);
% Initialize attention matrix
firm.aMat{1} = zeros(params.n,params.n);
% Initialize utility matrix
firm.uMat = zeros(params.n,params.maxT);

%% Populate outcome-level variables
firm.diffM=zeros(1,params.maxT);
firm.SM=ones(1,params.maxT);
firm.varSM=zeros(1,params.maxT);
firm.maxDiff=zeros(1,params.maxT);
firm.minDiff=zeros(1,params.maxT);
firm.skew=repmat(skewness(theta),1,params.maxT);
firm.CSRatio=repmat((params.gamma(1))/(params.gamma(1)+params.gamma(3)),1,params.maxT);
% Here we save the consolidation value
firm.startCons=params.cons;
% The actual consolidation achieved within the sample
actualCons=ConsolidationMeasure(theta,identity);
firm.cons=repmat(actualCons.fcons,1,params.maxT);

%% Time period the firm lives in
firm.T=2;

%% Pass along parameters
firm.minEqmT=params.minEqmT;
firm.maxEqmT=params.maxEqmT;
firm.globalsearch=params.globalsearch;
firm.e=params.e;
firm.n=params.n;
firm.rationality=params.rationality;
firm.maxT=params.maxT;
firm.thetaRepShockVar=params.thetaRepShockVar;
% Hash the firm to get the ID of starting values.
firm.firmID=DataHash(firm);

end