function firm=initFirm(params)
% initFirm.m - For given parameter ranges, initiate a firm
%%
% @param: params - struct of parameters
% @return: firm - a firm
%%
% Set the random stream
s = RandStream('mcg16807','Seed',params.m);
RandStream.setGlobalStream(s);
firm.m=params.m;


% First period is initialization
params.maxT=params.maxT+1;

%% Generation of types and quality vectors
% Matrix to hold identity and thetas
TIVec = zeros(params.n,2);
pd = makedist('Beta',params.thetaD(1),params.thetaD(2));

firm.Talpha=params.thetaD(1);
firm.Tbeta=params.thetaD(2);
% If exogenous scaling is desired use this and set fourth and third parameters
%TIVec(:,1) = (random(pd,n,1)-pd.mean).*thetaD(4)+thetaD(3);

% We scale automatically by fixing variance.
TIVec(:,1) = random(pd,params.n,1);
c=sqrt(1/var(TIVec(:,1)))*params.thetaVar;
% This fixes variance to thetaVar
TIVec(:,1) =TIVec(:,1).*c;
% De-mean theta, set mean to thetaMean
TIVec(:,1)=TIVec(:,1)-mean(TIVec(:,1))+params.thetaMean;

%% Generate identities
R = mnrnd(params.n,params.gamma);
TIVec(:,2) = [1.*ones(1,R(1)), 0.*ones(1,R(2)), -1.*ones(1,R(3))]';

%% Prepare consolidation
% We will sort a proportion of thetas by their size.
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

%% Shuffling of positions
firm.shufflePositions="Shuffled";
idx = randperm(size(TIVec,1));
TIVec=TIVec(idx,:);

% Prepare vectors
theta=TIVec(:,1);
identity=TIVec(:,2);

%% Create helper functions
% Calculate range of theta
firm.thetaRange = abs(max(theta)-min(theta));
% Perceived Benefit functions
gemA=firm.thetaRange;
firm.psiWatcher=@(theta_i,theta_j) -gemA.*(abs(theta_j-theta_i)-firm.thetaRange);
firm.psiClimber=@(theta_i,theta_j) -gemA.*(theta_i-theta_j);
firm.psiSlacker=@(theta_i,theta_j) -gemA.*(theta_j-theta_i);

%% Create G matrix DISABLED
firm.gMethod="Full";
firm.gMat=ones(params.n,params.n)-eye(params.n,params.n);

%% Generate CEO type DISABLED
firm.ceoAct="Off";
firm.startCeoAct="Off";

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
firm.sumX=zeros(1,params.maxT);
firm.avgTheta=zeros(1,params.maxT);
firm.avgX=zeros(1,params.maxT);
firm.maxX=zeros(1,params.maxT);
firm.minX=zeros(1,params.maxT);
firm.maxTheta=zeros(1,params.maxT);
firm.minTheta=zeros(1,params.maxT);
firm.varX=zeros(1,params.maxT);
firm.varTheta=zeros(1,params.maxT);
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
firm.gamma=params.gamma;
firm.ceoStartT=params.ceoStartT;
firm.learningRate=params.learningRate;
firm.globalsearch=params.globalsearch;
firm.e=params.e;
firm.n=params.n;
firm.rationality=params.rationality;
firm.maxT=params.maxT;
firm.thetaRepShockVar=params.thetaRepShockVar;
firm.ceoAct=params.ceoAct;
firm.conUtil=params.conUtil;
firm.conParam=params.conParam;
firm.maxDegree=params.maxDegree;
firm.probTypeSwitch=params.probTypeSwitch;
% Hash the firm to get the ID of starting values.
firm.firmID=DataHash(firm);

end