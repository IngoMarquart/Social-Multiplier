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
% @param: gMethod
% @param: pn,mn - Exogenous social matrix configuration
% @param: cons - Consolidation
% @return: firm - a firm
%%
function firm=initFirmSingle(theta,mu,params)

% Set the random stream
s = RandStream('mcg16807','Seed',params.m);
RandStream.setGlobalStream(s);
firm.m=params.m;
% First period is initialization
params.maxT=params.maxT+1;
firm.T=2;
% Prepare vectors
identity=mu;
firm.Talpha=0;
firm.Tbeta=0;
firm.shufflePositions=0;
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
firm.gMethod=params.gMethod;
if params.gMethod=="JR" % Jackson Rogers Social network
    p=min(0.1,rand);
    m=ceil(rand*min(10,params.n));
    Mnr=ceil(m/(p*params.n));
    firm.gMn=round(Mnr*rand);
    firm.gMr=Mnr-firm.gMn;
    firm.gPn=p;
    firm.gPr=p;
    firm.gMat=JacksonRogersNW(params.n,firm.gMn,firm.gPn,firm.gMr,firm.gPr, params.m);
elseif params.gMethod=="Task" % Task network
    % We need to intialize modularity, cluster and symmetry
    % Symmetry:
    firm.gSymmetry=0;
    if rand>0.5
        firm.gSymmetry=1;
    end
    firm.gAssembly=0;
    firm.gModularity=0.5*rand();
    firm.gCluster=randi([2 floor(params.n./3)],1,1);
    [firm.gMat,~,firm.gLinks]=TaskNetwork(params.n,firm.gCluster,firm.gModularity, firm.gSymmetry, params.m);
elseif params.gMethod=="TaskAssembly"
    % We need to intialize modularity, cluster and symmetry
    % Forward/Backward Symmetry:
    firm.gSymmetry=0;
    if rand>0.5
        firm.gSymmetry=1;
    end
    % Assemblyliness
    firm.gAssembly=rand;
    % Modularity - number of ties
    firm.gModularity=rand();
    firm.gCluster=randi([2 floor(params.n./3)],1,1);
    [firm.gMat,~,firm.gLinks]=TaskNetworkAss(params.n,firm.gCluster,firm.gModularity,firm.gAssembly, firm.gSymmetry, params.m);
else
    firm.gMat=ones(params.n,params.n)-eye(params.n,params.n);
end

%% Generate CEO type

if params.ceoAct=="Random" % Shuffle CEO Type
    type=unidrnd(3,1,1);
    switch type
        case 1 % Embed
            params.ceoAct="Embed";
        case 2 % Embed
            params.ceoAct="Decouple";
            %case 3 % Off
            %     params.ceoAct="Off";
        case 3 % lowEmbed
            params.ceoAct="LowEmbed";
    end
    firm.ceoAct=params.ceoAct;
    
end
firm.startCeoAct=params.ceoAct;

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
% Hash the firm to get the ID of starting values.
firm.firmID=DataHash(firm);


%% Calculate network measures on G
n=params.n;
G=firm.gMat;
if G==G'
    GPgraph=graph(G);
else
    GPgraph=digraph(G);
end

% Network density
PC=(n*(n-1))./2;
EC=height(GPgraph.Edges);
firm.gDensity=EC/PC;

% Average path length
% We calculate each path length separately for each
% connected component in the symmetric graph of A
% an then use a weighted average
bincell = conncomp(GPgraph, 'OutputForm', 'cell');
nrsubgraphs = length(bincell);
nrNodes=zeros(1,nrsubgraphs);
sumDist=zeros(1,nrsubgraphs);
nrEdges=zeros(1,nrsubgraphs);
weightfactor=zeros(1,nrsubgraphs);
for ii = 1:nrsubgraphs
    subg=subgraph(GPgraph, bincell{ii});
    DM=distances(subg);
    nrNodes(ii)=length(DM);
    sumDist(ii)=sum(sum(DM));
    weightfactor(ii)=nrNodes(ii)/n;
    nrEdges(ii)=(nrNodes(ii).*(nrNodes(ii)-1));
    % Zero weight for single nodes (Path Length not defined)
    if(nrNodes(ii) == 1)
        weightfactor(ii)=0;
        nrEdges(ii)=1;
    end
    % Renormalize weights to sum to one
    weightfactor=weightfactor./sum(weightfactor);
    
end
% Average path length is given by the sum of distances divided by the
firm.gAvgPathLength=((1./nrEdges).*weightfactor)*sumDist';
firm.gNrComponents=nrsubgraphs;

% Largest eigenvector
firm.gMaxEV=max(eig(G));

% Clustering & average Degree
deg = sum(G, 2); %Determine node degrees
cn = diag(G*triu(G)*G); %Number of triangles for each node
%The local clustering coefficient of each node
c = zeros(size(deg));
c(deg > 1) = 2 * cn(deg > 1) ./ (deg(deg > 1).*(deg(deg > 1) - 1));
firm.gAvgClustering=mean(c);
firm.gAvgDegree=mean(deg);


% Diameter
d=distances(GPgraph);
d(~isfinite(d))=1000;
% Eccentricity
ev=max(d,[],2);
% Radius and diameter
firm.gRadius=min(ev(ev>0));
firm.gDiameter=max(ev(ev>0));
end