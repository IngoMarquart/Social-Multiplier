%%
% % calcAggVars
% This function calculates the aggregate variables needed to be saved in firm
% @param: firm
%%

function firm= calcAggVars(firm)
% Unpack firm variables
PMat=firm.aMat{firm.T};
X=firm.xMat(:,firm.T);
theta=firm.thetaMat(:,firm.T);
identity=firm.muMat(:,firm.T);

NrClimbers=length(identity(identity==1));
NrWatchers=length(identity(identity==0));
NrSlackers=length(identity(identity==-1));


%% Populate outcome-level variables
firm.avgX(firm.T)=mean(X);
firm.sumX(firm.T)=sum(X);
firm.avgTheta(firm.T)=mean(theta);
firm.diffM(firm.T)=mean(X-theta);
firm.maxX(firm.T)=max(X);
firm.minX(firm.T)=min(X);
firm.maxTheta(firm.T)=max(theta);
firm.minTheta(firm.T)=min(theta);
firm.varX(firm.T)=var(X);
firm.varTheta(firm.T)=var(theta);
firm.maxDiff(firm.T)=max(X-theta);
firm.minDiff(firm.T)=min(X-theta);
firm.SM(firm.T)=mean(X)./mean(theta);
firm.skew(firm.T)=skewness(theta);
firm.CSRatio(firm.T)=NrClimbers/(NrClimbers+NrSlackers);
% The actual consolidation achieved within the sample
actualCons=ConsolidationMeasure(theta,identity);
firm.cons(firm.T)=actualCons.fcons;
firm.NrC=NrClimbers;
firm.NrW=NrWatchers;
firm.NrS=NrSlackers;
% Correlation measure across outputs and thetas
if firm.T>2
    XMatrix=firm.xMat(:,max(firm.T-12,1):firm.T);
    thetaMatrix=firm.thetaMat(:,max(firm.T-12,1):firm.T);
    diffX=diff(XMatrix')';
    diffTheta=diff(thetaMatrix')';
    firm.flucX=mean(mean(abs(diffX')));
    firm.flucT=mean(mean(abs(diffTheta')));
else
    firm.flucX=0;
    firm.flucT=0;
end

%% Network stuff
Pgraph=digraph(PMat);
firm.pgrank=centrality(Pgraph,'pagerank');
firm.pgrank=firm.pgrank(:);
firm.indegree=centrality(Pgraph,'indegree');
firm.indegree=firm.indegree(:);
firm.peerX=PMat*X;
firm.peerX=firm.peerX(:);
[firm.peer,~]=find(PMat');
firm.peer=firm.peer(:);
if length(firm.peer) > firm.n
    firm.peer=zeros(firm.n,1);
    firm.peerMu=zeros(firm.n,1);
else
   firm.peerMu=identity(firm.peer); 
   firm.peerMu=firm.peerMu(:);
end