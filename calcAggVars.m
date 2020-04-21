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


end