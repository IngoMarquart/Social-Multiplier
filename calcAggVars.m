%%
% % GraphNetwork
% This function graphs the G network and highlights the A network on top of it.
% @param: Gmat - Symmetric G matrix
% @param: SPmat - Symmetric P matrix
% @param: identity, theta, X: nx1 vectors of identity, theta and X
%%

function firm= calcAggVars(firm)
% Unpack firm variables
PMat=firm.aMat{firm.T};
X=firm.xMat(:,firm.T);
theta=firm.thetaMat(:,firm.T);
identity=firm.muMat(:,firm.T);

%% TODO fix here
NrClimbers=length(identity(identity==1));
NrWatchers=length(identity(identity==0));
NrSlackers=length(identity(identity==-1));


%% Populate outcome-level variables
firm.diffM(firm.T)=mean(X-theta);
firm.varSM(firm.T)=var(X./theta);
firm.SM(firm.T)=mean(X./theta);
firm.skew(firm.T)=skewness(theta);
firm.CSRatio(firm.T)=NrClimbers/(NrClimbers+NrSlackers);
% The actual consolidation achieved within the sample
actualCons=ConsolidationMeasure(theta,identity);
firm.cons(firm.T)=actualCons.fcons;


end