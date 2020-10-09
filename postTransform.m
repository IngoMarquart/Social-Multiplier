%%
% % postTransform
% This function calculates the transformation and learning in the following steps
% 1.
% @param: firm
%%

function firm=postTransform(firm)

% Set the random stream
s = RandStream('mcg16807','Seed',sum(int32([firm.firmID]))+firm.T);
RandStream.setGlobalStream(s);

%% Adaptation // Learning
firm.thetaMat(:,firm.T+1)=(1-firm.learningRate).*firm.thetaMat(:,firm.T)+firm.learningRate.*firm.xMat(:,firm.T);

%% Shock theta representation
if firm.thetaRepShockVar>0
    pd = makedist('Normal','mu',0,'sigma',firm.thetaRepShockVar);
    thetaShock=random(pd,firm.n,firm.n);
    thetaShock=thetaShock-eye(size(thetaShock)).*thetaShock;
    
    firm.thetaRep=repmat(firm.thetaMat(:,firm.T+1)',firm.n,1);
    firm.thetaRep=firm.thetaRep+thetaShock;
else
    firm.thetaRep=repmat(firm.thetaMat(:,firm.T+1)',firm.n,1);
end

%% Shock types
if firm.probTypeSwitch>0
    % Generate new identities
    R = mnrnd(firm.n,firm.gamma);
    new_identities = [1.*ones(1,R(1)), 0.*ones(1,R(2)), -1.*ones(1,R(3))]';
    old_identities = firm.muMat(:,firm.T);
    % Get Bernoulli (0/1) sample
    bernVec=binornd(1,firm.probTypeSwitch,firm.n,1);
    % 1-> new identity, 0-> old identity
    new_identities=bernVec.*new_identities+(1-bernVec).*old_identities;

    % update Identity
    firm.muMat(:,firm.T+1)=new_identities;

else
    firm.muMat(:,firm.T+1)=firm.muMat(:,firm.T);
end

end