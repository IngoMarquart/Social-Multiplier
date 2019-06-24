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
    firm.thetaMat(:,firm.T+1)=0.0.*firm.thetaMat(:,firm.T)+1.*firm.xMat(:,firm.T);
    
    %% Shock theta representation
    if firm.thetaRepShockVar>0
    pd = makedist('Normal','mu',0,'sigma',firm.thetaRepShockVar);
    thetaShock=random(pd,firm.n,firm.n);
    thetaShock=thetaShock-eye(size(thetaShock)).*thetaShock;
    
    firm.thetaRep=repmat(firm.thetaMat(:,1)',firm.n,1);
    firm.thetaRep=firm.thetaRep+thetaShock;
    end

end