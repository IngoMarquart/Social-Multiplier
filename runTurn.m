function firm = runTurn(firm)
% runTurn - runs one T on the given firm
%%
% @param: firm- a firm
% @return: firm - a firm
%%
%% Pre-turn Transform
%firm = preTransform(firm);
firm.e=firm.eMat(firm.T);

%% Confirm embedding size - Matlab limit on float size to invert matrices
if firm.e>=1.0000e+14
    firm.e=1.0000e+14;
end
%% Agents run
firm = agentsActionPFT(firm);

%% Post-turn transformation
firm= postTransform(firm);
firm.eMat(firm.T)=firm.e;
firm.eMat(firm.T+1)=firm.e;

%% Calculate aggregate Variables
firm= calcAggVars(firm);

end