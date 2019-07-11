%%
%   runTurn.m
%   runs one T on the given firm
%%
% @param: firm- a firm
% @return: firm - a firm
%%
function firm = runTurn(firm)

%% Pre-turn Transform
%firm = preTransform(firm);
firm.e=firm.eMat(firm.T);


%% CEO Action
if firm.ceoAct~="Off"
firm = ceoAction(firm);
end

%% Confirm embedding size - Matlab limit on float size to invert matrices
if firm.e>=1.0000e+14
    firm.e=1.0000e+14;
end
%% Agents run
firm = agentsAction(firm);

%% Post-turn transformation
firm= postTransform(firm);
firm.eMat(firm.T)=firm.e;


%% Calculate aggregate Variables
firm= calcAggVars(firm);

end