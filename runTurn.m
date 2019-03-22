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

%% CEO Action
%firm = ceoAction(firm);

%% Agents run
firm = agentsAction(firm);

%% Post-turn transformation
%firm= postTransform(firm);

%% Calculate aggregate Variables
firm= calcAggVars(firm);

end