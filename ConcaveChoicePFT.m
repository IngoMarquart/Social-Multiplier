%%
% DiscreteChoice
% Assuming that each actor picks only one peer, this function find the one
% giving the most utility
% @param: a_t_1 - Attention matrix of last time period
% @param: e,delta - Embedding and scaling variable
% @param: theta_i, theta - theta vector and theta of i
% @param: Psi - Expected benefit function handle
% @param: Choice - Vector given the indecies of valid peers (where connection in G exists)
% @param: nrChoice - Number of choices
% @return: util - Utility value
% @return: a_i_star - Attention choice in the space of all n actors
%%
function [util,a_i_star]=ConcaveChoicePFT(a_t_1,x_t_1,e, theta, Psi,i,Choice,nrChoice, rationality, conParam,maxDegree)

theta=theta(:);
n=length(theta);
% Define One vector
ez=ones(n,1);

% Default values
util=0;
a_i_star=zeros(n,1);

% Constants
PsiVec = Psi(theta(i),theta);

options = optimoptions('fmincon' ,'Algorithm', 'sqp', 'Display', 'none', 'UseParallel', false);
gs = GlobalSearch('Display', 'off', 'StartPointsToRun', 'all');


fnc = @(aiRest)valueFunction(aiRest,x_t_1,theta,e,PsiVec,i,maxDegree,conParam,Choice,nrChoice)

problem = createOptimProblem('fmincon', 'x0', curAiRest, 'objective', fnc, ...
    'lb', lb, 'ub', ub, 'Aineq', A, 'bineq', b, 'options', options);
[curAiRest, util] = run(gs, problem);

a_i_star = RecoverPi(curAiRest, Choice, n);

%% Regularize a_i
% Set extremely small values to zero
a_i_star(a_i_star<0.05)=0;
a_i_star = bsxfun(@times, a_i_star, 1 ./ (sum(a_i_star)));
a_i_star(isnan(a_i_star)) = 0;

end

function util=valueFunction(a_i_Rest,x_t_1,theta,e,PsiVec,i,maxDegree,conParam,Choice,nrChoice)
    % Norm vectors to correspond to column vectors for now
    theta=theta(:);
    aiRest=aiRest(:);
    % First we need to recover the true attention vector
    % from the restricted choices given G
    n = length(theta);
    a_i = RecoverPi(aiRest, Choice, n);
    util=utilityPFT(x(i),x,a_i,theta,e,PsiVec,i,maxDegree,conParam)
end