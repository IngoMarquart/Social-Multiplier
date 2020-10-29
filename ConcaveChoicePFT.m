function [util,a_i_star]=ConcaveChoicePFT(a_t_1,x_t_1,e, theta, Psi,i,Choice,nrChoice, rationality, conParam,maxDegree,ConA, Conb)
% ConcaveChoice - Global optimization routine if utility is concave wrt. monitoring weights
%%
% @param: a/x_t_1 - Attention/output matrix of last time period
% @param: e - Embedding and scaling variable
% @param: theta_i, theta - theta vector and theta of i
% @param: Psi - Expected benefit function handle
% @param: Choice - Vector given the indecies of valid peers (where connection in G exists)
% @param: nrChoice - Number of choices
% @param: conParam - v parameter of utility function
% @param: maxDegree - bound to degree
% @param: conA,ConB - Constraints for optimization
% @return: util - Utility value
% @return: a_i_star - Attention choice in the space of all n actors
%%
theta=theta(:);
n=length(theta);
% Define One vector
ez=ones(n,1);

% Default values
util=0;
a_i_star=zeros(n,1);

%% Generate Constraint matrices
curAiRest = a_i_star(Choice);

A = ConA{i};
b = Conb{i};
lb = zeros(size(curAiRest));
ub = ones(size(curAiRest));

% Constants
PsiVec = Psi(theta(i),theta);

options = optimoptions('fmincon' ,'Algorithm', 'sqp', 'Display', 'none', 'UseParallel', false);
gs = GlobalSearch('Display', 'off', 'StartPointsToRun', 'all');
ms = MultiStart('Display', 'off');


fnc = @(aiRest)valueFunction(aiRest,a_t_1,x_t_1,theta,e,PsiVec,i,maxDegree,conParam,Choice,nrChoice);

problem = createOptimProblem('fmincon', 'x0', curAiRest, 'objective', fnc, ...
    'lb', lb, 'ub', ub, 'Aineq', A, 'bineq', b, 'options', options);

% Global search
[curAiRest, util] = run(gs, problem);

% Multi search with startpoints given by discrete mixes
startpoints=eye(length(curAiRest));
startpoints=CustomStartPointSet(startpoints);
[curAiRest, util] = run(ms, problem,startpoints);

a_i_star = RecoverPi(curAiRest, Choice, n);

%% Regularize a_i
% Set extremely small values to zero
a_i_star(a_i_star<0.05)=0;
a_i_star = bsxfun(@times, a_i_star, 1 ./ (sum(a_i_star)));
a_i_star(isnan(a_i_star)) = 0;

end

function util=valueFunction(a_i_Rest,a_t_1,x_t_1,theta,e,PsiVec,i,maxDegree,conParam,Choice,nrChoice)
    % Norm vectors to correspond to column vectors for now
    theta=theta(:);
    aiRest=a_i_Rest(:);
    % First we need to recover the true attention vector
    % from the restricted choices given G
    n = length(theta);
    a_i = RecoverPi(aiRest, Choice, n);
    a=a_t_1;
    if (min(a_i < 0)) || max(a_i > 1)
        %disp(['Pstar out of range. Fixing. Occured with g: ', num2str(g), ', identity: ', num2str(identity(i)), ', theta: ', num2str(theta_i), ...
        %        ' in period ', num2str(t), ...
        %        ' with pstar min/max: ', num2str(min(curAi)), ', ', num2str(max(curAi))]);
        a_i(a_i < 0) = 0;
        a_i(a_i > 1) = 1;
    end

    if sum(a_i > 1)
        a_i = bsxfun(@times, a_i, 1 ./ (sum(a_i)));
        a_i(isnan(a_i)) = 0;
    end
    % Assemble new P
    if sum(isnan(a_i))>0
        util=100000;
    else
        a(i,:)=a_i';    
        x_new=XFOCPFT(x_t_1,a,theta,e);
        x=x_t_1;
        x(i)=x_new(i);
        util=-utilityPFT(x(i),x,a_i,theta,e,PsiVec,i,maxDegree,conParam);
        sx=size(util);
    end
    if length(util)>1
        display(x(i))
        display(x)
        display(a_i)
        display(util)
        display(theta)
        display(PsiVec)
        display(i)
        display(conParam)
        
        error("Whoops")
        util=util(1);
    end
end