function firm = agentsActionPFT(firm)
% agentsActionPFT - runs one T on the given firm
%%
% @param firm - a firm
% @return firm - a firm
%%
%% Unpack some common variables
n = firm.n;
% Row of thetas and prior thetas
theta = firm.thetaMat(:, firm.T);
%priorTheta = firm.thetaMat(:, firm.T-1);
% Row of  prior outputs
priorX = firm.xMat(:, firm.T-1);
identity = firm.muMat(:, firm.T);
% Realize current embedding
e=firm.e;
% Set end point
if firm.rationality==0
    endpoint=2;
else % We have disabled rationality for now, so one optimization run is sufficient
    endpoint=firm.maxEqmT;
end
%% Generate choice sets
ChoiceCell = {};
nrChoices = zeros(1, n);
[ChoiceCell, nrChoices] = GetChoiceSet(firm.gMat, n);

%% Generate Constraint matrices
ConA = {};
Conb = {};
[ConA, Conb] = CalcConstraints(n, nrChoices);

%% Create Convergence matrices
% Saving
attention = cell(1, endpoint);
% In this we will save the X-values
output = repmat(theta, 1, endpoint);
% Initialize p_matrix
attention{1} = zeros(n, n);
% Initialize u_matrix
utility = zeros(n, endpoint);
% Variable for convergence behavior
utilityVec = zeros(n, 1);



% Loop over time periods until convergence

for t = 2:(endpoint)
    
    % Initialize cell content for period t
    attention{t} = zeros(n, n);
    
    % Copy matrices such that parallel works better
    curAttention = attention{t};
    prevAttention = attention{t - 1};
    
    %%
    % If not using parfor outside of this function (ie. simulating one
    % company), you can enable parfor here to work through the employees
    % faster. For many companies, it is faster to do this outside this
    % function.
    % parfor i = 1:n
    % for i = 1:n
    for i = 1:n
        
        %% Initialize actor level stuff
        curUi = 0;
        % Theta
        % Get theta representation
        thetaRep = theta;
        % Initialize a_i vector
        prevAi = prevAttention(i, :)';
        curAi = prevAttention(i, :)';
        if firm.conUtil == 0
            %% Discrete optimization to find focal peer
            % Using the current representation of theta by agent i
            % Set up objective function
            if identity(i) == 1% Self-Improver
                [curUi, curAi] = ConvexDiscreteChoicePFT(prevAttention,priorX, firm.e, thetaRep, firm.psiClimber, i, ChoiceCell{i}, nrChoices(i), firm.rationality,1, firm.maxDegree);
            elseif identity(i) == 0% Self-Assessor
                [curUi, curAi] = ConvexDiscreteChoicePFT(prevAttention,priorX, firm.e,  thetaRep, firm.psiWatcher, i, ChoiceCell{i}, nrChoices(i), firm.rationality,1, firm.maxDegree);
            else % Self-Enhancer
                [curUi, curAi] = ConvexDiscreteChoicePFT(prevAttention,priorX, firm.e,  thetaRep, firm.psiSlacker, i, ChoiceCell{i}, nrChoices(i), firm.rationality,1, firm.maxDegree);
            end
        else
            %% Concave utility detected.
            % We run both discrete (fast) and global optimization (slow)
            % where the former checks border cases of the latter
            % We then pick the most beneficial for the employee!
            %% Discrete optimization to find focal peer
            % Using the current representation of theta by agent i
            % Set up objective function
            if identity(i) == 1% Self-Improver
                [curUiD, curAiD] = ConvexDiscreteChoicePFT(prevAttention,priorX, firm.e, thetaRep, firm.psiClimber, i, ChoiceCell{i}, nrChoices(i), firm.rationality,1, firm.maxDegree);
            elseif identity(i) == 0% Self-Assessor
                [curUiD, curAiD] = ConvexDiscreteChoicePFT(prevAttention,priorX, firm.e,  thetaRep, firm.psiWatcher, i, ChoiceCell{i}, nrChoices(i), firm.rationality,1, firm.maxDegree);
            else % Self-Enhancer
                [curUiD, curAiD] = ConvexDiscreteChoicePFT(prevAttention,priorX, firm.e,  thetaRep, firm.psiSlacker, i, ChoiceCell{i}, nrChoices(i), firm.rationality,1, firm.maxDegree);
            end
            %% Global optimization routine
            % Account for non-linear effects in monitoring, either concave or convex
            % Set up objective function
            if identity(i) == 1% Climber
                [curUi, curAi] = ConcaveChoicePFT(prevAttention,priorX, firm.e, thetaRep, firm.psiClimber, i, ChoiceCell{i}, nrChoices(i), firm.rationality,firm.conParam, firm.maxDegree,ConA, Conb);
            elseif identity(i) == 0% Watcher
                [curUi, curAi] = ConcaveChoicePFT(prevAttention,priorX, firm.e,  thetaRep, firm.psiWatcher, i, ChoiceCell{i}, nrChoices(i), firm.rationality,firm.conParam, firm.maxDegree,ConA, Conb);
            else % Slacker
                [curUi, curAi] = ConcaveChoicePFT(prevAttention,priorX, firm.e,  thetaRep, firm.psiSlacker, i, ChoiceCell{i}, nrChoices(i), firm.rationality,firm.conParam, firm.maxDegree,ConA, Conb);
            end 
            if (round(curUiD,5) < round(curUi,5)) && round(curUiD,5) <0
                %disp(['curUID with', num2str(-curUiD), ' higher than ', num2str(-curUi)])
                %disp(['curAID='])
                %disp(curAiD)
                %disp(['curAI='])
                %disp(curAi)
                %disp(['i=', num2str(i)])
                if round(curAi,5)==round(curAiD,5)
                    error("Uh oh, something is wrong with global opt!")
                end
                curAi=curAiD;
                curUi=curUiD;
            end
        end
        
        %% Post optimization checks
        %% Constraint check
        % If Fmincon doesn't satisfy the constraints
        % Which sometimes happens for tiny values, this keeps the simulation from spinning out
        % of control
        
        if (min(curAi < 0)) || max(curAi > 1)
            %disp(['Pstar out of range. Fixing. Occured with g: ', num2str(g), ', identity: ', num2str(identity(i)), ', theta: ', num2str(theta_i), ...
            %        ' in period ', num2str(t), ...
            %        ' with pstar min/max: ', num2str(min(curAi)), ', ', num2str(max(curAi))]);
            curAi(curAi < 0) = 0;
            curAi(curAi > 1) = 1;
        end
        
        if sum(curAi > 1)
            curAi = bsxfun(@times, curAi, 1 ./ (sum(curAi)));
            curAi(isnan(curAi)) = 0;
            % disp(['sum Pstar larger than 1. Fixing. Occured with g: ', num2str(g), ', identity: ', num2str(identity(i)), ', theta: ', num2str(theta_i), ...
            %         ' in period ', num2str(t), ...
            %         ' with sum pstar: ', num2str(sum(curAi))]);
        end
        
        %% Option to disconnect
        % Disable if disconnection is not part of possible solution
        % but allowed in model, so we allow for disconnection if no
        % positive utility can be found
        
        if 0 >- curUi
            curAi = zeros(length(curAi), 1);
        end
        
        %% Final calculations
        % Save utility vector
        prevUi = utility(i, t - 1);
        
        if prevUi == -curUi
            %disp(['Last period same for i=',num2str(i),' and ',num2str(-curUi)]);
            curAi = prevAi;
            utilityVec(i) = prevUi;
        elseif prevUi >- curUi
            %disp(['Last period better for i=',num2str(i),' diff is ',num2str(abs(-curUi-curUit_prev))]);
            curAi = prevAi;
            utilityVec(i) = prevUi;
        else
            utilityVec(i) = -curUi;
        end
        
        % Save P_i and X_i in the matrices
        curAttention(i, :) = curAi';
    end
    
    
    % curAttention now has all period best-replies
    % PFT Version: x_t_1 for j != i
    x_pft=XFOCPFT(priorX,curAttention,theta,e);
    x=x_pft;
    % Save values to matrices
    output(:, t) = x;
    attention{t} = curAttention;
    utility(:, t) = utilityVec(:);
    finalt = t;
end

%% Update Firm
% Attention
firm.aMat{firm.T} = attention{finalt};
% Output
firm.xMat(:, firm.T) = output(:, finalt);
% utility
firm.uMat(:, firm.T) = utility(:, finalt);

end
