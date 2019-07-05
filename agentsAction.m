%%
% agentsAction.m
% runs one T on the given firm
%%
% @param firm - a firm
% @return firm - a firm
% TODO: Fix optimization
%%
function firm = agentsAction(firm)

%% Unpack some common variables
n = firm.n;
theta = firm.thetaMat(:, firm.T);
identity = firm.muMat(:, firm.T);
% Realize current embedding
firm.e = firm.eMat(firm.T);

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
attention = cell(1, firm.maxEqmT);
% In this we will save the X-values
output = repmat(theta, 1, firm.maxEqmT);
% Initialize p_matrix
attention{1} = zeros(n, n);
% Initialize u_matrix
utility = zeros(n, firm.maxEqmT);
% Variable for convergence behavior
utilityVec = zeros(n, 1);
% Initital convergence difference
xdif = 100;

% Set up global maximizer if needed
options = optimset('Algorithm', 'sqp', 'Display', 'none', 'UseParallel', false);
options = optimoptions('fmincon' ,'Algorithm', 'sqp', 'Display', 'none', 'UseParallel', false);
%options = optimoptions('fmincon' ,'Algorithm', 'sqp', 'Display', 'none', 'UseParallel', true);
if firm.conUtil ~= 0
    gs = GlobalSearch('Display', 'off', 'StartPointsToRun', 'all');
else
    gs=0;
end
% Loop over time periods until convergence
for t = 2:(firm.maxEqmT)
    
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
        theta_i = theta(i);
        % Get theta representation
        thetaRep = firm.thetaRep(i, :)';
        % Initialize a_i vector
        prevAi = prevAttention(i, :)';
        curAi = prevAttention(i, :)';
        
        if firm.conUtil == 0
            %% Discrete optimization to find focal peer
            % Using the current representation of theta by agent i
            % Set up objective function
            if identity(i) == 1% Climber
                [curUi, curAi] = ConvexDiscreteChoice(prevAttention, firm.e, 1, theta(i), thetaRep, firm.psiClimber, i, ChoiceCell{i}, nrChoices(i), firm.rationality,1);
            elseif identity(i) == 0% Watcher
                [curUi, curAi] = ConvexDiscreteChoice(prevAttention, firm.e, 1, theta(i), thetaRep, firm.psiWatcher, i, ChoiceCell{i}, nrChoices(i), firm.rationality,1);
            else % Slacker
                [curUi, curAi] = ConvexDiscreteChoice(prevAttention, firm.e, 1, theta(i), thetaRep, firm.psiSlacker, i, ChoiceCell{i}, nrChoices(i), firm.rationality,1);
            end
            
            % Convention is utility is set to negative, same as fmincon
            curUi = -curUi;
        else
            %% Global optimization routine to find focal group
            % Initialize  restricted a_i taking into account GMat
            prevAiRest = prevAi(ChoiceCell{i});
            curAiRest = prevAiRest;
            % Constraints
            A = ConA{i};
            b = Conb{i};
            lb = zeros(size(curAiRest));
            ub = ones(size(curAiRest));
            %                 lb = zeros(size(curAiRest));
            % ub = ones(size(curAiRest));
            % A = ones(size(curAiRest));
            % b = [1];
            % Aeq = [];
            % beq = [];
            % Set up objective function
            if identity(i) == 1% Climber
                fnc = @(aiRest)concaveChoice(aiRest, prevAttention, theta, theta(i), firm.e, firm.psiClimber, firm.rationality, firm.conParam, ChoiceCell{i},i);
            elseif identity(i) == 0% Watcher
                fnc = @(aiRest)concaveChoice(aiRest, prevAttention, theta, theta(i), firm.e, firm.psiWatcher, firm.rationality, firm.conParam, ChoiceCell{i},i);
            else % Slacker
                fnc = @(aiRest)concaveChoice(aiRest, prevAttention, theta, theta(i), firm.e, firm.psiSlacker, firm.rationality, firm.conParam, ChoiceCell{i},i);
            end
            %
            %x = fmincon(fnc,aStart,A,b,Aeq,beq,lb,ub,[],options)';'nonlcon',@(x) nonLinConstraint(x,firm.maxDegree)
            problem = createOptimProblem('fmincon', 'x0', curAiRest, 'objective', fnc, ...
                'lb', lb, 'ub', ub, 'Aineq', A, 'bineq', b, 'options', options);
            [curAiRest, curUi] = run(gs, problem);
            
            curAi = RecoverPi(curAiRest, ChoiceCell{i}, n);
            
            %% Regularize a_i
            % Set extremely small values to zero
            curAi(curAi<0.05)=0;
            curAi = bsxfun(@times, curAi, 1 ./ (sum(curAi)));
            curAi(isnan(curAi)) = 0;
            
            if firm.conParam>1
                %% Discrete optimization to find focal peer under convexity
                % In case global optimization does not catch border solution
                % Using the current representation of theta by agent i
                % Set up objective function
                if identity(i) == 1% Climber
                    [curUi2, curAi2] = ConvexDiscreteChoice(prevAttention, firm.e, 1, theta(i), thetaRep, firm.psiClimber, i, ChoiceCell{i}, nrChoices(i), firm.rationality,firm.conParam);
                elseif identity(i) == 0% Watcher
                    [curUi2, curAi2] = ConvexDiscreteChoice(prevAttention, firm.e, 1, theta(i), thetaRep, firm.psiWatcher, i, ChoiceCell{i}, nrChoices(i), firm.rationality,firm.conParam);
                else % Slacker
                    [curUi2, curAi2] = ConvexDiscreteChoice(prevAttention, firm.e, 1, theta(i), thetaRep, firm.psiSlacker, i, ChoiceCell{i}, nrChoices(i), firm.rationality,firm.conParam);
                end
                
                % Convention is utility is set to negative, same as fmincon
                curUi2 = -curUi2;
                
                aiRest=curAi(ChoiceCell{i});
                aiRest2=curAi2(ChoiceCell{i});
                concaveChoice(aiRest, prevAttention, theta, theta(i), firm.e, firm.psiClimber, firm.rationality, firm.conParam, ChoiceCell{i},i);
                concaveChoice(aiRest2, prevAttention, theta, theta(i), firm.e, firm.psiClimber, firm.rationality, firm.conParam, ChoiceCell{i},i);
                
                if round(-curUi2,2) > round(-curUi,2)
                    %disp(['Convex Discrete utility for agent ',num2str(i),' is higher with ',num2str(-curUi2),' > ',num2str(-curUi)]);
                    curAi = curAi2;
                    curUi = curUi2;
                elseif round(-curUi2,2) < round(-curUi,2)
                    %disp(['Convex Continous utility for agent ',num2str(i),' is higher with ',num2str(-curUi),' > ',num2str(-curUi2)]);
                    
                end
                
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
        
        %utilityVec(i)=-curUi;
        
        % Save P_i and X_i in the matrices
        curAttention(i, :) = curAi';
    end
    
    % curAttention now has all period best-replies
    % calculate the  current x
    % based on real theta
    x = XFOCSPNE(curAttention, 1, theta, firm.e);
    
    % Save values to matrices
    output(:, t) = x;
    attention{t} = curAttention;
    utility(:, t) = utilityVec(:);
    
    % Check convergence. Fluctuations should be small enough, and the
    % Simulation has run for many periods already
    udif = abs(max(abs(utility(:, t) - utility(:, t - 1))));
    xdif = abs(max(abs(output(:, t) - output(:, t - 1))));
    pdif = abs(sum(sum(abs(curAttention - prevAttention))));
    finalt = t;
    
    %%% Uncomment "disp" to see convergence per iteration
    
    if (((pdif <= 1.0000e-5) && (udif <= 1.0000e-5) && (xdif <= 1.0000e-5)) && (t > firm.minEqmT))
        %disp(['Finish on t ',num2str(t),' with current convergence in utils ',num2str(udif),', in x ',num2str(xdif),' in p ',num2str(pdif)])
        break;
    end
    
    % If half-time is achieved and p values fluctuate only by tiny amounts,
    % end the simulation
    if ((abs(pdif) <= 1.0000e-6 * t) && (t >= firm.maxEqmT / 2))
        warning('Convergence achieved due to half-time');
        break;
    end
    
end

%% Update Firm
% Attention
firm.aMat{firm.T} = attention{finalt};
% Output
firm.xMat(:, firm.T) = output(:, finalt);
% utility
firm.uMat(:, firm.T) = utility(:, finalt);

end
