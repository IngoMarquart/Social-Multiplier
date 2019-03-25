%%
%   agentsAction.m
%   runs one T on the given firm
%%
% @param: firm- a firm
% @return: firm - a firm
%%
function firm = agentsAction(firm)

%% Unpack some common variables
n=firm.n;
theta=firm.thetaMat(:,firm.T);
identity=firm.muMat(:,firm.T);

%% Generate choice sets
ChoiceCell={};
nrChoices=zeros(1,n);
[ChoiceCell, nrChoices]=GetChoiceSet(firm.gMat,n);


%% Generate Constraint matrices
ConA={};
Conb={};
[ConA,Conb] = CalcConstraints(n, nrChoices);

%% Create Convergence matrices
% Saving
attention=cell(1,firm.maxEqmT);
% In this we will save the X-values
output=repmat(theta,1,firm.maxEqmT);
% Initialize p_matrix
attention{1} = zeros(n,n);
% Initialize u_matrix
utility = zeros(n,firm.maxEqmT);
% Variable for convergence behavior
utilityVec=zeros(n,1);
% Initital convergence difference
xdif=100;


% Loop over time periods until convergence
for t = 2:(firm.maxEqmT)
    
    % Initialize cell content for period t
    attention{t}=zeros(n,n);
    
    % Copy matrices such that parallel works better
    curAttention=attention{t};
    prevAttention=attention{t-1};
    
    %%
    % If not using parfor outside of this function (ie. simulating one
    % company), you can enable parfor here to work through the employees
    % faster. For many companies, it is faster to do this outside this
    % function.
    % parfor i = 1:n
    % for i = 1:n
    for i = 1:n
        
        %% Initialize actor level stuff
        curUi=0;
        % Constraints
        A=ConA{i};
        b=Conb{i};
        % Theta
        theta_i = theta(i);
        % Get theta representation
        thetaRep=firm.thetaRep(i,:)';
        % Initialize a_i vector
        prevAi=prevAttention(i,:)';
        curAi=prevAttention(i,:)';
        % Initialize  restricted a_i taking into account GMat
        prevAiRest=prevAi(ChoiceCell{i});
        curAiRest=prevAiRest;
        
        
        %% Discrete optimization
        % Using the current representation of theta by agent i
        % Set up objective function
        if identity(i)==1 % Climber
            [curUi,curAi]=DiscreteChoice(prevAttention,firm.e,1,theta(i), thetaRep, firm.psiClimber,i,ChoiceCell{i},nrChoices(i), firm.rationality);
        elseif identity(i)==0 % Watcher
            [curUi,curAi]=DiscreteChoice(prevAttention,firm.e,1,theta(i), thetaRep, firm.psiWatcher,i,ChoiceCell{i},nrChoices(i), firm.rationality);
        else % Slacker
            [curUi,curAi]=DiscreteChoice(prevAttention,firm.e,1,theta(i), thetaRep, firm.psiSlacker,i,ChoiceCell{i},nrChoices(i), firm.rationality);
        end
        % Convention is utility is negative for fmincon
        curUi=-curUi;
        
        
        
        %% Post optimization checks
        %% Constraint check
        % If Fmincon doesn't satisfy the constraints
        % Which sometimes happens for tiny values, this keeps the simulation from spinning out
        % of control
        
        if (min(curAi < 0)) || max(curAi >1)
            disp(['Pstar out of range. Fixing. Occured with g: ',num2str(g),', identity: ',num2str(identity(i)),', theta: ',num2str(theta_i), ...
                ' in period ',num2str(t), ...
                ' with pstar min/max: ',num2str(min(curAi)),', ',num2str(max(curAi))]);
            curAi(curAi<0)=0;
            curAi(curAi>1)=1;
        end
        if sum(curAi>1)
            curAi = bsxfun(@times, curAi, 1./(sum(curAi)));
            curAi(isnan(curAi))=0;
            disp(['sum Pstar larger than 1. Fixing. Occured with g: ',num2str(g),', identity: ',num2str(identity(i)),', theta: ',num2str(theta_i), ...
                ' in period ',num2str(t), ...
                ' with sum pstar: ',num2str(sum(curAi))]);
        end
        
        
        %% Option to disconnect
        % Disable if disconnection is not part of possible solution
        % but allowed in model, so we allow for disconnection if no
        % positive utility can be found
        
        if 0 > -curUi
            curAi=zeros(length(curAi),1);
        end
        
        %% Final calculations
        % Save utility vector
        prevUi=utility(i,t-1);
        
        
        if prevUi == -curUi
            %disp(['Last period same for i=',num2str(i),' and ',num2str(-curUi)]);
            %curAi=(pmat_prev(i,:)'+pmat_prev2(i,:)'+pmat_prev3(i,:)'+pmat_prev4(i,:)'+curAi)/5;
            % curAi=(pmat_prev(i,:)'+curAi)/2;
            % utilityVec(i)=(curUit_prev-curUi)/2;
            curAi=prevAi;
            utilityVec(i)=prevUi;
        elseif prevUi > -curUi
            %disp(['Last period better for i=',num2str(i),' diff is ',num2str(abs(-curUi-curUit_prev))]);
            % curAi=(pmat_prev(i,:)'+pmat_prev2(i,:)'+pmat_prev3(i,:)'+pmat_prev4(i,:)'+curAi)/5;
            % utilityVec(i)=(curUit_prev-curUi)/2;
            curAi=prevAi;
            utilityVec(i)=prevUi;
        else
            utilityVec(i)=-curUi;
        end
        %utilityVec(i)=-curUi;
        
        
        % Save P_i and X_i in the matrices
        curAttention(i,:)=curAi';
    end
    
    
    % curAttention now has all period best-replies
    % calculate the  current x
    % based on real theta
    x=XFOCSPNE(curAttention,1,theta,firm.e);
    
    % Save values to matrices
    output(:,t)=x;
    attention{t}=curAttention;
    utility(:,t)=utilityVec(:);
    
    % Check convergence. Fluctuations should be small enough, and the
    % Simulation has run for many periods already
    udif = max(abs(utility(:,t)-utility(:,t-1)));
    xdif = max(abs(output(:,t)-output(:,t-1)));
    pdif = sum(sum(abs(curAttention-prevAttention)));
    finalt=t;
    
    %%% Uncomment "disp" to see convergence per iteration
    %disp([' t ',num2str(t),' with current convergence in utils ',num2str(udif),', in x ',num2str(xdif),' in p ',num2str(pdif)])
    if (((pdif <= 1.0000e-4) && (udif <= 1.0000e-5) && (xdif<= 1.0000e-5)) && (t>2*firm.minEqmT))
        %disp(['Finish on t ',num2str(t),' with current convergence in utils ',num2str(udif),', in x ',num2str(xdif),' in p ',num2str(pdif)])
        break;
    end
    if (((pdif <= 1.0000e-6) && (udif <= 1.0000e-6) && (xdif<= 1.0000e-10)) && (t>firm.minEqmT))
        %disp(['Finish on t ',num2str(t),' with current convergence in utils ',num2str(udif),', in x ',num2str(xdif),' in p ',num2str(pdif)])
        break;
    end
    % If half-time is achieved and p values fluctuate only by tiny amounts,
    % end the simulation
    if  ((abs(pdif)<= 1.0000e-6*t) && (t>=firm.maxEqmT/2))
        warning('Convergence achieved due to half-time');
        break;
    end
end

%% Update Firm
% Attention
firm.aMat{firm.T}=attention{finalt};
% Output
firm.xMat(:,firm.T)=output(:,finalt);

end