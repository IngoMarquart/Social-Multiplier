function returndata=SimulateAttSubgame(n,T,gamma, thetaD,ConBen, gemA, gemL, delta, g,m, minT,graphit,globalsearch,convexp,Gmat)
%%
%
%   SimulateAttSubgame
%   Ingo Marquart, Nghi Truong
%   v. 0.9, 17.07.2018
%   For given parameter ranges, this function calculates the
%   SPNE as a steady state of a dynamic simulation
%%

% We used different distributions but have settled for Beta, hence forcing
% it here
BetaDist=1;


% Set Random Stream
s = RandStream('mcg16807','Seed',m);
RandStream.setGlobalStream(s);

%%% Setup
%% Population parameters
% Type distribution
gammaC = gamma(1);
gammaW = gamma(2);
gammaS = gamma(3);


%% Parameters of utility


%% Generation of types and quality vectors
% Matrix to hold identity and thetas
TIVec = zeros(n,2);

%% Generate thetas
if BetaDist==1
    pd = makedist('Beta',thetaD(1),thetaD(2));
    
    % If exogenous scaling is desired use this
    %TIVec(:,1) = (random(pd,n,1)-pd.mean).*thetaD(4)+thetaD(3);
    
    % We scale automatically by fixing variance.
    TIVec(:,1) = random(pd,n,1);
    c=sqrt(1/var(TIVec(:,1)));
    % This fixes variance to 1
    TIVec(:,1) =TIVec(:,1).*c;
    % De-mean theta
    TIVec(:,1)=TIVec(:,1)-mean(TIVec(:,1));
    % Ensure theta is positive
    TIVec(:,1)=TIVec(:,1)-min(TIVec(:,1))+1;
else
    % Lognormal does not allow to vary skew both ways, so we do not use it
    % anymore
    pd = makedist('Lognormal','mu',thetaD(1),'sigma',thetaD(2));
    TIVec(:,1) = random(pd,n,1);
end


%% Generate identities
p = [gammaC, gammaW, gammaS];
R = mnrnd(n,p);
TIVec(:,2) = [1.*ones(1,R(1)), 0.*ones(1,R(2)), -1.*ones(1,R(3))]';

% prepare vectors
theta=TIVec(:,1);
identity=TIVec(:,2);

%% Create helper functions
% calculate range of theta
thetaRange = abs(max(theta)-min(theta));
% Motivation functions
PsiL=@(theta_i,theta_j) -gemL.*abs(theta_j-theta_i)+thetaRange;
PsiA=@(theta_i,theta_j) -gemA.*(theta_i-theta_j);
PsiS=@(theta_i,theta_j) -gemA.*(theta_j-theta_i);



%% Populate Choice matrix

% pmat  saves the P_i choice for each agent, for each t
pmat=cell(1,T);

% In this we will save the X-values
xmat=repmat(theta,1,T);
% Initialize p_matrix
pmat{1} = zeros(n,n);
% Initialize u_matrix
umat = zeros(n,T);
% Variable for convergence behavior
uvec=zeros(n,1);
% Initital convergence difference
xdif=100;

% Loop over time periods until convergence
for t = 2:(T)
    
    % Initialize cell content for period t
    pmat{t}=zeros(n,n);
    
    % Copy matrices such that parallel works better
    pmatT=pmat{t};
    pmatPrev=pmat{t-1};
    
    xmatprev=xmat(:,t-1);
    xmatcur=xmat(:,t);
    
    %%
    % If not using parfor outside of this function (ie. simulating one
    % company), you can enable parfor here to work through the employees
    % faster. For many companies, it is faster to do this outside this
    % function.   
    % parfor i = 1:n
    % for i = 1:n
    for i = 1:n
        
        
        % Starting values for local search is last period choice
        p_prev=pmatPrev(i,:)';
        
        % Since we calculate the SPNE, we technically don't need this
        x_prev=xmatprev;
        
        % Initialize convenience variables
        theta_i = theta(i);
        
        %% Optimization Setup
        % Set up the global search
        if globalsearch==1 || n<10
            gs = GlobalSearch('Display','off','StartPointsToRun','bounds-ineqs');
            %gs = GlobalSearch('Display','off','StartPointsToRun','all');
        else
            gs=0;
        end
        
        % Set up options for the local search
        % We find that not parallizing here is much faster for some reason
        options  =  optimset('Algorithm','sqp','Display', 'none','UseParallel',false);
        
        %% Constraints
        
        % Inequality contraints, p_i must be a capacity and sum to less
        % than one
        A=[eye(length(p_prev));-1.*eye(length(p_prev));ones(1,length(p_prev));-1.*ones(1,length(p_prev))];
        b=[ones(length(p_prev),1);zeros(length(p_prev),1);1;0];
        
        % If needed, equality constraints that p_i must be a probability
        Aeq=eye(n,n);
        Aeq=Aeq(i,:);
        % Note that we do not end up using this constraint
        Aeq=[Aeq;ones(length(p_prev),1)'];
        beq=[0;1];
        
        %% Optimize according to type
        UCon2=0;
        
        % In case we want global optimization only during certain
        % intervalls
        fsintervall=max(1,round((T/2)/t));
        %%%%%%%% Climber
        if identity(i)==1
            
            fncC=@(p_i) -utilitySPNE(p_i,pmatPrev,g,delta,theta_i, theta, PsiA, ConBen,i,convexp,Gmat);
            if globalsearch==1
                % Equality constraints, p_i is probability (unused)
                problemCons = createOptimProblem('fmincon','x0',p_prev,'objective', fncC,...
                    'lb',zeros(length(p_prev),1),'ub',ones(length(p_prev),1), ...
                    'Aineq',A,'bineq',b,'Aeq',Aeq,'beq',beq,'options',options);
                
                % Only inequality constraints, p_i is capacity
                problem = createOptimProblem('fmincon','x0',p_prev,'objective', fncC,...
                    'lb',zeros(length(p_prev),1),'ub',ones(length(p_prev),1), ...
                    'Aineq',A,'bineq',b,'options',options);
            end
            
            % Global or local optimization?
            if ((~mod(t,fsintervall)) || (t<=100)||(t>=1000)) && globalsearch==1
                [pstar,UCon2] = run(gs,problem);
            else
                %[pstar,UCon2]= fmincon(problem);
                [pstar,UCon2]=fmincon(fncC,p_prev,A,b,[],[],zeros(length(p_prev),1),ones(length(p_prev),1),[],options);
                
                
            end
            % Get utility value
            %[util,~]=utilitySPNE(pstar,pmatPrev,g,delta,theta_i, theta, PsiA, ConBen,i,convexp);
            %[util_old,~]=utilitySPNE(p_prev,pmatPrev,g,delta,theta_i, theta, PsiA, ConBen,i,convexp);
            
            
            %%%%%%% Watcher
        elseif identity(i)==0
            
            fncW=@(p_i) -utilitySPNE(p_i,pmatPrev,g,delta,theta_i, theta, PsiL, ConBen,i,convexp,Gmat);
            
            if globalsearch==1
                problem = createOptimProblem('fmincon','x0',p_prev,'objective', fncW,...
                    'lb',zeros(length(p_prev),1),'ub',ones(length(p_prev),1), ...
                    'Aineq',A,'bineq',b,'options',options);
            end
            
            % Global or local optimization?
            if ((~mod(t,fsintervall)) || (t<=100)||(t>=1000)) && globalsearch==1
                [pstar,UCon2] = run(gs,problem);
            else
                %[pstar,UCon2]= fmincon(problem);
                [pstar,UCon2]=fmincon(fncW,p_prev,A,b,[],[],zeros(length(p_prev),1),ones(length(p_prev),1),[],options);
                
            end
            % Get utility value
            %[util,~]=utilitySPNE(pstar,pmatPrev,g,delta,theta_i, theta, PsiL, ConBen,i,convexp);
            %[util_old,~]=utilitySPNE(p_prev,pmatPrev,g,delta,theta_i, theta, PsiL, ConBen,i,convexp);
            %%%%%%% Slacker
        else
            fncS=@(p_i) -utilitySPNE(p_i,pmatPrev,g,delta,theta_i, theta, PsiS, ConBen,i,convexp,Gmat);
            
            
            if  globalsearch==1
                problem = createOptimProblem('fmincon','x0',p_prev,'objective', fncS,...
                    'lb',zeros(length(p_prev),1),'ub',ones(length(p_prev),1), ...
                    'Aineq',A,'bineq',b,'options',options);
            end
            
            % Global or local optimization?
            if ((~mod(t,fsintervall)) || (t<=100)||(t>=1000)) && globalsearch==1
                [pstar,UCon2] = run(gs,problem);
            else
                %[pstar,UCon2]= fmincon(problem);
                [pstar,UCon2]=fmincon(fncS,p_prev,A,b,[],[],zeros(length(p_prev),1),ones(length(p_prev),1),[],options);
                
                
            end
            % Get utility value
            %[util,~]=utilitySPNE(pstar,pmatPrev,g,delta,theta_i, theta, PsiS, ConBen,i,convexp);
            %[util_old,~]=utilitySPNE(p_prev,pmatPrev,g,delta,theta_i, theta, PsiS, ConBen,i,convexp);
        end
        
        
        %% Constraint check
        % If Fmincon doesn't satisfy the constraints
        % Which sometimes happens for tiny values, this keeps the simulation from spinning out
        % of control
        
        if (min(pstar < 0)) || max(pstar >1)
            %disp(['Pstar out of range. Fixing. Occured with g: ',num2str(g),', identity: ',num2str(identity(i)),', theta: ',num2str(theta_i), ...
            %    ' in period ',num2str(t), ...
            %    ' with pstar min/max: ',num2str(min(pstar)),', ',num2str(max(pstar))]);
            pstar(pstar<0)=0;
            pstar(pstar>1)=1;
        end
        if sum(pstar>1)
            pstar = bsxfun(@times, pstar, 1./(sum(pstar)));
            pstar(isnan(pstar))=0;
            warning('sum Pstar larger than 1. Fixing.');
            disp(['sum Pstar larger than 1. Fixing. Occured with g: ',num2str(g),', identity: ',num2str(identity(i)),', theta: ',num2str(theta_i), ...
                ' in period ',num2str(t), ...
                ' with sum pstar: ',num2str(sum(pstar))]);
        end
        
        
        %% Option to disconnect
        % Disable if disconnection is not part of possible solution
        % but allowed in model, so we allow for disconnection if no
        % positive utility can be found
        
        if 0 > -UCon2
            pstar=zeros(length(pstar),1);
        end
        
        %% Final calculations
        % Save utility vector
        UCont_1=umat(i,t-1);
        if UCont_1 == -UCon2
            %disp(['Last period same for i=',num2str(i),' and ',num2str(-UCon2)]);
            pstar=pmatPrev(i,:)';
            uvec(i)=UCont_1;
        elseif UCont_1 > -UCon2
            %disp(['Last period better for i=',num2str(i),' diff is ',num2str(abs(-UCon2-UCont_1))]);
            pstar=pmatPrev(i,:)';
            uvec(i)=UCont_1;
        else
            uvec(i)=-UCon2;
        end
        
        
        
        % Save P_i and X_i in the matrices
        %xmatcur(i)=x(i);
        pmatT(i,:)=pstar';
        
    end
    
    
    % pmatT now has all period best-replies
    % calculate the  current x
    x=XFOCSPNE(pmatT,delta,theta,g);
    
    % Save values to matrices
    xmat(:,t)=x;
    pmat{t}=pmatT;
    umat(:,t)=uvec(:);
    
    % Check convergence. Fluctuations should be small enough, and / or the
    % simulation has run for many periods already
    prevxdif=xdif;
    prevpmat=pmat{t-1};
    udif = max(abs(umat(:,t)-umat(:,t-1)));
    xdif = max(abs(xmat(:,t)-xmat(:,t-1)));
    pdif = max(max(abs(pmatT-prevpmat)));
    finalt=t;
    
    %%% Uncomment to see convergence per iteration
    disp([' t ',num2str(t),' with current convergence in utils ',num2str(udif),', in x ',num2str(xdif),' in p ',num2str(pdif)])
    if (((pdif <= 1.0000e-4) && (udif <= 1.0000e-5) && (xdif<= 1.0000e-5)) && (t>2*minT))
        disp(['Finish on t ',num2str(t),' with current convergence in utils ',num2str(udif),', in x ',num2str(xdif),' in p ',num2str(pdif)])
        break;
    end
    if (((pdif <= 1.0000e-6) && (udif <= 1.0000e-6) && (xdif<= 1.0000e-10)) && (t>minT))
        disp(['Finish on t ',num2str(t),' with current convergence in utils ',num2str(udif),', in x ',num2str(xdif),' in p ',num2str(pdif)])
        break;
    end
    % If half-time is achieved and p values fluctuate only by tiny amounts,
    % end the simulation
    if  ((abs(pdif)<= 1.0000e-6*t) && (t>=T/2))
        warning('Convergence achieved due to half-time');
        break;
    end
    
end

%% Calculate Return Variables
T=finalt;
Pmat=pmat{T};
X=xmat(:,T);
diff=X-theta;
utils=umat(:,T);

%% Calculate graph properties
type=string(ones(size(identity)));
type(identity==1)="C";
type(identity==0)="W";
type(identity==-1)="S";

% Rounding
Pmat=round(Pmat,4);
% Lowest Cutoff Value if desired for graphing
% Pmat(Pmat < 0.001)=0;

% We need a symmetric matrix to calculate cohesion
SPMat=(Pmat+Pmat')./2;
SPMat(SPMat>0)=1;
GPgraph=digraph(Pmat);
Sgraph=graph(SPMat);
%GPgraph=Sgraph;

% Calculate articulation points
[edgebins,iC] = biconncomp(Sgraph);
critNodes=length(iC);

% Delete Watchers from Graph
Sgraph.Nodes.Type=identity;
Sgraph.Nodes.X=X;
Sgraph.Nodes.Num=[1: length(identity)]';
Sgraph.Nodes.theta=round(theta,3);

del=Sgraph.Nodes.Num(identity==0)';
Sgraph2=Sgraph.rmnode(del);
% Save critical nodes for slackers and watchers only
[edgebins,iC] = biconncomp(Sgraph2);
critNodesCS=length(iC);

% Centrality we currently do not use
CentPR=centrality(GPgraph,"pagerank");
CentPR=CentPR./(max(CentPR));
CentPRSym=centrality(Sgraph,"pagerank");
CentPRSym=CentPRSym./(max(CentPRSym));

CentBet=centrality(GPgraph,"betweenness");
CentBet=CentBet./(max(CentBet));
CentBetSym=centrality(Sgraph,"betweenness");
CentBetSym=CentBetSym./(max(CentBetSym));


%% Deliver return values
% We will return two objects,
% RetStruct: A struct that can be turned into a table easily, including
% numerical values
% RetCell: A cell including non-numerical values such as the adjacency
% matrix, type vector and so forth.

%% Create Return Table

% Nr of types
RetStruct.NrClimbers=length(identity(identity==1));
RetStruct.NrWatchers=length(identity(identity==0));
RetStruct.NrSlackers=length(identity(identity==-1));

% Theta Data
RetStruct.ThetaMean=mean(theta);
RetStruct.ThetaMin=min(theta);
RetStruct.ThetaMax=max(theta);
RetStruct.TSkew=skewness(theta);
RetStruct.ThetaVar=var(theta);

if RetStruct.ThetaVar==0
    error('ThetaVar=0');
end

% X Data
RetStruct.XMean=mean(X);
RetStruct.XMin=min(X);
RetStruct.XMax=max(X);
RetStruct.XSkew=skewness(X);
RetStruct.DiffMean=mean(diff);
RetStruct.DiffMin=min(diff);
RetStruct.DiffMax=max(diff);
RetStruct.XVar=var(X);

% Social Multiplier
RetStruct.SM=mean(X./theta);
RetStruct.LSM=mean(log(X./theta));


% Calculate first order effects
P=pmat{T};
DP=eye(n,n).*sum(P,2);

for k =1:n
    prow=P(k,:);
    XP1(k)=(1/(1+sum(prow)*g))*theta(k)+(g/(1+sum(prow)*g))*prow*theta;
end
XP1=XP1';
% First order DiffMean, SM and log SM
RetStruct.PDiffMean=mean(XP1-theta);
RetStruct.PSM=mean(XP1./theta);
RetStruct.PLSM=real(mean(log(XP1./theta)));


% Information for the maximum theta actor
% His quality
RetStruct.Tmax=max(theta);
% His social multiplier (first order)
P1SM=XP1./theta;
RetStruct.TmaxX1=P1SM(theta==RetStruct.Tmax).*RetStruct.Tmax;
% His X
RetStruct.TmaxX=X(theta==RetStruct.Tmax);
% His SM
RetStruct.TmaxSM=RetStruct.TmaxX/RetStruct.Tmax;
% His identity
RetStruct.TmaxId=identity(theta==RetStruct.Tmax);

% Parameters used in distribution

RetStruct.ThetaD_a=thetaD(1);
RetStruct.ThetaD_b=thetaD(2);
RetStruct.ThetaD_mean=thetaD(3);
RetStruct.ThetaD_sd=thetaD(4);


RetStruct.ProbC=gammaC;
RetStruct.ProbW=gammaW;
RetStruct.ProbS=gammaS;

RetStruct.T=T;
RetStruct.g=g;
RetStruct.m=m;

RetStruct.n=n;

RetStruct.critNodes=critNodes;
RetStruct.critNodesCS=critNodesCS;

% Enable this to return centralities
% RetStruct.CentPRSymMean=mean(CentPRSym);
% RetStruct.CentPRSymMin=min(CentPRSym);
% RetStruct.CentPRSymMax=max(CentPRSym);
%
% RetStruct.CentPRMean=mean(CentPR);
% RetStruct.CentPRMin=min(CentPR);
% RetStruct.CentPRMax=max(CentPR);
%
% RetStruct.CentBetSymMean=mean(CentBetSym);
% RetStruct.CentBetSymMin=min(CentBetSym);
% RetStruct.CentBetSymMax=max(CentBetSym);
%
% RetStruct.CentBetMean=mean(CentBet);
% RetStruct.CentBetMin=min(CentBet);
% RetStruct.CentBetMax=max(CentBet);


RetCell={theta,diff,X,Pmat,identity,utils};


returndata.RetCell=RetCell;
returndata.RetStruct=RetStruct;

%% GRAPH ROUTINES
%%
if graphit==1
    GPgraph=digraph(Gmat);
    GPgraph=digraph(Pmat);
    GPgraph.Nodes.Type=identity;
    GPgraph.Nodes.X=X;
    GPgraph.Nodes.Num=[1: length(identity)]';
    GPgraph.Nodes.theta=round(theta,3);
    GPgraph.Nodes.diff=round(diff,3);
    %GPgraph.Nodes.theta=100.*GPgraph.Nodes.theta./(max(GPgraph.Nodes.theta));
    %GPgraph.Nodes.diff=100.*GPgraph.Nodes.diff./(max(GPgraph.Nodes.theta));
    
    
    %GPgraph.Nodes.NodeColors=CentBet;
    strtheta=string(GPgraph.Nodes.theta);
    strType=string(identity);
    strleer="|";
    strdiff=string(GPgraph.Nodes.diff);
    strX=string(GPgraph.Nodes.X);
    type(:) = strcat(strType,strleer,strtheta(:),strleer,strleer,strX(:),strleer,strdiff(:));
    type(:) = strcat(strtheta(:),strleer,strdiff(:));
    asdf= cellstr(type');
    h=plot(GPgraph,'Layout','force','NodeLabel',asdf);
    %h=plot(GPgraph,'Layout','force');
    h.NodeCData=identity;   
    h.MarkerSize=12;
    layout(h,'force')
    %highlight(h,Pgraph,'EdgeColor','r','LineWidth',1.5)
    
    txt = {'Attention network:',['N=',num2str(n)],['g=',num2str(g)], ... 
        ['NrC=',num2str(RetStruct.NrClimbers)],['NrW=',num2str(RetStruct.NrWatchers)],['NrS=',num2str(RetStruct.NrSlackers)],['Skew=',num2str(skewness(theta))], ... 
        ['SM=',num2str(RetStruct.SM)],['AvgTheta=',num2str(RetStruct.ThetaMean)],['AvgX=',num2str(RetStruct.XMean)]};
    annotation('textbox',...
    [0.14 0.9 0 0],...
    'String',txt);

    
    %h.LineWidth=2.2;
    %layout(h,'force','Iterations',5000)
    
    %h.EdgeCData = GPgraph.Edges.Weight;
    
    %Highlight articulation points
    %highlight(h, iC)
    
end
end