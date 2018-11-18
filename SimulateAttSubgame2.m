%%
%   SimulateAttSubgame
%   For given parameter ranges, this function calculates the
%   SPNE as a steady state of a dynamic simulation
%%
function returndata=SimulateAttSubgame2(n,T,gamma, thetaD,ConBen, gemA, gemL, delta, g,m, minT,graphit,globalsearch,convexp,Gmat,cons)


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

%% Prepare consolidation
% We will sort a number of thetas by their size.
nrSort=round(abs(cons)*n);
% Randomly select thetas
[~,idx]=datasample(TIVec(:,1),nrSort,'Replace',false);
% Sort this part of the sample
if cons < 0 % Consolidation < 0: We sort ascending, since slackers at the bottom
    vec=sort(TIVec(idx,1));
    sidx=sort(idx);
    TIVec(sidx,1)=vec;
else % Consolidation > 0: We sort descending, since climbers at the top
    vec=sort(TIVec(idx,1),'descend');
    sidx=sort(idx);
    TIVec(sidx,1)=vec;
end

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

%% Generate choice sets
ChoiceCell={};
nrChoices=zeros(1,n);
[ChoiceCell, nrChoices]=GetChoiceSet(Gmat,n);


%% Generate Constraint matrices
ConA={};
Conb={};
[ConA,Conb] = CalcConstraints(n, nrChoices);

% For small values of n we want global optimization
if n<15 && globalsearch >= 0
    globalsearch=1;
end

% Loop over time periods until convergence
for t = 2:(T)
    
    % Initialize cell content for period t
    pmat{t}=zeros(n,n);
    
    % Copy matrices such that parallel works better
    pmatT=pmat{t};
    pmat_prev=pmat{t-1};
    
    %
    %%
    % If not using parfor outside of this function (ie. simulating one
    % company), you can enable parfor here to work through the employees
    % faster. For many companies, it is faster to do this outside this
    % function.
    % parfor i = 1:n
    % for i = 1:n
    for i = 1:n
        
        % Initialize other stuff
        UCon=0;
        A=ConA{i};
        b=Conb{i};
        theta_i = theta(i);
        % Initialize previous p_i vector
        pi_prev=pmat_prev(i,:)';
        % Initialize previous restricted p_i
        prest_prev=pi_prev(ChoiceCell{i});
        prest_star=prest_prev;
        
        
        %% Discrete optimization or continous optimization
        if globalsearch==-1
            %% Discrete optimization
            
            % Set up objective function
            if identity(i)==1 % Climber
                [UCon,pi_star]=DiscreteChoice(pmat_prev,g,delta,theta(i), theta, PsiA,i,ChoiceCell{i},nrChoices(i));
            elseif identity(i)==0 % Watcher
                [UCon,pi_star]=DiscreteChoice(pmat_prev,g,delta,theta(i), theta, PsiL,i,ChoiceCell{i},nrChoices(i));
            else % Slacker
                [UCon,pi_star]=DiscreteChoice(pmat_prev,g,delta,theta(i), theta, PsiS,i,ChoiceCell{i},nrChoices(i));
            end
            % Convention is utility is negative for fmincon
            UCon=-UCon;           
            
        else
            %%  Continous Optimization Setup

            % Set up the global search
            if globalsearch==1
                gs = GlobalSearch('Display','off','StartPointsToRun','bounds-ineqs');
                %gs = GlobalSearch('Display','off','StartPointsToRun','all');
            else
                gs=0;
            end
            
            % Set up options for the local search
            % We find that not parallizing here is much faster for some reason
            options  =  optimset('Algorithm','sqp','Display', 'none','UseParallel',false);
            
            
            
            %% Set up objective function
            if identity(i)==1 % Climber
                fnc=@(prest) -utilitySPNE2(prest,pmat_prev,g,delta,theta_i, theta, PsiA, ConBen,i,convexp,ChoiceCell{i});
            elseif identity(i)==0 % Watcher
                fnc=@(prest) -utilitySPNE2(prest,pmat_prev,g,delta,theta_i, theta, PsiL, ConBen,i,convexp,ChoiceCell{i});
            else % Slacker
                fnc=@(prest) -utilitySPNE2(prest,pmat_prev,g,delta,theta_i, theta, PsiS, ConBen,i,convexp,ChoiceCell{i});
            end
            
            
            %% Set up optimization problems
            if globalsearch==1
                % Only inequality constraints, p_i is capacity
                problem = createOptimProblem('fmincon','x0',prest_prev,'objective', fnc,...
                    'lb',zeros(length(prest_prev),1),'ub',ones(length(prest_prev),1), ...
                    'Aineq',A,'bineq',b,'options',options);
            end
            
            %% Run optimization
            
            if nrChoices(i)>0
                % Global or local optimization?
                if (((t<=100)||(t>=1000)) && globalsearch==1)
                    [prest_star,UCon] = run(gs,problem);
                else
                    %[pstar,UCon2]= fmincon(problem);
                    [prest_star,UCon]=fmincon(fnc,prest_prev,A,b,[],[],zeros(length(prest_prev),1),ones(length(prest_prev),1),[],options);
                end
            end
            %% Recover pi_star
            % Add Code
            pi_star = RecoverPi(prest_star,ChoiceCell{i},n);
        end
        
        
        
        %% Post optimization checks
        %% Constraint check
        % If Fmincon doesn't satisfy the constraints
        % Which sometimes happens for tiny values, this keeps the simulation from spinning out
        % of control
        
        if (min(pi_star < 0)) || max(pi_star >1)
            disp(['Pstar out of range. Fixing. Occured with g: ',num2str(g),', identity: ',num2str(identity(i)),', theta: ',num2str(theta_i), ...
                ' in period ',num2str(t), ...
                ' with pstar min/max: ',num2str(min(pstar)),', ',num2str(max(pstar))]);
            pi_star(pi_star<0)=0;
            pi_star(pi_star>1)=1;
        end
        if sum(pi_star>1)
            pi_star = bsxfun(@times, pi_star, 1./(sum(pi_star)));
            pi_star(isnan(pi_star))=0;
            disp(['sum Pstar larger than 1. Fixing. Occured with g: ',num2str(g),', identity: ',num2str(identity(i)),', theta: ',num2str(theta_i), ...
                ' in period ',num2str(t), ...
                ' with sum pstar: ',num2str(sum(pi_star))]);
        end
        
        
        %% Option to disconnect
        % Disable if disconnection is not part of possible solution
        % but allowed in model, so we allow for disconnection if no
        % positive utility can be found
        
        if 0 > -UCon
            pi_star=zeros(length(pi_star),1);
        end
        
        %% Final calculations
        % Save utility vector
        UCont_prev=umat(i,t-1);
        if UCont_prev == -UCon
            %disp(['Last period same for i=',num2str(i),' and ',num2str(-UCon)]);
            pi_star=pmat_prev(i,:)';
            uvec(i)=UCont_prev;
        elseif UCont_prev > -UCon
            %disp(['Last period better for i=',num2str(i),' diff is ',num2str(abs(-UCon-UCont_prev))]);
            pi_star=pmat_prev(i,:)';
            uvec(i)=UCont_prev;
        else
            uvec(i)=-UCon;
        end
        
        
        
        % Save P_i and X_i in the matrices
        pmatT(i,:)=pi_star';
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
    prevpmat=pmat{t-1};
    udif = max(abs(umat(:,t)-umat(:,t-1)));
    xdif = max(abs(xmat(:,t)-xmat(:,t-1)));
    pdif = max(max(abs(pmatT-prevpmat)));
    finalt=t;
    
    %%% Uncomment to see convergence per iteration
    %disp([' t ',num2str(t),' with current convergence in utils ',num2str(udif),', in x ',num2str(xdif),' in p ',num2str(pdif)])
    if (((pdif <= 1.0000e-4) && (udif <= 1.0000e-5) && (xdif<= 1.0000e-5)) && (t>2*minT))
        %disp(['Finish on t ',num2str(t),' with current convergence in utils ',num2str(udif),', in x ',num2str(xdif),' in p ',num2str(pdif)])
        break;
    end
    if (((pdif <= 1.0000e-6) && (udif <= 1.0000e-6) && (xdif<= 1.0000e-10)) && (t>minT))
        %disp(['Finish on t ',num2str(t),' with current convergence in utils ',num2str(udif),', in x ',num2str(xdif),' in p ',num2str(pdif)])
        break;
    end
    % If half-time is achieved and p values fluctuate only by tiny amounts,
    % end the simulation
    if  ((abs(pdif)<= 1.0000e-6*t) && (t>=T/2))
        warning('Convergence achieved due to half-time');
        break;
    end
end
%% Deliver return values
% We will return two objects,
% RetStruct: A struct that can be turned into a table easily, including
% numerical values
% RetCell: A cell including non-numerical values such as the adjacency
% matrix, type vector and so forth.
%% Calculate Return Variables
T=finalt;
P=pmat{T};
X=xmat(:,T);
diff=X-theta;
utils=umat(:,T);


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

%% TODO CHECK XSKEW VS TSKEW in REULT

RetStruct.XSkew=skewness(X);
RetStruct.DiffMean=mean(diff);
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

%% Calculate graph properties
%% G
GPgraph=graph(Gmat);

% Network density
PC=(n*(n-1))./2;
EC=height(GPgraph.Edges);
RetStruct.GDensity=EC/PC;

% Average path length
% We calculate each path length separately for each
% connected component in the symmetric graph of A
% an then use a weighted average
bincell = conncomp(GPgraph, 'OutputForm', 'cell');
nrsubgraphs = length(bincell);
nrNodes=zeros(1,nrsubgraphs);
sumDist=zeros(1,nrsubgraphs);
nrEdges=zeros(1,nrsubgraphs);
weightfactor=zeros(1,nrsubgraphs);
for ii = 1:nrsubgraphs
    subg=subgraph(GPgraph, bincell{ii});
    DM=distances(subg);
    nrNodes(ii)=length(DM);
    sumDist(ii)=sum(sum(DM));
    weightfactor(ii)=nrNodes(ii)/n;
    nrEdges(ii)=(nrNodes(ii).*(nrNodes(ii)-1));
    % Zero weight for single nodes (Path Length not defined)
    if(nrNodes(ii) == 1)
        weightfactor(ii)=0;
        nrEdges(ii)=1;
    end
    % Renormalize weights to sum to one
    weightfactor=weightfactor./sum(weightfactor);
   
end
% Average path length is given by the sum of distances divided by the 
RetStruct.GAvgPathLength=((1./nrEdges).*weightfactor)*sumDist';


% Largest eigenvector
RetStruct.GMaxEV=max(eig(Gmat));


% Clustering & average Degree
deg = sum(Gmat, 2); %Determine node degrees 
cn = diag(Gmat*triu(Gmat)*Gmat); %Number of triangles for each node
%The local clustering coefficient of each node 
c = zeros(size(deg)); 
c(deg > 1) = 2 * cn(deg > 1) ./ (deg(deg > 1).*(deg(deg > 1) - 1));
RetStruct.GAvgClustering=mean(c);
RetStruct.GAvgDegree=mean(deg);


% Bonacich family centralities
B=BonacichCentrality(0,0,1,1,0,Gmat);
% Power centrality
BNeg=BonacichCentrality(0,0,1,1,1,Gmat);

% Create aggregate measures for Centralities
RetStruct.GmeanBonacich=mean(B);
RetStruct.GmaxBonacich=max(B);
RetStruct.GminBonacich=min(B);
RetStruct.GvarBonacich=var(B);

RetStruct.GmeanNegBonacich=mean(BNeg);
RetStruct.GmaxNegBonacich=max(BNeg);
RetStruct.GminNegBonacich=min(BNeg);
RetStruct.GvarNegBonacich=var(BNeg);


%% A
SPMat=(P+P')./2;
SPMat(SPMat>0)=1;
Agraph=digraph(P);
ASymgraph=graph(SPMat);

% Nr of vertices
ECSym=height(ASymgraph.Edges);
EC=height(Agraph.Edges);

% Network density
PC=(n*(n-1))./2;
RetStruct.ADensity=EC/PC;

% Average path length
% We calculate each path length separately for each
% connected component in the symmetric graph of A
% an then use a weighted average
bincell = conncomp(ASymgraph, 'OutputForm', 'cell');
nrsubgraphs = length(bincell);
nrNodes=zeros(1,nrsubgraphs);
sumDist=zeros(1,nrsubgraphs);
nrEdges=zeros(1,nrsubgraphs);
weightfactor=zeros(1,nrsubgraphs);
for ii = 1:nrsubgraphs
    subg=subgraph(GPgraph, bincell{ii});
    DM=distances(subg);
    nrNodes(ii)=length(DM);
    sumDist(ii)=sum(sum(DM));
    weightfactor(ii)=nrNodes(ii)/n;
    nrEdges(ii)=(nrNodes(ii).*(nrNodes(ii)-1));
    % Zero weight for single nodes (Path Length not defined)
    if(nrNodes(ii) == 1)
        weightfactor(ii)=0;
        nrEdges(ii)=1;
    end
    % Renormalize weights to sum to one
    weightfactor=weightfactor./sum(weightfactor);
   
end
% Average path length is given by the sum of distances divided by the 
RetStruct.AAvgPathLength=((1./nrEdges).*weightfactor)*sumDist';


% Largest eigenvector
RetStruct.AMaxEV=max(eig(P));
RetStruct.ASymMaxEV=max(eig(SPMat));

% Clustering & average Degree
deg = sum(SPMat, 2); %Determine node degrees
deg2 = sum(P, 2);
cn = diag(SPMat*triu(SPMat)*SPMat); %Number of triangles for each node
%The local clustering coefficient of each node 
c = zeros(size(deg)); 
c(deg > 1) = 2 * cn(deg > 1) ./ (deg(deg > 1).*(deg(deg > 1) - 1));
RetStruct.ASymAvgClustering=mean(c);
RetStruct.ASymAvgDegree=mean(deg);
RetStruct.AAvgDegree=mean(deg2);


% Bonacich family centralities
ABonacich=BonacichCentrality(0,0,1,1,0,P);
% In-Degree Bonacich family centralities
AInBonacich=BonacichCentrality(0,0,1,1,0,P');
% Power centrality
AnegBonacich=BonacichCentrality(0,0,1,1,1,P);


% Bonacich family centralities: P Symmetric
ASymBonacich=BonacichCentrality(0,0,1,1,0,SPMat);
% Power centrality
ASymnegBonacich=BonacichCentrality(0,0,1,1,1,SPMat);

% Create aggregate measures for Centralities
RetStruct.AmeanBonacich=mean(ABonacich);
RetStruct.AmaxBonacich=max(ABonacich);
RetStruct.AminBonacich=min(ABonacich);
RetStruct.AvarBonacich=var(ABonacich);

RetStruct.AInmeanBonacich=mean(AInBonacich);
RetStruct.AInmaxBonacich=max(AInBonacich);
RetStruct.AInminBonacich=min(AInBonacich);
RetStruct.AInvarBonacich=var(AInBonacich);

RetStruct.AmeanNegBonacich=mean(AnegBonacich);
RetStruct.AmaxNegBonacich=max(AnegBonacich);
RetStruct.AminNegBonacich=min(AnegBonacich);
RetStruct.AvarNegBonacich=var(AnegBonacich);

RetStruct.ASymmeanBonacich=mean(ASymBonacich);
RetStruct.ASymmaxBonacich=max(ASymBonacich);
RetStruct.ASymminBonacich=min(ASymBonacich);
RetStruct.ASymvarBonacich=var(ASymBonacich);

RetStruct.ASymmeanNegBonacich=mean(ASymnegBonacich);
RetStruct.ASymmaxNegBonacich=max(ASymnegBonacich);
RetStruct.ASymminNegBonacich=min(ASymnegBonacich);
RetStruct.ASymvarNegBonacich=var(ASymnegBonacich);


% Calculate articulation points
[~,iC] = biconncomp(ASymgraph);
critNodes=length(iC);

% Delete Watchers from Graph
ASymgraph.Nodes.Type=identity;
ASymgraph.Nodes.Num=[1: length(identity)]';
ASymgraph.Nodes.theta=round(theta,3);

del=ASymgraph.Nodes.Num(identity==0)';
Sgraph2=ASymgraph.rmnode(del);
% Save critical nodes for slackers and watchers only
[~,iC] = biconncomp(Sgraph2);
critNodesCS=length(iC);

RetStruct.critNodes=critNodes;
RetStruct.critNodesCS=critNodesCS;

RetCell={identity,theta,diff,X,utils,B,BNeg,ABonacich,AnegBonacich,P};


returndata.RetCell=RetCell;
returndata.RetStruct=RetStruct;

%% GRAPH ROUTINES
%%
if graphit==1
    GraphNetwork(Gmat,SPMat,identity,theta,X)
end
end