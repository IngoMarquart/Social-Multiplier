%%
% % MAINSPNE
% Authors: Ingo Marquart, Nghi Truong, Matthew Bothner, Richard Haynes
%
% This simulation generates a set of firms based on the paper
% "Shaping the social architecture of a start-up: What raises the social multiplier?"
% In particular, this calculates the subgame-perfect Nash equilibrium
%
% Note that if saveit=1, the simulation will save each firm advacency matrix as well
% as a table of aggregates in a folder ../Datasave/Timestamp
%% 


%% CONFIGURATION

% Restrict the maximum time periods for the simulation
T=500;
% Define a minimum of time periods each simulation runs
minT=4;
% Convexity parameter of Psi wrt. e, set to 1 (unused).
% Use this if you wish to examine very low values of g, as for
% very low values of g the differences become so small that
% Matlab otherwise has trouble finding the optimum (g<0.01).
% Set this value to 3 to examine low g.
convexp=1;
% Plot a network graph. Note that this should only be enabled for
% simulating a single firm!
graphit=0;
% Use of search algorithm for SPNE
% We usually use "-1" and use the other, continous numeric options to confirm correctness of our approach
% 1 - Use Matlab globalsearch - necessary for n<10, time intensive for high n
% 0 - Use Matlab minimizer - Faster than global search, almost always converges to same value
% -1 - Discrete optimizer - Custom, based on mathematical results in paper, fastest algorithm
globalsearch=-1;
% Save results. Creates a folder with timestamp and saves aggregate data
saveit=1;
% If saveadjmat = 1 in addition to saveit=1, then the simulation will create
% a subfolder "PMats" and save in it all adjacency matrices for the firm
saveadjmat=1;
% Force matlab to store table out of memory. Slower if enough RAM, but necessary
% if ResultTable exceeds RAM (leads to crash). 
longtable=0;


%% CURRENTLY UNUSED PARAMETERS
% Please do not change
% Exogenous Connection Benefit to force no isolates (unused).
ConBen=0;
% Different marginal Psi utility for climbers and slackers, set to 1 (unused).
gemA=1;
% Different marginal Psi utility for watchers, set to 1 (unused).
gemL=1;
% Exogenous scaling parameter for social effects, set to 1 (unused).
delta=1;


%% DEFINE PARAMETER RANGES OF FIRMS
% The following sets cell arrays for the parameters in use.
% The simulation will produce an observation for each combination.

%%
% gammaVec defines a cell array of sets of probabilities P(C),P(W),P(S)
% You may define a scale to check both for Climbers&Slackers, as well as
% Watchers.

% Probabilities of climbers relative to slackers. Simulation will check symmetrically for slackers
PCscale=0:0.1:0.5;
% Overall probability of watchers.
Wscale=[1/3,2/3];
% Inititalize
gammaVec={};
iC=1;
for watchP = Wscale
    for scale = PCscale
        normC=scale.*(1-watchP);
        normS=(1-scale).*(1-watchP);
        gammaVec{iC}=[normC,watchP,normS];
        if normC == normS
         iC=iC+1;    
        else
        gammaVec{iC+1}=[normS,watchP,normC]
        iC=iC+2;            
        end

    end
end

%% Previously generated vectors
% to replicate data.
% gammaVec={[1/3, 1/3, 1/3], ...
%     [1/2, 1/2, 0/3], ...
%     [0/3, 1/2, 1/2], ...
%     [2/9,1/3,4/9], ...
%     [2/6,1/3,2/6], ...
%     [4/9,1/3,2/9], ...
%     [1/9,2/3,2/9], ...
%     [1/6,2/3,1/6], ...
%     [2/9,2/3,1/9]};

% "Archetype examples"
% gammaVec={[1/3, 1/3, 1/3], ...
%     [2/9,2/9,5/9], ...
%     [5/9,2/9,2/9], ...
%     [2/9,5/9,2/9]};


%%
% thetaVec defines a cell array of parameters for the Beta distribution
% Each vector includes [a,b,c,d]
% a,b - Beta shape parameters
% b,c - Scaling of variance and mean - unused in the current version
% Note that the current version scales fixes variance to 1.
thetascale=2.5:0.5:7.5;
thetaVec={};
iz=1;
for scale = thetascale
    thetaVec{iz}=[2,scale,1,1];
    thetaVec{iz+1}=[scale,2,1,1];
    iz=iz+2;
end
thetaVec{iz}=[2,2,1,1];


%% "Archetype examples"
% thetascale=[5];
% thetastart=2;
% thetaVec={};

%%
% gVec is a cell array of all g (or e) values to run
% Note, we parallelize over this so its best to have it equal to the number
% of cores
% In this case, we use 8 values.
gVec={ 0.10,  0.30,  0.50,  0.70, 0.90, 1.00, 5.00, 50};

%%
% nVec is a cell array of all firm sizes to run
nVec={10,15,20,25,30,35,40,45,50,55,60};

%%
% mVec includes the random seeds to run for each configuration.
% We usually want many firms for each combination of paramters. 
% Each value is used to initialize the random number generator.
% Set manually as: mVec={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
iz=1;
for scale = 1:30
    mVec{iz}=scale;
    iz=iz+1;
end

%% Settings for G-Matrix
% If you set these values to non-zero, firms will
% start with a random network G at the beginning
% using Jackon&Rogers algorithm.

% 2xAverage degree
% Example:
% mnVec={1,3,5,10};
mnVec={0};
% Probability of connecting to peer once found
% Example:
% pnVec={0.1,0.3,0.5,1};
pnVec={0};

%% Settings for consolidation
consVec={-1,-0.8,-0.5,-0.3,0,0.3,0.5,0.8,1};
% Set to zero for no consolidation
% consVec={0};

%% SINGLE SIMULATION
% You can uncomment the following lines and run a single simulation

% thetaVec={[2,15 ,1,1]};
% nVec={50};
% mVec={1};
% gVec={100};
% gammaVec={[1/3, 1/3, 1/3]};
% mnVec={1};
% pnVec={0.3};
% consVec={0};



%% END OF CONFIGURATION


% Current runtime
runtime=posixtime(datetime);

% Define here the directories to save data
dirnameP=['../DataSave/',num2str(runtime),'/Pmats/'];
dirname=['../DataSave/',num2str(runtime),'/'];
if saveit==1
    mkdir(dirnameP);
end

%% Preallocation
NrSims = length(gammaVec)*length(thetaVec)*length(gVec)*length(nVec)*length(mVec)*length(mnVec)*length(pnVec)*length(consVec);
returncell={};

% Used cells for consistency. Parallelized dimensions we want as vector.
gVec=cell2mat(gVec);

% PreAllocation Test
ResultTable=table();

% We run one small simulation to test the system and create the table
% such that we can preallocate.
% Note that we only use this with longtable out of RAM.
if saveit==1
    GTest=ones(10,10)-eye(10,10);
    returndata=SimulateAttSubgame2(10,T,[1/3, 1/3, 1/3], [2,2,1,1],ConBen, gemA, gemL, delta, 1,1, 2,0,globalsearch,convexp,GTest,1);
    addrow=struct2table(returndata.RetStruct);
    addrow.timestamp=1;
    addrow.Archetype=strcat("Test");
    addrow.convexp=convexp;
    addrow.pn=1;
    addrow.mn=1;
    addrow.cons=0;
    ResultTable=addrow;
    
    if longtable==1
        TableName=strcat(dirname,'SPNEResAll.csv');
        writetable(ResultTable,TableName);
        storename=strcat(dirname,'*.csv');
        ds=datastore(storename);
        ResultTable=tall(ds);
    else
        % Preallocation without longtable does not seem to improve performance
    end
end

% timer
tic


%% START OF MAIN LOOP

count=1;

% Firms sizes
for cn = nVec
    n=cn{1};
    
    % Random seeds
    for cm = mVec
        m=cm{1};
        
        % For each new random seed we display how far along we are
        disp(['Simulating Iteration ',num2str(count),' of ',num2str(NrSims)])
        
        
        % Type distributions
        for cgamma = gammaVec
            gamma=cgamma{1};
            
            % Theta distributions
            for cthetaD=thetaVec
                thetaD=cthetaD{1};
                for ccons=consVec
                    cons=ccons{1};
                for cpn=pnVec
                    pn=cpn{1};
                    for cmn=mnVec
                        mn=cmn{1};
                        % Now we parallelize over g
                        % This requires each simulation to write independently into
                        % a prepared vector
                        tablecell={1,length(gVec)};
                        tablenamecell={1,length(gVec)};
                        
                        % As identifier of each simulation we just use count
                        timestamps=repmat(count,length(gVec),1)+(1:length(gVec))'-1;
                        
                        
                        % Pre Return Table to keep tall array out of parfor loop
                        % because the tall array in parfor really slows things down
                        PreReturnTable = table();
                        
                        % If you run only one simulation (for graphing etc), you
                        % can switch from parfor to for and parallelize in the
                        % simulation. Use
                        % for i = 1:length(gVec)
                        % and use your parfor elsewhere.
                        % If you want many firms, it is far more efficient to
                        % parallelize per firm! Use
                        % parfor i = 1:length(gVec)
                        parfor i = 1:length(gVec)
                            g=gVec(i);
                            
                            % Save current timestamp to associate matrix with table
                            % output
                            timestamp=timestamps(i);
                            
                            
                            %% Run simulation
                            
                            % Create random matrix
                            
                            if pn==0
                                Gmat=ones(n,n)-eye(n,n);
                            else
                                Gmat=JacksonRogersNW(n,mn,pn, m);
                            end
                            returndata=SimulateAttSubgame(n,T,gamma, thetaD,ConBen, gemA, gemL, delta, g,m, minT,graphit,globalsearch,convexp,Gmat,cons);
                            
                            % Row to add to table
                            addrow=struct2table(returndata.RetStruct);
                            addrow.timestamp=timestamp;
                            addrow.Archetype=strcat("CS:",num2str(gamma(1)./(gamma(3)+gamma(1))),"S:",num2str(thetaD(1)/thetaD(2)));
                            addrow.convexp=convexp;
                            addrow.pn=pn;
                            addrow.mn=mn;
                            addrow.cons=cons;
                            % Concat table and cell
                            PreReturnTable=[PreReturnTable;addrow];
                            returncell=[returncell;returndata.RetCell];
                            
                            
                            %% Save adjacency matrix of simulation
                            % Create table names in the form of P1...Pn
                            nName=arrayfun(@num2str,[1:length(returndata.RetCell{4})],'un',0);
                            nName=strcat('P',nName);
                            GName=arrayfun(@num2str,[1:length(returndata.RetCell{4})],'un',0);
                            GName=strcat('G',GName);
                            % Create other names of rows for each actor
                            Names=[{'Identity'},{'Theta'},{'X'},{'Utils'},{'GBonacich'},{'GNegBonacich'},{'ABonacich'},{'ANegBonacich'},nName,GName];
                            
                            % Create matrix (see "Names" for content of RetCell)
                            % RetCell={identity,theta,diff,X,utils,B,BNeg,ABonacich,AnegBonacich,P};
                            a=array2table([returndata.RetCell{1},returndata.RetCell{2},returndata.RetCell{4},returndata.RetCell{5},returndata.RetCell{6},returndata.RetCell{7},returndata.RetCell{8},returndata.RetCell{9},returndata.RetCell{10},Gmat]);
                            a.Properties.VariableNames=Names;
                            
                            % Save it to a cell
                            TableName=strcat(dirnameP,num2str(timestamp),'_','P_n',num2str(n),'_m',num2str(m),'_mn',num2str(mn),'_pn',num2str(pn),'.csv');
                            tablecell{1,i}=a;
                            tablenamecell{1,i}=TableName;
                            
                            % Increment count
                            count=count+1;
                        end
                        
                        % We now add all parallelized results for this run back
                        % into the original table, which may be an out-memory array
                        ResultTable=[ResultTable;PreReturnTable];
                        
                        % If saving is enabled, we save the adjacency matrices
                        % Saving is done outside the for loop because it seems
                        % faster
                        if saveit==1 && saveadjmat == 1
                            for i = 1:length(gVec)
                                writetable(tablecell{1,i},tablenamecell{1,i});
                            end
                        end
                    end
                end
                end
            end
        end
    end
    % For every n we also save the table separately, in case the computer
    % crashes midway through
    if saveit==1
        if longtable==1
            TableToWrite=gather(ResultTable(ResultTable.n==n,:));
            TableName=strcat(dirname,'SPNERes_n',num2str(n),'.csv');
            writetable(TableToWrite,TableName);
        else
            TableToWrite=ResultTable(ResultTable.n==n,:);
            TableName=strcat(dirname,'SPNERes_n',num2str(n),'.csv');
            writetable(TableToWrite,TableName);
        end
    end
end
toc

% Save final table.
if saveit==1
    if longtable==1
        TableToWrite=gather(ResultTable);
        TableName=strcat(dirname,'SPNEResAll.csv');
        writetable(TableToWrite,TableName);
    else
        TableToWrite=ResultTable;
        TableName=strcat(dirname,'SPNEResAll.csv');
        writetable(TableToWrite,TableName);
    end
end
toc
