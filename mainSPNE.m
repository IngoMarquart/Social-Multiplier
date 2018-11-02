%%
%   Matlab Simulation CWS Model - Simplified, SPNE
%   Ingo Marquart, Nghi Truong
%
%   This is the main file to set configuration and loop over
%   SimulateAttSubgame.m, which simulates on firm
%%



%% CONFIGURATION

% Restrict the maximum time periods for the simulation
T=100;
% Define a minimum of time periods each simulation runs
minT=4;
% Exogenous Connection Benefit to force no isolates. Unused.
ConBen=0;
% Different marginal Psi utility for climbers and slackers, set to 1 (unused).
gemA=1;
% Different marginal Psi utility for watchers, set to 1 (unused).
gemL=1;
% Exogenous scaling parameter for social effects, set to 1 (unused).
delta=1;
% Convexity parameter of Psi wrt. g, set to 1 (unused).
% Use this if you wish to examine very low values of g, as for
% very low values of g the differences become so small that
% Matlab otherwise has trouble finding the optimum (g<0.01).
% Set this value to 3 to examine low g.
convexp=1;
% Plot a network graph. Note that this should only be enabled for
% simulating a single firm!
graphit=1;
% Use globalsearch algorithm. Time intensive!
globalsearch=-1;
% Save results. Creates a folder with timestamp and saves aggregate data as
% well as ALL matrices
saveit=0;
% Force matlab to store table out of memory. Slower if enough RAM.
longtable=0;

%%
% Template for setting variables:
% Setting variables

%%
% gammaVec defines a cell array of sets of probabilities P(C),P(W),P(S)
gammaVec={[1/3, 1/3, 1/3], ...
    [1/2, 1/2, 0/3], ...
    [0/3, 1/2, 1/2], ...
    [2/9,1/3,4/9], ...
    [2/6,1/3,2/6], ...
    [4/9,1/3,2/9], ...
    [1/9,2/3,2/9], ...
    [1/6,2/3,1/6], ...
    [2/9,2/3,1/9]};

% gammaVec defines a cell array of sets of probabilities P(C),P(W),P(S)
gammaVec={[1/3, 1/3, 1/3], ...
    [2/9,2/9,5/9], ...
    [5/9,2/9,2/9], ...
    [2/9,5/9,2/9]};


%%
% thetaVec defines a cell array of parameters for the Beta distribution
% Each vector includes [a,b,c,d]
% a,b - Beta shape parameters
% b,c - Scaling of variance and mean - unused in the current version
% Note that the current version scales fixes variance to 1.
thetascale=2.5:0.5:8;
thetaVec={};
iz=1;
for scale = thetascale
    thetaVec{iz}=[2,scale,1,1];
    thetaVec{iz+1}=[scale,2,1,1];
    iz=iz+2;
end
thetaVec{iz}=[2,2,1,1];


%%
% thetaVec defines a cell array of parameters for the Beta distribution
% Each vector includes [a,b,c,d]
% a,b - Beta shape parameters
% b,c - Scaling of variance and mean - unused in the current version
% Note that the current version scales fixes variance to 1.
thetascale=[5];
thetastart=2;
thetaVec={};
iz=1;
for scale = thetascale
    thetaVec{iz}=[2,scale,1,1];
    thetaVec{iz+1}=[scale,2,1,1];
    iz=iz+2;
end
thetaVec{iz}=[2,2,1,1];

%%
% gVec is a cell array of all g values to run
% Note, we parallelize over this so its best to have it equal to the number
% of cores
gVec={ 0.10,  0.30,  0.50,  0.70, 0.90, 1.00, 5.00, 50};

%%
% nVec is a cell array of all firm sizes to run
nVec={75};

%%
% mVec includes the random seeds to run for each configuration.
% Set manually as: mVec={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
iz=1;
for scale = 1:30
    mVec{iz}=scale;
    iz=iz+1;
end

%% Settings for G-Matrix
% 2xAverage degree
%mnVec={1,3,5,10};
mnVec={0};
% Probability of connecting to peer once found
%pnVec={0.1,0.3,0.5,1};
pnVec={0};

%% Settings for consolidation
consVec={-1,-0.8,-0.5,-0.3,0,0.3,0.5,0.8,1};

%% Use this to run only one simulation
thetaVec={[2,15 ,1,1]};
nVec={100};
mVec={1};
gVec={0.1};
gammaVec={[1/3, 1/3, 1/3]};
mnVec={3};
pnVec={0};
consVec={0};
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
        % Preallocation here does not seem to improve performance
    end
end

% timer
tic


%% Running through loops

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
                        %
                        % Note, graphing doesn't work with parfor!
                        for i = 1:length(gVec)
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
                            returndata=SimulateAttSubgame2(n,T,gamma, thetaD,ConBen, gemA, gemL, delta, g,m, minT,graphit,globalsearch,convexp,Gmat,cons);
                            
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
                        if saveit==1
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
