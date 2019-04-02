
%% Simulation config
graphIt=0; % Graph the last firm
saveIt=1; % Save simulation results to a new folder
maxT=1;

%% Baseline parameters
paramsDefault.maxT=maxT; % Time periods for the company to run
paramsDefault.minEqmT=5; % Minimum periods to achieve equilibrium convergence
paramsDefault.maxEqmT=100; % Maximum period to achieve eqm convergence
paramsDefault.globalsearch=-1; % Unused
paramsDefault.thetaRepShockVar=0; % Standard deviation of shock to theta representations
paramsDefault.rationality=0; % First-stage rationality of agents
paramsDefault.pn=0; % P parameter for G network. Set to 0 for full G!
paramsDefault.mn=0; % M parameter for G network (Jackson&Rogers 2014 algorithm)


%% List of parameters to run
nList=[10:5:80];
mList=1:30;
eList=[0.10,  0.30,  0.50,  0.70, 0.90, 1.00, 5.00, 50];
consList=[-1,0,1];

%% Type settings
% Probabilities of climbers relative to slackers. 
% Simulation will check symmetrically for slackers
PCscale=0.25:0.25:0.5;
% Overall probability of watchers. Can also be a vector to check for several levels of watchers
Wscale=[2/9,2/3];

%% Theta settings
thetascale=[2,5];

%% Create a cell array of parameters to run
paramsCell=createParamCell(PCscale,Wscale,thetascale,mList,nList,eList,consList,paramsDefault);
% Nr of simulations
NrSims=length(paramsCell);

%% Initialize first row
% We need the outputs of a first row to set up tables
% So we run the first firm as a "test"
firm=initFirm(paramsCell{1});
firm=runTurn(firm);
startRow=returnRow(firm,2);


% A cell to keep resulting tables for each firm
% Sadly, matlab requires this to run a parfor loop
resultCell={1,NrSims};
% A cell to keep the last firm for graphing
graphFirm={1,NrSims};

%% Main loop over firms
disp(['Preparing for ',num2str(NrSims),' firms over ',num2str(maxT),' periods for a total of ',num2str(maxT*NrSims),' runs.'])
tic
parfor i = 1:NrSims
    params=paramsCell{i};
    
    % This is the temporary table to be filled for firm i
    tableToFill=repmat(startRow,maxT,1);
    
    % Initialize firm
    firm=initFirm(params);
    %% Main loop over T, where T=1 is initial setup
    for T=firm.T:maxT+1
        % Change time period of firm
        firm.T=T;
        % Run the current turn for the firm
        firm=runTurn(firm);
        % Return results as a row
        row=returnRow(firm,T);
        % Save the row into the table for the firm
        tableToFill(T-1,:)=row;
    end
    % Fill entirely table into the result cell
    resultCell{i}=tableToFill;
    
    % This saves the last firm entirely in a way compatible with parfor
    if graphIt==1 && i==NrSims
        graphFirm{i}=firm;
    end
    
    % Rough progress display
    if  mod(i,(floor(NrSims*0.1)))==0
        disp(['Iteration m=',num2str(params.m),', n=',num2str(params.n),', e=',num2str(params.e),', cons=',num2str(params.cons),', theta=',num2str(params.thetaD),', gamma=',num2str(params.gamma)])
    end
end
toc
%% Create resulting Table
% Due to matlab being matlab we had
% to save results in a cell
% Convert now to full table
mainTable=startRow;
for rows=resultCell
    mainTable=[mainTable;rows{:}];
end
% Delete first row and other objects
mainTable=mainTable(2:height(mainTable),:);
clear resultCell
% Graph the last firm if requested
if graphIt==1
    GraphNetwork(graphFirm{NrSims});
end

%% Saving mechanics
% Current runtime
runtime=posixtime(datetime);

% Define here the directories to save data
dirnameP=['../DataSave/',num2str(runtime),'/Pmats/'];
dirname=['../DataSave/',num2str(runtime),'/'];
if saveIt==1
    mkdir(dirnameP);

    TableName=strcat(dirname,'Results.csv');
    writetable(mainTable,TableName);
end