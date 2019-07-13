%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONFIGURATION %%%%%%%%%%%%%%%%%%%%
% run outsources config.m in path
config;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF CONFIG %%%%%%%%%%%%%%%%%%%%
%% Create a cell array of parameters to run
paramsCell = createParamCell(PCscale, Wscale, thetascale, mList, nList, eList, consList, paramsDefault, symmetric, conUtil, conParam);
% Nr of simulations
NrSims = length(paramsCell);

%% Initialize first row
% We need the outputs of a first row to set up tables
% So we run the first firm as a "test"
firm = initFirm(paramsCell{1});
firm = runTurn(firm);
startRow = returnRow(firm, 2);

% A cell to keep resulting tables for each firm
% Sadly, matlab requires this to run a parfor loop
if avgOverT == 1
    resultCell = {NrSims / maxT};
else
    resultCell = {NrSims};
end


%% Main loop over firms
disp(['Preparing for ', num2str(NrSims), ' firms over ', num2str(maxT), ' periods for a total of ', num2str(maxT * NrSims), ' runs.'])

%% Take into account T
NrSims=maxT * NrSims;

tic

% Pre-allocate main table
mainTable=repmat(startRow,NrSims,1);

% Overall counter
counter=0;

% The number of blocks required for a given maxCellLength
% Each Cell includes T rows
nrBlocks=ceil(NrSims/(maxCellLength*maxT));

for block=1:nrBlocks
    % Create result cell for block
    resultCell = cell(1,maxCellLength);
    % Maximal end point
    endi=counter+maxCellLength;
    cellLength=endi-counter;
    % Last block may be partial as we can only have
    % NrSims/maxT firms
    if endi> NrSims/maxT
        endi=ceil(NrSims/maxT);
        cellLength=endi-counter;
    end
    
    
    % Restrict param cell slice since it is broadcast variable
    % in parfor loop
    blockParamsCell=paramsCell(counter+1:endi);

% for or parfor
parfor i = 1:cellLength
    
    % Fill in parameters
    params = blockParamsCell{i};
    % This is the temporary table to be filled for firm i
    tableToFill = repmat(startRow, maxT, 1);
    
    % Initialize firm
    firm = initFirm(params);
    %% Main loop over T, where T=1 is initial setup
    for T = firm.T:maxT + 1
        % Change time period of firm
        firm.T = T;
        % Run the current turn for the firm
        firm = runTurn(firm);
        % Return results as a row
        row = returnRow(firm, T);
        % Save the row into the table for the firm
        tableToFill(T - 1, :) = row;
    end
    
    % Fill entirely table into the result cell
    if avgOverT == 1
        row = tableToFill(T - 1, :);
        meanData = mean(tableToFill{:, 3:end});
        row{:, 3:end} = meanData;
        row.varXT = var(tableToFill.avgX);
        row.varThetaRangeT = var(tableToFill.thetaRange);
        row.varDiffMeanT = var(tableToFill.diffM);
        
        resultCell{i} = row;
    else
        resultCell{i} = tableToFill;
    end
    
    % This saves the last firm entirely in a way compatible with parfor
    if graphIt == 1 && endi>=NrSims/maxT
        graphFirm = firm;
    end
    
    % Rough progress display
    %if  mod(i,(floor(NrSims*0.1)))==0
    %    disp(['Iteration m=',num2str(params.m),', n=',num2str(params.n),', e=',num2str(params.e),', cons=',num2str(params.cons),', theta=',num2str(params.thetaD),', gamma=',num2str(params.gamma)])
    %end
end
%% Create resulting Table
% Due to matlab being matlab we had
% to save results in a cell
% Convert now to full table

    addTable=vertcat(resultCell{:});
    % Derive Index for Table
    sStart=(counter)*maxT+1;
    sEnd=(endi)*maxT;
    mainTable(sStart:sEnd,:)=addTable;
    clear resultCell
    clear addTable
    counter=endi;
    
end
toc



%% Graphing functions
% Graph the last firm if requested
if graphIt == 1
    
    if toGraph == "NW"
        GraphNetwork(graphFirm);
    elseif toGraph == "SM"
        % Get firm
        firm = graphFirm;
        % Set first embedding to zero
        firm.eMat(1) = 0;
        % Normalize embedding
        firm.eMat = firm.eMat ./ (1 + firm.eMat);
        aggDiff=mean(firm.xMat)-mean(firm.thetaMat(:,1));
        tsMat = table(firm.avgTheta', firm.avgX');
        tsTheta=timeseries(firm.avgTheta(2:end)','name','AvgTheta');
        tsX=timeseries(firm.avgX(2:end)','name','AvgX');
        tseBar=timeseries(firm.eMat(2:end)','name','Ebar');
        %tsO = timeseries(tsMat, 'name', 'SM over time');
        plot(tsTheta,'-xb','Displayname','Average Theta')
        hold on
        plot(tsX,'-.xm','Displayname','Average X')
        plot(tseBar,'Displayname','E Bar')
        legend('show','Location','NorthEast')
        txt = {['N=', num2str(firm.n)], ['ebar=', num2str(firm.eMat(firm.T))], ...
            ['NrC=', num2str(firm.NrC)], ['NrW=', num2str(firm.NrW)], ['NrS=', num2str(firm.NrS)], ['Skew=', num2str(firm.skew(firm.T))], ...
            ['X(t)-Theta(1)=', num2str(aggDiff(firm.T))]};
        annotation('textbox', ...
            [0.14 0.9 0 0], ...
            'String', txt);
        hold off
    elseif toGraph == "SM-stacked" % Display all firms over time
        % Fill matrices for X and theta
        for firm=graphFirm
           firm=firm{1}; 
           Theta=[firm.avgTheta(2:end)'];
           X=[firm.avgX(2:end)'];
                   firm.eMat = firm.eMat ./ (1 + firm.eMat);
        tsTheta=timeseries(Theta,'name','AvgTheta');
        tsX=timeseries(X,'name','AvgX');
        tsSM=timeseries(X-Theta(1),'name','SM');
        plot(tsSM,'Displayname',['x(t)-t(1), e ',num2str(firm.eMat(firm.T))])
        hold on
                %plot(tsTheta,'Displayname',['AvgT, e ',num2str(firm.eMat(firm.T))])

        %plot(tsX,'Displayname',['AvgX, e ',num2str(firm.eMat(firm.T))])
        end
      
  
        legend('show','Location','NorthEast')
        txt = {['N=', num2str(firm.n)], ['ebar=', num2str(firm.eMat(firm.T))], ...
            ['NrC=', num2str(firm.NrC)], ['NrW=', num2str(firm.NrW)], ['NrS=', num2str(firm.NrS)], ['Skew=', num2str(firm.skew(firm.T))]};
        annotation('textbox', ...
            [0.14 0.9 0 0], ...
            'String', txt);
        hold off
    else % Plot average network
        % Get firm
        firm = graphFirm{NrSims};
        % Set first embedding to zero
        firm.eMat(1) = 0;
        % Create Average Network
        avgP = maxT - 1;
        aMat = firm.aMat{firm.T};
        
        for i = 1:(avgP - 1)
            aMat = aMat + firm.aMat{firm.T - i};
        end
        
        rowSum = sum(aMat, 2);
        rowSum = min(1, rowSum.^(-1));
        aMat = aMat .* rowSum;
        firm.aMat{firm.T} = aMat;
        
        GraphNetwork(firm);
    end
    
end

%% Saving mechanics
% Current runtime
runtime = posixtime(datetime);

% Define here the directories to save data
dirnameP = ['../DataSave/', num2str(runtime), '/Pmats/'];
dirname = ['../DataSave/', num2str(runtime), '/'];

if saveIt == 1
    mkdir(dirnameP);
    
    TableName = strcat(dirname, 'Results.csv');
    writetable(mainTable, TableName);
end
