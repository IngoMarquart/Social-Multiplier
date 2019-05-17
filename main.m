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

% A cell to keep the last firm for graphing
graphFirm = {NrSims};

%% Main loop over firms
disp(['Preparing for ', num2str(NrSims), ' firms over ', num2str(maxT), ' periods for a total of ', num2str(maxT * NrSims), ' runs.'])
tic

for i = 1:NrSims
    params = paramsCell{i};

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
    if graphIt == 1
        graphFirm{i} = firm;
    end

    % Rough progress display
    %if  mod(i,(floor(NrSims*0.1)))==0
    %    disp(['Iteration m=',num2str(params.m),', n=',num2str(params.n),', e=',num2str(params.e),', cons=',num2str(params.cons),', theta=',num2str(params.thetaD),', gamma=',num2str(params.gamma)])
    %end
end

toc
%% Create resulting Table
% Due to matlab being matlab we had
% to save results in a cell
% Convert now to full table
firstT = resultCell{1};
% Pre-Allocate
mainTable = repmat(firstT(1, :), NrSims, 1);

parfor i = 1:length(resultCell)
    mainTable(i, :) = resultCell{i}
end

% Delete first row and other objects
mainTable = mainTable(2:height(mainTable), :);
clear resultCell

%% Graphing functions
% Graph the last firm if requested
if graphIt == 1

    if toGraph == "NW"
        GraphNetwork(graphFirm{NrSims});
    elseif toGraph == "SM"
        % Get firm
        firm = graphFirm{NrSims};
        % Set first embedding to zero
        firm.eMat(1) = 0;
        % Normalize embedding
        firm.eMat = firm.eMat ./ (1 + firm.eMat);
        tsMat = [firm.eMat', firm.diffM'];
        tsO = timeseries(tsMat, 'name', 'SM vs. embedding');
        plot(tsO);
        txt = {'SM vs. e:', ['N=', num2str(firm.n)], ['e=', num2str(firm.eMat(firm.T))], ...
            ['NrC=', num2str(firm.NrC)], ['NrW=', num2str(firm.NrW)], ['NrS=', num2str(firm.NrS)], ['Skew=', num2str(firm.skew(firm.T))], ...
            ['AvgDiff=', num2str(firm.diffM(firm.T))]};
        annotation('textbox', ...
            [0.14 0.9 0 0], ...
            'String', txt);
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
