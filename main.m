%% Initial comments
% Use config.m to set parameters, then run main.m
%% FIGURE 1: Skew and embedding
% Please see config.m for further details
%% Plan of simulation
% Each firm we run saves its state as a "firm" struct, which is the
% argument of most of our model functions.
% Results of the simulation are saved in a table object for simple
% exporting.
% 1. Parameter Cell
%   With a list of variables passed from the config file, we create a Cell
%   object that includes parameters for each firm we plan to run, including
%   initial values, random networks and so on.
% 2. Execution plan
%   Since we do not want to hard-code which results are returned, but still
%   need to pre-allocate (see below) the result objects, we run a test firm
%   once to return a single "example" row, to see how many variables are
%   returned.
%   If not graphing, we will run several firms in parallel. Each firm runs
%   over maxT time periods.
%   Results from each time period t are saved as a row in tableToFill.
%   All these T rows for each firm are then saved in resultCell. This is
%   necessary to run firms in parallel.
%   Therefore, resultCell{i} holds all maxT
%   results for firm i.
%   After the simulation block, the cell is reduced (syncronously) to a
%   table called mainTable. This is the final return object.
% 3. Reduction into mainTable
%   Reduction proceeds in two steps. First, all cell items from resultCell
%   are vertically concatenated (across firms) 
%   into vertcat. vertcat is then filled into
%   mainTable at the appropriate location (see below why it is not just
%   appended).
% 4. Block execution
%   Firms need to be ran for maxT time periods which implies a large
%   resultCell. Since Cells take more memory than
%   a table, we may need to reduce resultCell into mainTable periodically
%   instead of once at the end. Set maxCellLength as the maximum size of
%   resultCell to hold in memory. When maxCellLength is larger than this
%   value, parallel executation is held and resultCell is reduced to
%   mainTable.
%   Thus, the main loop runs along
%   -> blocks
%   --> firms
%   ---> time periods
%   But mainTable saves a row for each time period. Thus, some index
%   juggling is necessary.
% 5. Working with parallel objects - code complexity.
%   Matlab's parallel programming model is at the time of
%   writing this program, limited.
%   You will note that we try to pre-allocate all tables and cells.
%   Generally, efficient and parallel MatLab code requires us to always
%   access tables or cells via a an index or range that is fixed prior to
%   executation. Essentially, this can only be the counter variable of a
%   loop. In contrast, if we could append results, the code would be much
%   simpler. However, when it comes to memory efficiency, appending is not
%   possible in MatLab2018.
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONFIGURATION %%%%%%%%%%%%%%%%%%%%
% run outsources config.m in path
config;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF CONFIG %%%%%%%%%%%%%%%%%%%%

nrMBlocks=ceil(size(mListTotal,2)/mIterations);
for bigSet = 1:nrMBlocks

% Restart parallel pool to avoid slowdown bug
% This is a Matlab Bug that occurs when looping over a parfor
poolobj = gcp('nocreate');
delete(poolobj);
parpool('local');
    
mList=[mListTotal((bigSet-1)*(mIterations)+1):mListTotal(bigSet*mIterations)];

disp(['Preparing for iteration set ', num2str(bigSet)])
%% Create a cell array of parameters to run
paramsCell = createParamCell(PCscale, Wscale, thetascale, mList, nList, eList, consList, paramsDefault, symmetric, conUtil, conParam,learningRates);
% Nr of simulations
NrFirms = length(paramsCell);

%% Initialize first row
% We need the outputs of a first row to set up tables
% So we run the first firm as a "test"
firm = initFirm(paramsCell{1});
firm = runTurn(firm);
% The startRow object holds our "test" row.
startRow = returnRow(firm, 2);

% A cell to keep resulting tables for each firm
% Sadly, matlab requires this to run a parfor loop
resultCell = {NrFirms};
% If graphing, we want to hold an entire firm object in memory.
graphFirm=cell(1,1);



% Nr of simulation runs amounts to T*NrFirms
NrSims=maxT * NrFirms;
%% Main loop over firms
disp(['Preparing for ', num2str(NrFirms), ' firms over ', num2str(maxT), ' periods for a total of ', num2str(NrSims), ' runs.'])

tic

% Pre-allocate main table: One row for each simulation run.
% Pre-allocation is essential for parallel processing.
mainTable=repmat(startRow,NrSims,1);



%% Block allocation
% The number of blocks required for a given maxCellLength
% Each Cell includes T rows
nrBlocks=ceil(NrFirms/maxCellLength);

% Overall counter of starting position in terms of firms
% The current "block" runs the firms between
% counterStartFirm and counterEndFirm.
% not including counterStartFirm - the first firm run is
% counterStartFirm+1.
% This is because of how we reduce to a table (see below)
counterStartFirm=0;

% We now run each block of firms separately.
% Our counting variable scheme is essentially
for block=1:nrBlocks
    % Create result cell for block
    resultCell = cell(1,maxCellLength);
    % Maximal end point if this block had maxCellLength firms
    counterEndFirm=counterStartFirm+maxCellLength;
    % In which case the length of this block is just
    % maxCellLength=endi-counter
    cellLength=counterEndFirm-counterStartFirm;
    % However, the last block may be partial. There may be less than
    % maxCellLength firms to run.
    if counterEndFirm> NrSims/maxT
        counterEndFirm=ceil(NrSims/maxT);
        cellLength=counterEndFirm-counterStartFirm;
    end
    
    % In addition to resultCell, we also need to worry about our parameter
    % cell we use to create the firm's initial state.
    % Since we access it from inside the parallel loop, we need to access
    % it with the loop counter only. Hence, here we slice off the
    % corresponding cellLength firms.
    blockParamsCell=paramsCell(counterStartFirm+1:counterEndFirm);
    % for or parfor
    % SET TO PARFOR FOR LARGE SAMPLE
    % Currently GRAPHING REQUIRES FOR
    parfor i = 1:cellLength
        
        % Fill in parameters
        params = blockParamsCell{i};
        % This is the temporary table to be filled for firm i
        % for all maxT time periods.
        tableToFill = repmat(startRow, maxT, 1);
        
        % Initialize firm
        firm = initFirm(params);
        %% Main loop over T, where T=1 is initial setup
        % note that we therefore run to firm.T:maxT + 1 time periods
        % because our firm starts at T=2.
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
        resultCell{i} = tableToFill;
        
        % This saves the last firm entirely in a way compatible with parfor
        if graphIt == 1 && counterEndFirm>=NrSims/maxT
            graphFirm = firm;
        end
        
        
    end
    
    % Rough progress display
    disp(['Finished Iteration Block start=',num2str(counterStartFirm),', end=',num2str(counterEndFirm),' in block ',num2str(block),'/',num2str(nrBlocks)])
    
    
    %% Create resulting Table
    % After finishing each block, we reduce resultCell into a table
    % First, concat all results from current block
    addTable=vertcat(resultCell{:});
    % This allows us to clear the resultCell
    clear resultCell
    % Derive Index for maintable (transform from firm-units to time-units)
    % counterStartFirm is the index of the last firm from the prior block.
    % Our starting position is hence just one index after its last time
    % period.
    sStart=(counterStartFirm)*maxT+1;
    % Consequently, our end point
    sEnd=(counterEndFirm)*maxT;
    % Add to main Table
    mainTable(sStart:sEnd,:)=addTable;
    % Now clear also temporary addTable
    clear addTable
    % Set our last ran firm as new starting point.
    counterStartFirm=counterEndFirm;
    
end
toc

clear blockParamsCell paramsCell


%% Graphing functions
% Graph the last firm if requested
if graphIt == 1
    %graphFirm=graphFirm{1};
    if toGraph == "NW"
        sumMat=graphFirm.aMat{1};
        for amat=graphFirm.aMat
            sumMat=sumMat+amat{1};
        end
        graphFirm.aMat{graphFirm.T}=sumMat/graphFirm.maxT;
        graphFirm.aMat{graphFirm.T}=graphFirm.aMat{graphFirm.T-1}
        GraphNetwork(graphFirm);
    elseif toGraph == "SM"
        % Get firm
        firm = graphFirm;
        % Set first embedding to zero
        firm.eMat(1) = 0;
        % Normalize embedding
        firm.eMat = firm.eMat./(1+firm.eMat);
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
            firm.eMat = firm.eMat;
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

%% Saving
% Current runtime
runtime = posixtime(datetime);

% Define here the directories to save data
if identifier=="time-of-run"
    dirname = ['../DataSave/', num2str(runtime), '-',num2str(mList(1)),'-',num2str(mList(end)),'/'];
else
    %dirname = ['../DataSave/', identifier, '-',num2str(mList(1)),'-',num2str(mList(end)),'/'];
    dirname = ['../DataSave/', identifier,'-',num2str(bigSet),'/'];

end


if saveIt == 1
    mkdir(dirname);
    
    TableName = strcat(dirname, 'Results.csv');
    writetable(mainTable, TableName);
    clear mainTable
end

end