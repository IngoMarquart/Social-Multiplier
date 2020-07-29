%% Initial comments
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONFIGURATION %%%%%%%%%%%%%%%%%%%%
% run outsources config.m in path
config;


toGraph="SM-firm";

maxT=100;
m=5;
e=0.95;
alpha=0.9;
theta=[1,2,3,4,10];
mu=[1,1,-1,1,1];


theta=theta(:);
mu=mu(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF CONFIG %%%%%%%%%%%%%%%%%%%%
%% Single firm setup to get parameters
n=size(theta,1);
nList = [n]; % Number of workers
mList = [m]; % Random Seeds / Iterations
eList = [e]; % Embedding levels
consList = [0]; % Consolidation levels
PCscale=[0.75]; Wscale=0; % Ratio of Improvers/Enhancers and P(Assessors)
thetascale=[2]; % Parameters of Beta Distribution
paramsDefault.maxT=maxT;
paramsDefault.shockTypes=1;
% The learning Rate delta
learningRates=[1-alpha];
graphIt=1;
%% Create a cell array of parameters to run
paramsCell = createParamCell(PCscale, Wscale, thetascale, mList, nList, eList, consList, paramsDefault, symmetric, conUtil, conParam,ceoActStartT,learningRates);
params=paramsCell{1};
%% Initialize first row
% We need the outputs of a first row to set up tables
% So we run the first firm as a "test"
firm = initFirm(params);
firm = runTurn(firm);
% The startRow object holds our "test" row.
startRow = returnRow(firm, 2);
resultTable = repmat(startRow, maxT, 1);

tic
% Initialize firm
firm = initFirmSingle(theta,mu,params);
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
    resultTable(T - 1, :) = row;
end

% This saves the last firm entirely in a way compatible with parfor
graphFirm = firm;
toc



%% Graphing functions
% Graph the last firm if requested
if graphIt == 1
    %graphFirm=graphFirm{1};
    if toGraph == "NW"
        GraphNetwork(graphFirm);
    elseif toGraph == "SM"
        % Get firm
        firm = graphFirm;
        % Set first embedding to zero
        firm.eMat(1) = 0;
        % Normalize embedding
        %firm.eMat = firm.eMat./(1+firm.eMat);
        aggDiff=mean(firm.xMat-firm.thetaMat(:,1));
        
        finalDiff=mean(firm.xMat(:,firm.T)-firm.thetaMat(:,T));
        
        tsMat = table(firm.avgTheta', firm.avgX');
        tsTheta=timeseries(firm.avgTheta(2:end)','name','AvgTheta');
        tsX=timeseries(firm.avgX(2:end)','name','AvgX');
        tseBar=timeseries(firm.eMat(2:end)','name','Ebar');
        %tsO = timeseries(tsMat, 'name', 'SM over time');
        plot(tsTheta,'-xb','Displayname','Average Theta')
        hold on
        plot(tsX,'-.xm','Displayname','Average X')
        %plot(tseBar,'Displayname','E Bar')
        legend('show','Location','NorthEast')
        txt = {['N=', num2str(firm.n)], ['ebar=', num2str(firm.eMat(firm.T))], ...
            ['NrC=', num2str(firm.NrC)], ['NrW=', num2str(firm.NrW)], ['NrS=', num2str(firm.NrS)], ['Skew=', num2str(firm.skew(firm.T))], ...
            ['avg(T)[X(t)-Theta(1)]=', num2str(mean(aggDiff))],['X(T)-Theta(1)=', num2str(aggDiff(firm.T))],['X(T)-Theta(T)=', num2str(finalDiff)]};
        annotation('textbox', ...
            [0.14 0.9 0 0], ...
            'String', txt);
        hold off
    elseif toGraph == "SM-firm"
        % Get firm
        firm = graphFirm;
        % Set first embedding to zero
        firm.eMat(1) = 0;
        % Normalize embedding
        %firm.eMat = firm.eMat./(1+firm.eMat);
        aggDiff=mean(firm.xMat-firm.thetaMat(:,1));
        
        finalDiff=mean(firm.xMat(:,firm.T)-firm.thetaMat(:,T));
        
        tsMat = table(firm.xMat);
        tsTheta=timeseries(firm.avgTheta(2:end)','name','AvgTheta');
        tsX=timeseries(firm.avgX(2:end)','name','AvgX');
        tseBar=timeseries(firm.eMat(2:end)','name','Ebar');
        tsO = timeseries(tsMat, 'name', 'X');
        xMat=firm.xMat';
        [a,nrN]=size(xMat)
        for k=1:nrN     
            type=firm.muMat(k,1);
            plot(xMat(:,k), 'color', rand(1,3), 'LineWidth', 2, 'DisplayName', sprintf('Worker: %i', type));
            hold on;          
            legend('-DynamicLegend');          
            legend('show');          
            drawnow;          
        end
        plot(tsTheta,'-xb','Displayname','Average Theta')
        plot(tsX,'-.xm','Displayname','Average X')
        legend('show','Location','NorthEast')
        txt = {['N=', num2str(firm.n)], ['ebar=', num2str(firm.eMat(firm.T))], ...
            ['NrC=', num2str(firm.NrC)], ['NrW=', num2str(firm.NrW)], ['NrS=', num2str(firm.NrS)], ['Skew=', num2str(firm.skew(firm.T))], ...
            ['avg(T)[X(t)-Theta(1)]=', num2str(mean(aggDiff))],['X(T)-Theta(1)=', num2str(aggDiff(firm.T))],['X(T)-Theta(T)=', num2str(finalDiff)]};
        annotation('textbox', ...
            [0.14 0.9 0 0], ...
            'String', txt);
        hold off
    else % Plot average network
        % Get firm
        firm = graphFirm;
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
