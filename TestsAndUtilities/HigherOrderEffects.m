%%
% Higher Order Effects
% Authors: Ingo Marquart, Nghi Truong, Matthew Bothner, Richard Haynes
%
% This is a utility to recalculate measured without re-running
% the simulations
%%



runtime='20181209 FullSpaceComplete';
ResName='20181209-ConsFullSpaceCOMPLETE.csv';
filenameTable=['../DataSave/',runtime,'/'];
filenamePMats=['../DataSave/',runtime,'/Pmats/'];

allTable=readtable(strcat(filenameTable,ResName));




allTable=allTable(2:end,:);
alllength=height(allTable);
%allTable=allTable(:,2:end);
newTable=allTable;

% Define new columns
newTable.TypeWeightOutCentMin=zeros(alllength,1);
newTable.TypeWeightOutCentMax=zeros(alllength,1);
newTable.TypeWeightOutCentMean=zeros(alllength,1);
newTable.TypeWeightInCentMin=zeros(alllength,1);
newTable.TypeWeightInCentMax=zeros(alllength,1);
newTable.TypeWeightInCentMean=zeros(alllength,1);
newTable.InfluenceCentMax=zeros(alllength,1);
newTable.InfluenceCentMin=zeros(alllength,1);
newTable.InfluenceCentMean=zeros(alllength,1);

newTable.InDegreeSimMax=zeros(alllength,1);
newTable.InDegreeSimMin=zeros(alllength,1);
newTable.InDegreeSimMean=zeros(alllength,1);


newTable.RolemodelDistance=zeros(alllength,1);

parfor i = 1:alllength
    row=newTable(i,:);
    g=row.g;
    
    filenamePMat=strcat(filenamePMats,num2str(row.timestamp),'_','P_n',num2str(row.n),'_m',num2str(row.m),'_mn',num2str(row.mn),'_pn',num2str(row.pn),'.csv');
    PMat=readtable(filenamePMat);
    if ~isempty(PMat)
        XMat=PMat(:,1:8);
        
        P=PMat{:,9:row.n+8};
        G=PMat{:,row.n+8:end};
        centralities=TypeCentralities(P,XMat.Identity);
        row.TypeWeightOutCentMin=min(centralities.outcentrality);
        row.TypeWeightOutCentMax=max(centralities.outcentrality);
        row.TypeWeightOutCentMean=mean(centralities.outcentrality);
        row.TypeWeightInCentMin=min(centralities.incentrality);
        row.TypeWeightInCentMax=max(centralities.incentrality);
        row.TypeWeightInCentMean=mean(centralities.incentrality);
        row.InfluenceCentMax=max(centralities.influencecentrality);
        row.InfluenceCentMin=min(centralities.influencecentrality);
        row.InfluenceCentMean=mean(centralities.influencecentrality);
        
        [~,indegreesim] = DegreeSimilarity(P);
        row.InDegreeSimMin=min(min(indegreesim));
        row.InDegreeSimMax=max(max(indegreesim));
        row.InDegreeSimMean=mean(mean(indegreesim));
        
        row.RolemodelDistance=RolemodelDistance(P);
        
        
        
        newTable(i,:)=row;
    else
        disp(['PMat was empty in i= ',num2str(i)])
    end
end


writetable(newTable,strcat(filenameTable,strcat('Update_',ResName)));