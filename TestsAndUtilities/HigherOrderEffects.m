%%
% Higher Order Effects
% Authors: Ingo Marquart, Nghi Truong, Matthew Bothner, Richard Haynes
%
% This is a utility to recalculate measured without re-running
% the simulations
%%


runtime='20181209 FullSpaceComplete';
ResName='Update_20181209-ConsFullSpaceCOMPLETE.csv';

% runtime='20181119 Con CentInDegree';
% ResName='Update_20181119-Cons-InDegFix.csv';

filenameTable=['../DataSave/',runtime,'/'];
filenamePMats=['../DataSave/',runtime,'/Pmats/'];

allTable=readtable(strcat(filenameTable,ResName));




allTable=allTable(2:end,:);
alllength=height(allTable);
%allTable=allTable(:,2:end);
newTable=allTable;

% Define new columns
% newTable.TypeWeightOutCentMin=zeros(alllength,1);
% newTable.TypeWeightOutCentMax=zeros(alllength,1);
% newTable.TypeWeightOutCentMean=zeros(alllength,1);
% newTable.TypeWeightInCentMin=zeros(alllength,1);
% newTable.TypeWeightInCentMax=zeros(alllength,1);
% newTable.TypeWeightInCentMean=zeros(alllength,1);
% newTable.InfluenceCentMax=zeros(alllength,1);
% newTable.InfluenceCentMin=zeros(alllength,1);
% newTable.InfluenceCentMean=zeros(alllength,1);
% newTable.InDegreeSimMax=zeros(alllength,1);
% newTable.InDegreeSimMin=zeros(alllength,1);
% newTable.InDegreeSimMean=zeros(alllength,1);
% newTable.RolemodelDistance=zeros(alllength,1);
newTable.AvgRgDiffC=zeros(alllength,1);
newTable.AvgRgDiffW=zeros(alllength,1);
newTable.AvgRgDiffS=zeros(alllength,1);
newTable.AvgRgDiffOverall=zeros(alllength,1);
newTable.samplecons = zeros(alllength,1);
newTable.sampleconsNeg = zeros(alllength,1);
newTable.CrossLinkPcC=zeros(alllength,1);
newTable.CrossLinkPcW=zeros(alllength,1);
newTable.CrossLinkPcS=zeros(alllength,1);
newTable.CrossLinkOverall=zeros(alllength,1);

parfor i = 1:alllength
    row=newTable(i,:);
    g=row.g;
    
    filenamePMat=strcat(filenamePMats,num2str(row.timestamp),'_','P_n',num2str(row.n),'_m',num2str(row.m),'_mn',num2str(row.mn),'_pn',num2str(row.pn),'.csv');
    PMat=readtable(filenamePMat);
    if ~isempty(PMat)
        XMat=PMat(:,1:8);
        
        P=PMat{:,9:row.n+8};
        G=PMat{:,row.n+8:end};
        
        % Average connection ranges
        ranges=RangeMeasure(XMat.Theta,XMat.Identity,P);
        row.AvgRgDiffC=ranges.rdiffC;
        row.AvgRgDiffW=ranges.rdiffW;
        row.AvgRgDiffS=ranges.rdiffS;
        row.AvgRgDiffOverall=ranges.rdiffOverall;
        
        
        % Actual consolidation
        samplecons = ConsolidationMeasure(XMat.Theta,XMat.Identity);
        row.samplecons=samplecons.cons;
        row.sampleconsNeg=samplecons.consNeg;
        % Cross linkage percentages
        crosslinkage = CrossLinkMeasure(XMat.Theta,XMat.Identity,P);
        row.CrossLinkPcC=crosslinkage.crosslinkageC;
        row.CrossLinkPcW=crosslinkage.crosslinkageW;
        row.CrossLinkPcS=crosslinkage.crosslinkageS;
        row.CrossLinkOverall=crosslinkage.overall;
        
        % Type centrality measures
%         centralities=TypeCentralities(P,XMat.Identity);
%         row.TypeWeightOutCentMin=min(centralities.outcentrality);
%         row.TypeWeightOutCentMax=max(centralities.outcentrality);
%         row.TypeWeightOutCentMean=mean(centralities.outcentrality);
%         row.TypeWeightInCentMin=min(centralities.incentrality);
%         row.TypeWeightInCentMax=max(centralities.incentrality);
%         row.TypeWeightInCentMean=mean(centralities.incentrality);
%         row.InfluenceCentMax=max(centralities.influencecentrality);
%         row.InfluenceCentMin=min(centralities.influencecentrality);
%         row.InfluenceCentMean=mean(centralities.influencecentrality);
%         
%         [~,indegreesim] = DegreeSimilarity(P);
%         row.InDegreeSimMin=min(min(indegreesim));
%         row.InDegreeSimMax=max(max(indegreesim));
%         row.InDegreeSimMean=mean(mean(indegreesim));
%         
%         row.RolemodelDistance=RolemodelDistance(P);
        
        
        
        newTable(i,:)=row;
    else
        disp(['PMat was empty in i= ',num2str(i)])
    end
end


writetable(newTable,strcat(filenameTable,strcat('Update_',ResName)));