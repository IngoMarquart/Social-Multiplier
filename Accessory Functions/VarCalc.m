%%
%   Distribution of CWS
%%
runtime='1525454272.045';
filenameTable=['DataSave/',runtime,'/'];
filenamePMats=['DataSave/',runtime,'/Pmats/'];

allTable=readtable(strcat(filenameTable,'SPNEResAll.csv'));


%Subset
allTable=allTable(allTable.ThetaMean>allTable.XMean,:);
%allTable=allTable((allTable.g<0.3),:);


alllength=height(allTable);

cdistW=zeros(alllength,4);
cdistS=zeros(alllength,4);
cdistC=zeros(alllength,4);
cdistAll=zeros(alllength,4);


parfor i = 1:alllength
row=allTable(i,:);


filenamePMat=strcat(filenamePMats,'P_n',num2str(row.n),'_m',num2str(row.m),'_',num2str(row.timestamp),'.csv');
PMat=readtable(filenamePMat);

P=PMat{:,end-(row.n-1):end};
% Recode to C=3, W=2, S=1 to separate from non connection 0
identity=PMat.X+2;
% Calculate the occurances of connections according to their target
nrC=P(identity==3,:)*identity;
nrS=P(identity==1,:)*identity;
nrW=P(identity==2,:)*identity;
nrAll=P*identity;

% Calculate the neighborhood distributions of C/W/S 
cdistW(i,:)=[sum(nrW==0)/length(nrW),sum(nrW==1)/length(nrW),sum(nrW==2)/length(nrW),sum(nrW==3)/length(nrW)];
cdistS(i,:)=[sum(nrS==0)/length(nrS),sum(nrS==1)/length(nrS),sum(nrS==2)/length(nrS),sum(nrS==3)/length(nrS)];
cdistC(i,:)=[sum(nrC==0)/length(nrC),sum(nrC==1)/length(nrC),sum(nrC==2)/length(nrC),sum(nrC==3)/length(nrC)];
cdistAll(i,:)=[sum(nrAll==0)/length(nrAll),sum(nrAll==1)/length(nrAll),sum(nrAll==2)/length(nrAll),sum(nrAll==3)/length(nrAll)];

end

cdistW(any(isnan(cdistW), 2), :) = [];
cdistS(any(isnan(cdistS), 2), :) = [];
cdistC(any(isnan(cdistC), 2), :) = [];
cdistAll(any(isnan(cdistAll), 2), :) = [];

mean(cdistS)
mean(cdistC)
mean(cdistW)
mean(cdistAll)