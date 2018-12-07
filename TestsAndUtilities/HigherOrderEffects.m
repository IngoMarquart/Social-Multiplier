%%
% Higher Order Effects
% Authors: Ingo Marquart, Nghi Truong, Matthew Bothner, Richard Haynes
%
% This is a utility to recalculate measured without re-running
% the simulations
%% 



runtime='TestFitty';
filenameTable=['../DataSave/',runtime,'/'];
filenamePMats=['../DataSave/',runtime,'/Pmats/'];

allTable=readtable(strcat(filenameTable,'SPNEResAll.csv'));




allTable=allTable(2:end,:)
alllength=height(allTable);
%allTable=allTable(:,2:end);
newTable=allTable;
newTable.critNodesCS=zeros(alllength,1);


for i = 1:alllength
    row=newTable(i,:);
    g=row.g;
    
    filenamePMat=strcat(filenamePMats,num2str(row.timestamp),'_','P_n',num2str(row.n),'_m',num2str(row.m),'_mn',num2str(row.mn),'_pn',num2str(row.pn),'.csv');
    PMat=readtable(filenamePMat);
    if ~isempty(PMat)
          XMat=PMat(:,1:8);
          P=PMat{:,9:end};
          
          
          
%         XMat.DiffMean=XMat.X-XMat.Theta;
%         XMat.SM=XMat.X./XMat.Theta;
%         XMat.LSM=log(XMat.X./XMat.Theta);
%         
%         
%         PMat.DiffMean=XMat.X-XMat.Theta;
%         PMat.SM=XMat.X./XMat.Theta;
%         PMat.LSM=log(XMat.X./XMat.Theta);
%         
%         n=length(XMat.Theta);
%         DP=eye(n,n).*sum(P,2);
%         L=DP-P;
%         
%         XMat.xfull=(eye(n,n)+g.*1.*L)\XMat.Theta;
%         
%         
%         for k =1:n
%             prow=P(k,:);
%             XMat.XP1(k)=(1/(1+sum(prow)*g))*XMat.Theta(k)+(g/(1+sum(prow)*g))*prow*XMat.Theta;
%         end
%         
%         
%         XMat.P1DiffMean=XMat.XP1-XMat.Theta;
%         PMat.P1DiffMean=XMat.XP1-XMat.Theta;
%         XMat.P1SM=XMat.XP1./XMat.Theta;
%         PMat.P1SM=XMat.XP1./XMat.Theta;
%         XMat.P1LSM=log(XMat.XP1./XMat.Theta);
%         PMat.P1LSM=log(XMat.XP1./XMat.Theta);
%         
%         row.LSM=mean(XMat.LSM);
%         row.SM=mean(XMat.SM);
%         row.PDiffMean=mean(XMat.P1DiffMean);
%         row.PSM=mean(XMat.P1SM);
%         row.PLSM=real(mean(XMat.P1LSM));
%         
%         row.TSkew=skewness(PMat.Theta);
%         row.XSkew=skewness(PMat.X);
%         
%         
%         row.Tmax=max(PMat.Theta);
%         row.TmaxX1=PMat.P1SM(PMat.Theta==row.Tmax).*row.Tmax;
%         row.TmaxX=PMat.X(PMat.Theta==row.Tmax);
%         row.TmaxSM=row.TmaxX/row.Tmax;
%         row.TmaxId=PMat.Identity(PMat.Theta==row.Tmax);
%         
%         
%         
%         
%         %% Critical Nodes for slackers and climbers
%         
%         
%         identity=PMat.Identity;
%         % Rounding
%         P=round(P,4);
%         % Symmetric Graph
%         SPMat=(P+P')./2;
%         SPMat(SPMat>0)=1;
%         GPgraph=graph(SPMat);
%         GPgraph.Nodes.Type=identity;
%         GPgraph.Nodes.X=PMat.X;
%         GPgraph.Nodes.Num=[1: length(identity)]';
%         GPgraph.Nodes.theta=round(PMat.Theta,3);
%         GPgraph.Nodes.diff=round(PMat.DiffMean,3);
%         % Delete Watchers from Graph
%         del=GPgraph.Nodes.Num(identity==0)';
%         GPgraph=GPgraph.rmnode(del);
%         % Save critical nodes
%         [edgebins,iC] = biconncomp(GPgraph);
%         row.critNodesCS=length(iC);
%         
%         % Plotting routine
%         %         strtheta=string(GPgraph.Nodes.theta);
%         %         strType=string(GPgraph.Nodes.Type);
%         %         strleer="|";
%         %         strdiff=string(GPgraph.Nodes.diff);
%         %         strX=string(GPgraph.Nodes.X);
%         %         type(:) = strcat(strType(:),strleer,strtheta(:),strleer,strX(:),strleer,strdiff(:));
%         %         asdf= cellstr(type');
%         %         h=plot(GPgraph,'Layout','force','NodeLabel',asdf);
%         %         h.NodeCData=GPgraph.Nodes.Type;
%         %         h.MarkerSize=8;
%         %         layout(h,'force')
%         %         h.LineWidth=2
%         
        
        
        
        newTable(i,:)=row;
    else
        disp(['PMat was empty in i= ',num2str(i)])
    end
end


%writetable(newTable,strcat(filenameTable,'SPNEResAll2.csv'));