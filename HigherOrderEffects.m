
runtime='HighN40-55';
filenameTable=['DataSave/',runtime,'/'];
filenamePMats=['DataSave/',runtime,'/Pmats/'];

allTable=readtable(strcat(filenameTable,'SPNEResAll.csv'));





alllength=height(allTable);
%allTable=allTable(:,2:end);
newTable=allTable;
newTable.LSM=zeros(alllength,1);
newTable.SM=zeros(alllength,1);
newTable.PDiffMean=zeros(alllength,1);
newTable.PSM=zeros(alllength,1);
newTable.PLSM=zeros(alllength,1);
newTable.TSkew=zeros(alllength,1);
newTable.XSkew=zeros(alllength,1);
newTable.Tmax=zeros(alllength,1);

newTable.TmaxX=zeros(alllength,1);
newTable.TmaxX1=zeros(alllength,1);
newTable.TmaxSM=zeros(alllength,1);
newTable.TmaxId=zeros(alllength,1);
newTable.critNodesCS=zeros(alllength,1);


parfor i = 1:alllength
    row=newTable(i,:);
    g=row.g;
    
    filenamePMat=strcat(filenamePMats,num2str(row.timestamp),'_P_n',num2str(row.n),'_m',num2str(row.m),'.csv');
    PMat=readtable(filenamePMat);
    if ~isempty(PMat)
        XMat=PMat(:,1:4);
        P=PMat{:,5:(4+row.n)};
        XMat.DiffMean=XMat.X-XMat.Theta;
        XMat.SM=XMat.X./XMat.Theta;
        XMat.LSM=log(XMat.X./XMat.Theta);
        
        
        PMat.DiffMean=XMat.X-XMat.Theta;
        PMat.SM=XMat.X./XMat.Theta;
        PMat.LSM=log(XMat.X./XMat.Theta);
        
        n=length(XMat.Theta);
        DP=eye(n,n).*sum(P,2);
        L=DP-P;
        
        XMat.xfull=(eye(n,n)+g.*1.*L)\XMat.Theta;
        
        
        for k =1:n
            prow=P(k,:);
            XMat.XP1(k)=(1/(1+sum(prow)*g))*XMat.Theta(k)+(g/(1+sum(prow)*g))*prow*XMat.Theta;
        end
        
        
        XMat.P1DiffMean=XMat.XP1-XMat.Theta;
        PMat.P1DiffMean=XMat.XP1-XMat.Theta;
        XMat.P1SM=XMat.XP1./XMat.Theta;
        PMat.P1SM=XMat.XP1./XMat.Theta;
        XMat.P1LSM=log(XMat.XP1./XMat.Theta);
        PMat.P1LSM=log(XMat.XP1./XMat.Theta);
        
        row.LSM=mean(XMat.LSM);
        row.SM=mean(XMat.SM);
        row.PDiffMean=mean(XMat.P1DiffMean);
        row.PSM=mean(XMat.P1SM);
        row.PLSM=real(mean(XMat.P1LSM));
        
        row.TSkew=skewness(PMat.Theta);
        row.XSkew=skewness(PMat.X);
        
        
        row.Tmax=max(PMat.Theta);
        row.TmaxX1=PMat.P1SM(PMat.Theta==row.Tmax).*row.Tmax;
        row.TmaxX=PMat.X(PMat.Theta==row.Tmax);
        row.TmaxSM=row.TmaxX/row.Tmax;
        row.TmaxId=PMat.Identity(PMat.Theta==row.Tmax);
        
        
        
        
        %% Critical Nodes for slackers and climbers
        
        
        identity=PMat.Identity;
        % Rounding
        P=round(P,4);
        % Symmetric Graph
        SPMat=(P+P')./2;
        SPMat(SPMat>0)=1;
        GPgraph=graph(SPMat);
        GPgraph.Nodes.Type=identity;
        GPgraph.Nodes.X=PMat.X;
        GPgraph.Nodes.Num=[1: length(identity)]';
        GPgraph.Nodes.theta=round(PMat.Theta,3);
        GPgraph.Nodes.diff=round(PMat.DiffMean,3);
        % Delete Watchers from Graph
        del=GPgraph.Nodes.Num(identity==0)';
        GPgraph=GPgraph.rmnode(del);
        % Save critical nodes
        [edgebins,iC] = biconncomp(GPgraph);
        row.critNodesCS=length(iC);
        
        % Plotting routine
        %         strtheta=string(GPgraph.Nodes.theta);
        %         strType=string(GPgraph.Nodes.Type);
        %         strleer="|";
        %         strdiff=string(GPgraph.Nodes.diff);
        %         strX=string(GPgraph.Nodes.X);
        %         type(:) = strcat(strType(:),strleer,strtheta(:),strleer,strX(:),strleer,strdiff(:));
        %         asdf= cellstr(type');
        %         h=plot(GPgraph,'Layout','force','NodeLabel',asdf);
        %         h.NodeCData=GPgraph.Nodes.Type;
        %         h.MarkerSize=8;
        %         layout(h,'force')
        %         h.LineWidth=2
        
        
        
        
        newTable(i,:)=row;
    else
        disp(['PMat was empty in i= ',num2str(i)])
    end
end


writetable(newTable,strcat(filenameTable,'SPNEResAll2.csv'));