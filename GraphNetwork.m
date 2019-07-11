%%
% % GraphNetwork
% This function graphs the G network and highlights the A network on top of it.
% @param: Gmat - Symmetric G matrix
% @param: SPmat - Symmetric P matrix
% @param: identity, theta, X: nx1 vectors of identity, theta and X
%%

function n=GraphNetwork(firm)

% Unpack firm variables
PMat=firm.aMat{firm.T};
X=firm.xMat(:,firm.T);
theta=firm.thetaMat(:,firm.T);
identity=firm.muMat(:,firm.T);

% Symmetrize network
SPMat=(PMat+PMat')./2;
SPMat(SPMat>0)=1;

PMat(PMat<0.1)=0;
SPMat(SPMat<0.1)=0;
GPgraph=digraph(firm.gMat);
Sgraph=graph(SPMat);

Pgraph=digraph(PMat);
diff=X-theta;
n=numel(theta);


%% TODO fix here
RetStruct.NrClimbers=length(identity(identity==1));
RetStruct.NrWatchers=length(identity(identity==0));
RetStruct.NrSlackers=length(identity(identity==-1));
RetStruct.TSkew=skewness(theta);
RetStruct.XMean=mean(X);
RetStruct.DiffMean=mean(diff);
RetStruct.ThetaMean=mean(theta);

if sum(sum(firm.gMat)) == ((n^2)-n) % Complete G network
    GPgraph=Pgraph;
    useG=0;
else
    useG=1;
end
type=string(ones(size(identity)));
type(identity==1)="C";
type(identity==0)="W";
type(identity==-1)="S";
GPgraph.Nodes.Type=identity;
GPgraph.Nodes.X=X;
GPgraph.Nodes.Num=[1: length(identity)]';
GPgraph.Nodes.theta=round(theta,3);
GPgraph.Nodes.diff=round(diff,3);

strtheta=string(GPgraph.Nodes.theta);
strType=string(identity);
strType(identity==1)="I";

strType(identity==0)="R";
strType(identity==-1)="H";
strleer="|";
strdiff=string(GPgraph.Nodes.diff);
strX=string(GPgraph.Nodes.X);
type(:) = strcat(strType,strleer,strtheta(:),strleer,strleer,strX(:),strleer,strdiff(:));
type(:) = strcat("    ",strtheta(:),strleer,strdiff(:));
%type(:)=strType;
asdf= cellstr(type');


if useG==0
LWidths = 6*GPgraph.Edges.Weight/max(GPgraph.Edges.Weight);
h=plot(GPgraph,'Layout','force','NodeLabel',asdf,'LineWidth',LWidths);
else
    h=plot(GPgraph,'Layout','force','NodeLabel',asdf);

end

%h.EdgeCData = GPgraph.Edges.Weight;

h.NodeCData=identity;
h.MarkerSize=15;
%h.Marker='<';
%h.NodeFontSize=10;
layout(h,'force','UseGravity',true,'WeightEffect','inverse');
%[~,b]=max(GPgraph.Nodes.theta);
%layout(h,'layered','Direction','right','Sources',[b]);
%layout(h,'circle','Center',[b]);

if useG==1
    highlight(h,Pgraph,'EdgeColor','r','LineWidth',4);
end
txt = {'Attention network:',['N=',num2str(n)],['e=',num2str(firm.e)], ...
    ['NrC=',num2str(RetStruct.NrClimbers)],['NrW=',num2str(RetStruct.NrWatchers)],['NrS=',num2str(RetStruct.NrSlackers)],['Skew=',num2str(skewness(theta))], ...
    ['AvgDiff=',num2str(RetStruct.DiffMean)],['AvgTheta=',num2str(RetStruct.ThetaMean)],['AvgX=',num2str(RetStruct.XMean)]};
annotation('textbox',...
    [0.14 0.9 0 0],...    
    'String',txt);



%h.LineWidth=2.2;
%layout(h,'force','Iterations',5000)


%Highlight articulation points
%highlight(h, iC)


end