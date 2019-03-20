%%
% % GraphNetwork
% This function graphs the G network and highlights the A network on top of it.
% @param: Gmat - Symmetric G matrix
% @param: SPmat - Symmetric P matrix
% @param: identity, theta, X: nx1 vectors of identity, theta and X
%% 

function n=GraphNetwork(Gmat, SPMat, PMat, identity, theta, X,g, RetStruct)

    PMat(PMat<0.1)=0;
    SPMat(SPMat<0.1)=0;
    GPgraph=graph(Gmat);
    Sgraph=graph(SPMat);
    Pgraph=digraph(PMat);
    diff=X-theta;
    n=numel(theta);

    if sum(sum(Gmat)) == ((n^2)-n) % Complete G network
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
    strleer="|";
    strdiff=string(GPgraph.Nodes.diff);
    strX=string(GPgraph.Nodes.X);
    type(:) = strcat(strType,strleer,strtheta(:),strleer,strleer,strX(:),strleer,strdiff(:));
    type(:) = strcat(strtheta(:),strleer,strdiff(:));
    asdf= cellstr(type');
    h=plot(GPgraph,'Layout','force','NodeLabel',asdf);
    %h=plot(GPgraph,'Layout','force');
    h.NodeCData=identity;
    h.MarkerSize=12;
    layout(h,'force')
    if useG==1
    highlight(h,Sgraph,'EdgeColor','r','LineWidth',1.5)
    end
    txt = {'Attention network:',['N=',num2str(n)],['g=',num2str(g)], ...
        ['NrC=',num2str(RetStruct.NrClimbers)],['NrW=',num2str(RetStruct.NrWatchers)],['NrS=',num2str(RetStruct.NrSlackers)],['Skew=',num2str(skewness(theta))], ...
        ['Diff=',num2str(RetStruct.DiffMean)],['AvgTheta=',num2str(RetStruct.ThetaMean)],['AvgX=',num2str(RetStruct.XMean)]};
    annotation('textbox',...
        [0.14 0.9 0 0],...
        'String',txt);
    
    
    %h.LineWidth=2.2;
    %layout(h,'force','Iterations',5000)
    
    %h.EdgeCData = GPgraph.Edges.Weight;
    
    %Highlight articulation points
    %highlight(h, iC)
 

 end