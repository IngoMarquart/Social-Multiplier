%%
%   returnRow.m
%   Returns the row of the current firm to save
%%
% @param: firm- a firm
% @return: row - a 1-row table
%%
function row = returnRow(firm, T)
firmID = string(firm.firmID);
n = firm.n;
e = firm.eMat(T);
m = firm.m;
diffM = firm.diffM(T);
SM = firm.SM(T);
avgX = firm.avgX(T);
avgTheta = firm.avgTheta(T);
rationality = firm.rationality;
maxDiff = firm.maxDiff(T);
minDiff = firm.minDiff(T);
sumX=firm.sumX(T);
maxX=firm.maxX(T);
minX=firm.minX(T);
maxTheta=firm.maxTheta(T);
minTheta=firm.minTheta(T);
varX=firm.varX(T);
varTheta=firm.varTheta(T);
thetaRange = firm.thetaRange;
realSkew = firm.skew(T);
realCons = firm.cons(T);
paramCons = firm.startCons;
NrC = firm.NrC;
NrW = firm.NrW;
NrS = firm.NrS;
Talpha = firm.Talpha;
Tbeta = firm.Tbeta;
ProbC=firm.gamma(1);
ProbW=firm.gamma(2);
ProbS=firm.gamma(3);


flucX=firm.flucX;
flucT=firm.flucT;

thetaStartStr=string(strcat(num2str(firm.thetaMat(:,1)',"%.2f,")));
thetaStr=string(strcat(num2str(firm.thetaMat(:,T)',"%.2f,")));
xStr=string(strcat(num2str(firm.xMat(:,T)',"%.2f,")));
muStr=string(strcat(num2str(firm.muMat(:,T)',"%.2f,")));


pgrank=string(strcat(num2str(firm.pgrank',"%.2f,")));
indegree=string(strcat(num2str(firm.indegree',"%.2f,")));
peerX=string(strcat(num2str(firm.peerX',"%.2f,")));
peer=string(strcat(num2str(firm.peer',"%.0f,")));
peerMu=string(strcat(num2str(firm.peerMu',"%.0f,")));


conUtil = firm.conUtil;
conParam = firm.conParam;
ShufflePositions=firm.shufflePositions;
ceoType=firm.ceoAct;
startCeoType=firm.startCeoAct;
gAvgPathLength=firm.gAvgPathLength;
gDensity=firm.gDensity;
gMaxEV=firm.gMaxEV;
gAvgClustering=firm.gAvgClustering;
gAvgDegree=firm.gAvgDegree;
gRadius=firm.gRadius;
gDiameter=firm.gDiameter;
ceoStartT=firm.ceoStartT;
learningRate=firm.learningRate;
gNrComponents=firm.gNrComponents;
if firm.gMethod=="Task" || firm.gMethod=="TaskAssembly"
    gMethod=firm.gMethod;
    gSymmetry=firm.gSymmetry;
    gClusters=firm.gCluster;
    gModularity=firm.gModularity;
    gLinks=firm.gLinks;
    gAssembly=firm.gAssembly;
    gMn=0;
    gMr=0;
    gPn=0;
    gPr=0;
else
    gMethod=firm.gMethod;
    gSymmetry=0;
    gClusters=0;
    gModularity=0;
    gLinks=0;
    gAssembly=0;
    gMn=0;
    gMr=0;
    gPn=0;
    gPr=0;
    if firm.gMethod=="JR"
        gMn=firm.gMn;
        gMr=firm.gMr;
        gPn=firm.gPn;
        gPr=firm.gPr;
    end
end

if firm.gMethod == "Full"
    row = table(firmID, T, n, m, e, rationality, NrC, NrW, NrS,ProbC,ProbW,ProbS, Talpha, Tbeta, realSkew, realCons, paramCons, avgX, avgTheta, sumX, maxX, minX, maxTheta,minTheta,varX,varTheta, SM, diffM, maxDiff, minDiff, thetaRange, conUtil, conParam, gMethod,ShufflePositions,thetaStartStr,thetaStr,xStr,muStr,ceoType,startCeoType,ceoStartT,learningRate,pgrank,indegree,peer,peerX,peerMu,flucX,flucT);
else
    row = table(firmID, T, n, m, e, rationality, NrC, NrW, NrS,ProbC,ProbW,ProbS, Talpha, Tbeta, realSkew, realCons, paramCons, avgX, avgTheta, sumX, maxX, minX, maxTheta,minTheta,varX,varTheta, SM, diffM, maxDiff, minDiff, thetaRange, conUtil, conParam, gMethod,gNrComponents,gRadius,gDiameter,gAvgPathLength,gDensity,gMaxEV,gAvgClustering,gAvgDegree, gSymmetry, gClusters, gModularity, gLinks,gAssembly,gMn,gMr,gPn,ShufflePositions,thetaStartStr,thetaStr,xStr,muStr,ceoType,startCeoType,ceoStartT,learningRate,pgrank,indegree,peer,peerX,peerMu,flucX,flucT);
end
