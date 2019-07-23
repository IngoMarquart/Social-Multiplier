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
expSM = firm.expSM(T);
avgX = firm.avgX(T);
avgTheta = firm.avgTheta(T);
varSM = firm.varSM(T);
rationality = firm.rationality;
maxDiff = firm.maxDiff(T);
minDiff = firm.minDiff(T);
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
thetaStartStr=string(strcat(num2str(sort(firm.thetaMat(:,1)'),"%.2f,")));
thetaStr=string(strcat(num2str(sort(firm.thetaMat(:,T)'),"%.2f,")));
xStr=string(strcat(num2str(sort(firm.xMat(:,T)'),"%.2f,")));
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
row = table(firmID, T, n, m, e, rationality, NrC, NrW, NrS,ProbC,ProbW,ProbS, Talpha, Tbeta, realSkew, realCons, paramCons, avgX, avgTheta, expSM, SM, varSM, diffM, maxDiff, minDiff, thetaRange, conUtil, conParam, gMethod,gNrComponents,gRadius,gDiameter,gAvgPathLength,gDensity,gMaxEV,gAvgClustering,gAvgDegree, gSymmetry, gClusters, gModularity, gLinks,gAssembly,gMn,gMr,gPn,ShufflePositions,thetaStartStr,thetaStr,xStr,ceoType,startCeoType);
end
