%%
%   returnRow.m
%   Returns the row of the current firm to save
%%
% @param: firm- a firm
% @return: row - a 1-row table
%%
function row = returnRow(firm,T)
firmID=string(firm.firmID);
n=firm.n;
e=firm.eMat(T);
m=firm.m;
diffM=firm.diffM(T);
SM=firm.SM(T);
varSM=firm.varSM(T);
rationality=firm.rationality;
maxDiff=firm.maxDiff(T);
minDiff=firm.minDiff(T);
thetaRange=firm.thetaRange;
realSkew=firm.skew(T);
realCons=firm.cons(T);
paramCons=firm.startCons;
NrC=firm.NrC;
NrW=firm.NrW;
NrS=firm.NrS;
Talpha=firm.Talpha;
Tbeta=firm.Tbeta;
row=table(firmID,T,n,m,e,rationality,NrC,NrW,NrS,Talpha,Tbeta,realSkew,realCons,paramCons,SM,varSM,diffM,maxDiff,minDiff,thetaRange);
end