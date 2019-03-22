
params.n=6;
params.m=1;
params.e=1;
params.gamma=[2/3, 1/6, 1/6];
params.thetaD=[2,5 ,1,1];
params.cons=0;
params.pn=0.2;
params.mn=2;
params.maxT=2;
params.minEqmT=5;
params.maxEqmT=100;
params.globalsearch=-1;

firm=initFirm(params);
for T=firm.T:firm.maxT
    firm.T=T;
    firm=runTurn(firm);
    
end



GraphNetwork(firm);