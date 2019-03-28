
params.n=30;
params.m=1;
params.e=5.5;
params.gamma=[2/3, 1/6, 1/6];
params.thetaD=[2,2 ,1,1];
params.cons=0;
params.pn=0.2;
params.mn=4;
params.maxT=1;
params.minEqmT=5;
params.maxEqmT=100;
params.globalsearch=-1;
params.thetaRepShockVar=0;
params.rationality=1;
firm=initFirm(params);
for T=firm.T:firm.maxT
    firm.T=T;
    firm=runTurn(firm);
    
end



GraphNetwork(firm);