
params.n=100;
params.m=1;
params.e=1;
params.gamma=[2/3, 1/6, 1/6];
params.thetaD=[2,5 ,1,1];
params.cons=0;
params.pn=0;
params.mn=4;
params.maxT=40;
params.minEqmT=5;
params.maxEqmT=100;
params.globalsearch=-1;
params.thetaRepShockVar=0.05;
params.rationality=1;
firm=initFirm(params);
for T=firm.T:firm.maxT
    firm.T=T;
    firm=runTurn(firm);
    
end



GraphNetwork(firm);