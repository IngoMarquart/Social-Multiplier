
params.n=20;
params.m=123;
params.e=1;
params.gamma=[1/3, 1/3, 1/3];
params.thetaD=[2,2 ,1,1];
params.cons=0;
params.pn=0;
params.mn=5;
params.maxT=1;
params.minEqmT=5;
params.maxEqmT=100;
params.globalsearch=-1;
params.thetaRepShockVar=0;
params.rationality=0;
firm=initFirm(params);
for T=firm.T:firm.maxT
    firm.T=T;
    firm=runTurn(firm);
    
end



GraphNetwork(firm);