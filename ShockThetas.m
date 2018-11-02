function [thetas, shocks] = ShockThetas(thetavec, boundpc,m)

% Set Random Stream
s = RandStream('mcg16807','Seed',m);
RandStream.setGlobalStream(s);


thetarange=max(thetavec)-min(thetavec);
dims=size(thetavec);
n=max(dims);
sd=thetarange*boundpc;
pd = makedist('Normal','mu',0,'sigma',sd);
shocks=random(pd,dims(1),dims(2));
thetas=thetavec+shocks;
end