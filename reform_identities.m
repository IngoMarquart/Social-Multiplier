%%
%   initFirm.m
%   For given parameter ranges, initiate a firm
%%
% @param: n - Number of actors
% @param: gamma - Type distribution
% @param: cons - Consolidation
% @param: shufflePositions - shuffle positions
% @param: theta - original theta
% @return: theta, identity - nx2 array corresponding to new identities
%%
function [theta,identity]=reform_identities(n, gamma, shufflePositions, theta)



%% Generation of types and quality vectors
% Matrix to hold identity and thetas
TIVec = zeros(n,2);
TIVec(:,1)=theta;


%% Generate identities
R = mnrnd(n,gamma);
TIVec(:,2) = [1.*ones(1,R(1)), 0.*ones(1,R(2)), -1.*ones(1,R(3))]';

%% Shuffling of positions
% For now, identities are ordered. Under consolidation, so are thetas.
% Here, we allow either shuffling or ordering by Mu

if shufflePositions=="Mu" % ordered by type
    [~,idx]=sort(TIVec(:,2),'descend');
    TIVec(:,2)=TIVec(idx,2);
else % Default: Shuffle
    idx = randperm(size(TIVec,1));
    TIVec(:,2)=TIVec(idx,2);
end


% Prepare vectors
theta=TIVec(:,1);
identity=TIVec(:,2);


end