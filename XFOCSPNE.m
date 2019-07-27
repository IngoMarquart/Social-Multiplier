%%
% % XFOCSPNE
% Solves the FOC of all actors, for given attention choices P and theta, for their effort
% levels X.
% @param: P - Attention matrix
% @param: g, delta - Embedding and exogenous weight
% @param: theta - Skill levels
% @return: x - vector of effort levels
%% 

function x=XFOCSPNE(P,delta,theta,g)
    theta=theta(:);
    % Create Laplacian matrix
    n=length(theta);
    DP=eye(n,n).*sum(P,2);
    L=DP-P;
    
    % Calculate X-vector given P
    %x=(eye(n,n)+delta.*L)\theta;
    % Use this BR function if g is in front of conf-cost
    x=(eye(n,n)+g.*delta.*L)\theta;
    
end