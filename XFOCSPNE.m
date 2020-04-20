%%
% % XFOCSPNE
% Solves the FOC of all actors, for given attention choices P and theta, for their effort
% levels X.
% @param: A - Attention matrix
% @param: e- Embedding
% @param: theta - theta levels
% @return: x - vector of effort levels
%% 

function x=XFOCSPNE(A,theta,e)
    theta=theta(:);
    % Create Laplacian matrix
    n=length(theta);
    D=eye(n,n).*sum(A,2);
    L=D-A;
    e_bar=e/(1-e);
    
    % Calculate X-vector given P
    %x=(eye(n,n)+delta.*L)\theta;
    % Use this BR function if g is in front of conf-cost
    x=(eye(n,n)+e_bar.*L)\theta;
    
end