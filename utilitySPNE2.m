%%
% % utilitySPNE
% Objective function to be optimized
% This version takes a restricted p vector, containing weights only for neighbors in G
% @param: prest - The restricted attention choice: Weights in the space of neighbor nodes
% @param: P_t_1 - Attention matrix of last time period
% @param: g,delta - Embedding and scaling variable
% @param: theta_i, theta - theta vector and theta of i
% @param: Psi - Expected benefit function handle
% @param: ConBen - Exogenous connection benefit
% @param: convexp - Exponent on the expected benefit of connection
% @param: ChoiceCell - Vector given the indecies of valid peers (where connection in G exists)
%% 

function [util,x]=utilitySPNE(prest,P_t_1,g,delta,theta_i, theta, Psi, ConBen,i,convexp,ChoiceCell)

    
    n=length(theta);
    % Define One vector
    ez=ones(n,1);
    
    
    % Create new P matrix with p_i
    P=P_t_1;
    % Recover p_i from p restriction
    p_i = RecoverPi(prest, ChoiceCell, n);
    % Assemble new P
    P(i,:)=p_i';
    

    
    % Get expected x
    x=XFOCSPNE(P,delta,theta,g);
    
    x_i=x(i);
    
    % Private utility
    PrivUtil=(x_i-theta_i)^2;
    
    
    PsiVec = Psi(theta_i,theta)+ConBen;
    PsiVec(i)=-100;
    CBenUtil=p_i'*(PsiVec);
    
    ConfUtil=p_i'*((x_i.*ez-x).*(x_i.*ez-x));
    
    
    util=-PrivUtil+g^(convexp).*delta.*CBenUtil-g.*delta.*ConfUtil;
end