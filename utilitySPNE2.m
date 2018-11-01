%%
% utilitySPNE
% Given a set of P-vectors from t-1
% this function gives the utility of p_i
% assuming that the SPNE will be played
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