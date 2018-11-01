%%
% utilitySPNE
% Given a set of P-vectors from t-1
% this function gives the utility of p_i
% assuming that the SPNE will be played
%%
function [util,p_i_star]=DiscreteChoice(P_t_1,g,delta,theta_i, theta, Psi,i,Choice,nrChoice)


n=length(theta);
% Define One vector
ez=ones(n,1);

% Default values
util=0;
p_i_star=zeros(n,1);

% Iterate over choice sets
for z = 1:nrChoice
    % Create a prest
    prest=zeros(nrChoice,1);
    prest(z)=1;
    % Recover p_i from p restriction
    p_i = RecoverPi(prest, Choice, n);
    
    % Create new P matrix
    P=P_t_1;
    % Assemble new P
    P(i,:)=p_i';
    % Get expected x
    x=XFOCSPNE(P,delta,theta,g);
    x_i=x(i);
    % Private utility, positive part
    PrivUtil=(x_i-theta_i)^2;
    % Calculate a benefit vector for each potential peer
    PsiVec = Psi(theta_i,theta);
    % Make sure no connection to oneself
    PsiVec(i)=-100;
    % Expected benefit
    CBenUtil=p_i'*(PsiVec);
    % Expected non-alignment cost
    ConfUtil=p_i'*((x_i.*ez-x).*(x_i.*ez-x));
    % Full utility
    NewUtil=-PrivUtil+g^(1).*delta.*CBenUtil-g.*delta.*ConfUtil;
    % If this choice is better than the previous one, do this
    if NewUtil > util
       p_i_star=p_i;
       util=NewUtil;
    end
end






end