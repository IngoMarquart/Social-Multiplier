function [util] =  concaveChoice(aiRest,A,theta,theta_i,e,Psi,rationality, conParam)

% First we need to recover the true attention vector
% from the restricted choices given G
a_i = RecoverPi(aiRest, ChoiceCell, n)

%% Second stage: Recover X
ebar=e/(e+1);
if rationality>0
    % rational x
    ratx=XFOCSPNE(A,1,theta,e);
else
    ratx=theta;
end
% All other x_j
x=(1-rationality).*theta+rationality.*ratx;
% x_i based on this
x(i)=(1-ebar).*theta_i+ebar.*a_i*x;


%% First stage: Recover utility for the given a_i

%% Private cost
PrivUtil=(x-theta_i).^2;
%% Social Cost
% Expected non-alignment cost
ConfUtil=a*((x-theta').*(x-theta'));
%% Concave Social Benefit function
% Calculate a benefit vector for each potential peer
PsiVec = Psi(theta_i,theta);
% Expected benefit
% CBenUtil=(a.^(1/3))*(PsiVec.^(1/1)).^(3/2);
PsiVec = max(0,theta'-theta_i);


isoutil=@(x,rho) (x.^(1-rho)-1)/(1-rho);
exputility=@(x,rho) (1-exp(-(rho.*x)));
rho=20;
CBenUtil=exputility(a,rho)*PsiVec;
%CBenUtil=(a.^(1/3))*(PsiVec.^(1/3)).^(3);

%% Full utility
NewUtil=-PrivUtil+e^(1).*1.*CBenUtil-e.*ConfUtil;
% return negative value for minimizer
util=-mean(NewUtil);

end