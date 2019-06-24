function [util] =  concaveChoice(aiRest,a,theta,theta_i,e,Psi,rationality, conParam, ChoiceCell,i)

% Norm vectors to correspond to column vectors for now
theta=theta(:);
aiRest=aiRest(:);
% First we need to recover the true attention vector
% from the restricted choices given G
n = length(theta);
a_i = RecoverPi(aiRest, ChoiceCell, n);
sumai=sum(a_i);
%% Second stage: Recover X
ebar=e*sumai/(1+e*sumai);
if rationality>0
    % rational x
    ratx=XFOCSPNE(a,1,theta,e);
else
    ratx=theta;
end
% All other x_j
x=(1-rationality).*theta+rationality.*ratx;
% x_i based on this
x(i)=1./(1+e*sumai).*theta_i+(e./(1+e*sumai)).*a_i'*x;


%% First stage: Recover utility for the given a_i

%% Private cost
PrivUtil=(x(i)-theta_i).^2;
%% Social Cost
% Expected non-alignment cost
ConfUtil=a_i'*((theta-x(i)).*(theta-x(i)));
%% Concave Social Benefit function
% Calculate a benefit vector for each potential peer
PsiVec = Psi(theta_i,theta);
% Expected benefit
% CBenUtil=(a.^(1/3))*(PsiVec.^(1/1)).^(3/2);


isoutil=@(r,rho) (r.^(1-rho))./(1-rho)*PsiVec;
exputility=@(r,rho) (1-exp(-(rho.*r)))*PsiVec;
cesutil=@(r,rho) ((r'.^(1/rho))*(PsiVec));


if sum(a_i>0) > 10
    CBenUtil=0;
else
CBenUtil=cesutil(a_i,conParam);
end
%% Full utility
NewUtil=-PrivUtil+e^(1).*1.*CBenUtil-e.*ConfUtil;
% return negative value for minimizer
util=-NewUtil;
end