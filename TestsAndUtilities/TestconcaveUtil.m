%% This was a failed attempt to write some automatic tests...

G5=[ 0 1 1 0;
    1 0 0 1;
    1 0 0 1;
    0 1 0 0];

G4=[ 0 1 1 0;
    1 0 0 0;
    1 0 0 1;
    0 0 1 0];

G2=[ 0 1 0 0;
    1 0 0 0;
    0 0 0 1;
    0 0 1 0];

G3=zeros(4,4);
G=ones(4,4)-eye(4,4);
theta=[1,2,3,4];
identity=[1,0,-1,0];

e=1;
n = 4;


gemA=1;
gemL=1;
thetaRange = abs(max(theta)-min(theta));
% Motivation functions
PsiL=@(theta_i,theta_j) -gemL.*abs(theta_j-theta_i)+thetaRange;
PsiA=@(theta_i,theta_j) -gemA.*(theta_i-theta_j);
PsiS=@(theta_i,theta_j) -gemA.*(theta_j-theta_i);



[ChoiceCell, NrChoices]=GetChoiceSet(G,n);
Choice=ChoiceCell{1};
prevAttention=zeros(n,n);
nrChoice=NrChoices(1);
i=1;
aiRest=[1,0,0]';

[curUi, curAi] = DiscreteChoice(prevAttention, e, 1, theta(i), theta, PsiA, i, ChoiceCell{i}, nrChoice, 0)

concaveChoice(aiRest, prevAttention, theta, theta(i),e, PsiA, 0, 1, ChoiceCell{i},i)
