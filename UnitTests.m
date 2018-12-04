G=[ 0 1 1 0;
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

theta=[1,2,3,4];
identity=[1,0,-1,0];

g=1;
n = 4;


gemA=1;
gemL=1;
thetaRange = abs(max(theta)-min(theta));
% Motivation functions
PsiL=@(theta_i,theta_j) -gemL.*abs(theta_j-theta_i)+thetaRange;
PsiA=@(theta_i,theta_j) -gemA.*(theta_i-theta_j);
PsiS=@(theta_i,theta_j) -gemA.*(theta_j-theta_i);

%% Get Choice Sets

[ChoiceCell, NrChoices]=GetChoiceSet(G,n);
if isequal(ChoiceCell{1}, [2 3]) && (NrChoices(1)==2)
    disp("GetChoiceSet: Success");
else
    error("GetChoieSet: Error");
end



%% RecoverPi.m

[ChoiceCell, NrChoices]=GetChoiceSet(G,n);
prest=[1, 0]';
p_i = RecoverPi(prest,ChoiceCell{1}, n);
prest=[1, 1]';
p_i2 = RecoverPi(prest,ChoiceCell{2}, n);
p_i3 = RecoverPi(prest,ChoiceCell{3}, n);

if isequal(p_i, [0;1;0;0]) && isequal(p_i2, [1;0;0;1]) && isequal(p_i3, [1;0;0;1])
    disp("RecoverPi: Success");
else
    error("RecoverPi: Error");
end


%% Constraint Calculation
ConA={};
Conb={};
[ChoiceCell, NrChoices]=GetChoiceSet(G,n);

[ConA,Conb] = CalcConstraints(n, NrChoices);

prest=[1, 0]';
if ConA{1}*prest-Conb{1} <= 0
    disp("CalcConstraints: Success");
else
    error("CalcConstraints: Error");
end


%% Discrete Choice

[ChoiceCell, NrChoices]=GetChoiceSet(G,n);
Choice=ChoiceCell{1};
P_t_1=zeros(n,n);
nrChoice=NrChoices(1);
delta=1;
i=1;
[util,p_i_star]=DiscreteChoice(P_t_1,g,delta,theta(1), theta', PsiA,i,Choice,nrChoice)


%% Bonacich centralities
normAlpha=1;
normBeta=1;
alpha=1;
beta=1;
Power=0;

centralities = BonacichCentrality(alpha,beta,normAlpha,normBeta,Power,G4)