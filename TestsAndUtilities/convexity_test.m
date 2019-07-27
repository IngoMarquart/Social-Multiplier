s = RandStream('mcg16807','Seed',100);
RandStream.setGlobalStream(s);

sigma=0;
trials=100;
start_theta=[5,8,12];
a=[0,0,1];
aStart=[0,0,1];
theta_i=0;
e=0.3;
n=length(start_theta);


sigmavec=0:0.5:5;
hvec=[];
pd = makedist('Normal','mu',0,'sigma',sigma);
thetaShock=random(pd,trials,length(start_theta));
theta=repmat(start_theta,trials,1);%+thetaShock;


evec=[0.1:0.1:10];
% y=[0;0;0];
% for e=evec
%    ay=[-utility([1,0,0],theta,theta_i,e);-utility([0,1,0],theta,theta_i,e);-utility([0,0,1],theta,theta_i,e)];
%    y=[y,ay]; 
% end
% y=y(:,2:end)';
% plot(evec',y);



options  =  optimset('Algorithm','sqp','Display', 'none','UseParallel',false);
lb = zeros(size(aStart));
ub = ones(size(aStart));
A = ones(size(aStart));
b = [1];
Aeq = [];
beq = [];
gs = GlobalSearch('Display','off','StartPointsToRun','all');




y=[0;0;0];
for e=evec
   fnc=@(a)utility(a,theta,theta_i,e);
   %x = fmincon(fnc,aStart,A,b,Aeq,beq,lb,ub,[],options)';
   problem = createOptimProblem('fmincon','x0',aStart,'objective', fnc,... 
                  'lb',lb,'ub',ub, ... 
                   'Aineq',A,'bineq',b,'options',options); 
   x=run(gs,problem)';
   y=[y,x]; 
end
y=y(:,2:end)';
%plot(evec',y);
hh=y./sum(y,2);
h=(sum(hh.^2,2)-(1/n))/(1-1/n);
vy=var(y')';
y=[y,h];
stackedplot(y);
hvec=[hvec,h];

%end




function [util] =  utility(a,theta,theta_i,e)
ebar=e/(e+1);
x=(1-ebar).*theta_i+ebar.*a*theta';

% Private utility, positive part
PrivUtil=(x-theta_i).^2;
% Calculate a benefit vector for each potential peer
PsiVec = max(0.0001,theta'-theta_i);
% Expected benefit
% CBenUtil=(a.^(1/3))*(PsiVec.^(1/1)).^(3/2);


isoutil=@(x,rho) (x.^(1-rho))./(1-rho)*PsiVec;
exputility=@(x,rho) (1-exp(-(rho.*x)))*PsiVec;
cesutil=@(x,rho) (x.^(1/rho))*(PsiVec).^(rho);
rho=2;
CBenUtil=cesutil(a,rho);
%CBenUtil=(a.^(1/3))*(PsiVec.^(1/3)).^(3);
% Expected non-alignment cost
ConfUtil=a*((x-theta').*(x-theta'));
% Full utility
NewUtil=-PrivUtil+e.*CBenUtil-e.*ConfUtil;
%NewUtil=sqrt(NewUtil);
util=-mean(NewUtil);

end