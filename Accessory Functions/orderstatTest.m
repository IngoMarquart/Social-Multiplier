nvec=[2,5,8,10,20,30,40,50,70,100];
nvec=100
for n=nvec
%     thetaDR=[2,15,5.4375,33.5];
%     thetaDL=[15,2,5.4375,33.5];
    thetaDR=[2,4,5.4375,13.8];
    thetaDL=[4,2,5.4375,13.8];
    
    maxIter=100000;
    minLeft=zeros(maxIter,1);
    minRight=zeros(maxIter,1);
    maxLeft=zeros(maxIter,1);
    maxRight=zeros(maxIter,1);
    meanLeft=zeros(maxIter,1);
    meanRight=zeros(maxIter,1);
    varLeft=zeros(maxIter,1);
    varRight=zeros(maxIter,1);   
    
    for Iter=1:maxIter
    pd = makedist('Beta',thetaDR(1),thetaDR(2));
    
    thetaRight = (random(pd,n,1)-pd.mean).*thetaDR(4)+thetaDR(3);
    pd = makedist('Beta',thetaDL(1),thetaDL(2));
    thetaLeft = (random(pd,n,1)-pd.mean).*thetaDL(4)+thetaDL(3);
    
    minLeft(Iter)=min(thetaLeft);
    minRight(Iter)=min(thetaRight);
    maxLeft(Iter)=max(thetaLeft);
    maxRight(Iter)=max(thetaRight);
    meanLeft(Iter)=mean(thetaLeft);
    meanRight(Iter)=mean(thetaRight);
    varLeft(Iter)=var(thetaLeft);
    varRight(Iter)=var(thetaRight);
    end

clear xvec yvec


[yvec(:,2),xvec(:,1)]=hist(maxRight,50);
[yvec(:,1),xvec(:,1)]=hist(maxLeft,xvec(:,1));

%set(0,'DefaultFigureVisible','off');
% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to bar
bar1 = bar(xvec,yvec,'BarWidth',1);
set(bar1(1),'DisplayName','Left Skewed Distribution');
set(bar1(2),'DisplayName','Right Skewed Distribution','FaceAlpha',0.7);


line1=line([median(maxLeft), median(maxLeft)], [0, max(yvec(:,1))], 'Color', 'r', 'LineWidth', 3);
line2=line([median(maxRight), median(maxRight)], [0, max(yvec(:,1))], 'Color', 'b', 'LineWidth', 3);

line3=line([mean(maxLeft), mean(maxLeft)], [0, max(yvec(:,1))], 'Color', 'm', 'LineWidth', 3);
line4=line([mean(maxRight), mean(maxRight)], [0, max(yvec(:,1))], 'Color', 'c', 'LineWidth', 3);

set(line1,'DisplayName','Median maxLeft');
set(line2,'DisplayName','Median maxRight');

set(line3,'DisplayName','Mean maxLeft');
set(line4,'DisplayName','Mean maxRight');

% Create plot
%plot(X1,Y1,'LineWidth',2,'Color',[1 0 0]);

% Create plot
%plot(X2,Y2,'LineWidth',2,'Color',[1 0 0]);

tit="F_n, n="+num2str(n)+", a="+num2str(thetaDL(1))+", b="+num2str(thetaDL(2));
tit2="avg mR="+num2str(mean(meanRight)) +", avg mL="+num2str(mean(meanLeft))+", avg vL="+num2str(mean(varLeft))+", avg vR="+num2str(mean(varRight));
% Create title
title({tit,tit2});

box(axes1,'on');
% Create legend
legend(axes1,'show');
filename="F_n-Beta"+num2str(n)+"-"+num2str(thetaDL(1))+"-"+num2str(thetaDL(2))+".png";
saveas(figure1, filename);
hold off


clear xvec yvec


[yvec(:,1),xvec(:,1)]=hist(minLeft,50);
[yvec(:,2),xvec(:,1)]=hist(minRight,xvec(:,1));


% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to bar
bar1 = bar(xvec,yvec,'BarWidth',1);
set(bar1(1),'DisplayName','Left Skewed Distribution');
set(bar1(2),'DisplayName','Right Skewed Distribution','FaceAlpha',0.7);


line1=line([median(minLeft), median(minLeft)], [0, max(yvec(:,2))], 'Color', 'r', 'LineWidth', 3);
line2=line([median(minRight), median(minRight)], [0, max(yvec(:,2))], 'Color', 'b', 'LineWidth', 3);

line3=line([mean(minLeft), mean(minLeft)], [0, max(yvec(:,2))], 'Color', 'm', 'LineWidth', 3);
line4=line([mean(minRight), mean(minRight)], [0, max(yvec(:,2))], 'Color', 'c', 'LineWidth', 3);

set(line1,'DisplayName','Median minLeft');
set(line2,'DisplayName','Median minRight');

set(line3,'DisplayName','Mean minLeft');
set(line4,'DisplayName','Mean minRight');
% Create plot
%plot(X1,Y1,'LineWidth',2,'Color',[1 0 0]);

% Create plot
%plot(X2,Y2,'LineWidth',2,'Color',[1 0 0]);

tit="F_1, n="+num2str(n)+", a="+num2str(thetaDL(1))+", b="+num2str(thetaDL(2));
tit2="avg mR="+num2str(mean(meanRight)) +", avg mL="+num2str(mean(meanLeft))+", avg vL="+num2str(mean(varLeft))+", avg vR="+num2str(mean(varRight));
% Create title
title({tit,tit2});

box(axes1,'on');
% Create legend
legend(axes1,'show');
filename="F_1-Beta"+num2str(n)+"-"+num2str(thetaDL(1))+"-"+num2str(thetaDL(2))+".png";
saveas(figure1, filename);
hold off

end

set(0,'DefaultFigureVisible','on');