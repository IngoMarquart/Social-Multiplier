SkewLevels = [2:150];
SkewLevels = [0.1:0.01:4];
n=1000;
firmsize=100;
meanRight=zeros(length(SkewLevels),2);
meanLeft=zeros(length(SkewLevels),2);

meanRightEmp=zeros(length(SkewLevels),2);
meanLeftEmp=zeros(length(SkewLevels),2);


meanRightEmpSD=zeros(length(SkewLevels),2);
meanLeftEmpSD=zeros(length(SkewLevels),2);

meanVecOne=zeros(length(SkewLevels),4);

iVec=[1:length(SkewLevels)];
parfor i=iVec
    skew=SkewLevels(i)
    WPvar=(2*skew)/((1+2+skew)*(2+skew)^2);
    WPsd=sqrt(WPvar);
    % Create distributions 
    pdR = makedist('Beta',1,1/skew); 
    pdL = makedist('Beta',1/skew,1); 
     thetaRight = (random(pdR,n,n));
     thetaLeft = (random(pdL,n,n));
    % Calculate nxn sample values from distribution for empirical order
    % statistics
   % thetaRight = ((random(pdR,n,firmsize))-pdR.mean)./WPsd;
   % thetaLeft = ((random(pdL,n,firmsize))-pdL.mean)./WPsd;
    % Calculate empirical moments
    minR=min(thetaRight,[],2);
    maxR=max(thetaRight,[],2);
    GmeanL=mean(thetaLeft,2);
    GmeanR=mean(thetaRight,2);
    avgminR=mean(minR);
    avgmaxR=mean(maxR);
    minL=min(thetaLeft,[],2);
    maxL=max(thetaLeft,[],2);
    avgminL=mean(minL);
    avgmaxL=mean(maxL);
%     % Define PDFs of orderstatistis
%     fminR = @(x) (1-cdf(pdR,x)).^(firmsize-1).*pdf(pdR,x).*firmsize;
%     fmaxR = @(x) firmsize.*(cdf(pdR,x)).^(firmsize-1).*pdf(pdR,x);
%     fminL = @(x) firmsize.*(1-cdf(pdL,x)).^(firmsize-1).*pdf(pdL,x);
%     fmaxL = @(x) firmsize.*(cdf(pdL,x)).^(firmsize-1).*pdf(pdL,x);
%     % Numerically integrate
%     exminR = integral(@(x) x.*fminR(x),-1,1)-pdR.mean;
%     exmaxR = integral(@(x) x.*fmaxR(x),-1,1)-pdR.mean;
%     exminL = integral(@(x) x.*fminL(x),-1,1)-pdL.mean;
%     exmaxL = integral(@(x) x.*fmaxL(x),-1,1)-pdL.mean;   
    % save results
    meanRight(i,:) = [pdR.mean, pdR.var];
    meanLeft(i,:) = [pdL.mean, pdL.var];
    meanRightEmp(i,:) = [avgminR; avgmaxR];
    meanLeftEmp(i,:) = [avgminL; avgmaxL];
   
    meanRightEmpSD(i,:) = [avgminR-pdR.mean,; avgmaxR-pdR.mean];
    meanLeftEmpSD(i,:) = [avgminL-pdL.mean; avgmaxL-pdL.mean];
    
    meanVecOne(i,:)=[avgminL; avgmaxL;pdL.mean];
end

Vec=[flipud(meanLeftEmp);meanRightEmp];
Vec2=[flipud(meanLeftEmpSD);meanRightEmpSD];

Vec=meanVecOne;
Vec2=[meanLeftEmpSD];
Vec3=Vec(:,1)-Vec(:,2);

% Create figure
figure1 = figure('InvertHardcopy','off','PaperUnits','points',...
    'Color',[1 1 1]);
colormap(winter);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(Vec,'LineWidth',2);
set(plot1(1),'DisplayName','E(min)');
set(plot1(2),'DisplayName','E(max)');

% Create ylabel
ylabel({'Exp. Value'},'FontSize',11.9144327105334);

% Create xlabel
xlabel('Skew of Beta Distribution','FontSize',11.9144327105334);

% Create title
title({'Changes of E(max) and E(min) quality for skew levels'},...
    'FontSize',11.9144327105334);

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 160]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[-5 5]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(axes1,[-1 1]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',10.8313024641213,'XTick',[0 40 80 120 160],...
    'XTickLabel',{'10/2','5/2','2/2','2/5','2/10'});
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'FontSize',9.74817221770918);