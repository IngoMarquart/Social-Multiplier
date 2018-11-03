n=10000;
thetaDR=[0.4,1,0.5,15];
thetaDL=[0.4,1,0.5,9];
%thetaDR=[0.5,0.5,5.4375,7];
%thetaDL=[0.5,0.5,5.4375,7];
pd = makedist('Beta',thetaDR(1),thetaDR(2)); 
thetaRight = (random(pd,n,1)-pd.mean).*thetaDR(4)+thetaDR(3);
pd = makedist('Beta',thetaDL(1),thetaDL(2));
thetaLeft = (random(pd,n,1)-pd.mean).*thetaDL(4)+thetaDL(3);


% 
hist(thetaRight,100,'kernel')
hist(thetaLeft,100,'kernel')
hax=axes; 
histogram(thetaRight,100,'Normalization','pdf')
h = findobj(gca,'Type','patch')
set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)

hold on
histogram(thetaLeft,100,'Normalization','pdf')
h = findobj(gca,'Type','patch');
set(h,'facealpha',0.75);
% 
% 
line1=line([median(thetaLeft), median(thetaLeft)], get(hax,'YLim'), 'Color', 'r', 'LineWidth', 3);
line2=line([median(thetaRight), median(thetaRight)],get(hax,'YLim'), 'Color', 'b', 'LineWidth', 3);

colormap winter
hold off