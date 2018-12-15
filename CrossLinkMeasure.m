%%%
% % Cross Linkage Measure
% This function calculates cross linkages between groups of types
% @param: theta - Theta vector
% @param: identity - Identity vector
% @param: A - Attention matrix
% @return: crosslinkage - vector of cross linkage percentages
%%%
function crosslinkage = CrossLinkMeasure(theta,identity,A)
	n=numel(identity);
	for k =1:n
    prow=A(k,:);
    if sum(prow)~=0
    linkid(k)=(1/sum(prow)).*prow*identity;
else
	linkid(k)=identity(k);
	end
	linkdiff(k)=abs(linkid(k)-identity(k));
end
	linkvec=zeros(n,1);
	linkvec=linkdiff>=0.5;
	

	crosslinkage.crosslinkageC=sum(linkvec(identity==1));
	crosslinkage.crosslinkageW=sum(linkvec(identity==0));
	crosslinkage.crosslinkageS=sum(linkvec(identity==-1));
	crosslinkage.overall=sum(linkvec)./n;
end