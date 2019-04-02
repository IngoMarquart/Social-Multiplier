%%%
% % RangeMeasure
% Calculates the range of connection of groups
% @param: theta - Theta vector
% @param: identity - Identity vector
% @param: A - Attention matrix
% @return: rangevec - consolidation value
%%%
function rangevec = RangeMeasure(theta,identity,A)
% Make sure we have column vectors
theta=theta(:);
identity=identity(:);
% Get number of actors
n=numel(identity);

% Get the average peer quality of k
for k =1:n
    prow=A(k,:);
    if sum(prow)~=0
        rdiff(k)=abs(theta(k)-(1/sum(prow)).*prow*theta);
    else
        rdiff(k)=0;
    end
end
rangevec.rdiffC=sum(rdiff(identity==1))./sum((identity==1));
rangevec.rdiffS=sum(rdiff(identity==0))./sum((identity==0));
rangevec.rdiffW=sum(rdiff(identity==-1))./sum((identity==-1));
rangevec.rdiffOverall=sum(rdiff)./n;
end