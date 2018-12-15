%%%
% % ConsolidationMeasure
% This function calculates the actual consolidation within the company.
% @param: theta - Theta vector
% @param: identity - Identity vector
% @return: samplecons - consolidation value
%%%
function samplecons = ConsolidationMeasure(theta,identity)
[~,sortindex]=sort(theta);
n=numel(theta);
same=round(sortindex(:)==[1:n]');
samplecons=sum(same)./n;
end