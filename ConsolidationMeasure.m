%%%
% % ConsolidationMeasure
% This function calculates the actual consolidation within the company.
% @param: theta - Theta vector
% @param: identity - Identity vector
% @return: samplecons - consolidation value
%%%
function samplecons = ConsolidationMeasure(theta,identity)
[~,sortindex]=sort(theta,'descend');
n=numel(theta);
same=(sortindex(:)==[1:n]');
sameNeg=(sortindex(:)==[n:-1:1]');
samplecons.cons=sum(same)./n;
samplecons.consNeg=-sum(sameNeg)./n;
if samplecons.cons >= abs(samplecons.consNeg)
    samplecons.fcons=samplecons.cons;
else
    samplecons.fcons=samplecons.consNeg;
end
end