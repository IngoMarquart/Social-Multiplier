function [c,ceq]=nonLinConstraint(x,cutoff)
c = sum(x>0)-cutoff;
ceq = [];
end