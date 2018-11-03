

function cost=ConfCost(p_i,Gz,x_t_1,delta,theta_i, x_i)
%sumpi = sum(p_i);
% Vector of ones for each z
%ez = ones(size(Gz,1),1).*1/(1+delta.*sumpi);
ez = ones(size(Gz,2),1);
% Delete g
Gz=Gz>0;
% Calculate cost for a given p_i
cost = p_i'*(Gz*((ez.*x_i-x_t_1).*(ez.*x_i-x_t_1)));
%cost = p_i'*((-Gz*x_t_1+ez.*(p_i'*Gz*x_t_1.*delta+theta_i)).*(-Gz*x_t_1+ez.*(p_i'*Gz*x_t_1.*delta+theta_i)));
end