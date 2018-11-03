function ben=ConfBenefit_L(p_i,Gz,theta,delta,theta_i,v,PsiL, ConBen,IndBen)
PsiVec = PsiL(theta_i,theta);
%sumpi = sum(p_i);
% Calculate cost for a given p_i
%ben = p_i'*(((Gz*PsiVec)+(Gz*(IndBen.*ones(size(Gz,2),1))).^v)+ConBen);

% Delete g
%Gz=Gz>0;

ben = p_i'*((Gz*PsiVec)+ConBen);
end