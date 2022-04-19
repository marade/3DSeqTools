%% [lambda,Ihat] = makeMLL( I, du, u0, s)
% this function compute the minus-log-likelihood
function [lambda,Ihat] = makeMLL( I, du, u0, s)

l1 = sum(du.^2);

l2 = sum((I-u0).*du);

Ihat = l2/l1;

if isnan(Ihat)
     Ihat = 1;
 end
 
lambda = (l1*Ihat^2-2*Ihat*l2)/(2*s^2);

end