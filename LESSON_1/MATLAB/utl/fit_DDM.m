function [loglik] = fit_DDM(x,rt,resp)
% LIKELIHOOD FUNCTION for RDM

% Transform parameters to their native space
a = exp(x(1));
v = exp(x(2));
w = 1/(1+exp(-x(3)));

% Compute likelihood of "upper" responses 
P1 = log(utl_wfpt(rt(resp)', -v, a, 1-w));
% Compute likelihood of "lower" responses 
P2 = log(utl_wfpt(rt(~resp)', v, a, w));

% Sum response likelihood
loglik=sum([P1 P2]);

end