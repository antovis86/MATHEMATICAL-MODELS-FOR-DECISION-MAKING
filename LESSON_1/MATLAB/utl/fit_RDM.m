function [loglik] = fit_RDM(x,rt,resp)
% LIKELIHOOD FUNCTION for RDM

% Transform parameters to their native space
a = exp(x(1));
v1 = exp(x(2));
v0 = exp(x(3));

% Compute likelihood of responses of Accumulator 1 
P1 = log(utl_inverse_gaussian_defective(rt(resp)',v1,v0,a,a));
% Compute likelihood of responses of Accumulator 2
P2 = log(utl_inverse_gaussian_defective(rt(~resp)',v0,v1,a,a));

% Sum response likelihood
loglik=sum([P1 P2]);

end