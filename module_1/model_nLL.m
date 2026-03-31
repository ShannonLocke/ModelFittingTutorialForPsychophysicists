function [nLL] = model_nLL(mu,sigma,lambda,X,resp)
% Computes the negative log-likelihood for a cumulative normal psychometric function.
% 
%
% Created by SML Nov 2020

p = lambda/2 + (1-lambda) * normcdf(X,mu,sigma); % probability respond clockwise according to proposed psychometric function
nLL = resp*log(p) + (1-resp)*log(1-p); % log of the bernoulli distribution with resp as row vector, p as column vector
nLL = -nLL; % actually make it the NEGATIVE log-likelihood...

end