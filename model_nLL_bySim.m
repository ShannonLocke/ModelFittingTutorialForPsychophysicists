function [nLL] = model_nLL_bySim(mu,sigma,lambda,X,resp)
% Computes the negative log-likelihood for a cumulative normal psychometric 
% function by simulating an observer. DO NOT use this for actual fitting,
% this is just for teaching purposes.
%
% Created by SML Nov 2020

% Simulate responses to get p(clockwise):
nSim = 100; % number of simulated trials per acutal trial
pSim = lambda/2 + (1-lambda) * normcdf(X,mu,sigma); 
pSim = repmat(pSim,[1,nSim]);
respSim = rand(size(pSim)) < pSim; % simulated responses

% Tally the responses to determine the proportion of the time the simulated
% observer selected clockwise, and use this to estimate the probability our
% observer selected clockwise for this parameter pair:
p = lambda/2 + (1-lambda) * mean(respSim,2); 

% Compute negative log-likelihood:
nLL = resp*log(p) + (1-resp)*log(1-p); % <== make sure resp is row vector, and p is column vector
nLL = -nLL;

end