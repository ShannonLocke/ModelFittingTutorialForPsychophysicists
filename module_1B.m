%% MODULE 1B
%
% Goal: Demonstrate different model comparison methods.
%
% Example: Ten participants perform the task described in part A. Three 
% (30%) are unbiased and seven (70%) are biased. Some of the biased observers 
% are biased towards reporting ?clockwise? and others ?counter-clockwise?. 
% For simplicity, we have simulated the biased observers to all have a bias 
% magnitude of 1 deg. We perform the model comparison to understand if the 
% mu parameter is necessary and if a cumulative normal is a good psychometric 
% function for this task. 
%
% Created by SML Nov 2020

%% Simulate data:

% fix random seed for reproducability:
rng(20202020) 

% true underlying psychometric function of 10 observers:
mu_true = [0 0 0 1 1 1 1 -1 -1 -1]; % PSE
nObs = length(mu_true);
sigma_true = repelem(1.5,nObs) + 0.25 * randn([1,nObs]); % psychometric function slope, with Gaussian noise jitter
lambda = 0.05; % fixed lapse rate

% Task parameters (method of constant stimuli):
stimLvl = (-6:2:6)'; % stimulus orientations
nLvls = length(stimLvl);
repeats = 20; % number of repeats per stimulus level
nTrials = repeats * nLvls; % number of trials

% Simulate responses:
X = repmat(stimLvl,[repeats,nObs]); % stimulus level in each trial (not shuffled because only simulation)
pCW = normcdf(X,repmat(mu_true,[nTrials,1]),repmat(sigma_true,[nTrials,1])); % probability of clockwise jugdement
pCW = lambda/2 + (1 - lambda) * pCW; % apply lapse rate
resp = rand(size(X)) < pCW; % respond clockwise=1, counter-clockwise=0

% Visualise:
figure; hold on
xvals = -6:0.1:6;
nx = length(xvals);
xvals = repmat(xvals',[1,nObs]);
yvals = lambda/2 + (1 - lambda) * normcdf(xvals,repmat(mu_true,[nx,1]),repmat(sigma_true,[nx,1]));
plot(xvals,yvals,'LineWidth',2)
xlabel('Orientation (deg)'); xlim([-6,6])
ylabel('Proportion Clockwise'); ylim([0,1])
set(findall(gcf,'-property','FontSize'),'FontSize',14)

% Save simulated responses to mat and text files:
saveYN = true;
if saveYN
    sidx = repelem((1:nObs),nTrials)';
    T = table(sidx,X(:),resp(:));
    T.Properties.VariableNames = {'subject','stimLvl','resp'};
    writetable(T,'simObs_module1B.txt','Delimiter',' ')
    save('simObs_module1B.mat','X','resp')
end

%% Find MLEs with BADS:

% Preallocate model-fit vectors for speed:
MLE = NaN([nObs,2,4]);
best_nll = NaN([nObs,4]);

% Starting location and boundaries for search:
sloc = [0 2]; % starting location
lb = [-3 0]; % lower bound
ub = [3 3];  % upper bound
plb = [-2 0.5]; % plausible lower bound
pub = [2 2.5]; % plausible upper bound

% Fit models:
% 1) fixed at 0, 2) fixed mu at 1, 3) fixed mu at -1, 4) free mu
for sidx = 1:nObs
    getX = X(:,sidx); % get subject-specific stimlulus levels
    getResp = resp(:,sidx)'; % get subject-specific responses
    [MLE(sidx,2,1), best_nll(sidx,1)] = bads(@(x) model_nLL(0,x(1),lambda,getX,getResp), sloc(2),lb(2),ub(2),plb(2),pub(2));
    [MLE(sidx,2,2), best_nll(sidx,2)] = bads(@(x) model_nLL(1,x(1),lambda,getX,getResp), sloc(2),lb(2),ub(2),plb(2),pub(2));
    [MLE(sidx,2,3), best_nll(sidx,3)] = bads(@(x) model_nLL(-1,x(1),lambda,getX,getResp), sloc(2),lb(2),ub(2),plb(2),pub(2));
    [MLE(sidx,:,4), best_nll(sidx,4)] = bads(@(x) model_nLL(x(1),x(2),lambda,getX,getResp), sloc,lb,ub,plb,pub);
end

% Save model fits to mat file:
saveYN = true;
if saveYN; save('modelFits_module1B.mat','MLE','best_nll'); end

%% Find model evidence by brute-force grid method:

% Define grid:
d_mu = 0.1; % step size for mu
muVals = -3:d_mu:3; % grid points for mu
sigmaVals = logspace(log(0.5),log(2.5),40); % grid points for sigma
d_sig = diff(log(sigmaVals)); d_sig = d_sig(1); % step size for sigma in logspace
[MU,SIGMA] = meshgrid(muVals,sigmaVals); % create the grid. Use ndgrid if number of parameters >2.
n_gp1 = numel(sigmaVals); % % number of grid-points (fixed mu)
n_gp2 = numel(MU); % number of grid-points (free mu)

% Fit for each subject & calculate model evidence:
logEv = NaN([4,nObs]); % preallocate 
for sidx = 1:nObs
    % get subject-specific data:
    getX = X(:,sidx);
    getResp = resp(:,sidx)';
    % preallocate LL grids for speed
    NLL_1 = NaN(size(sigmaVals));
    NLL_2 = NaN(size(sigmaVals));
    NLL_3 = NaN(size(sigmaVals));
    NLL_4 = NaN(size(MU));
    % fit fixed-mu models (unbiased, biased-CW, biased-CCW):
    for ii = 1:n_gp1
        NLL_1(ii) = model_nLL(0,sigmaVals(ii),lambda,getX,getResp);
        NLL_2(ii) = model_nLL(1,sigmaVals(ii),lambda,getX,getResp);
        NLL_3(ii) = model_nLL(-1,sigmaVals(ii),lambda,getX,getResp);
    end
    % fit free-mu models:
    for ii = 1:n_gp2
        NLL_4(ii) = model_nLL(MU(ii),SIGMA(ii),lambda,getX,getResp);
    end
    % compute model evidence:
    logEv(1,sidx) = getModelEvidence(-NLL_1', d_sig);
    logEv(2,sidx) = getModelEvidence(-NLL_2', d_sig);
    logEv(3,sidx) = getModelEvidence(-NLL_3', d_sig);
    logEv(4,sidx) = getModelEvidence(-NLL_4, [d_mu, d_sig]);
end

% Save model evidence to mat file:
saveYN = true;
if saveYN; save('modelEvidence_module1B.mat','logEv'); end

%% Load simulation data and computed values:

% load('modelFits_module1B.mat','MLE','best_nll');
% load('modelEvidence_module1B.mat','logEv')

%% 7a) Bootstrapping non-parametric:

sel_idx = 1; % select an observer to model
nSims = 100; % set this much higher for real analysis (>1000, <10000)
simMLE = NaN([nSims,2]); % preallocate for speed
sloc = [0 2]; % starting location
for ii = 1:nSims
    trialidx = randi(nTrials,[nTrials,1]); % sample with replacement from trials
    getX = X(trialidx,sel_idx); % get sampled stimulus levels
    getResp = resp(trialidx,sel_idx)'; % get sampled responses
    simMLE(ii,:) = fminsearch(@(x) model_nLL(x(1),x(2),lambda,getX,getResp), sloc); % Find maximum likelihood estimates of new dataset
end

% Compute 95% confidence intervals:
mu_95CI = prctile(simMLE(:,1),[2.5 97.5]);
sigma_95CI = prctile(simMLE(:,2),[2.5 97.5]);

% Visualise:
figure
subplot(2,2,1); hold on
hist(simMLE(:,1)); % mu
plot([mu_true(sel_idx),mu_true(sel_idx)],[0,10],'r-','LineWidth',2)
xlabel('\mu (deg)')
subplot(2,2,3); hold on
plot(mu_95CI,[1 1],'b-','LineWidth',2) % mu
plot(mu_true(sel_idx),1,'r*')
xlabel('\mu (deg)')
xlim([floor(mu_95CI(1)), ceil(mu_95CI(2))])
subplot(2,2,2); hold on
hist(simMLE(:,2)); % sigma
plot([sigma_true(sel_idx),sigma_true(sel_idx)],[0,10],'r-','LineWidth',2)
xlabel('\sigma (deg)')
subplot(2,2,4); hold on
plot(sigma_95CI,[1 1],'b-','LineWidth',2) % mu
plot(sigma_true(sel_idx),1,'r*')
xlabel('\sigma (deg)')
xlim([floor(sigma_95CI(1)), ceil(sigma_95CI(2))])
set(findall(gcf,'-property','FontSize'),'FontSize',14)


%% 7b) Bootstrapping: parametric:

sel_idx = 1; % select an observer to model
sim_pCW = lambda/2 + (1 - lambda) * normcdf(X(:,sel_idx),MLE(sel_idx,1,4),MLE(sel_idx,2,4)); % probability of clockwise jugdement according to observer fit
nSims = 100; % set this much higher for real analysis (>1000, <10000)
simMLE = NaN([nSims,2]); % preallocate for speed
sloc = [0 2]; % starting location
for ii = 1:nSims
    getResp = rand([1,nTrials]) < sim_pCW'; % get sampled responses
    simMLE(ii,:) = fminsearch(@(x) model_nLL(x(1),x(2),lambda,X(:,sel_idx),getResp), sloc); % Find maximum likelihood estimates of new dataset
end

% Compute 95% confidence intervals:
mu_95CI = prctile(simMLE(:,1),[2.5 97.5]);
sigma_95CI = prctile(simMLE(:,2),[2.5 97.5]);

% Visualise:
figure
subplot(2,2,1); hold on
hist(simMLE(:,1)); % mu
plot([mu_true(sel_idx),mu_true(sel_idx)],[0,10],'r-','LineWidth',2)
xlabel('\mu (deg)')
subplot(2,2,3); hold on
plot(mu_95CI,[1 1],'b-','LineWidth',2) % mu
plot(mu_true(sel_idx),1,'r*')
xlabel('\mu (deg)')
xlim([floor(mu_95CI(1)), ceil(mu_95CI(2))])
subplot(2,2,2); hold on
hist(simMLE(:,2)); % sigma
plot([sigma_true(sel_idx),sigma_true(sel_idx)],[0,10],'r-','LineWidth',2)
xlabel('\sigma (deg)')
subplot(2,2,4); hold on
plot(sigma_95CI,[1 1],'b-','LineWidth',2) % mu
plot(sigma_true(sel_idx),1,'r*')
xlabel('\sigma (deg)')
xlim([floor(sigma_95CI(1)), ceil(sigma_95CI(2))])
set(findall(gcf,'-property','FontSize'),'FontSize',14)

%% 8a) AIC scores:

k = [1 1 1 2]; % number of free parameters in each model
AIC = 2*best_nll + repmat(2*k,[nObs,1]); % plain AIC score
correction = (2.*k.*(k+1))./(nTrials-k-1); % correction for sample size
AICc = AIC + repmat(correction,[nObs,1]); % corrected AIC score
relAICc = AICc - repmat(AICc(:,4),[1,4]); % corrected AIC score relative to any-bias model

% Visualise individual relative AIC scores:
figure;
subplot(1,2,1)
bar(relAICc(:,1:3))
xlabel('Subject')
ylabel('Relative AICc score (>0 = worse!)')
legend({'Unbiased','Biased-CW','Biased-CCW'})
set(findall(gcf,'-property','FontSize'),'FontSize',14)

%% 8b) AIC: Visualise group-average relative AIC scores:

meanVals = mean(relAICc);
semVals = std(relAICc)/sqrt(nObs);
subplot(1,2,2); hold on
bar(meanVals)
errorbar(1:4,meanVals',semVals','k.','HandleVisibility','off')
xlabel('Model')
ylabel('Mean Relative AICc Score (>0 = worse!)')
xticks(1:3)
xticklabels({'Unbiased','Biased-CW','Biased-CCW'})
xlim([0.5,3.5])
set(findall(gcf,'-property','FontSize'),'FontSize',14)

%% 9a) BIC scores:

k = [1 1 1 2]; % number of free parameters in each model
BIC = 2*best_nll + repmat(k*log(nTrials),[nObs,1]); % BIC score
relBIC = BIC - repmat(BIC(:,4),[1,4]); % BIC score relative to any-bias model

% Visualise individual relative BIC scores:
figure;
subplot(1,2,1)
bar(relBIC(:,1:3))
xlabel('Subject')
ylabel('Relative BIC score (>0 = worse!)')
legend({'Unbiased','Biased-CW','Biased-CCW'})
set(findall(gcf,'-property','FontSize'),'FontSize',14)

% Visualise group-average relative BIC scores:
meanVals = mean(relBIC);
semVals = std(relBIC)/sqrt(nObs);
subplot(1,2,2); hold on
bar(meanVals)
errorbar(1:4,meanVals',semVals','k.','HandleVisibility','off')
xlabel('Model')
ylabel('Mean Relative BICc Score (>0 = worse!)')
xticks(1:3)
xticklabels({'Unbiased','Biased-CW','Biased-CCW'})
xlim([0.5,3.5])
set(findall(gcf,'-property','FontSize'),'FontSize',14)

%% 9b) BIC and approx. Bayes Factor:

approxBF = exp(-0.5 * (BIC(:,1:3) - BIC(:,4))); % calculate the approximate Bayes Factor

% Visualise individual relative Bayes factors:
figure;
subplot(1,2,1)
bar(approxBF)
xlabel('Subject')
ylabel('Approximate Bayes Factor (>1 = better!)')
legend({'Unbiased','Biased-CW','Biased-CCW'})
set(findall(gcf,'-property','FontSize'),'FontSize',14)

% Visualise group-average relative Bayes factors:
meanVals = mean(approxBF);
semVals = std(approxBF)/sqrt(nObs);
subplot(1,2,2); hold on
bar(meanVals)
errorbar(1:3,meanVals',semVals','k.','HandleVisibility','off')
xlabel('Model')
ylabel('Mean Approximate Bayes Factors (>1 = better!)')
xticks(1:3)
xticklabels({'Unbiased','Biased-CW','Biased-CCW'})
xlim([0.5,3.5])
set(findall(gcf,'-property','FontSize'),'FontSize',14)

%% 10) Bayes Factor:
% Note to self: find out why this is sooooo much larger than the
% approximate values!!!

BF = exp(logEv');  % calculate Bayes factor using model evidence
BF = BF(:,1:3)./repmat(BF(:,4),[1,3]);

% Visualise individual relative Bayes factors:
figure;
subplot(1,2,1)
bar(BF)
xlabel('Subject')
ylabel('Approximate Bayes Factor (>1 = better!)')
legend({'Unbiased','Biased-CW','Biased-CCW'})
set(findall(gcf,'-property','FontSize'),'FontSize',14)

% Visualise group-average relative Bayes factors:
meanVals = mean(BF);
semVals = std(BF)/sqrt(nObs);
subplot(1,2,2); hold on
bar(meanVals)
errorbar(1:3,meanVals',semVals','k.','HandleVisibility','off')
xlabel('Model')
ylabel('Mean Approximate Bayes Factors (>1 = better!)')
xticks(1:3)
xticklabels({'Unbiased','Biased-CW','Biased-CCW'})
xlim([0.5,3.5])
set(findall(gcf,'-property','FontSize'),'FontSize',14)

%% 11a) Protected exceedance probabilities:

% Use the VBA package to get EPs and calc PEPs:
[posterior,out] = VBA_groupBMC(logEv);
PEP = (1-out.bor)*out.ep + out.bor/length(out.ep); 
disp(['The BOR is: ' num2str(out.bor)])

% Visualise model EPs and PEPs:
h = figure;
subplot(1,2,1)
bar([out.ep' PEP'])
title('Model Comparison Results')
xlabel('Model')
ylabel('Protected Exceedance Probability')
xticks(1:4)
modelName = {'Unbiased', 'Biased-CW', 'Biased-CCW', 'Any bias'};
xticklabels(modelName)
ylim([0,1])
xlim([0.5,4.5])
legend({'Exceedance Prob.','Protected Exceedance Prob.'})
set(findall(gcf,'-property','FontSize'),'FontSize',14)

%% 11b) Model partitioning:

% Use the VBA package to get EPs and calc PEPs:
options.families = {[1], [2,3,4]}; % unbiased model versus biased models
[posterior,out] = VBA_groupBMC(logEv, options);
PEP_partition = (1-out.bor)*out.families.ep + out.bor/length(out.families.ep); 

% Visualise model partition unbiased/biased:
figure(h)
subplot(1,2,2)
bar(PEP_partition)
title('Model Partition Results')
xlabel('Model Type')
ylabel('Protected Exceedance Probability')
xticks(1:2)
modelName = {'Unbiased','Biased'};
xticklabels(modelName)
ylim([0,1])
xlim([0.5,2.5])
set(findall(gcf,'-property','FontSize'),'FontSize',14)
