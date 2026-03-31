%% MODULE 1A
%
% Goal: Demonstrate different methods of optimization (fminsearch, fmincon, 
% BADS, and grid brute-force). 
%
% Example: Fitting the psychometric function of a biased observer in an 
% orientation discrimination task measured with the method of constant stimuli. 
%
% Created by SML Nov 2020

%% 1) Simulate data:

% fix random seed for reproducability:
rng(2020) 

% true underlying psychometric function:
mu_true = 1; % <=== biased PSE
sigma_true = 1.5; % psychometric function slope
lambda = 0.05; % fixed lapse rate

% Task parameters (method of constant stimuli):
stimLvl = (-6:2:6)'; % stimulus orientations
nLvls = length(stimLvl);
repeats = 20; % number of repeats per stimulus level

% Simulate responses:
X = repelem(stimLvl,repeats); % stimulus level in each trial (not shuffled because only simulation)
pCW = lambda/2 + (1 - lambda) * normcdf(X,mu_true,sigma_true); % probability of clockwise jugdement
resp = rand(size(X)) < pCW; % respond clockwise=1, counter-clockwise=0

% Save simulated responses to text file:
saveYN = false;
if saveYN
    T = table(X,resp);
    writetable(T,'simObs_module1A.txt','Delimiter',' ')
end

% Tally responses:
nCW = NaN([1,nLvls]); % preallocate for speed
for ii = 1:nLvls
   nCW(ii) = sum(resp(X==stimLvl(ii))); % number judged clockwise
end
propCW = nCW/repeats; % proportion judged clockwise
resp = resp'; % transpose matrix here for fast computation

% Visualise:
figure; hold on
xvals = -6:0.1:6;
yvals = lambda/2 + (1 - lambda) * normcdf(xvals,mu_true,sigma_true);
plot(xvals,yvals,'LineWidth',2)
plot(stimLvl,propCW,'ro','LineWidth',2)
xlabel('Orientation (deg)'); xlim([-6,6])
ylabel('Proportion Clockwise'); ylim([0,1])
set(findall(gcf,'-property','FontSize'),'FontSize',14)

% Preallocate model-fit vectors for speed:
MLE = NaN([5,2]);
best_nll = NaN([5,1]);
timeCheck = NaN([5,1]);

%% 2) Fit with fminsearch:

tic;
sloc = [0 2]; % starting location
[MLE(1,:), best_nll(1)] = fminsearch(@(x) model_nLL(x(1),x(2),lambda,X,resp), sloc);
timeCheck(1) = toc; % time to run

%% 3) Fit with fmincon:

tic;
sloc = [0 2]; % starting location
lb = [-3 0]; % lower bound
ub = [3 3]; % upper bound
[MLE(2,:), best_nll(2)] = fmincon(@(x) model_nLL(x(1),x(2),lambda,X,resp), sloc,[],[],[],[],lb,ub);
timeCheck(2) = toc; % time to run

%% 4) Fit with BADS:

tic
sloc = [0 2]; % starting location
lb = [-3 0]; % lower bound
ub = [3 3];  % upper bound
plb = [-2 0.5]; % plausible lower bound
pub = [2 2.5]; % plausible upper bound
[MLE(3,:), best_nll(3)] = bads(@(x) model_nLL(x(1),x(2),lambda,X,resp), sloc,lb,ub,plb,pub);
timeCheck(3) = toc; % time to run

%% 5) Fitting by simulation (BADS):
% Also see Matlab's function "patternsearch" for an alternative algorithm
% suitable for noisy objective functions.

tic
OPTIONS.UncertaintyHandling = 1; % Tell BADS that you are calculating the posterior from simulations
OPTIONS.NoiseFinalSamples = 100; % How many noise samples are collected and averaged for the final nLL value
sloc = [0 2]; % starting location
lb = [-3 0]; % lower bound
ub = [3 3]; % upper bound
plb = [-2 0.5]; % plausible lower bound
pub = [2 2.5]; % plausible upper bound
[MLE(4,:), best_nll(4)] = bads(@(x) model_nLL_bySim(x(1),x(2),lambda,X,resp), sloc,lb,ub,plb,pub,OPTIONS);
timeCheck(4) = toc; % time to run

%% 6) Fit with brute-force grid:

tic
muVals = -3:0.1:3; % grid points for mu
sigmaVals = logspace(log(0.5),log(2.5),30); % grid points for sigma
[MU,SIGMA] = meshgrid(muVals,sigmaVals); % create the grid. Use ndgrid if number of parameters >2.
NLL = NaN(size(MU)); % preallocate for speed
n_gp = numel(MU); % number of grid-points
for ii = 1:n_gp
    NLL(ii) = model_nLL(MU(ii),SIGMA(ii),lambda,X,resp);
end
idx = find(NLL==min(NLL(:))); % index of minimum
MLE(5,:) = [MU(idx), SIGMA(idx)]; % parameters at minimum
best_nll(5) = NLL(idx); % negative log-likelihood at minimum
timeCheck(5) = toc; % time to run

%% Compare the parameter estimates:

figure; hold on
surf(MU,SIGMA,-NLL)
shading interp
view(2)
colorbar
xlabel('\mu (deg)')
ylabel('\sigma (deg)')
ylim([sigmaVals(1) sigmaVals(end)])
plot(mu_true,sigma_true,'r*')
plot(MLE(:,1),MLE(:,2),'bo')
set(findall(gcf,'-property','FontSize'),'FontSize',14)

%% Compare negative log-likelihoods:

disp('Negative log-likelihoods: ')
disp(['fminsearch = ' num2str(best_nll(1))])
disp(['fmincon = ' num2str(best_nll(2))])
disp(['BADS = ' num2str(best_nll(3))])
disp(['simulation and BADS = ' num2str(best_nll(4))])
disp(['brute-force grid = ' num2str(best_nll(5))])

%% Compare timing:
% Compute time of gradient decent methods as a fraction of the grid method:

relTime = timeCheck(1:4)/timeCheck(5); % relative time to complete compared to grid brute-force
figure
bar(relTime)
xticklabels({'fminsearch','fmincon','BADS','sim. BADS'})
xlabel('Method')
ylabel('% of Time taken relative to grid method')
set(findall(gcf,'-property','FontSize'),'FontSize',14)
