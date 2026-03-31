%% MODULE 3
%
% Goal: Demonstrate how to fit basic and hierarchical d-prime models in 
% RStan. Also, to understand the advantages and disadvantages of using the 
% hierarchical method.
%
% Example: Ten participants perform a simplified version of the orientation 
% discrimination  task described in Module 1. This time, they only see one 
% of two stimuli on every trial (either counterclockwise or clockwise). 
% Prior to the main experiment, each observer completed a staircase 
% calibration procedure to equate stimulus difficulty to d?=1. There is 
% also no reason to suspect an observer will be biased in this task. 
%
% Created by SML Dec 2020

%% Simulate single unbiased observer with d'=1:

% fix random seed for reproducability:
rng(2020)

% true underlying SDT parameters:
dprime = 1; % sensitivity
k = 0; % unbiased criterion

% Task parameters (method of constant stimuli):
repeats = 50; % number of repeats per stimulus (CCW/CW tilt)
stimidx = repelem([1;2],repeats); % stimulus level on each trial (1=CCW, 2=CW)

% Simulate responses:
pCW_givenCCW = 1 - normcdf(k,-0.5*dprime,1); % probability of responding CW, when stim CCW
pCW_givenCW = 1 - normcdf(k,0.5*dprime,1); % probability of responding CW, when stim CW
pCW = stimidx; % fill temporarily with stimulus index
pCW(pCW==1) = pCW_givenCCW; % fill CCW trials
pCW(pCW==2) = pCW_givenCW; % fill CW trials
resp = rand(size(stimidx)) < pCW; % respond clockwise=1, counter-clockwise=0

% Compute and report the empirical values:
HR = mean(resp(stimidx==2)); % reporting CW when stimulus CW
FA = mean(resp(stimidx==1)); % reporting CW when stimulus CCW
dprime_emp = norminv(HR) - norminv(FA);
k_emp = -0.5*(norminv(HR) + norminv(FA));
disp(['The empirical d-prime is : ' num2str(dprime_emp) ' (true value is ' num2str(dprime), ')'])
disp(['The empirical criterion is : ' num2str(k_emp) ' (true value is ' num2str(k), ')'])

% Save simulated responses to text file:
T = table(stimidx,resp);
writetable(T,'simIndivObs_module3.txt','Delimiter',' ')

%% Simulate 10 observers in hierarchical fashion:

% true underlying SDT parameters:
nObs = 10; % number of observers
dprime_pop_mu = 1; % average sensitivity, set by staircasing procedure
dprime_pop_sd = 0.25; % sd in sensitivity, noise from observer or staircase measurement
k_pop_mu = 0; % average criterion, set to be unbiased
k_pop_sd = 0.1; % sd in criterion
dprime = dprime_pop_mu + dprime_pop_sd * randn([1,nObs]); % sample d'
k = k_pop_mu + k_pop_sd * randn([1,nObs]); % sample k

% Task parameters (method of constant stimuli):
repeats = 50; % number of repeats per stimulus (CCW/CW tilt)
stimidx = repelem([1;2],repeats); % stimulus level on each trial (1=CCW, 2=CW)

% Simulate responses:
pCW_givenCCW = 1 - normcdf(k,-0.5*dprime,1); % probability of responding CW, when stim CCW
pCW_givenCW = 1 - normcdf(k,0.5*dprime,1); % probability of responding CW, when stim CW
pCW = NaN([repeats*2,nObs]); % fill temporarily with stimulus index
for ii = 1:nObs
    pCW(stimidx==1,ii) = pCW_givenCCW(ii); % fill CCW trials
    pCW(stimidx==2,ii) = pCW_givenCW(ii); % fill CW trials
end
resp = rand(size(pCW)) < pCW; % respond clockwise=1, counter-clockwise=0

% Compute and report the empirical values:
HR = mean(resp(stimidx==2,:)); % reporting CW when stimulus CW
FA = mean(resp(stimidx==1,:)); % reporting CCW when stimulus CW
dprime_emp = norminv(HR) - norminv(FA);
k_emp = -0.5*(norminv(HR) + norminv(FA));
disp('Estimated versus true d-prime for each subject:')
disp([dprime_emp' dprime'])
disp('Estimated versus true criterion for each subject:')
disp([k_emp' k'])

% Save simulated responses to text file:
T = table(stimidx,resp);
writetable(T,'simMultObs_module3.txt','Delimiter',' ')

%% Visualise how the hierarchical model fit the subject data compared to the standard method and the basic fits:
% MAKE SURE TO RUN THIS AFTER FITTING THE MODELS IN R!!!

% Load the RStan fits:
fitResults = readtable('data_compareFits.txt','Delimiter',' ');

% 
figure; 
subplot(1,2,1); hold on
plot([0,2],[0,2],'k--','HandleVisibility','off')
plot(dprime,dprime_emp,'bo')
plot(dprime,fitResults.dBasic,'ro')
plot(dprime,fitResults.dHierarchical,'go')
legend({'Formula','Basic','Hierarchical'},'Location','northwest')
xlabel('True d-prime')
ylabel('Estimated d-prime')
title('D-Prime')
subplot(1,2,2); hold on
plot([-0.5,0.5],[-0.5,0.5],'k--','HandleVisibility','off')
plot(k,k_emp,'bo')
plot(k,fitResults.kBasic,'ro')
plot(k,fitResults.kHierarchical,'go')
legend({'Formula','Basic','Hierarchical'},'Location','northwest')
xlabel('True criterion')
ylabel('Estimated criterion')
xlim([-0.2,0.2])
title('Criterion')
set(findall(gcf,'-property','FontSize'),'FontSize',14)