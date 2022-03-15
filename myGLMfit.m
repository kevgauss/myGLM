tic
fprintf('Loading data: ');

stim = load('Stimulus.mat');
stim = stim.Stim;
spC = load('spcounts.mat');
spC = spC.sps;
th_k = load('kfilter.mat');
th_k = th_k.k;
th_h = load('hfilter.mat');
th_h = th_h.H;

toc
%% 1. Data preparation %%-------
tic
fprintf('Preparation: ');
fs = 30000;

nT = size(stim,1);

% tsp2 = tsp(tsp<nT);

N=1;
xbins = (N/2:N:nT);
% spC = hist(tsp2,xbins)';

NumberOfSpikes = numel(find(spC));

dt = 1/fs;
window = 1:nT;
twindow = window*dt;
Window = 1:round(nT/1);
tWindow = Window*dt;

stim = stim(Window)/10;
spC = spC(Window);

[sLen,sWid] = size(stim);
[rLen,rWid] = size(spC);


kLen = 2000;  % select the k filter length
toc

%% 2. Construct STA  %%--------
fprintf('------------------------\n');
fprintf('Building STA: ');
tic
iisp = find(spC);
nsp = sum(spC);
Xpos = makeStimRows(stim, kLen, iisp); 
sta = (Xpos'*spC(iisp))/nsp;
toc
subplot(211);
plot(sta);
title('STA');
subplot(212);
plot(th_k);
title('theoretical k filter');
%% 3. Initiate parameters  %%---------
h = rand([kLen,1]); % random weights for history filter

prs0.k = sta;
prs0.h  = h;
prs0.dc = 0; % direct current = baseline spikerate
prs0.dt = dt;
prs0.nlfun = @expf;
prs0.spc = spC;

%% initial -log L
fprintf('Init -log L: ');
tic

Struct.stim = stim;
Struct.spc = spC;
Struct.dt = dt;
Struct.nlfun = prs0.nlfun;

Prs0 = [prs0.k(:) prs0.h(:)];
neglogL0 = neglogL(Prs0,Struct);

toc
fprintf('------------------------\n');
fprintf('Initial negative log-likelihood: %.5f\n', neglogL0);

%% fitting prep


algopts = getFminOptsForVersion(version);
Opts = {'display', 'iter', 'maxiter', 100};
opts = optimset(algopts{:}, Opts{:});


lfun = @(prs)neglogL(prs,Struct);
%% minimization

prs1 = fminunc(lfun,Prs0,opts);  %gradient and Hessian are needed


function [f,df,ddf] = expf(x)

f = exp(x);
df = f;
ddf = df;
end