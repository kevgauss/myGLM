%% preparation
%set parameters
fprintf('------------------------------------- \n');
fprintf('PREPARATIONS \n');
nT = 100000;
fprintf('Stimulus length: %d\n',nT);
% N=1;
% tbins = (N/2:N:nT);
swid = 1;
fs = 3e+04;
fprintf('Sampling frequency fs: %d\n',fs);
dt = 1/fs;
dtStim = dt;
dtSp = dt;
nt = nT*dt*1000;
fprintf('Time window: %.5f',nt);
fprintf('ms \n');


%% create noisy stimulus
fprintf('________________________________________ \n');
fprintf('Creating noisy stimulus vector \n');
Stim = rand(nT,swid)*2-1;  


%% construct a stimulus filter
fprintf('________________________________________ \n');
coeff = 2;
nkt = 1000*coeff;
fprintf('Simulus filter length nkt: %d\n',nkt);
T = 1:nkt;
Tau = 60*coeff; %ms
Tau2 = 70*coeff;
Tau3 = 100*coeff; %must be > Tau,Tau2, increase to have a deeper and longer refractory period
Tau4 = 30*coeff;

fprintf('Constructing the theoretical stimulus filter k \n');
filter = s(T,Tau,Tau2,Tau3,Tau4);
gg.k=[];
kcoef = 4; %will increase the probability of spiking
gg.k = flip(filter')*kcoef;

%% construct the design matrix
ntfilt = nkt;
fprintf('________________________________________ \n');
fprintf('Constructing the Stimulus Design matrix Xs \n');
paddedStim = [zeros(ntfilt-1,1); Stim];
Xstim = hankel(paddedStim(1:end-ntfilt+1), Stim(end-ntfilt+1:end));
% Xdsgn = [Xstim, Xsp];

%% construct theoretical history filter

Xmax=10;
X = 1:(Xmax-1)/nkt:Xmax-(Xmax-1)/nkt;
c = 1/10;
cc = 2;
cx = 0;
pc = 2;
px = 2;
py = 3;

H = h(X,c,cc,cx,pc,px,py)';
dx = nkt/Xmax;
xplot = X*dx;



%% compute rate
fprintf('Computing firing rate = exp(k*Xs) \n');
% rate = exp(Xstim*gg1.k);
rate = Xstim*gg.k;
k = gg.k;

norm_k = k/norm(k);

windowstart = 1;
windowend = 100000;
window = windowstart:windowend;

rate1 = rate(window);

% norm_rate1 = rate1/norm(rate1);
% check = sum(norm_rate1);

deltaT = window*dt*1000; %in ms

%% simulation
fprintf('------------------------------------- \n');
fprintf('PREDICTIONS \n');
fprintf('prediction window: %d',windowstart);
fprintf(' : %d\n',windowend);

reps = 10;
fprintf('Number of repetitions: %d\n',reps);
sim = exprnd(10,[reps,numel(window)]);

Sim = zeros(reps,numel(window));
for i=1:reps
   for j=window
      if sim(i,j) <= rate1(j)
         Sim(i,j)=1; 
      end
   end
end

%% selecting sps in sim
% 
repchoice = 3;
sps1 = Sim(repchoice,:)';
NumberOfSpikes = numel(sps1(sps1>0));

%% history design matrix and rate
fprintf('________________________________________ \n');
fprintf('Constructing the History Design matrix Xh \n');


nthist = nkt;
paddedSts = [zeros(nthist,1); sps1(1:end-1)];
Xsp = hankel(paddedSts(1:end-nthist+1), paddedSts(end-nthist+1:end));

%% computing rate with history filter
fprintf('Computing firing rate = exp(k*Xs+h*Xh) \n');
Xstim2 = Xstim(window,:);
rate2 = Xstim2*gg.k+Xsp*H;
lincomb = Xsp*H;
fprintf('------------------------------------- \n');
fprintf('PREDICTIONS 2nd round\n');
fprintf('prediction window: %d',windowstart);
fprintf(' : %d\n',windowend);

%% simulation 2nd round
reps = 10;
fprintf('Number of repetitions: %d\n',reps);
sim2 = exprnd(10,[reps,numel(window)]);

Sim2 = zeros(reps,numel(window));
for i=1:reps
   for j=window
      if sim2(i,j) <= rate2(j)
         Sim2(i,j)=1; 
      end
   end
end

%% number of spike
repchoice2 = repchoice;
sps2 = Sim2(repchoice2,:)';
NumberOfSpikes2 = numel(sps2(sps2>0));

%% test GLM fit
fprintf('---------------------------------------- \n');
fprintf('FITTING GLM \n');
fprintf('Repetition choice: %d\n',repchoice);
fprintf('Number of spikes: %.5f\n', NumberOfSpikes);
fprintf('_______________________________ \n');

sps = sps2;
Stimulus = Stim(window);
sta = simpleSTC(Stimulus,sps,nkt); 
sta = reshape(sta,nkt,[]); 

exptmask= [];  

nkbasis = 5;  
nhbasis = 5;  
hpeakFinal = .2;   
gg0 = makeFittingStruct_GLM(dtStim,dtSp,nkt,nkbasis,sta,nhbasis,hpeakFinal);
gg0.sps = sps; 
gg0.mask = exptmask; 
gg0.ihw = randn(size(gg0.ihw))*1;  
[negloglival0,rr] = neglogli_GLM(gg0,Stimulus);
opts = {'display', 'iter', 'maxiter', 200};
[gg1, negloglival] = MLfit_GLM(gg0,Stimulus,opts);
fprintf('Number of spikes: %d\n', NumberOfSpikes);
fprintf('Initial negative log-likelihood: %.5f\n', negloglival0);
fprintf('Final negative log-likelihood: %.5f\n',negloglival);   
gg1k = gg1.k;
norm_gg1k  = gg1k/norm(gg1k);
gg1h = gg1.ih(1:nkt);
norm_gg1h  = gg1h/norm(gg1h);
% 
% %% ISI histogram
% 
% t = window*dt*1000;
% tsp = [];
% for i = window
%     if sps(i) == 1
%         Tsp = t(i);
%         tsp = [tsp Tsp];
%     end
% end
% 
% ISI = diff(tsp);

%% savings

% save('simulation1.mat','Sim');
% save('Stimulus.mat','Stim');
% % save('gg1k1.mat','gg1k');
% save('kfilter.mat','k');
% save('hfilter.mat','H');
% save('spcounts.mat','sps');


%% plotting

% plotwindow = 1:1000;
plotwindow = window;
twindow = window*dt*1000;

iiplot = 1:nkt;
iplot = -nkt+1:0;
tplot = iplot*dt*1000;

% 
% plot(1:90000,exp(lincomb)/100,1:90000,exp(rate1)/10000,1:90000,exp(rate2)/100000);
% legend('hist rate exp(h*Xh)','stim rate exp(k*Xs)','full rate exp(k*Xs+h*Xh)');

plot(iiplot,flip(norm_k),iiplot,flip(norm_gg1k),iiplot,H,iiplot,norm_gg1h,iiplot,H*0,'k--');
legend('k_{th}','k','h_{th}','h');

% figure;
% subplot(211);
% plot(iiplot,flip(gg.k),iiplot,H,iiplot,H*0,'k--');
% subplot(212);
% % plot(plotwindow,exp(rate(plotwindow)));
% plot(plotwindow,exp(rate(plotwindow)),plotwindow,exp(rate2(plotwindow)));
% title('firing rate');
% 
% figure;
% subplot(711);
% plot(plotwindow,exp(rate(plotwindow)));
% title('firing rate');
% subplot(7,1,2:7);
% imagesc(window,1:reps,Sim2);
% title('simulation');
% ylabel('Trials');

%%
% figure;
% subplot(8,1,1:2);
% plot(twindow,sps,'k',twindow,norm_rate(plotwindow)*100,'r');
% title(['Firing rate and spikes: ',num2str(NumberOfSpikes), ' spikes']);
% % xlabel('t (ms)');
% subplot(8,1,3:4);
% histogram(ISI,100);
% set(gca,'Yscale','log');
% title('ISI');
% subplot(8,1,5:8);
% plot(tplot,norm_k,tplot,norm_gg1k,tplot,k*0,'k--');
% legend('normalized theoretical filter','normalized GLM fitted filter');
% title('Stimulus filter reconstruction');
% xlabel('t (ms)');

% figure;
% histogram(ISI,100);
% set(gca,'Yscale','log');

% figure;
% plot(norm_k);

function[hist] = h(x,coef,ccoef,cxshift,pcoef,pxshift,pyshift)
    hist = coef*(ccoef./(cosh(x-cxshift)).*(pcoef.*(x-pxshift).^2-pyshift));

end

function[alpha]=s(t,tau1,tau2,tau3,tau4)
    alpha = exp(-t/tau1)+exp(-t/tau2)-exp(-t/tau3)-exp(-t/tau4);
end
