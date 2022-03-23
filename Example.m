% For k filter only
% Gradient:  dl = (rate - spc)*Xstim


stim = [1 2 3 4 5]';
delta = [10 20 30 40 50]'; % delta = rate - spc
n = numel(stim);
m = 3;

%with design matrix (method used in Pillow's Lab tutorial)
ntfilt = m;
paddedStim = [zeros(ntfilt-1,1); stim];
Xstim = hankel(paddedStim(1:end-ntfilt+1), stim(end-ntfilt+1:end));
dldk3 = (delta'*Xstim)';

%with convolution product
conv = convolve(delta,stim);
dldk = conv(m:2*m-1);
