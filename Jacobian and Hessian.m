% For stimulus and k filter only


rate = [10 0 20 0 30]';
stim = [1 2 3 4 5]';
delta = [10 -1 20 1 30]'; % delta = rate - spc
n = numel(stim);
m = 3; % size of k filter


%% ----------------------------------------
% Gradient:  dl = (rate - spc)*Xstim
%with design matrix (method used in Pillow's Lab tutorial):
ntfilt = m;
paddedStim = [zeros(ntfilt-1,1); stim];
Xstim = hankel(paddedStim(1:end-ntfilt+1), stim(end-ntfilt+1:end));
dldk_X = (delta'*Xstim)';

%with convolution product:
conv = convolve(delta,stim);
dldk = flip(conv(n:n+m-1));

%% ---------------------------------------
%Hessian: -H = rate X^2
%with design matrix:
irate = rate'.*eye(n,n);
H_X = Xstim'*irate*Xstim;

%with convolution product:
Xrate = [];

for i=1:n
    select = irate(:,i);
    xrate = convolve(select,stim);
    xrate = flip(xrate(n:n+m-1));
    Xrate = [Xrate xrate];

end

H = [];

for i=1:m
    select2 = Xrate(i,:)';
    h = convolve(select2,stim);
    h = flip(h(n:n+m-1));
    H = [H  h];
end


function fcg = convolve(f,g)

    [fm,~] = size(f);
    [gm,~] = size(g);
    fcgn = fm+gm-1;
    
    fcg = ifft(sum(fft(f,fcgn).*fft(flipud(g),fcgn),2));
    
end

