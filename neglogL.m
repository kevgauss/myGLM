function [neglogL, dL, H] = neglogL(prs,struct);

    stim = struct.stim;
    spc = struct.spc;
    dt = struct.dt;
    nlfun = struct.nlfun;
    
    k = prs(:,1);
    h = prs(:,2);
    
    rLen = size(spc,1);
    
    xk = convolve(stim,k);

    if ~isempty(h)  % ~isempty = isnotempty
        yh = convolve(spc,h);
    end

    xTheta = xk + yh;
    rate = nlfun(xTheta);

    bmask = true(rLen,1);
    rrmask = rate(bmask); 
    trm1 = sum(rrmask)*dt; 
    trm2 = -sum(log(rrmask(find(spc(bmask))))); 
    neglogL = trm1 + trm2;
    
    