
function fcg = convolve(f,g);

% sameconv.m from Dr. Alison I. Weber
% https://github.com/aiweber/GLM_and_Izhikevich/blob/master/sameconv.m

% convolution theorem: convolve f with g with inverse Fourier transform
% of the product of Fourier transform


    [fm,fn] = size(f);
    [gm,gn] = size(g);
    fcgn = fm+gm-1;
    
    fcg = ifft(sum(fft(f,fcgn).*fft(flipud(g),fcgn),2));
    
    fcg = fcg(1:fm,:);
end