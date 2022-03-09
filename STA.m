function sta = STA(stimulus,spikes,klen)


% Part of 2010 Pillow Lab simpleSTC.m: 
% https://github.com/pillowlab/GLMspiketools/blob/master/glmtools_misc/simpleSTC.m

iisp = find(spikes);
nsp = sum(spikes);
Xpos = makeStimRows(stimulus, klen, iisp); 
%design matrix with only the stimulus rows for spC = 1
sta = (Xpos'*spikes(iisp))/nsp;

end