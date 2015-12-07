function erpamp = erpamp_phasesorted_fft( freq, tlock, chan, phasesubtract)
% function erpamp = erpamp_phasesorted( freq, tlock, chan)

% each time point should be optimal for given freq
if size(freq.fourierspctrm,3)~=size(freq.fourierspctrm,4)
  error('freq should be computed so there is optimal time point for each freq bin')
end

erpamp=nan(4,size(freq.fourierspctrm,3),2);
if size(freq.fourierspctrm,1)>=40 % we want at least 10 trials per phase bin to make a sensible average
  for ff=1:size(freq.fourierspctrm,3)
    [angsort,angind]=sort(wrapToPi(angle(squeeze(freq.fourierspctrm(:,match_str(freq.label,chan),ff,ff)))-phasesubtract));
    phasebin{ff,1} =angind(angsort < -pi*3/4 | angsort > pi*3/4);
    phasebin{ff,2} =angind((angsort > -pi*3/4 & angsort < -pi*2/4)  | (angsort > pi*2/4 & angsort < pi*3/4) );
    phasebin{ff,3} =angind((angsort > -pi*2/4 & angsort < -pi*1/4)  | (angsort > pi*1/4 & angsort < pi*2/4) );
    phasebin{ff,4} =angind(angsort > -pi*1/4 & angsort < pi*1/4);    
    for bb=1:4
      erpamp(bb,ff,1)=median(tlock.trial(phasebin{ff,bb},match_str(tlock.label,chan),dsearchn(tlock.time',0.1)));
      erpamp(bb,ff,2)=median(tlock.trial(phasebin{ff,bb},match_str(tlock.label,chan),dsearchn(tlock.time',0.2)));
    end
  end
end
end
