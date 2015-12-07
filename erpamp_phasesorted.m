function erpamp = erpamp_phasesorted( freq, tlock, chan, timpt )
% function erpamp = erpamp_phasesorted( freq, tlock, chan, timpt )

[angsort,angind]=sort(angle(squeeze(freq.fourierspctrm(:,match_str(freq.label,chan),:,dsearchn(freq.time',timpt)))),1);
erpamp=nan(4,size(angsort,2),2);
if size(angsort,1)>=40 % we want at least 10 trials per phase bin to make a sensible average
  for ff=1:size(angsort,2)
    phasebin{ff,1} =angind(angsort(:,ff) < -pi*3/4 | angsort(:,ff) > pi*3/4,ff);
    phasebin{ff,2} =angind([angsort(:,ff) > -pi*3/4 & angsort(:,ff) < -pi*2/4]  | [angsort(:,ff) > pi*2/4 & angsort(:,ff) < pi*3/4],ff);
    phasebin{ff,3} =angind([angsort(:,ff) > -pi*2/4 & angsort(:,ff) < -pi*1/4]  | [angsort(:,ff) > pi*1/4 & angsort(:,ff) < pi*2/4],ff);
    phasebin{ff,4} =angind(angsort(:,ff) > -pi*1/4 & angsort(:,ff) < pi*1/4,ff);
    
    for bb=1:4
%       subplot(4,size(angsort,2),(bb-1)*size(angsort,2)+ff);
%       hist(tlock.trial(phasebin{ff,bb},match_str(tlock.label,chan),dsearchn(tlock.time',0.1)))
%       axis([-40 40 0 15])
      erpamp(bb,ff,1)=median(tlock.trial(phasebin{ff,bb},match_str(tlock.label,chan),dsearchn(tlock.time',0.1)));
      erpamp(bb,ff,2)=median(tlock.trial(phasebin{ff,bb},match_str(tlock.label,chan),dsearchn(tlock.time',0.2)));
    end
  end
end
end
