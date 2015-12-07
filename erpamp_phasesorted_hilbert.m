function erpamp = erpamp_phasesorted_hilbert( tlock_hilbert, tlock, chan, timpt )
% function erpamp = erpamp_phasesorted( tlock_hilbert, tlock, chan, timpt )

erpamp=nan(4,length(tlock_hilbert),2);
for ff=1:length(tlock_hilbert)
  [angsort,angind]=sort(angle(tlock_hilbert{ff}.trial(:,match_str(tlock_hilbert{ff}.label,chan),dsearchn(tlock_hilbert{ff}.time',timpt))));
  if size(angsort,1)>=40 % we want at least 10 trials per phase bin to make a sensible average
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
