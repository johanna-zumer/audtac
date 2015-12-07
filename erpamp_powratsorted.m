function [erpamp, powrat] = erpamp_powratsorted( freq, tlock, erpchan, timwin )
% function [erpamp, powrat] = erpamp_powratsorted( freq, tlock, erpchan, timwin )
%
% freq:  providing power for theta and alpha

tt=3;

% get time points of freq power to use
freqtimwin=dsearchn(freq{tt,12}.time',timwin(1)):dsearchn(freq{tt,12}.time',timwin(2));

% get freq range of power to use
freqtheta=dsearchn(freq{tt,12}.freq',4):dsearchn(freq{tt,12}.freq',7);
freqalpha=dsearchn(freq{tt,12}.freq',9):dsearchn(freq{tt,12}.freq',12);

% average over time, frequency, and all channels to get power per trial per band
for ss=10:12
  %   abspow_theta{ss}=squeeze(nanmean(nanmean(abs(freq{tt,ss}.fourierspctrm(:,match_str(freq{tt,ss}.label,chan),freqtheta,freqtimwin)),4),3));
  %   abspow_alpha{ss}=squeeze(nanmean(nanmean(abs(freq{tt,ss}.fourierspctrm(:,match_str(freq{tt,ss}.label,chan),freqalpha,freqtimwin)),4),3));
  if ~isempty(freq{tt,ss})
    abspow_theta=squeeze(nanmean(nanmean(nanmean(abs(freq{tt,ss}.fourierspctrm(:,:,freqtheta,freqtimwin)),4),3),2));
    abspow_alpha=squeeze(nanmean(nanmean(nanmean(abs(freq{tt,ss}.fourierspctrm(:,:,freqalpha,freqtimwin)),4),3),2));
    powrat{ss}=abspow_theta./abspow_alpha;
  else
    powrat{ss}=[];
  end
end

try nt0=size(tlock{tt,10}.trial,1); catch nt0=0; end
nt1=size(tlock{tt,11}.trial,1);
nt2=size(tlock{tt,12}.trial,1);

% does combining stages make sense?  (depends on question being asked of data)
if ~isempty(powrat{10}) && ~isempty(powrat{11})
  powrat{24}=[powrat{10}; powrat{11}];
else
  powrat{24}=[];
end
if ~isempty(powrat{11}) && ~isempty(powrat{12})
  powrat{25}=[powrat{11}; powrat{12}];
else
  powrat{25}=[];
end
if ~isempty(powrat{10}) && ~isempty(powrat{11}) && ~isempty(powrat{12})
  powrat{26}=[powrat{10}; powrat{11}; powrat{12}];
else
  powrat{26}=[];
end

erpamp=nan(26,3,2);
for ss=[10 11 12 24 25 26]
  if ~isempty(powrat{ss}) && length(powrat{ss})>=30
    [powsort{ss},powind{ss}]=sort(powrat{ss});
    powbot=powind{ss}(1:floor(length(powsort{ss})/3));
    powmid=powind{ss}(ceil(length(powsort{ss})/3):floor(2*length(powsort{ss})/3));
    powtop=powind{ss}(ceil(2*length(powsort{ss})/3):end);
    
%     figure;subplot(3,1,1);hist(powbot);subplot(3,1,2);hist(powmid);subplot(3,1,3);hist(powtop)
    
    if ss<13
      erpamp(ss,1,1)=median(tlock{tt,ss}.trial(powbot,match_str(tlock{tt,ss}.label,erpchan),dsearchn(tlock{tt,ss}.time',0.1)));
      erpamp(ss,2,1)=median(tlock{tt,ss}.trial(powmid,match_str(tlock{tt,ss}.label,erpchan),dsearchn(tlock{tt,ss}.time',0.1)));
      erpamp(ss,3,1)=median(tlock{tt,ss}.trial(powtop,match_str(tlock{tt,ss}.label,erpchan),dsearchn(tlock{tt,ss}.time',0.1)));
      erpamp(ss,1,2)=median(tlock{tt,ss}.trial(powbot,match_str(tlock{tt,ss}.label,erpchan),dsearchn(tlock{tt,ss}.time',0.2)));
      erpamp(ss,2,2)=median(tlock{tt,ss}.trial(powmid,match_str(tlock{tt,ss}.label,erpchan),dsearchn(tlock{tt,ss}.time',0.2)));
      erpamp(ss,3,2)=median(tlock{tt,ss}.trial(powtop,match_str(tlock{tt,ss}.label,erpchan),dsearchn(tlock{tt,ss}.time',0.2)));
    elseif ss==24
      erpamp(ss,1,1)=median([tlock{tt,10}.trial(powbot(powbot<=nt0),match_str(tlock{tt,10}.label,erpchan),dsearchn(tlock{tt,10}.time',0.1)); tlock{tt,11}.trial(powbot(powbot>nt0)-nt0,match_str(tlock{tt,11}.label,erpchan),dsearchn(tlock{tt,11}.time',0.1)) ]);
      erpamp(ss,2,1)=median([tlock{tt,10}.trial(powmid(powmid<=nt0),match_str(tlock{tt,10}.label,erpchan),dsearchn(tlock{tt,10}.time',0.1)); tlock{tt,11}.trial(powmid(powmid>nt0)-nt0,match_str(tlock{tt,11}.label,erpchan),dsearchn(tlock{tt,11}.time',0.1)) ]);
      erpamp(ss,3,1)=median([tlock{tt,10}.trial(powtop(powtop<=nt0),match_str(tlock{tt,10}.label,erpchan),dsearchn(tlock{tt,10}.time',0.1)); tlock{tt,11}.trial(powtop(powtop>nt0)-nt0,match_str(tlock{tt,11}.label,erpchan),dsearchn(tlock{tt,11}.time',0.1)) ]);
      erpamp(ss,1,2)=median([tlock{tt,10}.trial(powbot(powbot<=nt0),match_str(tlock{tt,10}.label,erpchan),dsearchn(tlock{tt,10}.time',0.2)); tlock{tt,11}.trial(powbot(powbot>nt0)-nt0,match_str(tlock{tt,11}.label,erpchan),dsearchn(tlock{tt,11}.time',0.2)) ]);
      erpamp(ss,2,2)=median([tlock{tt,10}.trial(powmid(powmid<=nt0),match_str(tlock{tt,10}.label,erpchan),dsearchn(tlock{tt,10}.time',0.2)); tlock{tt,11}.trial(powmid(powmid>nt0)-nt0,match_str(tlock{tt,11}.label,erpchan),dsearchn(tlock{tt,11}.time',0.2)) ]);
      erpamp(ss,3,2)=median([tlock{tt,10}.trial(powtop(powtop<=nt0),match_str(tlock{tt,10}.label,erpchan),dsearchn(tlock{tt,10}.time',0.2)); tlock{tt,11}.trial(powtop(powtop>nt0)-nt0,match_str(tlock{tt,11}.label,erpchan),dsearchn(tlock{tt,11}.time',0.2)) ]);
    elseif ss==25
      erpamp(ss,1,1)=median([tlock{tt,11}.trial(powbot(powbot<=nt1),match_str(tlock{tt,11}.label,erpchan),dsearchn(tlock{tt,11}.time',0.1)); tlock{tt,12}.trial(powbot(powbot>nt1)-nt1,match_str(tlock{tt,12}.label,erpchan),dsearchn(tlock{tt,12}.time',0.1)) ]);
      erpamp(ss,2,1)=median([tlock{tt,11}.trial(powmid(powmid<=nt1),match_str(tlock{tt,11}.label,erpchan),dsearchn(tlock{tt,11}.time',0.1)); tlock{tt,12}.trial(powmid(powmid>nt1)-nt1,match_str(tlock{tt,12}.label,erpchan),dsearchn(tlock{tt,12}.time',0.1)) ]);
      erpamp(ss,3,1)=median([tlock{tt,11}.trial(powtop(powtop<=nt1),match_str(tlock{tt,11}.label,erpchan),dsearchn(tlock{tt,11}.time',0.1)); tlock{tt,12}.trial(powtop(powtop>nt1)-nt1,match_str(tlock{tt,12}.label,erpchan),dsearchn(tlock{tt,12}.time',0.1)) ]);
      erpamp(ss,1,2)=median([tlock{tt,11}.trial(powbot(powbot<=nt1),match_str(tlock{tt,11}.label,erpchan),dsearchn(tlock{tt,11}.time',0.2)); tlock{tt,12}.trial(powbot(powbot>nt1)-nt1,match_str(tlock{tt,12}.label,erpchan),dsearchn(tlock{tt,12}.time',0.2)) ]);
      erpamp(ss,2,2)=median([tlock{tt,11}.trial(powmid(powmid<=nt1),match_str(tlock{tt,11}.label,erpchan),dsearchn(tlock{tt,11}.time',0.2)); tlock{tt,12}.trial(powmid(powmid>nt1)-nt1,match_str(tlock{tt,12}.label,erpchan),dsearchn(tlock{tt,12}.time',0.2)) ]);
      erpamp(ss,3,2)=median([tlock{tt,11}.trial(powtop(powtop<=nt1),match_str(tlock{tt,11}.label,erpchan),dsearchn(tlock{tt,11}.time',0.2)); tlock{tt,12}.trial(powtop(powtop>nt1)-nt1,match_str(tlock{tt,12}.label,erpchan),dsearchn(tlock{tt,12}.time',0.2)) ]);
    elseif ss==26
      erpamp(ss,1,1)=median([tlock{tt,10}.trial(powbot(powbot<=nt0),match_str(tlock{tt,10}.label,erpchan),dsearchn(tlock{tt,10}.time',0.1)); tlock{tt,11}.trial(powbot(powbot>nt0 & powbot<=nt1+nt0)-nt0,match_str(tlock{tt,11}.label,erpchan),dsearchn(tlock{tt,11}.time',0.1)); tlock{tt,12}.trial(powbot(powbot>nt1+nt0)-nt1-nt0,match_str(tlock{tt,12}.label,erpchan),dsearchn(tlock{tt,12}.time',0.1)) ]);
      erpamp(ss,2,1)=median([tlock{tt,10}.trial(powmid(powmid<=nt0),match_str(tlock{tt,10}.label,erpchan),dsearchn(tlock{tt,10}.time',0.1)); tlock{tt,11}.trial(powmid(powmid>nt0 & powmid<=nt1+nt0)-nt0,match_str(tlock{tt,11}.label,erpchan),dsearchn(tlock{tt,11}.time',0.1)); tlock{tt,12}.trial(powmid(powmid>nt1+nt0)-nt1-nt0,match_str(tlock{tt,12}.label,erpchan),dsearchn(tlock{tt,12}.time',0.1)) ]);
      erpamp(ss,3,1)=median([tlock{tt,10}.trial(powtop(powtop<=nt0),match_str(tlock{tt,10}.label,erpchan),dsearchn(tlock{tt,10}.time',0.1)); tlock{tt,11}.trial(powtop(powtop>nt0 & powtop<=nt1+nt0)-nt0,match_str(tlock{tt,11}.label,erpchan),dsearchn(tlock{tt,11}.time',0.1)); tlock{tt,12}.trial(powtop(powtop>nt1+nt0)-nt1-nt0,match_str(tlock{tt,12}.label,erpchan),dsearchn(tlock{tt,12}.time',0.1)) ]);
      erpamp(ss,1,2)=median([tlock{tt,10}.trial(powbot(powbot<=nt0),match_str(tlock{tt,10}.label,erpchan),dsearchn(tlock{tt,10}.time',0.2)); tlock{tt,11}.trial(powbot(powbot>nt0 & powbot<=nt1+nt0)-nt0,match_str(tlock{tt,11}.label,erpchan),dsearchn(tlock{tt,11}.time',0.2)); tlock{tt,12}.trial(powbot(powbot>nt1+nt0)-nt1-nt0,match_str(tlock{tt,12}.label,erpchan),dsearchn(tlock{tt,12}.time',0.2)) ]);
      erpamp(ss,2,2)=median([tlock{tt,10}.trial(powmid(powmid<=nt0),match_str(tlock{tt,10}.label,erpchan),dsearchn(tlock{tt,10}.time',0.2)); tlock{tt,11}.trial(powmid(powmid>nt0 & powmid<=nt1+nt0)-nt0,match_str(tlock{tt,11}.label,erpchan),dsearchn(tlock{tt,11}.time',0.2)); tlock{tt,12}.trial(powmid(powmid>nt1+nt0)-nt1-nt0,match_str(tlock{tt,12}.label,erpchan),dsearchn(tlock{tt,12}.time',0.2)) ]);
      erpamp(ss,3,2)=median([tlock{tt,10}.trial(powtop(powtop<=nt0),match_str(tlock{tt,10}.label,erpchan),dsearchn(tlock{tt,10}.time',0.2)); tlock{tt,11}.trial(powtop(powtop>nt0 & powtop<=nt1+nt0)-nt0,match_str(tlock{tt,11}.label,erpchan),dsearchn(tlock{tt,11}.time',0.2)); tlock{tt,12}.trial(powtop(powtop>nt1+nt0)-nt1-nt0,match_str(tlock{tt,12}.label,erpchan),dsearchn(tlock{tt,12}.time',0.2)) ]);
    end
    
  end
end
