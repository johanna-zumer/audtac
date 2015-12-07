function tlock=trim_nans(tlock)

begsample=1;
endsample=size(tlock.trial,3);
zerosample=dsearchn(tlock.time',0);
for rr=1:size(tlock.trial,1)
  tmp=find(isnan(squeeze(tlock.trial(rr,1,:))));
  if ~isempty(tmp)
    tmpr=max(tmp(tmp<zerosample));
    tmpo=min(tmp(tmp>zerosample));
    if begsample<tmpr
      begsample=tmpr;
    end
    if endsample>tmpo
      endsample=tmpo;
    end
  end
end
cfg=[];
cfg.latency=[tlock.time(begsample+1) tlock.time(endsample-1)];
tlock=ft_selectdata(cfg,tlock);
