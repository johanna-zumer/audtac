function tlock_fake = addbeforeFFT(data1,data2,trlkeep1,trlkeep2)

for at=1:length(trlkeep1)
  cfg=[];
  cfg.trials=trlkeep1(at);
  tmp1=ft_selectdata(cfg,data1);
  cfg.trials=trlkeep2(at);
  tmp2=ft_selectdata(cfg,data2);
  cfg=[];
  cfg.operation='add';
  cfg.parameter='trial';
  tmpsum=ft_math(cfg,tmp1,tmp2);
  if at==1
    tlock_fake=tmpsum;
  end
  tlock_fake.trial(at,:,:)=tmpsum.trial(1,:,:);
end
