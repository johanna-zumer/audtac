function tlock_fake = addbeforeFFT(data1,data2,trlkeep)

for at=1:length(trlkeep)
  cfg=[];
  cfg.trials=trlkeep(at);
  tmpt=ft_selectdata(cfg,data1);
  tmpa=ft_selectdata(cfg,data2);
  cfg=[];
  cfg.operation='add';
  cfg.parameter='trial';
  tmpsum=ft_math(cfg,tmpt,tmpa);
  if at==1
    tlock_fake=tmpsum;
  end
  tlock_fake.trial(at,:,:)=tmpsum.trial(1,:,:);
end
