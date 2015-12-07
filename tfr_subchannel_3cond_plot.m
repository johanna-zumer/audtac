function tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
% function tfr_subchannel_3cond_plot(tmp,tmpu,tmpm,chansel,zlim,figinds,figstrings)

if size(zlim)~=[4 2]
  error('wrong zlim size')
end
if length(figinds)~=2
  error('wrong size figinds')
end
if length(figstrings)~=2
  error('wrong size figstrings')
end

cfg=[];
cfg.avgoverchan='yes';
cfg.channel=chansel;

for pp=1:length(pow)
  tmp{pp}=ft_selectdata(cfg,pow{pp});
  tmp{pp}.mask=logical(ceil(tmp{pp}.mask));
end

% tmp1=ft_selectdata(cfg,tmp);
% tmp1.mask=logical(ceil(tmp1.mask));
% tmpu1=ft_selectdata(cfg,tmpu);
% tmpu1.mask=logical(ceil(tmpu1.mask));
% tmpm1=ft_selectdata(cfg,tmpm);
% tmpm1.mask=logical(ceil(tmpm1.mask));

cfg=[];
cfg.parameter='plvavgabs';
cfg.layout='elec1010.lay';
cfg.maskparameter='mask';
cfg.maskalpha=0.5;
cfg.xlim=[-.6 1.0];

figure(figinds(1));
for pp=1:length(pow)
  subplot(length(pow),1,pp);
  if pp<4
    cfg.zlim=zlim(1,:);
  elseif pp==4
    cfg.zlim=zlim(2,:);
  end
  ft_singleplotTFR(cfg,tmp{pp});
  axis([cfg.xlim(1) cfg.xlim(2) tmp{1}.freq(1) tmp{1}.freq(end)]);
end
% print(figinds(1),figstrings{1},'-dpng')
print(figinds(1),figstrings{1},'-depsc2')

% subplot(3,1,3);
% ft_singleplotTFR(cfg,tmp1);
% subplot(3,1,1);
% cfg.zlim=zlim(2,:);
% ft_singleplotTFR(cfg,tmpu1);
% subplot(3,1,2);
% ft_singleplotTFR(cfg,tmpm1);


figure(figinds(2));
cfg.parameter='plvavgang';
for pp=1:length(pow)
  subplot(length(pow),1,pp);
  if pp<4
    cfg.zlim=zlim(3,:);
  elseif pp==4
    cfg.zlim=zlim(4,:);
  end
  ft_singleplotTFR(cfg,tmp{pp});
  axis([cfg.xlim(1) cfg.xlim(2) tmp{1}.freq(1) tmp{1}.freq(end)]);
end
% print(figinds(2),figstrings{2},'-dpng')
print(figinds(2),figstrings{2},'-depsc2')

% subplot(3,1,3);
% cfg.zlim=zlim(3,:);
% ft_singleplotTFR(cfg,tmp1);
% subplot(3,1,1);
% cfg.zlim=zlim(4,:);
% ft_singleplotTFR(cfg,tmpu1);
% subplot(3,1,2);
% ft_singleplotTFR(cfg,tmpm1);




