function tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline)
% function tfr_subchannel_3cond_plot_pow(tmp,tmpu,tmpm,chansel,zlim,figinds,figstrings,baseline)

% pow is struct with multiple indexed 'freq' data structures

if size(zlim)~=[1 2]
  error('wrong zlim size')
end
if all(size(baseline)==[1 2]) || strcmp(baseline,'no')
else
  error('wrong baseline size')
end
if length(figinds)~=1
  error('wrong size figinds')
end
if length(figstrings)~=1
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

figure(figinds(1));
cfg=[];
cfg.parameter='powspctrm';
cfg.layout='elec1010.lay';
cfg.maskparameter='mask';
cfg.maskalpha=0.5;
cfg.zlim=zlim;
cfg.xlim=[-.6 1.6];
cfg.baseline=baseline;
cfg.baselinetype='absolute';
for pp=1:length(pow)
  subplot(length(pow),1,pp);
%   subplot(1,length(pow),pp);
  ft_singleplotTFR(cfg,tmp{pp});
  axis([cfg.xlim(1) cfg.xlim(2) tmp{1}.freq(1) tmp{1}.freq(end)]);
end  
tmpfig=get(figinds(1));
for ff=2:2:(2*length(tmp)), set(tmpfig.Children(ff),'FontSize',18);end  % this controls the x,y, and colorbar fontsize
for ff=2:2:(2*length(tmp)), set(tmpfig.Children(ff).Title,'Visible','off'); end
% subplot(3,1,3);
% ft_singleplotTFR(cfg,tmp1);
% subplot(3,1,1);
% ft_singleplotTFR(cfg,tmpu1);
% subplot(3,1,2);
% ft_singleplotTFR(cfg,tmpm1);
% print(figinds(1),figstrings{1},'-dpng')
print(figinds(1),figstrings{1},'-depsc2')
