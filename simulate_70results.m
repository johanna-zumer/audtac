% script to show possibility of different response profiles to two
% differnet stimulus conditions: same ERP result but opposing ITC

% In addition to needing FieldTrip installed on path, this also requires
% 'getplv' also available in the same Github directory.

%%
% sets up path
eeg_legomagic_preamble

%%

dataATpN_AT70=[];
dataApT_AT70=[];
dataATpN_TA70=[];
dataApT_TA70=[];

dataATpN_AT70.label{1}='chan';
dataApT_AT70.label{1}='chan';
dataATpN_TA70.label{1}='chan';
dataApT_TA70.label{1}='chan';

time=-1:.001:2;
for nn=1:100
  dataATpN_AT70.time{nn}=time;
  dataApT_AT70.time{nn}=time;
  dataATpN_TA70.time{nn}=time;
  dataApT_TA70.time{nn}=time;

  dataATpN_AT70.trial{nn}= [cos(4*2*pi*time(1:1000)+rand*2*pi), 1.1*cos(4*2*pi*time(1001:2000)+0.5*rand*2*pi), cos(4*2*pi*time(2001:3001)+rand*2*pi)];
  dataApT_AT70.trial{nn}=[cos(4*2*pi*time(1:1000)+rand*2*pi), 3.0*cos(4*2*pi*time(1001:2000)+0.6*rand*2*pi), cos(4*2*pi*time(2001:3001)+rand*2*pi)];

  dataATpN_TA70.trial{nn}= [cos(4*2*pi*time(1:1000)+rand*2*pi), 3.0*cos(4*2*pi*time(1001:2000)+0.6*rand*2*pi), cos(4*2*pi*time(2001:3001)+rand*2*pi)];
  dataApT_TA70.trial{nn}=[cos(4*2*pi*time(1:1000)+rand*2*pi), 5.0*cos(4*2*pi*time(1001:2000)+0.5*rand*2*pi), cos(4*2*pi*time(2001:3001)+rand*2*pi)];
end

cfg=[];
tlockATpN_AT70=ft_timelockanalysis(cfg,dataATpN_AT70);
tlockApT_AT70=ft_timelockanalysis(cfg,dataApT_AT70);
tlockATpN_TA70=ft_timelockanalysis(cfg,dataATpN_TA70);
tlockApT_TA70=ft_timelockanalysis(cfg,dataApT_TA70);


cfg=[];
cfg.method='mtmconvol';
cfg.foi=4;
cfg.taper='hanning';
cfg.t_ftimwin=0.5; % s
cfg.toi=-.5:.1:1.5;
cfg.output='fourier';
freqATpN_AT70=ft_freqanalysis(cfg,dataATpN_AT70);
freqApT_AT70=ft_freqanalysis(cfg,dataApT_AT70);
freqATpN_TA70=ft_freqanalysis(cfg,dataATpN_TA70);
freqApT_TA70=ft_freqanalysis(cfg,dataApT_TA70);

%%
figure(1);
subplot(3,2,1);plot(tlockATpN_AT70.time,tlockATpN_AT70.avg,'b');hold on;plot(tlockApT_AT70.time,tlockApT_AT70.avg,'g');ylim([-2 2])
ylabel('ERP')
title('AT70')
subplot(3,2,2);plot(tlockATpN_TA70.time,tlockATpN_TA70.avg,'b');hold on;plot(tlockApT_TA70.time,tlockApT_TA70.avg,'g');ylim([-2 2])
title('TA70')

subplot(3,2,3);plot(freqATpN_AT70.time,abs(squeeze(getplv(freqATpN_AT70.fourierspctrm))),'b');hold on;plot(freqApT_AT70.time,abs(squeeze(getplv(freqApT_AT70.fourierspctrm))),'g');ylim([0 1])
ylabel('ITC')
subplot(3,2,4);plot(freqATpN_AT70.time,abs(squeeze(getplv(freqATpN_TA70.fourierspctrm))),'b');hold on;plot(freqApT_AT70.time,abs(squeeze(getplv(freqApT_TA70.fourierspctrm))),'g');ylim([0 1])

subplot(3,2,5);plot(freqATpN_AT70.time,squeeze(mean(abs(freqATpN_AT70.fourierspctrm),1)),'b');hold on;plot(freqApT_AT70.time,squeeze(mean(abs(freqApT_AT70.fourierspctrm),1)),'g');ylim([0 3])
ylabel('TFP')
legend({'AT+N' 'A+T'})
subplot(3,2,6);plot(freqATpN_AT70.time,squeeze(mean(abs(freqATpN_TA70.fourierspctrm),1)),'b');hold on;plot(freqApT_AT70.time,squeeze(mean(abs(freqApT_TA70.fourierspctrm),1)),'g');ylim([0 3])

print(1,[fdir 'simulation70ms.png'],'-dpng')
print(1,[fdir 'simulation70ms.eps'],'-painters','-depsc')



