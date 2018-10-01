

dataMS_AT70=[];
dataApT_AT70=[];
dataMS_TA70=[];
dataApT_TA70=[];

dataMS_AT70.label{1}='chan';
dataApT_AT70.label{1}='chan';
dataMS_TA70.label{1}='chan';
dataApT_TA70.label{1}='chan';

time=-1:.001:2;
for nn=1:100
  dataMS_AT70.time{nn}=time;
  dataApT_AT70.time{nn}=time;
  dataMS_TA70.time{nn}=time;
  dataApT_TA70.time{nn}=time;

  dataMS_AT70.trial{nn}= [cos(4*2*pi*time(1:1000)+rand*2*pi), 1.1*cos(4*2*pi*time(1001:2000)+0.5*rand*2*pi), cos(4*2*pi*time(2001:3001)+rand*2*pi)];
  dataApT_AT70.trial{nn}=[cos(4*2*pi*time(1:1000)+rand*2*pi), 3.0*cos(4*2*pi*time(1001:2000)+0.6*rand*2*pi), cos(4*2*pi*time(2001:3001)+rand*2*pi)];

  dataMS_TA70.trial{nn}= [cos(4*2*pi*time(1:1000)+rand*2*pi), 3.0*cos(4*2*pi*time(1001:2000)+0.6*rand*2*pi), cos(4*2*pi*time(2001:3001)+rand*2*pi)];
  dataApT_TA70.trial{nn}=[cos(4*2*pi*time(1:1000)+rand*2*pi), 5.0*cos(4*2*pi*time(1001:2000)+0.5*rand*2*pi), cos(4*2*pi*time(2001:3001)+rand*2*pi)];
end

cfg=[];
tlockMS_AT70=ft_timelockanalysis(cfg,dataMS_AT70);
tlockApT_AT70=ft_timelockanalysis(cfg,dataApT_AT70);
tlockMS_TA70=ft_timelockanalysis(cfg,dataMS_TA70);
tlockApT_TA70=ft_timelockanalysis(cfg,dataApT_TA70);


cfg=[];
cfg.method='mtmconvol';
cfg.foi=4;
cfg.taper='hanning';
cfg.t_ftimwin=0.5; % s
cfg.toi=-.5:.1:1.5;
cfg.output='fourier';
freqMS_AT70=ft_freqanalysis(cfg,dataMS_AT70);
freqApT_AT70=ft_freqanalysis(cfg,dataApT_AT70);
freqMS_TA70=ft_freqanalysis(cfg,dataMS_TA70);
freqApT_TA70=ft_freqanalysis(cfg,dataApT_TA70);

%%
figure;
subplot(3,2,1);plot(tlockMS_AT70.time,tlockMS_AT70.avg,'b');hold on;plot(tlockApT_AT70.time,tlockApT_AT70.avg,'g');ylim([-2 2])
ylabel('ERP')
subplot(3,2,2);plot(tlockMS_TA70.time,tlockMS_TA70.avg,'b');hold on;plot(tlockApT_TA70.time,tlockApT_TA70.avg,'g');ylim([-2 2])

subplot(3,2,3);plot(freqMS_AT70.time,abs(squeeze(getplv(freqMS_AT70.fourierspctrm))),'b');hold on;plot(freqApT_AT70.time,abs(squeeze(getplv(freqApT_AT70.fourierspctrm))),'g');ylim([0 1])
ylabel('ITC')
subplot(3,2,4);plot(freqMS_AT70.time,abs(squeeze(getplv(freqMS_TA70.fourierspctrm))),'b');hold on;plot(freqApT_AT70.time,abs(squeeze(getplv(freqApT_TA70.fourierspctrm))),'g');ylim([0 1])

subplot(3,2,5);plot(freqMS_AT70.time,squeeze(mean(abs(freqMS_AT70.fourierspctrm),1)),'b');hold on;plot(freqApT_AT70.time,squeeze(mean(abs(freqApT_AT70.fourierspctrm),1)),'g');ylim([0 3])
ylabel('TFP')
subplot(3,2,6);plot(freqMS_AT70.time,squeeze(mean(abs(freqMS_TA70.fourierspctrm),1)),'b');hold on;plot(freqApT_AT70.time,squeeze(mean(abs(freqApT_TA70.fourierspctrm),1)),'g');ylim([0 3])


