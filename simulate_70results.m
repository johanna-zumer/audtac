% script to show possibility of different response profiles to two
% differnet stimulus conditions: same ERP result but opposing ITC

% In addition to needing FieldTrip installed on path, this also requires
% 'getplv' also available in the same Github directory.

% The mean and standard error for ERP and TFP can be found for a given
% iteration; however to get the standard error for ITC, the simulation is
% run multiple times.

%%
% sets up path
eeg_legomagic_preamble

%%

for iter=1:100
  
  dataATpN_AT70=[];
  dataApT_AT70=[];
  dataATpN_TA70=[];
  dataApT_TA70=[];
  
  dataATpN_AT70.label{1}='chan';
  dataApT_AT70.label{1}='chan';
  dataATpN_TA70.label{1}='chan';
  dataApT_TA70.label{1}='chan';
  
  time=-1:.001:2;
  for nn=1:100  % trials
    dataATpN_AT70.time{nn}=time;
    dataApT_AT70.time{nn}=time;
    dataATpN_TA70.time{nn}=time;
    dataApT_TA70.time{nn}=time;
    
    dataATpN_AT70.trial{nn}= [cos(4*2*pi*time(1:1000)+rand*2*pi), 1.1*cos(4*2*pi*time(1001:2000)+0.5*rand*2*pi), cos(4*2*pi*time(2001:3001)+rand*2*pi)];
    dataApT_AT70.trial{nn}=[cos(4*2*pi*time(1:1000)+rand*2*pi), 3.0*cos(4*2*pi*time(1001:2000)+0.6*rand*2*pi), cos(4*2*pi*time(2001:3001)+rand*2*pi)];
    
    dataATpN_TA70.trial{nn}= [cos(4*2*pi*time(1:1000)+rand*2*pi), 3.0*cos(4*2*pi*time(1001:2000)+0.6*rand*2*pi), cos(4*2*pi*time(2001:3001)+rand*2*pi)];
    dataApT_TA70.trial{nn}=[cos(4*2*pi*time(1:1000)+rand*2*pi), 5.0*cos(4*2*pi*time(1001:2000)+0.5*rand*2*pi), cos(4*2*pi*time(2001:3001)+rand*2*pi)];
  end % nn
  
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
  
  itc_ATpN_AT70(:,iter)=getplv(freqATpN_AT70.fourierspctrm);
  itc_ApT_AT70(:,iter) =getplv(freqApT_AT70.fourierspctrm);
  itc_ATpN_TA70(:,iter)=getplv(freqATpN_TA70.fourierspctrm);
  itc_ApT_TA70(:,iter) =getplv(freqApT_TA70.fourierspctrm);
  
end % iterations

save([edir 'sim70.mat'],'itc*','freq*','tlock*')


%%
% This will plot the mean for ERP, ITC, TFP for the last iteration
% and standard error of the mean for ERP and TFP.
% and then use the mean/stE over all iterations for final plot of ITC.

try
  close 1
end
figure(1);
subplot(4,2,1);
plot(tlockATpN_AT70.time,tlockATpN_AT70.avg,'b');hold on;
plot(tlockATpN_AT70.time,tlockATpN_AT70.avg-sqrt(tlockATpN_AT70.var/99),'b--');hold on;
plot(tlockATpN_AT70.time,tlockATpN_AT70.avg+sqrt(tlockATpN_AT70.var/99),'b--');hold on;
plot(tlockApT_AT70.time,tlockApT_AT70.avg,'g');
plot(tlockApT_AT70.time,tlockApT_AT70.avg-sqrt(tlockApT_AT70.var/99),'g--');
plot(tlockApT_AT70.time,tlockApT_AT70.avg+sqrt(tlockApT_AT70.var/99),'g--');
ylim([-2 2])
ylabel('ERP')
title('AT70')

subplot(4,2,2);
plot(tlockATpN_TA70.time,tlockATpN_TA70.avg,'b');hold on;
plot(tlockATpN_TA70.time,tlockATpN_TA70.avg-sqrt(tlockATpN_TA70.var/99),'b--');hold on;
plot(tlockATpN_TA70.time,tlockATpN_TA70.avg+sqrt(tlockATpN_TA70.var/99),'b--');hold on;
plot(tlockApT_TA70.time,tlockApT_TA70.avg,'g');
plot(tlockApT_TA70.time,tlockApT_TA70.avg-sqrt(tlockApT_TA70.var/99),'g');
plot(tlockApT_TA70.time,tlockApT_TA70.avg+sqrt(tlockApT_TA70.var/99),'g');
ylim([-2 2])
title('TA70')

subplot(4,2,3);plot(freqATpN_AT70.time,abs(squeeze(getplv(freqATpN_AT70.fourierspctrm))),'b');hold on;
plot(freqApT_AT70.time,abs(squeeze(getplv(freqApT_AT70.fourierspctrm))),'g');ylim([0 1])
ylabel('ITC')

subplot(4,2,4);plot(freqATpN_AT70.time,abs(squeeze(getplv(freqATpN_TA70.fourierspctrm))),'b');hold on;
plot(freqApT_AT70.time,abs(squeeze(getplv(freqApT_TA70.fourierspctrm))),'g');ylim([0 1])

subplot(4,2,5);
plot(freqATpN_AT70.time,squeeze(mean(abs(freqATpN_AT70.fourierspctrm),1)),'b');hold on;
plot(freqATpN_AT70.time,squeeze(mean(abs(freqATpN_AT70.fourierspctrm),1))-squeeze(std(abs(freqATpN_AT70.fourierspctrm),[],1)/sqrt(99)),'b--');hold on;
plot(freqATpN_AT70.time,squeeze(mean(abs(freqATpN_AT70.fourierspctrm),1))+squeeze(std(abs(freqATpN_AT70.fourierspctrm),[],1)/sqrt(99)),'b--');hold on;
plot(freqApT_AT70.time,squeeze(mean(abs(freqApT_AT70.fourierspctrm),1)),'g');
plot(freqApT_AT70.time,squeeze(mean(abs(freqApT_AT70.fourierspctrm),1))-squeeze(std(abs(freqApT_AT70.fourierspctrm),[],1)/sqrt(99)),'g--');
plot(freqApT_AT70.time,squeeze(mean(abs(freqApT_AT70.fourierspctrm),1))+squeeze(std(abs(freqApT_AT70.fourierspctrm),[],1)/sqrt(99)),'g--');
ylim([0 3])
ylabel('TFP')
legend({'AT+N' '' '' 'A+T'})

subplot(4,2,6);
plot(freqATpN_AT70.time,squeeze(mean(abs(freqATpN_TA70.fourierspctrm),1)),'b');hold on;
plot(freqATpN_AT70.time,squeeze(mean(abs(freqATpN_TA70.fourierspctrm),1))-squeeze(std(abs(freqATpN_TA70.fourierspctrm),[],1)/sqrt(99)),'b--');hold on;
plot(freqATpN_AT70.time,squeeze(mean(abs(freqATpN_TA70.fourierspctrm),1))+squeeze(std(abs(freqATpN_TA70.fourierspctrm),[],1)/sqrt(99)),'b--');hold on;
plot(freqApT_AT70.time,squeeze(mean(abs(freqApT_TA70.fourierspctrm),1)),'g');
plot(freqApT_AT70.time,squeeze(mean(abs(freqApT_TA70.fourierspctrm),1))-squeeze(std(abs(freqApT_TA70.fourierspctrm),[],1)/sqrt(99)),'g--');
plot(freqApT_AT70.time,squeeze(mean(abs(freqApT_TA70.fourierspctrm),1))+squeeze(std(abs(freqApT_TA70.fourierspctrm),[],1)/sqrt(99)),'g--');
ylim([0 3])

subplot(4,2,7);
plot(freqATpN_AT70.time,mean(abs(itc_ATpN_AT70),2),'b');hold on;
plot(freqATpN_AT70.time,mean(abs(itc_ATpN_AT70),2)-std(abs(itc_ATpN_AT70),[],2)/sqrt(99),'b--');hold on;
plot(freqATpN_AT70.time,mean(abs(itc_ATpN_AT70),2)+std(abs(itc_ATpN_AT70),[],2)/sqrt(99),'b--');hold on;
plot(freqApT_AT70.time, mean(abs(itc_ApT_AT70),2),'g');ylim([0 1])
plot(freqApT_AT70.time, mean(abs(itc_ApT_AT70),2)-std(abs(itc_ApT_AT70),[],2)/sqrt(99),'g--');ylim([0 1])
plot(freqApT_AT70.time, mean(abs(itc_ApT_AT70),2)+std(abs(itc_ApT_AT70),[],2)/sqrt(99),'g--');ylim([0 1])
ylabel('ITC (100 iterations)')

subplot(4,2,8);
plot(freqATpN_AT70.time,mean(abs(itc_ATpN_TA70),2),'b');hold on;
plot(freqATpN_AT70.time,mean(abs(itc_ATpN_TA70),2)-std(abs(itc_ATpN_TA70),[],2)/sqrt(99),'b--');
plot(freqATpN_AT70.time,mean(abs(itc_ATpN_TA70),2)+std(abs(itc_ATpN_TA70),[],2)/sqrt(99),'b--');
plot(freqApT_AT70.time, mean(abs(itc_ApT_TA70),2),'g');
plot(freqApT_AT70.time, mean(abs(itc_ApT_TA70),2)-std(abs(itc_ApT_TA70),[],2)/sqrt(99),'g--');
plot(freqApT_AT70.time, mean(abs(itc_ApT_TA70),2)+std(abs(itc_ApT_TA70),[],2)/sqrt(99),'g--');
ylim([0 1])



print(1,[fdir 'simulation70ms.png'],'-dpng')
print(1,[fdir 'simulation70ms.eps'],'-painters','-depsc')



