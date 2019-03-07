eeg_legomagic_preamble

%%  Stats on phase-sorting of unisensory ERP: attempt 1 (ignore now)

% related to phaset0=1 and phaset0use>0
plotflag=1;
printflag=0;

soalist=[1 3 4 5 6 7 9];
tt=3;
sleep=0;
tacaud=1;
if sleep
  ss=12;
  iter=11;
  iiuse=setdiff(iiBuse,1:7);
  trialkc=0;
else
  ss=10;
  iter=27;
  iiuse=setdiff(iiSuse,1:7);
  trialkc=-1;
end

subuseind=1;
for ii=iiuse
  % for ii=8:18
  cd([edir sub{ii}])
  load(['tlock_ERPphasebin_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat']);
  for bb=1:4
    tlockavg_tac5_phasebin{bb}{subuseind}=tlock_tactlock_phasebin{5,tt,ss,bb};
    tlockavg_aud5_phasebin{bb}{subuseind}=tlock_audtlock_phasebin{5,tt,ss,bb};
    tlockavg_nul5_phasebin{bb}{subuseind}=tlock_nultlock_phasebin{5,tt,ss,bb};
    tlockavg_ms15_phasebin{bb}{subuseind}=tlock_ms1tlock_phasebin{5,tt,ss,bb};
    
    tac5_dof(bb,subuseind)=tlockavg_tac5_phasebin{bb}{subuseind}.dof(18,1000);
    aud5_dof(bb,subuseind)=tlockavg_aud5_phasebin{bb}{subuseind}.dof(18,1000);
    nul5_dof(bb,subuseind)=tlockavg_nul5_phasebin{bb}{subuseind}.dof(18,1000);
    ms15_dof(bb,subuseind)=tlockavg_ms15_phasebin{bb}{subuseind}.dof(18,1000);
  end
  
  subuseind=subuseind+1;
  clear tlock*tlock_phasebin
end % ii

for bb=1:4,
  cfg=[];
  grave_tac_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_tac5_phasebin{bb}{:})
  grave_aud_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_aud5_phasebin{bb}{:})
  grave_nul_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_nul5_phasebin{bb}{:})
  grave_ms1_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_ms15_phasebin{bb}{:})
  
  cfg=[];
  cfg.keepindividual='yes';
  grind_tac_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_tac5_phasebin{bb}{:})
  grind_aud_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_aud5_phasebin{bb}{:})
  grind_nul_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_nul5_phasebin{bb}{:})
  grind_ms1_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_ms15_phasebin{bb}{:})
end

if plotflag
  if sleep
  else
    for bb=1:4
      topoplot_highlight(100+bb,grave_tac_pb{bb},[.0 .37],[]);
      if printflag
        %       print(111,['D:\audtac\figs\gravetacPaud_topoOverSTime_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      end
    end
  end
  
  chanplot{1}={'Fz' 'FC1' 'FC2' 'F1' 'F2' 'C1' 'C2' 'Cz'};
  %   subplot(1,7,figind);
  cfg=[];
  cfg.channel=chanplot{1};
  cfg.xlim=[-0.5 1.1];
  cfg.ylim=[-7 7];
  figure;
  ft_singleplotER(cfg, grave_tac_pb{:})
  legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
  %   xlabel(['Tactile at time 0, ' sleepcond])
  %   ylabel(chanlabel{cc})
  title('Tactile Alone')
  figure;
  ft_singleplotER(cfg, grave_aud_pb{:})
  legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
  title('Auditory Alone')
  figure;
  ft_singleplotER(cfg, grave_nul_pb{:})
  legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
  title('Null Alone')
  figure;
  ft_singleplotER(cfg, grave_ms1_pb{:})
  legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
  title('Multisensory AT0')
  
  
end


load eeg1010_neighb
nsub=subuseind-1;
cfg=[];
cfg.latency=[.0 .4];
% cfg.channel=chanuse;
cfg.neighbours=neighbours;
cfg.parameter='individual';
cfg.method='montecarlo';
cfg.numrandomization=1000;
cfg.correctm='cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.statistic='depsamplesFunivariate';
cfg.design=zeros(2,4*nsub);
cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub) 3*ones(1,nsub) 4*ones(1,nsub)];
cfg.design(2,:)=[1:nsub 1:nsub 1:nsub 1:nsub];
cfg.ivar=1;
cfg.uvar=2;
cfg.tail=1;
stat_tacpb=ft_timelockstatistics(cfg, grind_tac_pb{:});
stat_audpb=ft_timelockstatistics(cfg, grind_aud_pb{:});
stat_nulpb=ft_timelockstatistics(cfg, grind_nul_pb{:});
stat_ms1pb=ft_timelockstatistics(cfg, grind_ms1_pb{:});

save([edir 'stat_pb_uni.mat'],'stat*','grave*','*dof');
save([edir 'grind_pb_uni.mat'],'grind*');

figure;imagesc(stat_nulpb.time,1:63,stat_nulpb.mask);title('Significance F-test mask: Null');xlabel('Time (s)');ylabel('Channel')
figure;imagesc(stat_nulpb.time,1:63,stat_tacpb.mask);title('Significance F-test mask: Tactile');xlabel('Time (s)');ylabel('Channel')
figure;imagesc(stat_nulpb.time,1:63,stat_audpb.mask);title('Significance F-test mask: Auditory');xlabel('Time (s)');ylabel('Channel')
figure;imagesc(stat_nulpb.time,1:63,stat_ms1pb.mask);title('Significance F-test mask: AT0');xlabel('Time (s)');ylabel('Channel')

% Nul is significant up to 120ms.   Do post-hoc t-test on anything showing
% diff after 120ms (namely Tac and AT0).
cfg=[];
cfg.latency=[.12 .4];
cfg.neighbours=neighbours;
cfg.parameter='individual';
cfg.method='montecarlo';
cfg.numrandomization=1000;
cfg.correctm='cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.statistic='depsamplesT';
cfg.design=zeros(2,2*nsub);
cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
cfg.design(2,:)=[1:nsub 1:nsub];
cfg.ivar=1;
cfg.uvar=2;
stat_tacpb_posthocT=ft_timelockstatistics(cfg, grind_tac_pb{1}, grind_tac_pb{4});
stat_ms1pb_posthocT=ft_timelockstatistics(cfg, grind_ms1_pb{1}, grind_ms1_pb{4});
% Tactile does show sig diff: (but not AT0)
figure;imagesc(stat_tacpb_posthocT.time,1:63,stat_tacpb_posthocT.mask);title('Significance T-test mask: Tactile');xlabel('Time (s)');ylabel('Channel')


% chanplot{2}=stat_tacpb.label(find(stat_tacpb.mask(:,145)));
%   cfg=[];
%   cfg.channel=chanplot{2};
%   cfg.xlim=[-0.5 1.1];
%   cfg.ylim=[-7 7];
%   figure;
%   ft_singleplotER(cfg, grave_tac_pb{:})
%   legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
% %   xlabel(['Tactile at time 0, ' sleepcond])
% %   ylabel(chanlabel{cc})
%   title('Tactile Alone')
%   figure;
%   ft_singleplotER(cfg, grave_aud_pb{:})
%   legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
%   title('Auditory Alone')
%   figure;
%   ft_singleplotER(cfg, grave_nul_pb{:})
%   legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
%   title('Null Alone')
%   figure;
%   ft_singleplotER(cfg, grave_ms1_pb{:})
%   legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
%   title('Multisensory AT0')
%
%   chanplot{3}=stat_tacpb.label(find(stat_ms1pb.mask(:,145)));
%   cfg=[];
%   cfg.channel=chanplot{3};
%   cfg.xlim=[-0.5 1.1];
%   cfg.ylim=[-7 7];
%   figure;
%   ft_singleplotER(cfg, grave_tac_pb{:})
%   legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
% %   xlabel(['Tactile at time 0, ' sleepcond])
% %   ylabel(chanlabel{cc})
%   title('Tactile Alone')
%   figure;
%   ft_singleplotER(cfg, grave_aud_pb{:})
%   legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
%   title('Auditory Alone')
%   figure;
%   ft_singleplotER(cfg, grave_nul_pb{:})
%   legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
%   title('Null Alone')
%   figure;
%   ft_singleplotER(cfg, grave_ms1_pb{:})
%   legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
%   title('Multisensory AT0')


%%  Stats on phase-sorting of unisensory and multisensry ERP versus Nul (use)

% related to phaset0=1 and phaset0use>0
plotflag=0;
statsflag=1;
soalist=[1 3 4 5 6 7 9];
tt=3;
sleep=1;
tacaud=1;
mcseed=13;

if sleep
  ss=12;
  iter=11;
  iiuse=setdiff(iiBuse,1:7);
  trialkc=-1;
else
  ss=10;
  iter=27;
  iiuse=setdiff(iiSuse,1:7);
  trialkc=-1;
end

subuseind=1;
for ii=iiuse
  % for ii=8:18
  cd([edir sub{ii}])
  load(['tlock_ERPphasebin_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat']);
  for bb=1:4
    
    tlockavg_tac5con1_phasebin{bb}{subuseind}=tlock_tac4mscon1_phasebin{5,3,ss,bb};
    tlockavg_aud5con1_phasebin{bb}{subuseind}=tlock_aud4mscon1_phasebin{5,3,ss,bb};
    tlockavg_nul5con1_phasebin{bb}{subuseind}=tlock_nul4mscon1_phasebin{5,3,ss,bb};
    
    if sleep
      tlockavg_tac5con1_absbin{bb}{subuseind}=tlock_tac4mscon1_absbin{5,3,ss,bb};
      tlockavg_aud5con1_absbin{bb}{subuseind}=tlock_aud4mscon1_absbin{5,3,ss,bb};
      tlockavg_nul5con1_absbin{bb}{subuseind}=tlock_nul4mscon1_absbin{5,3,ss,bb};
    end
    
    for ll=soalist
      tlockavg_msllcon1_phasebin{ll,bb}{subuseind}=tlock_ms4mscon1_phasebin{ll,3,ss,bb};
      tlockavg_nulllcon1_phasebin{ll,bb}{subuseind}=tlock_nul4mscon1_phasebin{ll,3,ss,bb};
      tlockavg_msllcon2_phasebin{ll,bb}{subuseind}=tlock_ms4mscon2_phasebin{ll,3,ss,bb};
      tlockavg_nulllcon2_phasebin{ll,bb}{subuseind}=tlock_nul4mscon2_phasebin{ll,3,ss,bb};
      tlockavg_audllcon1_phasebin{ll,bb}{subuseind}=tlock_aud4mscon1_phasebin{ll,3,ss,bb};
      if sleep
        tlockavg_msllcon1_absbin{ll,bb}{subuseind}=tlock_ms4mscon1_absbin{ll,3,ss,bb};
        tlockavg_nulllcon1_absbin{ll,bb}{subuseind}=tlock_nul4mscon1_absbin{ll,3,ss,bb};
        tlockavg_msllcon2_absbin{ll,bb}{subuseind}=tlock_ms4mscon2_absbin{ll,3,ss,bb};
        tlockavg_nulllcon2_absbin{ll,bb}{subuseind}=tlock_nul4mscon2_absbin{ll,3,ss,bb};
        tlockavg_audllcon1_absbin{ll,bb}{subuseind}=tlock_aud4mscon1_absbin{ll,3,ss,bb};
      end
    end % ll
    
  end
  
  subuseind=subuseind+1;
  clear tlock*tlock_phasebin
end % ii

for bb=1:4,
  cfg=[];
  grave_tac_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_tac5con1_phasebin{bb}{:});
  grave_aud_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_aud5con1_phasebin{bb}{:});
  grave_nul_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_nul5con1_phasebin{bb}{:});
  
  cfg=[];
  cfg.keepindividual='yes';
  grind_tac_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_tac5con1_phasebin{bb}{:});
  grind_aud_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_aud5con1_phasebin{bb}{:});
  grind_nul_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_nul5con1_phasebin{bb}{:});
  
  cfg=[];
  cfg.parameter='individual';
  cfg.operation='subtract';
  grind_tacMnul_pb{bb}=ft_math(cfg,grind_tac_pb{bb},grind_nul_pb{bb});
  grind_audMnul_pb{bb}=ft_math(cfg,grind_aud_pb{bb},grind_nul_pb{bb});
  
  grave_tacMnul_pb{bb}=grind_tacMnul_pb{bb};
  grave_tacMnul_pb{bb}.avg=squeeze(mean(grave_tacMnul_pb{bb}.individual,1));
  grave_tacMnul_pb{bb}=rmfield(grave_tacMnul_pb{bb},'individual');
  grave_tacMnul_pb{bb}.dimord='chan_time';
  
  grave_audMnul_pb{bb}=grind_audMnul_pb{bb};
  grave_audMnul_pb{bb}.avg=squeeze(mean(grave_audMnul_pb{bb}.individual,1));
  grave_audMnul_pb{bb}=rmfield(grave_audMnul_pb{bb},'individual');
  grave_audMnul_pb{bb}.dimord='chan_time';
  
  if sleep
    cfg=[];
    grave_tac_ab{bb}=ft_timelockgrandaverage(cfg,tlockavg_tac5con1_absbin{bb}{:});
    grave_aud_ab{bb}=ft_timelockgrandaverage(cfg,tlockavg_aud5con1_absbin{bb}{:});
    grave_nul_ab{bb}=ft_timelockgrandaverage(cfg,tlockavg_nul5con1_absbin{bb}{:});
    
    cfg=[];
    cfg.keepindividual='yes';
    grind_tac_ab{bb}=ft_timelockgrandaverage(cfg,tlockavg_tac5con1_absbin{bb}{:});
    grind_aud_ab{bb}=ft_timelockgrandaverage(cfg,tlockavg_aud5con1_absbin{bb}{:});
    grind_nul_ab{bb}=ft_timelockgrandaverage(cfg,tlockavg_nul5con1_absbin{bb}{:});
    
    cfg=[];
    cfg.parameter='individual';
    cfg.operation='subtract';
    grind_tacMnul_ab{bb}=ft_math(cfg,grind_tac_ab{bb},grind_nul_ab{bb});
    grind_audMnul_ab{bb}=ft_math(cfg,grind_aud_ab{bb},grind_nul_ab{bb});
    
    grave_tacMnul_ab{bb}=grind_tacMnul_ab{bb};
    grave_tacMnul_ab{bb}.avg=squeeze(mean(grave_tacMnul_ab{bb}.individual,1));
    grave_tacMnul_ab{bb}=rmfield(grave_tacMnul_ab{bb},'individual');
    grave_tacMnul_ab{bb}.dimord='chan_time';
    
    grave_audMnul_ab{bb}=grind_audMnul_ab{bb};
    grave_audMnul_ab{bb}.avg=squeeze(mean(grave_audMnul_ab{bb}.individual,1));
    grave_audMnul_ab{bb}=rmfield(grave_audMnul_ab{bb},'individual');
    grave_audMnul_ab{bb}.dimord='chan_time';
  end
  
  for ll=soalist
    cfg=[];
    grave_msllcon1_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_msllcon1_phasebin{ll,bb}{:});
    grave_nulllcon1_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_nulllcon1_phasebin{ll,bb}{:});
    grave_msllcon2_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_msllcon2_phasebin{ll,bb}{:});
    grave_nulllcon2_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_nulllcon2_phasebin{ll,bb}{:});
    grave_audllcon1_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_audllcon1_phasebin{ll,bb}{:});
    cfg=[];
    cfg.keepindividual='yes';
    grind_msllcon1_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_msllcon1_phasebin{ll,bb}{:});
    grind_nulllcon1_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_nulllcon1_phasebin{ll,bb}{:});
    grind_msllcon2_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_msllcon2_phasebin{ll,bb}{:});
    grind_nulllcon2_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_nulllcon2_phasebin{ll,bb}{:});
    grind_audllcon1_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_audllcon1_phasebin{ll,bb}{:});
    
    cfg=[];
    cfg.parameter='individual';
    cfg.operation='subtract';
    grind_msMnul1_pb{ll,bb}=ft_math(cfg,grind_msllcon1_pb{ll,bb},grind_nulllcon1_pb{ll,bb});
    grind_msMnul2_pb{ll,bb}=ft_math(cfg,grind_msllcon2_pb{ll,bb},grind_nulllcon2_pb{ll,bb});
    grind_msMaud1_pb{ll,bb}=ft_math(cfg,grind_msllcon1_pb{ll,bb},grind_audllcon1_pb{ll,bb});
    
    grave_msMnul1_pb{ll,bb}=grind_msMnul1_pb{ll,bb};
    grave_msMnul1_pb{ll,bb}.avg=squeeze(mean(grave_msMnul1_pb{ll,bb}.individual,1));
    grave_msMnul1_pb{ll,bb}=rmfield(grave_msMnul1_pb{ll,bb},'individual');
    grave_msMnul1_pb{ll,bb}.dimord='chan_time';
    
    grave_msMnul2_pb{ll,bb}=grind_msMnul2_pb{ll,bb};
    grave_msMnul2_pb{ll,bb}.avg=squeeze(mean(grave_msMnul2_pb{ll,bb}.individual,1));
    grave_msMnul2_pb{ll,bb}=rmfield(grave_msMnul2_pb{ll,bb},'individual');
    grave_msMnul2_pb{ll,bb}.dimord='chan_time';
    
    grave_msMaud1_pb{ll,bb}=grind_msMaud1_pb{ll,bb};
    grave_msMaud1_pb{ll,bb}.avg=squeeze(mean(grave_msMaud1_pb{ll,bb}.individual,1));
    grave_msMaud1_pb{ll,bb}=rmfield(grave_msMaud1_pb{ll,bb},'individual');
    grave_msMaud1_pb{ll,bb}.dimord='chan_time';
    
    if sleep
      cfg=[];
      grave_msllcon1_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_msllcon1_absbin{ll,bb}{:});
      grave_nulllcon1_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_nulllcon1_absbin{ll,bb}{:});
      grave_msllcon2_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_msllcon2_absbin{ll,bb}{:});
      grave_nulllcon2_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_nulllcon2_absbin{ll,bb}{:});
      grave_audllcon1_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_audllcon1_absbin{ll,bb}{:});
      cfg=[];
      cfg.keepindividual='yes';
      grind_msllcon1_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_msllcon1_absbin{ll,bb}{:});
      grind_nulllcon1_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_nulllcon1_absbin{ll,bb}{:});
      grind_msllcon2_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_msllcon2_absbin{ll,bb}{:});
      grind_nulllcon2_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_nulllcon2_absbin{ll,bb}{:});
      grind_audllcon1_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_audllcon1_absbin{ll,bb}{:});
      
      cfg=[];
      cfg.parameter='individual';
      cfg.operation='subtract';
      grind_msMnul1_ab{ll,bb}=ft_math(cfg,grind_msllcon1_ab{ll,bb},grind_nulllcon1_ab{ll,bb});
      grind_msMnul2_ab{ll,bb}=ft_math(cfg,grind_msllcon2_ab{ll,bb},grind_nulllcon2_ab{ll,bb});
      grind_msMaud1_ab{ll,bb}=ft_math(cfg,grind_msllcon1_ab{ll,bb},grind_audllcon1_ab{ll,bb});
      
      grave_msMnul1_ab{ll,bb}=grind_msMnul1_ab{ll,bb};
      grave_msMnul1_ab{ll,bb}.avg=squeeze(mean(grave_msMnul1_ab{ll,bb}.individual,1));
      grave_msMnul1_ab{ll,bb}=rmfield(grave_msMnul1_ab{ll,bb},'individual');
      grave_msMnul1_ab{ll,bb}.dimord='chan_time';
      
      grave_msMnul2_ab{ll,bb}=grind_msMnul2_ab{ll,bb};
      grave_msMnul2_ab{ll,bb}.avg=squeeze(mean(grave_msMnul2_ab{ll,bb}.individual,1));
      grave_msMnul2_ab{ll,bb}=rmfield(grave_msMnul2_ab{ll,bb},'individual');
      grave_msMnul2_ab{ll,bb}.dimord='chan_time';
      
      grave_msMaud1_ab{ll,bb}=grind_msMaud1_ab{ll,bb};
      grave_msMaud1_ab{ll,bb}.avg=squeeze(mean(grave_msMaud1_ab{ll,bb}.individual,1));
      grave_msMaud1_ab{ll,bb}=rmfield(grave_msMaud1_ab{ll,bb},'individual');
      grave_msMaud1_ab{ll,bb}.dimord='chan_time';
      
    end
  end % ll
end % bb
% save([edir 'grind_pb_uninul.mat'],'grind*')
save([edir 'grind_pb_uninul_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'grind*');

if statsflag
  load eeg1010_neighb
  nsub=subuseind-1;
  cfg=[];
  cfg.neighbours=neighbours;
  cfg.parameter='individual';
  cfg.method='montecarlo';
  cfg.numrandomization=1000;
  cfg.correctm='cluster';
  cfg.clusteralpha = 0.05;
  cfg.clusterstatistic = 'maxsum';
  cfg.minnbchan = 2;
  cfg.ivar=1;
  cfg.uvar=2;
  cfg.statistic='depsamplesFunivariate';
  cfg.design=zeros(2,4*nsub);
  cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub) 3*ones(1,nsub) 4*ones(1,nsub)];
  cfg.design(2,:)=[1:nsub 1:nsub 1:nsub 1:nsub];
  cfg.tail=1;
  cfg.randomseed=mcseed;
  
  if sleep
    cfg.latency=[.1 .75];
  else
    cfg.latency=[.1 .45];
  end
  
  stat_tacMnul_pb=ft_timelockstatistics(cfg, grind_tacMnul_pb{:});
  stat_audMnul_pb=ft_timelockstatistics(cfg, grind_audMnul_pb{:});
  
  if sleep
    stat_tacMnul_ab=ft_timelockstatistics(cfg, grind_tacMnul_ab{:});
    stat_audMnul_ab=ft_timelockstatistics(cfg, grind_audMnul_ab{:});
  end
  
  for ll=soalist
    if sleep
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[.1 .75];
      elseif ll==6
        cfg.latency=[.12 .77];
      elseif ll==7
        cfg.latency=[.17 .82];
      elseif ll==9
        cfg.latency=[.6 1.25];
      end
    else
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[.1 .45];
      elseif ll==6
        cfg.latency=[.12 .47];
      elseif ll==7
        cfg.latency=[.17 .52];
      elseif ll==9
        cfg.latency=[.6 .95];
      end
    end
    stat_msMnul1_pb{ll}=ft_timelockstatistics(cfg, grind_msMnul1_pb{ll,:});
    stat_msMnul2_pb{ll}=ft_timelockstatistics(cfg, grind_msMnul2_pb{ll,:});
    stat_msMaud1_pb{ll}=ft_timelockstatistics(cfg, grind_msMaud1_pb{ll,:});
    if sleep
      stat_msMnul1_ab{ll}=ft_timelockstatistics(cfg, grind_msMnul1_ab{ll,:});
      stat_msMnul2_ab{ll}=ft_timelockstatistics(cfg, grind_msMnul2_ab{ll,:});
      stat_msMaud1_ab{ll}=ft_timelockstatistics(cfg, grind_msMaud1_ab{ll,:});
    end
    save([edir 'stat_pb_uninul_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '_mcseed' num2str(mcseed) '.mat'],'stat*')
  end % ll
  
end % statsflag




if plotflag
  if ~exist('grave_tacMnul_pb','var')
    if ~exist('grind_tacMnul_pb','var')
      load([edir 'grind_pb_uninul_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat']);
    end
    
    for bb=1:4
      
      grave_tacMnul_pb{bb}=grind_tacMnul_pb{bb};
      grave_tacMnul_pb{bb}.avg=squeeze(mean(grave_tacMnul_pb{bb}.individual,1))
      grave_tacMnul_pb{bb}=rmfield(grave_tacMnul_pb{bb},'individual');
      grave_tacMnul_pb{bb}.dimord='chan_time';
      
      grave_audMnul_pb{bb}=grind_audMnul_pb{bb};
      grave_audMnul_pb{bb}.avg=squeeze(mean(grave_audMnul_pb{bb}.individual,1))
      grave_audMnul_pb{bb}=rmfield(grave_audMnul_pb{bb},'individual');
      grave_audMnul_pb{bb}.dimord='chan_time';
      
      if sleep
        grave_tacMnul_ab{bb}=grind_tacMnul_ab{bb};
        grave_tacMnul_ab{bb}.avg=squeeze(mean(grave_tacMnul_ab{bb}.individual,1))
        grave_tacMnul_ab{bb}=rmfield(grave_tacMnul_ab{bb},'individual');
        grave_tacMnul_ab{bb}.dimord='chan_time';
        
        grave_audMnul_ab{bb}=grind_audMnul_ab{bb};
        grave_audMnul_ab{bb}.avg=squeeze(mean(grave_audMnul_ab{bb}.individual,1))
        grave_audMnul_ab{bb}=rmfield(grave_audMnul_ab{bb},'individual');
        grave_audMnul_ab{bb}.dimord='chan_time';
        
      end
      
      for ll=soalist
        grave_msMnul1_pb{ll,bb}=grind_msMnul1_pb{ll,bb};
        grave_msMnul1_pb{ll,bb}.avg=squeeze(mean(grave_msMnul1_pb{ll,bb}.individual,1))
        grave_msMnul1_pb{ll,bb}=rmfield(grave_msMnul1_pb{ll,bb},'individual');
        grave_msMnul1_pb{ll,bb}.dimord='chan_time';
        
        grave_msMnul2_pb{ll,bb}=grind_msMnul2_pb{ll,bb};
        grave_msMnul2_pb{ll,bb}.avg=squeeze(mean(grave_msMnul2_pb{ll,bb}.individual,1))
        grave_msMnul2_pb{ll,bb}=rmfield(grave_msMnul2_pb{ll,bb},'individual');
        grave_msMnul2_pb{ll,bb}.dimord='chan_time';
        
        grave_msMaud1_pb{ll,bb}=grind_msMaud1_pb{ll,bb};
        grave_msMaud1_pb{ll,bb}.avg=squeeze(mean(grave_msMaud1_pb{ll,bb}.individual,1))
        grave_msMaud1_pb{ll,bb}=rmfield(grave_msMaud1_pb{ll,bb},'individual');
        grave_msMaud1_pb{ll,bb}.dimord='chan_time';
        
        if sleep
          grave_msMnul1_ab{ll,bb}=grind_msMnul1_ab{ll,bb};
          grave_msMnul1_ab{ll,bb}.avg=squeeze(mean(grave_msMnul1_ab{ll,bb}.individual,1))
          grave_msMnul1_ab{ll,bb}=rmfield(grave_msMnul1_ab{ll,bb},'individual');
          grave_msMnul1_ab{ll,bb}.dimord='chan_time';
          
          grave_msMnul2_ab{ll,bb}=grind_msMnul2_ab{ll,bb};
          grave_msMnul2_ab{ll,bb}.avg=squeeze(mean(grave_msMnul2_ab{ll,bb}.individual,1))
          grave_msMnul2_ab{ll,bb}=rmfield(grave_msMnul2_ab{ll,bb},'individual');
          grave_msMnul2_ab{ll,bb}.dimord='chan_time';
          
          grave_msMaud1_ab{ll,bb}=grind_msMaud1_ab{ll,bb};
          grave_msMaud1_ab{ll,bb}.avg=squeeze(mean(grave_msMaud1_ab{ll,bb}.individual,1))
          grave_msMaud1_ab{ll,bb}=rmfield(grave_msMaud1_ab{ll,bb},'individual');
          grave_msMaud1_ab{ll,bb}.dimord='chan_time';
        end
      end
      
    end  % bb
  end
  
  chanplot{1}={'Fz' 'Cz' 'F1' 'F2' 'FC1' 'FC2' 'C1' 'C2'}; % frontocentral
  chanplot{2}={'CP5' 'POz' 'Pz' 'P3' 'P4' 'C4' 'O1' 'O2' 'P7' 'PO7'}; % occipital
  
  
  cfg=[];
  cfg.channel=chanplot{1};
  cfg.xlim=[-0.5 1.1];
  cfg.ylim=[-7 7];
  figure(2);
  ft_singleplotER(cfg, grave_tacMnul_pb{:})
  legend({'Peak' 'P to T' 'Trough' 'T to P'})
  title('T-N')
  
  cfg=[];
  cfg.channel=chanplot{1};
  cfg.xlim=[-0.5 1.1];
  cfg.ylim=[-7 7];
  figure(8);
  ft_singleplotER(cfg, grave_audMnul_pb{:})
  legend({'Peak' 'P to T' 'Trough' 'T to P'})
  title('A-N')
  
  if sleep
    cfg=[];
    cfg.channel=chanplot{1};
    cfg.xlim=[-0.5 1.1];
    cfg.ylim=[-7 7];
    figure(12);
    ft_singleplotER(cfg, grave_tacMnul_ab{:})
    legend({'Low' 'Mid-low' 'Mid-high' 'High'})
    title('T-N')
    
    cfg=[];
    cfg.channel=chanplot{1};
    cfg.xlim=[-0.5 1.1];
    cfg.ylim=[-7 7];
    figure(18);
    ft_singleplotER(cfg, grave_audMnul_ab{:})
    legend({'Low' 'Mid-low' 'Mid-high' 'High'})
    title('A-N')
  end
  
  for ll=soalist
    
    cfg=[];
    cfg.channel=chanplot{1};
    cfg.xlim=[-0.5 1.1];
    cfg.ylim=[-7 7];
    figure(ll);
    ft_singleplotER(cfg, grave_msMnul1_pb{ll,:})
    legend({'Peak' 'P to T' 'Trough' 'T to P'})
    title('AT*-N')
    
    %     cfg=[];
    %     cfg.channel=chanplot{1};
    %     cfg.xlim=[-0.5 1.1];
    %     cfg.ylim=[-7 7];
    %     figure(100+ll);
    %     ft_singleplotER(cfg, grave_msMnul2_pb{ll,:})
    %     legend({'Peak' 'P to T' 'Trough' 'T to P'})
    %     title('AT*-N')
    
    cfg=[];
    cfg.channel=chanplot{1};
    cfg.xlim=[-0.5 1.1];
    cfg.ylim=[-7 7];
    figure(20+ll);
    ft_singleplotER(cfg, grave_msMaud1_pb{ll,:})
    legend({'Peak' 'P to T' 'Trough' 'T to P'})
    title('AT*-N')
    
    if sleep
      cfg=[];
      cfg.channel=chanplot{1};
      cfg.xlim=[-0.5 1.1];
      cfg.ylim=[-7 7];
      figure(10+ll);
      ft_singleplotER(cfg, grave_msMnul1_ab{ll,:})
      legend({'Low' 'Mid-low' 'Mid-high' 'High'})
      title('AT*-N')
      
      %       cfg=[];
      %       cfg.channel=chanplot{1};
      %       cfg.xlim=[-0.5 1.1];
      %       cfg.ylim=[-7 7];
      %       figure(100+10+ll);
      %       ft_singleplotER(cfg, grave_msMnul2_ab{ll,:})
      %       legend({'Low' 'Mid-low' 'Mid-high' 'High'})
      %       title('AT*-N')
      
      cfg=[];
      cfg.channel=chanplot{1};
      cfg.xlim=[-0.5 1.1];
      cfg.ylim=[-7 7];
      figure(30+ll);
      ft_singleplotER(cfg, grave_msMaud1_ab{ll,:})
      legend({'Low' 'Mid-low' 'Mid-high' 'High'})
      title('AT*-N')
    end
  end
  
end

% followup stats





%% Stats on phase-sorting MS ERP

plotflag=1;
printflag=1;
statsflag=0;
stats1flag=0;
sleepshortstatflag=1;

soalist=[1 3 4 5 6 7 9];
tt=3;
sleep=1;
tacaud=1;
if sleep
  ss=12;
  iter=11;
  iiuse=setdiff(iiBuse,1:7);
  trialkc=0;
  if trialkc==0
    iiuse=setdiff(iiuse,[18 24 ]);
  end
else
  ss=10;
  iter=27;
  iiuse=setdiff(iiSuse,1:7);
  trialkc=-1;
end

subuseind=1;
for ii=iiuse
  cd([edir sub{ii}])
  load(['tlock_ERPphasebin_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat']);
  for ll=soalist
    for bb=1:4
      tlockavg_tacPaud1_phasebin{ll,bb}{subuseind}=tlock_tacPaud_4mscon1_phasebin{ll,tt,ss,bb};
      tlockavg_tacMSpN1_phasebin{ll,bb}{subuseind}=tlock_tacMSpN_4mscon1_phasebin{ll,tt,ss,bb};
      if sleep
        tlockavg_tacPaud1_absbin{ll,bb}{subuseind}=tlock_tacPaud_4mscon1_absbin{ll,tt,ss,bb};
        tlockavg_tacMSpN1_absbin{ll,bb}{subuseind}=tlock_tacMSpN_4mscon1_absbin{ll,tt,ss,bb};
      else
        tlockavg_tacPaud1IAF_phasebin{ll,bb}{subuseind}=tlock_tacPaud_4mscon1IAF_phasebin{ll,tt,ss,bb};
        tlockavg_tacMSpN1IAF_phasebin{ll,bb}{subuseind}=tlock_tacMSpN_4mscon1IAF_phasebin{ll,tt,ss,bb};
      end
    end % bb
  end % ll
  freqsubfinal(subuseind)=freqsub(ii);
  try freqsleepfinal(subuseind)=freqsleep(ii); catch end
  
  subuseind=subuseind+1;
  clear tlock*tlock_phasebin
end  % ii
save([edir 'freqfinal_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'freq*final');

% subuseind=1;
% for ii=iiuse
%   cd([edir sub{ii}])
%   load(['tlock_ERPphasebin_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'freqsub');
%   freqsubfinal(subuseind)=freqsub(ii);
%
%   subuseind=subuseind+1;
%   clear tlock*tlock_phasebin
% end  % ii
figure(600);hist(freqsubfinal,7);

for ll=soalist
  for bb=1:4,
    %     cfg=[];
    %     cfg.parameter='avg';
    %     cfg.operation='subtract';
    %     tlockavg_TPAmMSPN1_pb{ll,bb}=ft_math(cfg,tlockavg_tacPaud1_phasebin{ll,bb}{:},tlockavg_tacMSpN1_phasebin{ll,bb}{:});
    
    cfg=[];
    grave_tacPaud1_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacPaud1_phasebin{ll,bb}{:});
    grave_tacMSpN1_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacMSpN1_phasebin{ll,bb}{:});
    if sleep
      grave_tacPaud1_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacPaud1_absbin{ll,bb}{:});
      grave_tacMSpN1_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacMSpN1_absbin{ll,bb}{:});
    else
      grave_tacPaud1IAF_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacPaud1IAF_phasebin{ll,bb}{:});
      grave_tacMSpN1IAF_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacMSpN1IAF_phasebin{ll,bb}{:});
    end
    %     grave_TPAmMSPN1_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_TPAmMSPN1_phasebin{ll,bb}{:});
    
    cfg=[];
    cfg.keepindividual='yes';
    grind_tacPaud1_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacPaud1_phasebin{ll,bb}{:});
    grind_tacMSpN1_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacMSpN1_phasebin{ll,bb}{:});
    if sleep
      grind_tacPaud1_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacPaud1_absbin{ll,bb}{:});
      grind_tacMSpN1_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacMSpN1_absbin{ll,bb}{:});
    else
      grind_tacPaud1IAF_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacPaud1IAF_phasebin{ll,bb}{:});
      grind_tacMSpN1IAF_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacMSpN1IAF_phasebin{ll,bb}{:});
    end
    %     grind_TPAmMSPN1_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_TPAmMSPN1_phasebin{ll,bb}{:});
    
    cfg=[];
    cfg.parameter='individual';
    cfg.operation='subtract';
    grind_TPAmMSPN1_pb{ll,bb}=ft_math(cfg,grind_tacPaud1_pb{ll,bb},grind_tacMSpN1_pb{ll,bb})
    if sleep
      grind_TPAmMSPN1_ab{ll,bb}=ft_math(cfg,grind_tacPaud1_ab{ll,bb},grind_tacMSpN1_ab{ll,bb})
    else
      grind_TPAmMSPN1IAF_pb{ll,bb}=ft_math(cfg,grind_tacPaud1IAF_pb{ll,bb},grind_tacMSpN1IAF_pb{ll,bb})
    end
    
    grave_TPAmMSPN1_pb{ll,bb}=grind_TPAmMSPN1_pb{ll,bb};
    grave_TPAmMSPN1_pb{ll,bb}.avg=squeeze(mean(grave_TPAmMSPN1_pb{ll,bb}.individual,1))
    grave_TPAmMSPN1_pb{ll,bb}=rmfield(grave_TPAmMSPN1_pb{ll,bb},'individual');
    grave_TPAmMSPN1_pb{ll,bb}.dimord='chan_time';
    
    if sleep
      grave_TPAmMSPN1_ab{ll,bb}=grind_TPAmMSPN1_ab{ll,bb};
      grave_TPAmMSPN1_ab{ll,bb}.avg=squeeze(mean(grave_TPAmMSPN1_ab{ll,bb}.individual,1))
      grave_TPAmMSPN1_ab{ll,bb}=rmfield(grave_TPAmMSPN1_ab{ll,bb},'individual');
      grave_TPAmMSPN1_ab{ll,bb}.dimord='chan_time';
    else
      grave_TPAmMSPN1IAF_pb{ll,bb}=grind_TPAmMSPN1IAF_pb{ll,bb};
      grave_TPAmMSPN1IAF_pb{ll,bb}.avg=squeeze(mean(grave_TPAmMSPN1IAF_pb{ll,bb}.individual,1))
      grave_TPAmMSPN1IAF_pb{ll,bb}=rmfield(grave_TPAmMSPN1IAF_pb{ll,bb},'individual');
      grave_TPAmMSPN1IAF_pb{ll,bb}.dimord='chan_time';
    end
  end
  cfg=[];
  cfg.parameter='individual';
  cfg.operation='subtract';
  grind_TPAmMSPN_4m2_pb{ll}=ft_math(cfg,grind_TPAmMSPN1_pb{ll,4},grind_TPAmMSPN1_pb{ll,2});
  grave_TPAmMSPN_4m2_pb{ll}=grind_TPAmMSPN_4m2_pb{ll};
  grave_TPAmMSPN_4m2_pb{ll}.avg=squeeze(mean(grave_TPAmMSPN_4m2_pb{ll}.individual,1));
  grave_TPAmMSPN_4m2_pb{ll}=rmfield(grave_TPAmMSPN_4m2_pb{ll},'individual');
  grave_TPAmMSPN_4m2_pb{ll}.dimord='chan_time';
  
  grind_TPAmMSPN_3m1_pb{ll}=ft_math(cfg,grind_TPAmMSPN1_pb{ll,3},grind_TPAmMSPN1_pb{ll,1});
  grave_TPAmMSPN_3m1_pb{ll}=grind_TPAmMSPN_3m1_pb{ll};
  grave_TPAmMSPN_3m1_pb{ll}.avg=squeeze(mean(grave_TPAmMSPN_3m1_pb{ll}.individual,1));
  grave_TPAmMSPN_3m1_pb{ll}=rmfield(grave_TPAmMSPN_3m1_pb{ll},'individual');
  grave_TPAmMSPN_3m1_pb{ll}.dimord='chan_time';
  
  if sleep
    grind_TPAmMSPN_4m1_ab{ll}=ft_math(cfg,grind_TPAmMSPN1_ab{ll,4},grind_TPAmMSPN1_ab{ll,1});
    grave_TPAmMSPN_4m1_ab{ll}=grind_TPAmMSPN_4m1_ab{ll};
    grave_TPAmMSPN_4m1_ab{ll}.avg=squeeze(mean(grave_TPAmMSPN_4m1_ab{ll}.individual,1));
    grave_TPAmMSPN_4m1_ab{ll}=rmfield(grave_TPAmMSPN_4m1_ab{ll},'individual');
    grave_TPAmMSPN_4m1_ab{ll}.dimord='chan_time';
  else
    grind_TPAmMSPNIAF_4m2_pb{ll}=ft_math(cfg,grind_TPAmMSPN1IAF_pb{ll,4},grind_TPAmMSPN1IAF_pb{ll,2});
    grave_TPAmMSPNIAF_4m2_pb{ll}=grind_TPAmMSPNIAF_4m2_pb{ll};
    grave_TPAmMSPNIAF_4m2_pb{ll}.avg=squeeze(mean(grave_TPAmMSPNIAF_4m2_pb{ll}.individual,1));
    grave_TPAmMSPNIAF_4m2_pb{ll}=rmfield(grave_TPAmMSPNIAF_4m2_pb{ll},'individual');
    grave_TPAmMSPNIAF_4m2_pb{ll}.dimord='chan_time';
    
    grind_TPAmMSPNIAF_3m1_pb{ll}=ft_math(cfg,grind_TPAmMSPN1IAF_pb{ll,3},grind_TPAmMSPN1IAF_pb{ll,1});
    grave_TPAmMSPNIAF_3m1_pb{ll}=grind_TPAmMSPNIAF_3m1_pb{ll};
    grave_TPAmMSPNIAF_3m1_pb{ll}.avg=squeeze(mean(grave_TPAmMSPNIAF_3m1_pb{ll}.individual,1));
    grave_TPAmMSPNIAF_3m1_pb{ll}=rmfield(grave_TPAmMSPNIAF_3m1_pb{ll},'individual');
    grave_TPAmMSPNIAF_3m1_pb{ll}.dimord='chan_time';
  end
end % ll


if plotflag
  for ll=soalist
    
    %     if sleep
    %     else
    %       for bb=1:4
    %         topoplot_highlight(100+bb,grave_tacMSpN1_pb{ll,bb},[.0 .37],[]);
    %         if printflag
    %           %       print(111,['D:\audtac\figs\gravetacPaud_topoOverSTime_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
    %         end
    %       end
    %     end
    
    chanplot{1}={'Fz' 'FC1' 'FC2' 'F1' 'F2' 'C1' 'C2' 'Cz'};
    %   subplot(1,7,figind);
    cfg=[];
    cfg.channel=chanplot{1};
    cfg.xlim=[-0.5 1.1];
    cfg.ylim=[-7 7];
    figure(ll);
    ft_singleplotER(cfg, grave_tacPaud1_pb{ll,:})
    legend({'Peak' 'P to T' 'Trough' 'T to P'})
    %     legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
    %   xlabel(['Tactile at time 0, ' sleepcond])
    %   ylabel(chanlabel{cc})
    title('T + A')
    
    figure(10+ll);
    ft_singleplotER(cfg, grave_tacMSpN1_pb{ll,:})
    legend({'Peak' 'P to T' 'Trough' 'T to P'})
    title('MS + N')
    
    figure(20+ll);
    ft_singleplotER(cfg, grave_TPAmMSPN1_pb{ll,:})
    legend({'Peak' 'P to T' 'Trough' 'T to P'})
    title('(T + A)- (MS + N)')
    
    figure(30+ll);
    ft_singleplotER(cfg, grave_TPAmMSPN_4m2_pb{ll}, grave_TPAmMSPN_3m1_pb{ll})
    title('(T + A)- (MS + N)')
    legend({'TtoP - PtoT' 'T - P'})
    
    if sleep
      figure(100+ll);
      ft_singleplotER(cfg, grave_tacPaud1_ab{ll,:})
      legend({'Low' 'Mid-low' 'Mid-high' 'High'})
      title('T + A')
      
      figure(100+10+ll);
      ft_singleplotER(cfg, grave_tacMSpN1_ab{ll,:})
      legend({'Low' 'Mid-low' 'Mid-high' 'High'})
      title('MS + N')
      
      figure(100+20+ll);
      ft_singleplotER(cfg, grave_TPAmMSPN1_ab{ll,:})
      legend({'Low' 'Mid-low' 'Mid-high' 'High'})
      title('(T + A)- (MS + N)')
      
      figure(100+30+ll);
      ft_singleplotER(cfg, grave_TPAmMSPN_4m1_ab{ll})
      title('(T + A)- (MS + N)')
      legend({'High vs Low'})
    else
      figure(100+ll);
      ft_singleplotER(cfg, grave_tacPaud1IAF_pb{ll,:})
      legend({'Peak' 'P to T' 'Trough' 'T to P'})
      title('T + A')
      
      figure(100+10+ll);
      ft_singleplotER(cfg, grave_tacMSpN1IAF_pb{ll,:})
      legend({'Peak' 'P to T' 'Trough' 'T to P'})
      title('MS + N')
      
      figure(100+20+ll);
      ft_singleplotER(cfg, grave_TPAmMSPN1IAF_pb{ll,:})
      legend({'Peak' 'P to T' 'Trough' 'T to P'})
      title('(T + A)- (MS + N)')
      
      figure(100+30+ll);
      ft_singleplotER(cfg, grave_TPAmMSPNIAF_4m2_pb{ll}, grave_TPAmMSPNIAF_3m1_pb{ll})
      title('(T + A)- (MS + N)')
      legend({'TtoP - PtoT' 'T - P'})
    end
    
    
  end
end


if statsflag
  load eeg1010_neighb
  nsub=subuseind-1;
  cfg=[];
  % cfg.channel=chanuse;
  cfg.neighbours=neighbours;
  cfg.parameter='individual';
  cfg.method='montecarlo';
  cfg.numrandomization=1000;
  cfg.correctm='cluster';
  cfg.clusteralpha = 0.05;
  cfg.clusterstatistic = 'maxsum';
  cfg.minnbchan = 2;
  cfg.ivar=1;
  cfg.uvar=2;
  cfg.statistic='depsamplesFunivariate';
  cfg.design=zeros(2,4*nsub);
  cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub) 3*ones(1,nsub) 4*ones(1,nsub)];
  cfg.design(2,:)=[1:nsub 1:nsub 1:nsub 1:nsub];
  cfg.tail=1;
  for ll=soalist
    if sleep && ~sleepshortstatflag
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[0 1];
      elseif ll==6
        cfg.latency=[0 1]+.02;
      elseif ll==7
        cfg.latency=[0 1]+.07;
      elseif ll==9
        cfg.latency=[0 1]+.5;
      end
    elseif sleep && sleepshortstatflag
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[0 .45];
      elseif ll==6
        cfg.latency=[0 .45]+.02;
      elseif ll==7
        cfg.latency=[0 .45]+.07;
      elseif ll==9
        cfg.latency=[0 .45]+.5;
      end
    else
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[.1 .45];
      elseif ll==6
        cfg.latency=[.12 .47];
      elseif ll==7
        cfg.latency=[.17 .52];
      elseif ll==9
        cfg.latency=[.6 .95];
      end
    end
    
    % these test F-test (1-way ANOVA) for any difference of the 4 bins against each other.
    stat_tacPaud1_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1_pb{ll,:});
    %   stat_tacPaud2_pb=ft_timelockstatistics(cfg, grind_tacPaud2_pb{ll,:});
    stat_tacMSpN1_pb{ll}=ft_timelockstatistics(cfg, grind_tacMSpN1_pb{ll,:});
    %   stat_tacMSpN2_pb=ft_timelockstatistics(cfg, grind_tacMSpN2_pb{ll,:});
    stat_TPAmMSPN1_pb{ll}=ft_timelockstatistics(cfg, grind_TPAmMSPN1_pb{ll,:});
    
    if sleep
      stat_tacPaud1_ab{ll}=ft_timelockstatistics(cfg, grind_tacPaud1_ab{ll,:});
      stat_tacMSpN1_ab{ll}=ft_timelockstatistics(cfg, grind_tacMSpN1_ab{ll,:});
      stat_TPAmMSPN1_ab{ll}=ft_timelockstatistics(cfg, grind_TPAmMSPN1_ab{ll,:});
    else
      stat_tacPaud1IAF_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1IAF_pb{ll,:});
      %   stat_tacPaud2_pb=ft_timelockstatistics(cfg, grind_tacPaud2_pb{ll,:});
      stat_tacMSpN1IAF_pb{ll}=ft_timelockstatistics(cfg, grind_tacMSpN1IAF_pb{ll,:});
      %   stat_tacMSpN2_pb=ft_timelockstatistics(cfg, grind_tacMSpN2_pb{ll,:});
      stat_TPAmMSPN1IAF_pb{ll}=ft_timelockstatistics(cfg, grind_TPAmMSPN1IAF_pb{ll,:});
    end
    %   save([edir 'stat_pb_mult.mat'],'stat*')
    %     save([edir 'stat_pb_mult_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'stat*')
    save([edir 'stat_pb_mult_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag) '.mat'],'stat*');
  end
  
  cfg.statistic='depsamplesT';
  cfg.design=zeros(2,2*nsub);
  cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
  cfg.design(2,:)=[1:nsub 1:nsub];
  cfg.tail=0;
  for ll=soalist
    if sleep && ~sleepshortstatflag
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[0 1];
      elseif ll==6
        cfg.latency=[0 1]+.02;
      elseif ll==7
        cfg.latency=[0 1]+.07;
      elseif ll==9
        cfg.latency=[0 1]+.5;
      end
    elseif sleep && sleepshortstatflag
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[0 .45];
      elseif ll==6
        cfg.latency=[0 .45]+.02;
      elseif ll==7
        cfg.latency=[0 .45]+.07;
      elseif ll==9
        cfg.latency=[0 .45]+.5;
      end
    else
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[.1 .45];
      elseif ll==6
        cfg.latency=[.12 .47];
      elseif ll==7
        cfg.latency=[.17 .52];
      elseif ll==9
        cfg.latency=[.6 .95];
      end
    end
    % these test t-test for any difference of the A+T vs MS+N for each bin separately
    stat_TPAmMSPN1_peak_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1_pb{ll,1}, grind_tacMSpN1_pb{ll,1});
    stat_TPAmMSPN1_ptot_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1_pb{ll,2}, grind_tacMSpN1_pb{ll,2});
    stat_TPAmMSPN1_trgh_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1_pb{ll,3}, grind_tacMSpN1_pb{ll,3});
    stat_TPAmMSPN1_ttop_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1_pb{ll,4}, grind_tacMSpN1_pb{ll,4});
    
    if sleep
      stat_TPAmMSPN1_low_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1_ab{ll,1}, grind_tacMSpN1_ab{ll,1});
      stat_TPAmMSPN1_midl_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1_ab{ll,2}, grind_tacMSpN1_ab{ll,2});
      stat_TPAmMSPN1_midh_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1_ab{ll,3}, grind_tacMSpN1_ab{ll,3});
      stat_TPAmMSPN1_high_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1_ab{ll,4}, grind_tacMSpN1_ab{ll,4});
    else
      stat_TPAmMSPN1IAF_peak_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1IAF_pb{ll,1}, grind_tacMSpN1IAF_pb{ll,1});
      stat_TPAmMSPN1IAF_ptot_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1IAF_pb{ll,2}, grind_tacMSpN1IAF_pb{ll,2});
      stat_TPAmMSPN1IAF_trgh_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1IAF_pb{ll,3}, grind_tacMSpN1IAF_pb{ll,3});
      stat_TPAmMSPN1IAF_ttop_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1IAF_pb{ll,4}, grind_tacMSpN1IAF_pb{ll,4});
    end
  end
  
  %   save([edir 'stat_pb_mult.mat'],'stat*')
  %   save([edir 'stat_pb_mult_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'stat*')
  save([edir 'stat_pb_mult_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag) '.mat'],'stat*');
  clear grind*IAF* grind_*_*m*
  %   save([edir 'grind_pb_mult.mat'],'grind*')
  save([edir 'grind_pb_mult_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'grind*')
end

if stats1flag  % specific opposite phase or power contrasts (rather than full ANOVA)
  load eeg1010_neighb
  nsub=size(grind_TPAmMSPN1_pb{1,4}.individual,1);
  cfg=[];
  cfg.neighbours=neighbours;
  cfg.parameter='individual';
  cfg.method='montecarlo';
  cfg.numrandomization=1000;
  cfg.correctm='cluster';
  cfg.clusteralpha = 0.05;
  cfg.clusterstatistic = 'maxsum';
  cfg.minnbchan = 2;
  cfg.ivar=1;
  cfg.uvar=2;
  cfg.statistic='depsamplesT';
  cfg.design=zeros(2,2*nsub);
  cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
  cfg.design(2,:)=[1:nsub 1:nsub];
  cfg.tail=0;
  for ll=soalist
    if sleep==0
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[.1 .45];
      elseif ll==6
        cfg.latency=[.12 .47];
      elseif ll==7
        cfg.latency=[.17 .52];
      elseif ll==9
        cfg.latency=[.6 .95];
      end
    elseif sleep==1 && ~sleepshortstatflag
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[0 1];
      elseif ll==6
        cfg.latency=[0 1]+.02;
      elseif ll==7
        cfg.latency=[0 1]+.07;
      elseif ll==9
        cfg.latency=[0 1]+.5;
      end
    elseif sleep==1 && sleepshortstatflag
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[0 .45];
      elseif ll==6
        cfg.latency=[0 .45]+.02;
      elseif ll==7
        cfg.latency=[0 .45]+.07;
      elseif ll==9
        cfg.latency=[0 .45]+.5;
      end
    end
    stat_TPAmMSPN_ptotMttop_pb{ll}=ft_timelockstatistics(cfg, grind_TPAmMSPN1_pb{ll,2}, grind_TPAmMSPN1_pb{ll,4});
    stat_TPAmMSPN_peakMtrgh_pb{ll}=ft_timelockstatistics(cfg, grind_TPAmMSPN1_pb{ll,1}, grind_TPAmMSPN1_pb{ll,3});
    
    if sleep
      stat_TPAmMSPN_highMlow_ab{ll}=ft_timelockstatistics(cfg, grind_TPAmMSPN1_ab{ll,1}, grind_TPAmMSPN1_ab{ll,4});
    end
    
    if sleep==0 && (ll==3 || ll==5 || ll==7)  % conditions where previously there was significance
      aa=nan(2,2);
      try  aa(1,:)=[min(stat_TPAmMSPN1_peak_pb{ll}.time(find(mean(stat_TPAmMSPN1_peak_pb{ll}.mask,1)))) max(stat_TPAmMSPN1_peak_pb{ll}.time(find(mean(stat_TPAmMSPN1_peak_pb{ll}.mask,1))))]; catch end
      try  aa(2,:)=[min(stat_TPAmMSPN1_trgh_pb{ll}.time(find(mean(stat_TPAmMSPN1_trgh_pb{ll}.mask,1)))) max(stat_TPAmMSPN1_trgh_pb{ll}.time(find(mean(stat_TPAmMSPN1_trgh_pb{ll}.mask,1))))]; catch end
      cfg.latency=[min(aa(:,1)) max(aa(:,2))];
      stat_TPAmMSPN_peakMtrgh_pb_posthoctime{ll}=ft_timelockstatistics(cfg, grind_TPAmMSPN1_pb{ll,1}, grind_TPAmMSPN1_pb{ll,3});
      
      if ll==3 || ll==5
        aa=nan(2,2);
        try  aa(1,:)=[min(stat_TPAmMSPN1_ptot_pb{ll}.time(find(mean(stat_TPAmMSPN1_ptot_pb{ll}.mask,1)))) max(stat_TPAmMSPN1_ptot_pb{ll}.time(find(mean(stat_TPAmMSPN1_ptot_pb{ll}.mask,1))))]; catch end
        try  aa(2,:)=[min(stat_TPAmMSPN1_ttop_pb{ll}.time(find(mean(stat_TPAmMSPN1_ttop_pb{ll}.mask,1)))) max(stat_TPAmMSPN1_ttop_pb{ll}.time(find(mean(stat_TPAmMSPN1_ttop_pb{ll}.mask,1))))]; catch end
        cfg.latency=[min(aa(:,1)) max(aa(:,2))];
        stat_TPAmMSPN_ptotMttop_pb_posthoctime{ll}=ft_timelockstatistics(cfg, grind_TPAmMSPN1_pb{ll,2}, grind_TPAmMSPN1_pb{ll,4});
      end
    end
    
    
  end
  
  %   save([edir 'stat_pb_mult.mat'],'stat_TPAmMSPN_*M*','-append')
  %   save([edir 'stat_pb_mult_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'stat_TPAmM*M*','-append')
  save([edir 'stat_pb_mult_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag) '.mat'],'stat_TPAmM*M*','-append');
end

% going through results
for ll=soalist
  anysig(ll,1)=length(find(stat_tacPaud1_pb{ll}.mask));
  anysig(ll,2)=length(find(stat_tacMSpN1_pb{ll}.mask));
  anysig(ll,3)=length(find(stat_TPAmMSPN1_pb{ll}.mask));
  anysig(ll,4)=length(find(stat_tacPaud1IAF_pb{ll}.mask));
  anysig(ll,5)=length(find(stat_tacMSpN1IAF_pb{ll}.mask));
  anysig(ll,6)=length(find(stat_TPAmMSPN1IAF_pb{ll}.mask));
  
  anysig1(ll,1)=length(find(stat_TPAmMSPN1_peak_pb{ll}.mask));
  anysig1(ll,2)=length(find(stat_TPAmMSPN1_ptot_pb{ll}.mask));
  anysig1(ll,3)=length(find(stat_TPAmMSPN1_trgh_pb{ll}.mask));
  anysig1(ll,4)=length(find(stat_TPAmMSPN1_ttop_pb{ll}.mask));
  anysig1(ll,5)=length(find(stat_TPAmMSPN1IAF_peak_pb{ll}.mask));
  anysig1(ll,6)=length(find(stat_TPAmMSPN1IAF_ptot_pb{ll}.mask));
  anysig1(ll,7)=length(find(stat_TPAmMSPN1IAF_trgh_pb{ll}.mask));
  anysig1(ll,8)=length(find(stat_TPAmMSPN1IAF_ttop_pb{ll}.mask));
  
end

%% Plotting phase-dependent results (mult)

% clearvars -except sub *dir
chanplot{1}={'Fz' 'Cz' 'F1' 'F2' 'FC1' 'FC2' 'C1' 'C2'}; % frontocentral
chanplot{2}={'CP5' 'POz' 'Pz' 'P3' 'P4' 'C4' 'O1' 'O2' 'P7' 'PO7'}; % occipital
chanlabel{1}='Frontocentral electrodes';
chanlabel{2}='Occipital-parietal electrodes';
iterflag=1;
tt=3;

sleep=0;
sleepshortstatflag=0;
if sleep
  iter=11;
  ss=12;
  trialkc=-1;
else
  iter=27;
  ss=10;
  trialkc=-1;
end

if iterflag
  if sleep
    try
      load([edir 'stat_pb_mult_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag) '.mat']);
    catch
      load([edir 'stat_pb_mult_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat']);
    end
    load([edir 'grind_pb_mult_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat']);
  else
    load([edir 'stat_pb_mult.mat']);
    load([edir 'grind_pb_mult.mat']);
  end
end

soalist=[1 3 4 5 6 7 9];
timwinstatflag=0; % =1 for using only stat.time, =0 for using [-0.5 1.0];
if sleepshortstatflag
  timwin=[-0.5 1];
else
  timwin=[-0.5 1.5];
end
soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];

scalediff=1;

bluered=varycolor(10);
coloruse=parula(8);

% For ANOVA results
for pb=1:2 % pb=1 means phase results, pb=2 means power results
  close all
  clear tmp*
  for ll=soalist
    
    cfg=[];
    if timwinstatflag==1
      cfg.latency=[stat_TPAmMSPN1_peak_pb{ll}.time(1) stat_TPAmMSPN1_peak_pb{ll}.time(end)];
    elseif timwinstatflag==0
      cfg.latency=timwin;
      stattimwin=[stat_TPAmMSPN1_peak_pb{ll}.time(1) stat_TPAmMSPN1_peak_pb{ll}.time(end)];
    end
    cfg.channel=stat_TPAmMSPN1_peak_pb{ll}.label;
    
    switch pb
      case 1
        tmp1=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,1});
        tmp1.dimord='chan_time';
        tmp1.avg=squeeze(mean(tmp1.individual,1));
        tmp1=rmfield(tmp1,'individual');
        
        tmp5=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,2});
        tmp5.dimord='chan_time';
        tmp5.avg=squeeze(mean(tmp5.individual,1));
        tmp5=rmfield(tmp5,'individual');
        
        tmp8=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,3});
        tmp8.dimord='chan_time';
        tmp8.avg=squeeze(mean(tmp8.individual,1));
        tmp8=rmfield(tmp8,'individual');
        
        tmp12=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,4});
        tmp12.dimord='chan_time';
        tmp12.avg=squeeze(mean(tmp12.individual,1));
        tmp12=rmfield(tmp12,'individual');
        
        tmpmask=stat_TPAmMSPN1_pb{ll}.mask;
      case 2
        tmp1=ft_selectdata(cfg,grind_TPAmMSPN1_ab{ll,1});
        tmp1.dimord='chan_time';
        tmp1.avg=squeeze(mean(tmp1.individual,1));
        tmp1=rmfield(tmp1,'individual');
        
        tmp5=ft_selectdata(cfg,grind_TPAmMSPN1_ab{ll,2});
        tmp5.dimord='chan_time';
        tmp5.avg=squeeze(mean(tmp5.individual,1));
        tmp5=rmfield(tmp5,'individual');
        
        tmp8=ft_selectdata(cfg,grind_TPAmMSPN1_ab{ll,3});
        tmp8.dimord='chan_time';
        tmp8.avg=squeeze(mean(tmp8.individual,1));
        tmp8=rmfield(tmp8,'individual');
        
        tmp12=ft_selectdata(cfg,grind_TPAmMSPN1_ab{ll,4});
        tmp12.dimord='chan_time';
        tmp12.avg=squeeze(mean(tmp12.individual,1));
        tmp12=rmfield(tmp12,'individual');
        
        tmpmask=stat_TPAmMSPN1_ab{ll}.mask;
    end
    
    
    if timwinstatflag==0
      tmp1.mask=zeros(size(tmp1.avg,1),length(tmp1.time));
      tmp1.mask(:,dsearchn(tmp1.time',stattimwin(1)):dsearchn(tmp1.time',stattimwin(end)))=tmpmask;
      tmp5.mask=zeros(size(tmp5.avg,1),length(tmp5.time));
      tmp5.mask(:,dsearchn(tmp5.time',stattimwin(1)):dsearchn(tmp5.time',stattimwin(end)))=tmpmask;
      tmp8.mask=zeros(size(tmp8.avg,1),length(tmp8.time));
      tmp8.mask(:,dsearchn(tmp8.time',stattimwin(1)):dsearchn(tmp8.time',stattimwin(end)))=tmpmask;
      tmp12.mask=zeros(size(tmp12.avg,1),length(tmp12.time));
      tmp12.mask(:,dsearchn(tmp12.time',stattimwin(1)):dsearchn(tmp12.time',stattimwin(end)))=tmpmask;
    else
      tmpu1{bb}.mask=tmpmask
      tmpm10{bb}.mask=tmpmask;
      tmpd5{bb}.mask=tmpmask;
    end
    
    for cg=1:length(chanplot)
      cfg=[];
      cfg.parameter='avg';
      cfg.layout='elec1010.lay';
      cfg.ylim=[-5 8];
      cfg.linewidth=3;
      cfg.xlim=timwin;
      if cg>length(chanplot)
        cfg.channel=tmp5.label(any(tmp5.mask,2));
      else
        cfg.channel=chanplot{cg};
      end
      cfg.graphcolor=coloruse([pb:2:8],:);
      cfg.interactive='no';
      cfg.maskparameter='mask';
      cfg.maskstyle='box'; % default
      figure(ll+10*(cg+1))
      ft_singleplotER(cfg,tmp1,tmp5,tmp8,tmp12);
      hold on;plot(tmp1.time,0,'k');
      set(gca,'XTick',[-.5:.1:1])
      set(gca,'XTickLabel',{'-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'})
      set(gca,'FontSize',30)
      title([])
      plot([0 0],cfg.ylim,'Color',bluered(4,:),'LineWidth',6)
      plot([soades(ll) soades(ll)],cfg.ylim,'Color',bluered(9,:),'LineWidth',6)
      plot([stattimwin(1) stattimwin(1)],cfg.ylim,'k--','LineWidth',6)
      plot([stattimwin(2) stattimwin(2)],cfg.ylim,'k--','LineWidth',6)
      axis([-0.55 1.5 cfg.ylim(1) cfg.ylim(2)])
      if cg==3
        legend('Sum Unisensory','MultSens + Null','SumUnisens - MultsensNull')
      end
    end
    switch pb
      case 1
        print(ll+20,[fdir 'erp_tacPaud_MSpN_diff_FC_pb_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag) '.png'],'-dpng')
        print(ll+30,[fdir 'erp_tacPaud_MSpN_diff_OP_pb_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag) '.png'],'-dpng')
      case 2
        print(ll+20,[fdir 'erp_tacPaud_MSpN_diff_FC_ab_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag) '.png'],'-dpng')
        print(ll+30,[fdir 'erp_tacPaud_MSpN_diff_OP_ab_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag) '.png'],'-dpng')
    end
    
    
    if any(tmpmask(:))
      masktime=find(any(tmp5.mask,1));
      cfg=[];
      cfg.parameter='avg';
      cfg.layout='elec1010.lay';
      cfg.maskalpha=0.5;
      cfg.zlim=[-5 5];
      cfg.highlight='on';
      cfg.highlightsize=12;
      cfg.xlim=[tmp5.time(masktime(1)) tmp5.time(masktime(end))];
      cfg.comment='no';
      sigchannels=tmp5.label(find(ceil(mean(tmp5.mask(:,dsearchn(tmp5.time',cfg.xlim(1)):dsearchn(tmp5.time',cfg.xlim(2))),2))));
      cfg.highlightchannel=sigchannels;
      figure(1000*bb+ll);
      ft_topoplotER(cfg,tmp1);
      switch pb
        case 1
          print(1000*bb+ll,[fdir 'erp_topoBin1_pb_' num2str(ll) num2str(tt) num2str(ss) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        case 2
          print(1000*bb+ll,[fdir 'erp_topoBin1_ab_' num2str(ll) num2str(tt) num2str(ss) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
      end
      figure(1000*bb+10+ll);
      ft_topoplotER(cfg,tmp5);
      switch pb
        case 1
          print(1000*bb+10+ll,[fdir 'erp_topoBin2_pb_' num2str(ll) num2str(tt) num2str(ss) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        case 2
          print(1000*bb+10+ll,[fdir 'erp_topoBin2_ab_' num2str(ll) num2str(tt) num2str(ss) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
      end
      figure(1000*bb+20+ll);
      ft_topoplotER(cfg,tmp8);
      switch pb
        case 1
          print(1000*bb+20+ll,[fdir 'erp_topoBin3_pb_' num2str(ll) num2str(tt) num2str(ss) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        case 2
          print(1000*bb+20+ll,[fdir 'erp_topoBin3_ab_' num2str(ll) num2str(tt) num2str(ss) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
      end
      figure(1000*bb+30+ll);
      ft_topoplotER(cfg,tmp12);
      switch pb
        case 1
          print(1000*bb+20+ll,[fdir 'erp_topoBin4_pb_' num2str(ll) num2str(tt) num2str(ss) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        case 2
          print(1000*bb+20+ll,[fdir 'erp_topoBin4_ab_' num2str(ll) num2str(tt) num2str(ss) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
      end
    end
    
    
  end % ll
end




% for individual bin results
for pb=1:2
  close all
  clear tmp*
  for ll=soalist
    for bb=1:4
      
      cfg=[];
      if timwinstatflag==1
        cfg.latency=[stat_TPAmMSPN1_peak_pb{ll}.time(1) stat_TPAmMSPN1_peak_pb{ll}.time(end)];
      elseif timwinstatflag==0
        cfg.latency=timwin;
        stattimwin=[stat_TPAmMSPN1_peak_pb{ll}.time(1) stat_TPAmMSPN1_peak_pb{ll}.time(end)];
      end
      cfg.channel=stat_TPAmMSPN1_peak_pb{ll}.label;
      
      switch pb
        case 1
          tmpu1{bb}=ft_selectdata(cfg,grind_tacPaud1_pb{ll,bb});
          tmpu1{bb}.dimord='chan_time';
          tmpu1{bb}.avg=squeeze(mean(tmpu1{bb}.individual,1));
          tmpu1{bb}=rmfield(tmpu1{bb},'individual');
          
          tmpm10{bb}=ft_selectdata(cfg,grind_tacMSpN1_pb{ll,bb});
          tmpm10{bb}.dimord='chan_time';
          tmpm10{bb}.avg=squeeze(mean(tmpm10{bb}.individual,1));
          tmpm10{bb}=rmfield(tmpm10{bb},'individual');
          
          tmpd5{bb}=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,bb});
          tmpd5{bb}.dimord='chan_time';
          tmpd5{bb}.avg=scalediff*squeeze(mean(tmpd5{bb}.individual,1));
          tmpd5{bb}=rmfield(tmpd5{bb},'individual');
          
          switch bb
            case 1
              tmpmask=stat_TPAmMSPN1_peak_pb{ll}.mask;
            case 2
              tmpmask=stat_TPAmMSPN1_ptot_pb{ll}.mask;
            case 3
              tmpmask=stat_TPAmMSPN1_trgh_pb{ll}.mask;
            case 4
              tmpmask=stat_TPAmMSPN1_ttop_pb{ll}.mask;
          end
        case 2
          tmpu1{bb}=ft_selectdata(cfg,grind_tacPaud1_ab{ll,bb});
          tmpu1{bb}.dimord='chan_time';
          tmpu1{bb}.avg=squeeze(mean(tmpu1{bb}.individual,1));
          tmpu1{bb}=rmfield(tmpu1{bb},'individual');
          
          tmpm10{bb}=ft_selectdata(cfg,grind_tacMSpN1_ab{ll,bb});
          tmpm10{bb}.dimord='chan_time';
          tmpm10{bb}.avg=squeeze(mean(tmpm10{bb}.individual,1));
          tmpm10{bb}=rmfield(tmpm10{bb},'individual');
          
          tmpd5{bb}=ft_selectdata(cfg,grind_TPAmMSPN1_ab{ll,bb});
          tmpd5{bb}.dimord='chan_time';
          tmpd5{bb}.avg=scalediff*squeeze(mean(tmpd5{bb}.individual,1));
          tmpd5{bb}=rmfield(tmpd5{bb},'individual');
          
          switch bb
            case 1
              tmpmask=stat_TPAmMSPN1_low_pb{ll}.mask;
            case 2
              tmpmask=stat_TPAmMSPN1_midl_pb{ll}.mask;
            case 3
              tmpmask=stat_TPAmMSPN1_midh_pb{ll}.mask;
            case 4
              tmpmask=stat_TPAmMSPN1_high_pb{ll}.mask;
          end
      end
      
      if timwinstatflag==0
        tmpu1{bb}.mask=zeros(size(tmpu1{bb}.avg,1),length(tmpu1{bb}.time));
        tmpu1{bb}.mask(:,dsearchn(tmpu1{bb}.time',stattimwin(1)):dsearchn(tmpu1{bb}.time',stattimwin(end)))=tmpmask;
        tmpm10{bb}.mask=zeros(size(tmpm10{bb}.avg,1),length(tmpm10{bb}.time));
        tmpm10{bb}.mask(:,dsearchn(tmpm10{bb}.time',stattimwin(1)):dsearchn(tmpm10{bb}.time',stattimwin(end)))=tmpmask;
        tmpd5{bb}.mask=zeros(size(tmpd5{bb}.avg,1),length(tmpd5{bb}.time));
        tmpd5{bb}.mask(:,dsearchn(tmpd5{bb}.time',stattimwin(1)):dsearchn(tmpd5{bb}.time',stattimwin(end)))=tmpmask;
      else
        tmpu1{bb}.mask=tmpmask
        tmpm10{bb}.mask=tmpmask;
        tmpd5{bb}.mask=tmpmask;
      end
      
      for cg=1:length(chanplot)
        cfg=[];
        cfg.parameter='avg';
        cfg.layout='elec1010.lay';
        cfg.ylim=[-5 8];
        cfg.linewidth=3;
        cfg.xlim=timwin;
        if cg>length(chanplot)
          cfg.channel=tmpd5{bb}.label(any(tmpd5{bb}.mask,2));
        else
          cfg.channel=chanplot{cg};
        end
        %         cfg.graphcolor=[bluered(1,:); bluered(10,:); coloruse(2*(bb-1)+pb,:)]; % color for different phase bins
        cfg.graphcolor=[bluered(1,:); bluered(10,:); coloruse(pb,:)];
        cfg.interactive='no';
        cfg.maskparameter='mask';
        cfg.maskstyle='box'; % default
        figure(100*bb+ll+10*(cg+1))
        ft_singleplotER(cfg,tmpu1{bb},tmpm10{bb},tmpd5{bb});
        hold on;plot(tmpu1{bb}.time,0,'k');
        set(gca,'XTick',[-.5:.1:1])
        set(gca,'XTickLabel',{'-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'})
        set(gca,'FontSize',30)
        title([])
        plot([0 0],cfg.ylim,'Color',bluered(4,:),'LineWidth',6)
        plot([soades(ll) soades(ll)],cfg.ylim,'Color',bluered(9,:),'LineWidth',6)
        plot([stattimwin(1) stattimwin(1)],cfg.ylim,'k--','LineWidth',6)
        plot([stattimwin(2) stattimwin(2)],cfg.ylim,'k--','LineWidth',6)
        axis([-0.55 1.5 cfg.ylim(1) cfg.ylim(2)])
        if cg==3
          legend('Sum Unisensory','MultSens + Null','SumUnisens - MultsensNull')
        end
      end
      switch pb
        case 1
          print(100*bb+ll+20,[fdir 'erp_tacPaud_MSpN_diff_FC_pb_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_bin' num2str(bb)  '_shortstat' num2str(sleepshortstatflag) '.png'],'-dpng')
          print(100*bb+ll+30,[fdir 'erp_tacPaud_MSpN_diff_OP_pb_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_bin' num2str(bb)  '_shortstat' num2str(sleepshortstatflag) '.png'],'-dpng')
        case 2
          print(100*bb+ll+20,[fdir 'erp_tacPaud_MSpN_diff_FC_ab_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_bin' num2str(bb)  '_shortstat' num2str(sleepshortstatflag) '.png'],'-dpng')
          print(100*bb+ll+30,[fdir 'erp_tacPaud_MSpN_diff_OP_ab_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_bin' num2str(bb)  '_shortstat' num2str(sleepshortstatflag) '.png'],'-dpng')
      end
      
      if any(tmpmask(:))
        masktime=find(any(tmpd5{bb}.mask,1));
        cfg=[];
        cfg.parameter='avg';
        cfg.layout='elec1010.lay';
        cfg.maskalpha=0.5;
        cfg.zlim=[-5 5];
        cfg.highlight='on';
        cfg.highlightsize=12;
        cfg.xlim=[tmpd5{bb}.time(masktime(1)) tmpd5{bb}.time(masktime(end))];
        cfg.comment='no';
        sigchannels=tmpd5{bb}.label(find(ceil(mean(tmpd5{bb}.mask(:,dsearchn(tmpd5{bb}.time',cfg.xlim(1)):dsearchn(tmpd5{bb}.time',cfg.xlim(2))),2))));
        cfg.highlightchannel=sigchannels;
        figure(1000*bb+ll);
        ft_topoplotER(cfg,tmpu1{bb});
        switch pb
          case 1
            print(1000*bb+ll,[fdir 'erp_topoU_pb_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
          case 2
            print(1000*bb+ll,[fdir 'erp_topoU_ab_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        end
        figure(1000*bb+10+ll);
        ft_topoplotER(cfg,tmpm10{bb});
        switch pb
          case 1
            print(1000*bb+10+ll,[fdir 'erp_topoM_pb_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
          case 2
            print(1000*bb+10+ll,[fdir 'erp_topoM_ab_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        end
        figure(1000*bb+20+ll);
        ft_topoplotER(cfg,tmpd5{bb});
        switch pb
          case 1
            print(1000*bb+20+ll,[fdir 'erp_topoDiff_pb_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
          case 2
            print(1000*bb+20+ll,[fdir 'erp_topoDiff_ab_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        end
      end
      
    end % bb
  end % ll
end % pb


close all
clear tmp*
% Now, for specific ptotMttop and peakMtrgh and highMlow

for ll=soalist
  
  cfg=[];
  cfg.parameter='individual';
  cfg.operation='subtract';
  %   grind_TPAmMSPN_4m2_pb{ll}=ft_math(cfg,grind_TPAmMSPN1_pb{ll,4},grind_TPAmMSPN1_pb{ll,2});
  %   grave_TPAmMSPN_4m2_pb{ll}=grind_TPAmMSPN_4m2_pb{ll};
  %   grave_TPAmMSPN_4m2_pb{ll}.avg=squeeze(mean(grave_TPAmMSPN_4m2_pb{ll}.individual,1));
  %   grave_TPAmMSPN_4m2_pb{ll}=rmfield(grave_TPAmMSPN_4m2_pb{ll},'individual');
  %   grave_TPAmMSPN_4m2_pb{ll}.dimord='chan_time';
  %
  %   grind_TPAmMSPN_3m1_pb{ll}=ft_math(cfg,grind_TPAmMSPN1_pb{ll,3},grind_TPAmMSPN1_pb{ll,1});
  %   grave_TPAmMSPN_3m1_pb{ll}=grind_TPAmMSPN_3m1_pb{ll};
  %   grave_TPAmMSPN_3m1_pb{ll}.avg=squeeze(mean(grave_TPAmMSPN_3m1_pb{ll}.individual,1));
  %   grave_TPAmMSPN_3m1_pb{ll}=rmfield(grave_TPAmMSPN_3m1_pb{ll},'individual');
  %   grave_TPAmMSPN_3m1_pb{ll}.dimord='chan_time';
  
  grind_TPAmMSPN_2m4_pb{ll}=ft_math(cfg,grind_TPAmMSPN1_pb{ll,2},grind_TPAmMSPN1_pb{ll,4});
  grave_TPAmMSPN_2m4_pb{ll}=grind_TPAmMSPN_2m4_pb{ll};
  grave_TPAmMSPN_2m4_pb{ll}.avg=squeeze(mean(grave_TPAmMSPN_2m4_pb{ll}.individual,1));
  grave_TPAmMSPN_2m4_pb{ll}=rmfield(grave_TPAmMSPN_2m4_pb{ll},'individual');
  grave_TPAmMSPN_2m4_pb{ll}.dimord='chan_time';
  
  grind_TPAmMSPN_1m3_pb{ll}=ft_math(cfg,grind_TPAmMSPN1_pb{ll,1},grind_TPAmMSPN1_pb{ll,3});
  grave_TPAmMSPN_1m3_pb{ll}=grind_TPAmMSPN_1m3_pb{ll};
  grave_TPAmMSPN_1m3_pb{ll}.avg=squeeze(mean(grave_TPAmMSPN_1m3_pb{ll}.individual,1));
  grave_TPAmMSPN_1m3_pb{ll}=rmfield(grave_TPAmMSPN_1m3_pb{ll},'individual');
  grave_TPAmMSPN_1m3_pb{ll}.dimord='chan_time';
  
  if exist('grind_TPAmMSPN1_ab','var')
    grind_TPAmMSPN_4m1_ab{ll}=ft_math(cfg,grind_TPAmMSPN1_ab{ll,4},grind_TPAmMSPN1_ab{ll,1});
    grave_TPAmMSPN_4m1_ab{ll}=grind_TPAmMSPN_4m1_ab{ll};
    grave_TPAmMSPN_4m1_ab{ll}.avg=squeeze(mean(grave_TPAmMSPN_4m1_ab{ll}.individual,1));
    grave_TPAmMSPN_4m1_ab{ll}=rmfield(grave_TPAmMSPN_4m1_ab{ll},'individual');
    grave_TPAmMSPN_4m1_ab{ll}.dimord='chan_time';
    kkmax=3;
  else
    kkmax=2;
  end
  
  for kk=1:kkmax
    cfg=[];
    if timwinstatflag==1
      cfg.latency=[stat_TPAmMSPN1_peak_pb{ll}.time(1) stat_TPAmMSPN1_peak_pb{ll}.time(end)];
    elseif timwinstatflag==0
      cfg.latency=timwin;
      stattimwin=[stat_TPAmMSPN1_peak_pb{ll}.time(1) stat_TPAmMSPN1_peak_pb{ll}.time(end)];
    end
    cfg.channel=stat_TPAmMSPN1_peak_pb{ll}.label;
    
    switch kk
      case 1
        tmpu1=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,2});
        tmpu1.dimord='chan_time';
        tmpu1.avg=squeeze(mean(tmpu1.individual,1));
        tmpu1=rmfield(tmpu1,'individual');
        
        tmpm10=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,4});
        tmpm10.dimord='chan_time';
        tmpm10.avg=squeeze(mean(tmpm10.individual,1));
        tmpm10=rmfield(tmpm10,'individual');
        
        tmpd5=ft_selectdata(cfg,grind_TPAmMSPN_2m4_pb{ll});
        tmpd5.dimord='chan_time';
        tmpd5.avg=scalediff*squeeze(mean(tmpd5.individual,1));
        tmpd5=rmfield(tmpd5,'individual');
      case 2
        tmpu1=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,1});
        tmpu1.dimord='chan_time';
        tmpu1.avg=squeeze(mean(tmpu1.individual,1));
        tmpu1=rmfield(tmpu1,'individual');
        
        tmpm10=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,3});
        tmpm10.dimord='chan_time';
        tmpm10.avg=squeeze(mean(tmpm10.individual,1));
        tmpm10=rmfield(tmpm10,'individual');
        
        tmpd5=ft_selectdata(cfg,grind_TPAmMSPN_1m3_pb{ll});
        tmpd5.dimord='chan_time';
        tmpd5.avg=scalediff*squeeze(mean(tmpd5.individual,1));
        tmpd5=rmfield(tmpd5,'individual');
      case 3
        tmpu1=ft_selectdata(cfg,grind_TPAmMSPN1_ab{ll,1});
        tmpu1.dimord='chan_time';
        tmpu1.avg=squeeze(mean(tmpu1.individual,1));
        tmpu1=rmfield(tmpu1,'individual');
        
        tmpm10=ft_selectdata(cfg,grind_TPAmMSPN1_ab{ll,4});
        tmpm10.dimord='chan_time';
        tmpm10.avg=squeeze(mean(tmpm10.individual,1));
        tmpm10=rmfield(tmpm10,'individual');
        
        tmpd5=ft_selectdata(cfg,grind_TPAmMSPN_4m1_ab{ll});
        tmpd5.dimord='chan_time';
        tmpd5.avg=scalediff*squeeze(mean(tmpd5.individual,1));
        tmpd5=rmfield(tmpd5,'individual');
    end
    
    switch kk
      case 1
        tmpmask=stat_TPAmMSPN_ptotMttop_pb{ll}.mask;
      case 2
        tmpmask=stat_TPAmMSPN_peakMtrgh_pb{ll}.mask;
      case 3
        if trialkc==-1
          tmpmask=stat_TPAmMSPN_lowMhigh_ab{ll}.mask;
        elseif trialkc==0
          tmpmask=stat_TPAmMSPN_highMlow_ab{ll}.mask;
        end
    end
    
    if timwinstatflag==0
      tmpu1.mask=zeros(size(tmpu1.avg,1),length(tmpu1.time));
      tmpu1.mask(:,dsearchn(tmpu1.time',stattimwin(1)):dsearchn(tmpu1.time',stattimwin(end)))=tmpmask;
      tmpm10.mask=zeros(size(tmpm10.avg,1),length(tmpm10.time));
      tmpm10.mask(:,dsearchn(tmpm10.time',stattimwin(1)):dsearchn(tmpm10.time',stattimwin(end)))=tmpmask;
      tmpd5.mask=zeros(size(tmpd5.avg,1),length(tmpd5.time));
      tmpd5.mask(:,dsearchn(tmpd5.time',stattimwin(1)):dsearchn(tmpd5.time',stattimwin(end)))=tmpmask;
    else
      tmpu1.mask=tmpmask;
      tmpm10.mask=tmpmask;
      tmpd5.mask=tmpmask;
    end
    
    for cg=1:length(chanplot)
      cfg=[];
      cfg.parameter='avg';
      cfg.layout='elec1010.lay';
      cfg.ylim=[-5 8];
      cfg.linewidth=3;
      cfg.xlim=timwin;
      if cg>length(chanplot)
        cfg.channel=tmpd5.label(any(tmpd5.mask,2));
      else
        cfg.channel=chanplot{cg};
      end
      switch kk
        case 1
          cfg.graphcolor=[coloruse(3,:); coloruse(7,:);  bluered(10,:)]
        case 2
          cfg.graphcolor=[coloruse(1,:); coloruse(5,:);  bluered(10,:)]
        case 3
          cfg.graphcolor=[coloruse(2,:); coloruse(8,:);  bluered(10,:)]
      end
      cfg.interactive='no';
      cfg.maskparameter='mask';
      cfg.maskstyle='box'; % default
      figure(500+ll+10*(cg+1))
      ft_singleplotER(cfg,tmpu1,tmpm10,tmpd5);
      hold on;plot(tmpu1.time,0,'k');
      set(gca,'XTick',[-.5:.1:1])
      set(gca,'XTickLabel',{'-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'})
      set(gca,'FontSize',30)
      title([])
      plot([0 0],cfg.ylim,'Color',bluered(4,:),'LineWidth',6)
      plot([soades(ll) soades(ll)],cfg.ylim,'Color',bluered(9,:),'LineWidth',6)
      plot([stattimwin(1) stattimwin(1)],cfg.ylim,'k--','LineWidth',6)
      plot([stattimwin(2) stattimwin(2)],cfg.ylim,'k--','LineWidth',6)
      axis([-0.55 1.5 cfg.ylim(1) cfg.ylim(2)])
      if cg==3
        legend('Sum Unisensory','MultSens + Null','SumUnisens - MultsensNull')
      end
    end
    %   print(30+ll,[fdir 'erp_tacPaud_MSpN_diff_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    switch kk
      case 1
        print(500+ll+20,[fdir 'erp_tacPaud_MSpN_diff_FC_pb_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag)  '_bindiff24.png'],'-dpng')
        print(500+ll+30,[fdir 'erp_tacPaud_MSpN_diff_OP_pb_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag)  '_bindiff24.png'],'-dpng')
      case 2
        print(500+ll+20,[fdir 'erp_tacPaud_MSpN_diff_FC_pb_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag)  '_bindiff13.png'],'-dpng')
        print(500+ll+30,[fdir 'erp_tacPaud_MSpN_diff_OP_pb_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag)  '_bindiff13.png'],'-dpng')
      case 3
        print(500+ll+20,[fdir 'erp_tacPaud_MSpN_diff_FC_ab_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag)  '_bindiff41.png'],'-dpng')
        print(500+ll+30,[fdir 'erp_tacPaud_MSpN_diff_OP_ab_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag)  '_bindiff41.png'],'-dpng')
    end
    
    if any(tmpmask(:))
      masktime=find(any(tmpm10.mask,1));
      cfg=[];
      cfg.parameter='avg';
      cfg.layout='elec1010.lay';
      cfg.maskalpha=0.5;
      cfg.zlim=[-5 5];
      cfg.highlight='on';
      cfg.highlightsize=12;
      cfg.xlim=[tmpm10.time(masktime(1)) tmpm10.time(masktime(end))];
      cfg.comment='no';
      sigchannels=tmpm10.label(find(ceil(mean(tmpm10.mask(:,dsearchn(tmpm10.time',cfg.xlim(1)):dsearchn(tmpm10.time',cfg.xlim(2))),2))));
      cfg.highlightchannel=sigchannels;
      figure(1000*bb+ll);
      ft_topoplotER(cfg,tmpu1);
      switch kk
        case 1
          print(5000+ll,[fdir 'erp_topoU_pb_2m4_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        case 2
          print(5000+ll,[fdir 'erp_topoU_pb_1m3_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        case 3
          print(5000+ll,[fdir 'erp_topoU_ab_4m1_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
      end
      figure(1000*bb+10+ll);
      ft_topoplotER(cfg,tmpm10);
      switch kk
        case 1
          print(5000+10+ll,[fdir 'erp_topoM_pb_2m4_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        case 2
          print(5000+10+ll,[fdir 'erp_topoM_pb_1m3_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        case 3
          print(5000+10+ll,[fdir 'erp_topoM_ab_4m1_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
      end
      figure(1000*bb+20+ll);
      ft_topoplotER(cfg,tmpd5);
      switch kk
        case 1
          print(5000+20+ll,[fdir 'erp_topoDiff_pb_2m4_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        case 2
          print(5000+20+ll,[fdir 'erp_topoDiff_pb_1m3_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        case 3
          print(5000+20+ll,[fdir 'erp_topoDiff_ab_4m1_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
      end
    end
    
  end % bb
end % ll


if sleep==0
  pb=1; ll=5; % 1vs3 and 2vs4 (ll=5 only significant finding in phase binning)
  cfg=[];
  if timwinstatflag==1
    cfg.latency=[stat_TPAmMSPN_peakMtrgh_pb_posthoctime{ll}.time(1) stat_TPAmMSPN_peakMtrgh_pb_posthoctime{ll}.time(end)];
  elseif timwinstatflag==0
    cfg.latency=timwin;
    stattimwin=[stat_TPAmMSPN_peakMtrgh_pb_posthoctime{ll}.time(1) stat_TPAmMSPN_peakMtrgh_pb_posthoctime{ll}.time(end)];
  end
  cfg.channel=stat_TPAmMSPN_peakMtrgh_pb_posthoctime{ll}.label;
  
  
  tmpu1=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,1});
  tmpu1.dimord='chan_time';
  tmpu1.avg=squeeze(mean(tmpu1.individual,1));
  tmpu1=rmfield(tmpu1,'individual');
  
  tmpm10=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,3});
  tmpm10.dimord='chan_time';
  tmpm10.avg=squeeze(mean(tmpm10.individual,1));
  tmpm10=rmfield(tmpm10,'individual');
  
  tmpd5=ft_selectdata(cfg,grind_TPAmMSPN_1m3_pb{ll});
  tmpd5.dimord='chan_time';
  tmpd5.avg=scalediff*squeeze(mean(tmpd5.individual,1));
  tmpd5=rmfield(tmpd5,'individual');
  
  
  tmpmask=stat_TPAmMSPN_peakMtrgh_pb_posthoctime{ll}.mask;
  if timwinstatflag==0
    tmpu1.mask=zeros(size(tmpu1.avg,1),length(tmpu1.time));
    tmpu1.mask(:,dsearchn(tmpu1.time',stattimwin(1)):dsearchn(tmpu1.time',stattimwin(end)))=tmpmask;
    tmpm10.mask=zeros(size(tmpm10.avg,1),length(tmpm10.time));
    tmpm10.mask(:,dsearchn(tmpm10.time',stattimwin(1)):dsearchn(tmpm10.time',stattimwin(end)))=tmpmask;
    tmpd5.mask=zeros(size(tmpd5.avg,1),length(tmpd5.time));
    tmpd5.mask(:,dsearchn(tmpd5.time',stattimwin(1)):dsearchn(tmpd5.time',stattimwin(end)))=tmpmask;
  else
    tmpu1{bb}.mask=tmpmask
    tmpm10{bb}.mask=tmpmask;
    tmpd5{bb}.mask=tmpmask;
  end
  
  for cg=1:length(chanplot)
    cfg=[];
    cfg.parameter='avg';
    cfg.layout='elec1010.lay';
    cfg.ylim=[-5 8];
    cfg.linewidth=3;
    cfg.xlim=timwin;
    if cg>length(chanplot)
      cfg.channel=tmpd5.label(any(tmpd5.mask,2));
    else
      cfg.channel=chanplot{cg};
    end
    cfg.graphcolor=coloruse([1 1 4],:);
    cfg.linestyle={'-', '--', '-'};
    cfg.interactive='no';
    cfg.maskparameter='mask';
    cfg.maskstyle='box'; % default
    figure(ll+10*(cg+1))
    ft_singleplotER(cfg,tmpu1,tmpm10,tmpd5);
    hold on;plot(tmpu1.time,0,'k');
    set(gca,'XTick',[-.5:.1:1])
    set(gca,'XTickLabel',{'-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'})
    set(gca,'FontSize',30)
    title([])
    plot([stattimwin(1) stattimwin(1)],cfg.ylim,'k--','LineWidth',6)
    plot([stattimwin(2) stattimwin(2)],cfg.ylim,'k--','LineWidth',6)
    plot([0 0],cfg.ylim,'Color',coloruse(4,:),'LineWidth',6)
    plot([soades(ll) soades(ll)],cfg.ylim,'Color',bluered(9,:),'LineWidth',6)
    axis([-0.55 1 cfg.ylim(1) cfg.ylim(2)])
    if cg==3
      legend('Sum Unisensory','MultSens + Null','SumUnisens - MultsensNull')
    end
  end
  %   print(ll+20,[fdir 'erp_tacPaud_MSpN_diff_PeakVsTrgh_FC_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  %   print(ll+30,[fdir 'erp_tacPaud_MSpN_diff_PeakVsTrgh_OP_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  print(ll+20,[fdir 'erp_tacPaud_MSpN_diff_FC_pb_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) 'bindiff13_posthoctime.png'],'-dpng')
  print(ll+30,[fdir 'erp_tacPaud_MSpN_diff_OP_pb_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) 'bindiff13_posthoctime.png'],'-dpng')
  
  if any(tmpmask(:))
    masktime=find(any(tmpd5.mask,1));
    cfg=[];
    cfg.parameter='avg';
    cfg.layout='elec1010.lay';
    cfg.maskalpha=0.5;
    cfg.zlim=[-5 5];
    cfg.highlight='on';
    cfg.highlightsize=12;
    cfg.xlim=[tmpd5.time(masktime(1)) tmpd5.time(masktime(end))];
    cfg.comment='no';
    sigchannels=tmpd5.label(find(ceil(mean(tmpd5.mask(:,dsearchn(tmpd5.time',cfg.xlim(1)):dsearchn(tmpd5.time',cfg.xlim(2))),2))));
    cfg.highlightchannel=sigchannels;
    figure(100+ll);
    ft_topoplotER(cfg,tmpu1);
    print(100+ll,[fdir 'erp_topoTPAmMSpN_bin1_' num2str(ll) num2str(tt) num2str(ss) '_posthoctime.png'],'-dpng')
    figure(100+10+ll);
    ft_topoplotER(cfg,tmpm10);
    print(100+10+ll,[fdir 'erp_topoTPAmMSpN_bin3_' num2str(ll) num2str(tt) num2str(ss) '_posthoctime.png'],'-dpng')
    figure(100+20+ll);
    ft_topoplotER(cfg,tmpd5);
    print(100+20+ll,[fdir 'erp_topoTPAmMSpN_bindiff13_' num2str(ll) num2str(tt) num2str(ss) '_posthoctime.png'],'-dpng')
  end
end % sleep

%% Plotting phase-dependent results (Uni vs Nul)

% clearvars -except sub *dir
chanplot{1}={'Fz' 'Cz' 'F1' 'F2' 'FC1' 'FC2' 'C1' 'C2'}; % frontocentral
chanplot{2}={'CP5' 'POz' 'Pz' 'P3' 'P4' 'C4' 'O1' 'O2' 'P7' 'PO7'}; % occipital
chanlabel{1}='Frontocentral electrodes';
chanlabel{2}='Occipital-parietal electrodes';
iterflag=1;
sleep=0;
tt=3;
if sleep
  iter=11;
  ss=12;
  trialkc=0;
else
  iter=27;
  ss=10;
  trialkc=-1;
end

if iterflag
  if sleep
  else
    load([edir 'stat_pb_uninul.mat']);
    load([edir 'grind_pb_uninul.mat']);
  end
end

soalist=[1 3 4 5 6 7 9];
timwinstatflag=0; % =1 for using only stat.time, =0 for using [-0.5 1.0];
timwin=[-0.5 1];
soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];

scalediff=2;

coloruse=varycolor(10);
% optimsied to be maximally apart and distinct
% 1  TacPAud
% 2  Nul
% 3  MS_synch
% 4  Tac
% 5  diff_TPAMSPN
% 6  diff_synchAsynch
% 7  MS
% 8  MS_asynch
% 9  Aud
% 10 MSpN

close all
clear tmp*

condname={'AT500' '' 'AT70' 'AT20' 'AT0' 'TA20' 'TA70' '' 'TA500' '' 'AUD' 'TAC'};
for ll=[soalist 11 12]
  
  cfg=[];
  cfg.channel=grind_tacMnul_pb{1}.label;
  
  switch ll
    case {1, 3, 4, 5, 6, 7, 9}
      if timwinstatflag==1
        cfg.latency=[stat_msMnul1_pb{ll}.time(1) stat_msMnul1_pb{ll}.time(end)];
      elseif timwinstatflag==0
        cfg.latency=timwin;
        stattimwin=[stat_msMnul1_pb{ll}.time(1) stat_msMnul1_pb{ll}.time(end)];
      end
      tmpu1=ft_selectdata(cfg,grind_msMnul1_pb{ll,1});
      tmpu1.dimord='chan_time';
      tmpu1.avg=squeeze(mean(tmpu1.individual,1));
      tmpu1=rmfield(tmpu1,'individual');
      
      tmpm10=ft_selectdata(cfg,grind_msMnul1_pb{ll,3});
      tmpm10.dimord='chan_time';
      tmpm10.avg=squeeze(mean(tmpm10.individual,1));
      tmpm10=rmfield(tmpm10,'individual');
      
      tmpd4=ft_selectdata(cfg,grind_msMnul1_pb{ll,2});
      tmpd4.dimord='chan_time';
      tmpd4.avg=squeeze(mean(tmpd4.individual,1));
      tmpd4=rmfield(tmpd4,'individual');
      
      tmpd7=ft_selectdata(cfg,grind_msMnul1_pb{ll,4});
      tmpd7.dimord='chan_time';
      tmpd7.avg=squeeze(mean(tmpd7.individual,1));
      tmpd7=rmfield(tmpd7,'individual');
      
      tmpmask=stat_msMnul1_pb{ll}.mask;
      
    case 11
      if timwinstatflag==1
        cfg.latency=[stat_audMnul_pb.time(1) stat_audMnul_pb.time(end)];
      elseif timwinstatflag==0
        cfg.latency=timwin;
        stattimwin=[stat_audMnul_pb.time(1) stat_audMnul_pb.time(end)];
      end
      tmpu1=ft_selectdata(cfg,grind_audMnul_pb{1});
      tmpu1.dimord='chan_time';
      tmpu1.avg=squeeze(mean(tmpu1.individual,1));
      tmpu1=rmfield(tmpu1,'individual');
      
      tmpm10=ft_selectdata(cfg,grind_audMnul_pb{3});
      tmpm10.dimord='chan_time';
      tmpm10.avg=squeeze(mean(tmpm10.individual,1));
      tmpm10=rmfield(tmpm10,'individual');
      
      tmpd4=ft_selectdata(cfg,grind_audMnul_pb{2});
      tmpd4.dimord='chan_time';
      tmpd4.avg=squeeze(mean(tmpd4.individual,1));
      tmpd4=rmfield(tmpd4,'individual');
      
      tmpd7=ft_selectdata(cfg,grind_audMnul_pb{4});
      tmpd7.dimord='chan_time';
      tmpd7.avg=squeeze(mean(tmpd7.individual,1));
      tmpd7=rmfield(tmpd7,'individual');
      
      tmpmask=stat_audMnul_pb.mask;
      
    case 12
      if timwinstatflag==1
        cfg.latency=[stat_tacMnul_pb.time(1) stat_tacMnul_pb.time(end)];
      elseif timwinstatflag==0
        cfg.latency=timwin;
        stattimwin=[stat_tacMnul_pb.time(1) stat_tacMnul_pb.time(end)];
      end
      tmpu1=ft_selectdata(cfg,grind_tacMnul_pb{1});
      tmpu1.dimord='chan_time';
      tmpu1.avg=squeeze(mean(tmpu1.individual,1));
      tmpu1=rmfield(tmpu1,'individual');
      
      tmpm10=ft_selectdata(cfg,grind_tacMnul_pb{3});
      tmpm10.dimord='chan_time';
      tmpm10.avg=squeeze(mean(tmpm10.individual,1));
      tmpm10=rmfield(tmpm10,'individual');
      
      tmpd4=ft_selectdata(cfg,grind_tacMnul_pb{2});
      tmpd4.dimord='chan_time';
      tmpd4.avg=squeeze(mean(tmpd4.individual,1));
      tmpd4=rmfield(tmpd4,'individual');
      
      tmpd7=ft_selectdata(cfg,grind_tacMnul_pb{4});
      tmpd7.dimord='chan_time';
      tmpd7.avg=squeeze(mean(tmpd7.individual,1));
      tmpd7=rmfield(tmpd7,'individual');
      
      tmpmask=stat_tacMnul_pb.mask;
      
  end
  
  
  
  
  if timwinstatflag==0
    tmpu1.mask=zeros(size(tmpu1.avg,1),length(tmpu1.time));
    tmpu1.mask(:,dsearchn(tmpu1.time',stattimwin(1)):dsearchn(tmpu1.time',stattimwin(end)))=tmpmask;
    tmpm10.mask=zeros(size(tmpm10.avg,1),length(tmpm10.time));
    tmpm10.mask(:,dsearchn(tmpm10.time',stattimwin(1)):dsearchn(tmpm10.time',stattimwin(end)))=tmpmask;
    tmpd4.mask=zeros(size(tmpd4.avg,1),length(tmpd4.time));
    tmpd4.mask(:,dsearchn(tmpd4.time',stattimwin(1)):dsearchn(tmpd4.time',stattimwin(end)))=tmpmask;
    tmpd7.mask=zeros(size(tmpd7.avg,1),length(tmpd7.time));
    tmpd7.mask(:,dsearchn(tmpd7.time',stattimwin(1)):dsearchn(tmpd7.time',stattimwin(end)))=tmpmask;
  else
    tmpu1.mask=tmpmask;
    tmpm10.mask=tmpmask;
    tmpd4.mask=tmpmask;
    tmpd7.mask=tmpmask;
  end
  
  
  for cg=1:length(chanplot)
    cfg=[];
    cfg.parameter='avg';
    cfg.layout='elec1010.lay';
    cfg.ylim=[-5 8];
    cfg.linewidth=3;
    cfg.xlim=timwin;
    if cg>length(chanplot)
      cfg.channel=tmpd4.label(any(tmpd4.mask,2));
    else
      cfg.channel=chanplot{cg};
    end
    cfg.graphcolor=coloruse([1 4 10 7],:);
    cfg.interactive='no';
    cfg.maskparameter='mask';
    cfg.maskstyle='box'; % default
    figure(ll*100+10*(cg+1))
    ft_singleplotER(cfg,tmpu1,tmpd7,tmpm10,tmpd4); % index 1, 4, 3, 2, corresponding to
    % see bin_my_angle for what each of 4 bins means
    % 1: peak, 2: peak to trough, 3: trough, 4: trough to peak
    hold on;plot(tmpu1.time,0,'k');
    set(gca,'XTick',[-.5:.1:1])
    set(gca,'XTickLabel',{'-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'})
    set(gca,'FontSize',30)
    title([])
    plot([stattimwin(1) stattimwin(1)],cfg.ylim,'k--','LineWidth',6)
    plot([stattimwin(2) stattimwin(2)],cfg.ylim,'k--','LineWidth',6)
    if ll<10 || ll==12
      plot([0 0],cfg.ylim,'Color',coloruse(4,:),'LineWidth',6)
    end
    if ll<10
      plot([soades(ll) soades(ll)],cfg.ylim,'Color',coloruse(9,:),'LineWidth',6)
    elseif ll==11
      plot([soades(5) soades(5)],cfg.ylim,'Color',coloruse(9,:),'LineWidth',6)
    end
    axis([-0.55 1 cfg.ylim(1) cfg.ylim(2)])
    if cg==3
      legend('Sum Unisensory','MultSens + Null','SumUnisens - MultsensNull')
    end
  end
  print(ll*100+20,[fdir 'erp_' condname{ll} 'vsNul_diff_PeakVsTrgh_FC_sleep' num2str(sleep) '_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  print(ll*100+30,[fdir 'erp_' condname{ll} 'vsNul_diff_PeakVsTrgh_OP_sleep' num2str(sleep) '_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  
  % topos of significant chunks only
  if ~isempty(find(mean(tmpmask,1))) %in other words, if there is some significant time point(s)
    masktime=find(any(tmpmask,1));
    cfg=[];
    cfg.parameter='avg';
    cfg.layout='elec1010.lay';
    cfg.maskalpha=0.5;
    cfg.zlim=[-5 5];
    cfg.highlight='on';
    cfg.highlightsize=12;
    cfg.xlim=[tmpu1.time(masktime(1)) tmpu1.time(masktime(end))];
    cfg.comment='no';
    sigchannels=tmpu1.label(find(ceil(mean(tmpu1.mask(:,dsearchn(tmpu1.time',cfg.xlim(1)):dsearchn(tmpu1.time',cfg.xlim(2))),2))));
    cfg.highlightchannel=sigchannels;
    figure(100+ll);
    ft_topoplotER(cfg,tmpu1);
    print(100+ll,[fdir 'erp_topo_Peak_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    figure(100+10+ll);
    ft_topoplotER(cfg,tmpd7);
    print(100+10+ll,[fdir 'erp_topo_Trgh2Peak_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    figure(100+20+ll);
    ft_topoplotER(cfg,tmpm10);
    print(100+20+ll,[fdir 'erp_topo_Trgh_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    figure(100+30+ll);
    ft_topoplotER(cfg,tmpd4);
    print(100+30+ll,[fdir 'erp_topo_Peak2Trgh_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  end
end


