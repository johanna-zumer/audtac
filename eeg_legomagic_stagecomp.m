eeg_legomagic_preamble

%% Compare Stages during 'sleep' data (W vs N1 vs N2)
% This had been written and then lost during a file transfer.
% Thus only bits of original code remain or are gussed at as to how it was set up.
% At least this was discovered only a few days after it was written, so still very fresh in my mind!

% To generate
% tlock_statmc_STAGECOMP_sleep1_ss12_iter11_trialkc-1.mat
% tlock_grind_STAGECOMP_sleep1_iter11.mat

plotflag=0;
printflag=0;
statsflag=1;
soalist=[1 3 4 5 6 7 9];
tophalfflag=0;


chanuse_sleep0={'all' '-F4'};
chanuse_sleep1={'all' '-AF7' '-AF3' '-Fp1'};
tt=3;

sleep=1;
chanuse=chanuse_sleep1;
iteruse=11;
iter=iteruse;
trialkc=-1;
usetr=1;

tacaud=1;
if sleep
  if tophalfflag
    load sortTacN2.mat
    subuseall=sort(sortTacN2(end-8:end)');
  else
    subuseall=setdiff(iiBuse,[3:7]);
  end
else
  subuseall=iiSuse;
end



for ll=soalist
% for ll=[7 9]
  clearvars -except ll tt sub edir ddir ii* sleep *flag soa* chanuse* iter* usetr trial* synch* tacaud  subuseall stat* grind*
  submin=subuseall(1)-1;
  subuseind=0;
  %         subuse=subuseall;
  subuse=nan(12,length(subuseall));
  subuse(10:12,:)=repmat(subuseall,[3 1]);
  for ii=subuseall
    subuseind=subuseind+1;
    cd([edir sub{ii} ])
    
    load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'])
    tk1=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter)  '_usetr' num2str(usetr) '_trialkc' num2str(trialkc) '.mat']);
    
    
    for ss=10:12
      
      numtrt(ll,tt,ss,subuseind)=numt_trials(ll,tt,ss); % does this need to be per 'sleep01' as well?
      if numtrt(ll,tt,ss,subuseind)<20
        subuse(ss,subuseind)=nan;
      end
      
      if ~isnan(subuse(ss,subuseind))
        
        tlock_tacPaud_each{ss}{subuseind}=tlock_tacPaud{ll,tt,ss};
        tlock_tacMSpN_each{ss}{subuseind}=tlock_tacMSpN{ll,tt,ss};
        tlock_MStlock_each{ss}{subuseind}=tlock_tMSAlone{ll,tt,ss}; % ll is tac plus auditory multisensory, tac-locked
        tlock_tactlock_each{ss}{subuseind}=tlock_tTacAlone{ll,tt,ss}; % ll+20 is tac alone, tac-locked
        tlock_audtlock_each{ss}{subuseind}=tlock_tAudAlone{ll,tt,ss}; % ll+20 is aud alone, aud-locked
        tlock_nulttlock_each{ss}{subuseind}=tlock_tNulAlone{ll,tt,ss}; % ll+50 is nul, tac-locked
        
        cfg=[];
        cfg.operation='subtract';
        cfg.parameter='avg';
        tlock_tacVSnul_each{ss}{subuseind}=ft_math(cfg,tlock_tTacAlone{ll,tt,ss},tlock_tNulAlone{ll,tt,ss});
        tlock_audVSnul_each{ss}{subuseind}=ft_math(cfg,tlock_tAudAlone{ll,tt,ss},tlock_tNulAlone{ll,tt,ss});
        tlock_msVSnul_each{ss}{subuseind}=ft_math(cfg,tlock_tMSAlone{ll,tt,ss},tlock_tNulAlone{ll,tt,ss});
      end
    end % ss
    
    clear tlock*N tlock*tac tlock*aud
  end % ii
  subuseindfinal=subuseind
  
  
  useAll=~isnan(subuse(10,:)) & ~isnan(subuse(11,:));
  useN1N2=~isnan(subuse(11,:));
  
  for ss=10:12
    
    cfg=[];
    cfg.keepindividual='yes';
    cfg.channel=chanuse;
    %     grind_tacPaud_stageAll{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_tacPaud_each{ss}{useAll});
    %     grind_tacMSpN_stageAll{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_tacMSpN_each{ss}{useAll});
    grind_tactlock_stageAll_save{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_tactlock_each{ss}{useAll});
    grind_audtlock_stageAll_save{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_audtlock_each{ss}{useAll});
    grind_nultlock_stageAll_save{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_nulttlock_each{ss}{useAll});
    grind_MStlock_stageAll_save{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_MStlock_each{ss}{useAll});
    
    if ss>10
      %       grind_tacPaud_stageN1N2{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_tacPaud_each{ss}{useN1N2});
      %       grind_tacMSpN_stageN1N2{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_tacMSpN_each{ss}{useN1N2});
      grind_tactlock_stageN1N2_save{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_tactlock_each{ss}{useN1N2});
      grind_audtlock_stageN1N2_save{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_audtlock_each{ss}{useN1N2});
      grind_nultlock_stageN1N2_save{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_nulttlock_each{ss}{useN1N2});
      grind_MStlock_stageN1N2_save{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_MStlock_each{ss}{useN1N2});
    end
    
    cfg=[];
    cfg.channel=chanuse;
    grave_tacPaud_stageAll{ss}=ft_timelockgrandaverage(cfg,tlock_tacPaud_each{ss}{useAll});
    grave_tacMSpN_stageAll{ss}=ft_timelockgrandaverage(cfg,tlock_tacMSpN_each{ss}{useAll});
    grave_tactlock_stageAll{ss}=ft_timelockgrandaverage(cfg,tlock_tactlock_each{ss}{useAll});
    grave_audtlock_stageAll{ss}=ft_timelockgrandaverage(cfg,tlock_audtlock_each{ss}{useAll});
    grave_nultlock_stageAll{ss}=ft_timelockgrandaverage(cfg,tlock_nulttlock_each{ss}{useAll});
    grave_MStlock_stageAll{ss}=ft_timelockgrandaverage(cfg,tlock_MStlock_each{ss}{useAll});
    
    grave_tacVSnul_stageAll{ss}=ft_timelockgrandaverage(cfg,tlock_tacVSnul_each{ss}{useAll});
    grave_audVSnul_stageAll{ss}=ft_timelockgrandaverage(cfg,tlock_audVSnul_each{ss}{useAll});
    grave_msVSnul_stageAll{ss}=ft_timelockgrandaverage(cfg,tlock_msVSnul_each{ss}{useAll});
    
    if ss>10
      grave_tacPaud_stageN1N2{ss}=ft_timelockgrandaverage(cfg,tlock_tacPaud_each{ss}{useN1N2});
      grave_tacMSpN_stageN1N2{ss}=ft_timelockgrandaverage(cfg,tlock_tacMSpN_each{ss}{useN1N2});
      grave_tactlock_stageN1N2{ss}=ft_timelockgrandaverage(cfg,tlock_tactlock_each{ss}{useN1N2});
      grave_audtlock_stageN1N2{ss}=ft_timelockgrandaverage(cfg,tlock_audtlock_each{ss}{useN1N2});
      grave_nultlock_stageN1N2{ss}=ft_timelockgrandaverage(cfg,tlock_nulttlock_each{ss}{useN1N2});
      grave_MStlock_stageN1N2{ss}=ft_timelockgrandaverage(cfg,tlock_MStlock_each{ss}{useN1N2});
      
      grave_tacVSnul_stageN1N2{ss}=ft_timelockgrandaverage(cfg,tlock_tacVSnul_each{ss}{useN1N2});
      grave_audVSnul_stageN1N2{ss}=ft_timelockgrandaverage(cfg,tlock_audVSnul_each{ss}{useN1N2});
      grave_msVSnul_stageN1N2{ss}=ft_timelockgrandaverage(cfg,tlock_msVSnul_each{ss}{useN1N2});
    end
    
    finduseAll=find(useAll);
    finduseN1N2=find(useN1N2);
    cfg=[];
    cfg.operation='subtract';
    cfg.parameter='avg';
    cfg.channel=chanuse;
    for ii=1:length(find(useAll))
      tlock_TPA_MSPN_stageAll{ss}{ii}=ft_math(cfg,tlock_tacPaud_each{ss}{finduseAll(ii)},tlock_tacMSpN_each{ss}{finduseAll(ii)});
    end
    if ss>10
      for ii=1:length(find(useN1N2))
        tlock_TPA_MSPN_stageN1N2{ss}{ii}=ft_math(cfg,tlock_tacPaud_each{ss}{finduseN1N2(ii)},tlock_tacMSpN_each{ss}{finduseN1N2(ii)});
      end
    end
    cfg=[];
    cfg.channel=chanuse;
    grave_TPA_MSPN_stageAll{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_TPA_MSPN_stageAll{ss}{:});
    if ss>10
      grave_TPA_MSPN_stageN1N2{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_TPA_MSPN_stageN1N2{ss}{:});
    end
    cfg.keepindividual='yes';
    grind_TPA_MSPN_stageAll{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_TPA_MSPN_stageAll{ss}{:});
    if ss>10
      grind_TPA_MSPN_stageN1N2{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_TPA_MSPN_stageN1N2{ss}{:});
    end
    
    %     cfg=[];
    %     cfg.latency=[-.5 1];
    %     grind_tacPaud_save{ll,tt,ss}=ft_selectdata(cfg,grind_tacPaud{ss});
    %     grind_tacMSpN_save{ll,tt,ss}=ft_selectdata(cfg,grind_tacMSpN{ss});
    %     grind_tacPaud_save{ll,tt,ss}=ft_selectdata(cfg,grind_tacPaud{ss});
    %     grind_tacMSpN_save{ll,tt,ss}=ft_selectdata(cfg,grind_tacMSpN{ss});
    
  end %ss
  
  
  if statsflag
    load eeg1010_neighb
    nsuball=length(find(sum(isnan(subuse(10:12,:)))==0));
    nsubN1N2=size(grind_TPA_MSPN_stageN1N2{ll,tt,11}.individual,1);
    
    
    cfg=[];
    if ll==1 || ll==3 || ll==4 || ll==5
      cfg.latency=[.1 .45];
    elseif ll==6
      cfg.latency=[.12 .47];
    elseif ll==7
      cfg.latency=[.17 .52];
    elseif ll==9
      cfg.latency=[.6 .95];
    end
    cfg.channel=chanuse;
    cfg.neighbours=neighbours;
    cfg.parameter='individual';
    cfg.method='montecarlo';
    cfg.numrandomization=2000;
    cfg.correctm='cluster';
    cfg.clusteralpha = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan = 2;
    cfg.statistic='depsamplesT';
    cfg.ivar=1;
    cfg.uvar=2;
    
    cfg.design=zeros(2,2*nsuball);
    cfg.design(1,:)=[ones(1,nsuball) 2*ones(1,nsuball)];
    cfg.design(2,:)=[1:nsuball 1:nsuball];
    
    statt_mc_WN1{ll}=ft_timelockstatistics(cfg,grind_TPA_MSPN_stageAll{ll,tt,11},grind_TPA_MSPN_stageAll{ll,tt,10});
    statt_ms_early_WN1{ll}=ft_timelockstatistics(cfg,grind_MStlock_stageAll_save{ll,tt,11},grind_MStlock_stageAll_save{ll,tt,10});
    
    cfg.design=zeros(2,2*nsubN1N2);
    cfg.design(1,:)=[ones(1,nsubN1N2) 2*ones(1,nsubN1N2)];
    cfg.design(2,:)=[1:nsubN1N2 1:nsubN1N2];
    
    statt_mc_N1N2{ll}=ft_timelockstatistics(cfg,grind_TPA_MSPN_stageN1N2{ll,tt,11},grind_TPA_MSPN_stageN1N2{ll,tt,12});
    statt_ms_early_N1N2{ll}=ft_timelockstatistics(cfg,grind_MStlock_stageN1N2_save{ll,tt,11},grind_MStlock_stageN1N2_save{ll,tt,12});
    
    if ll==1
      cfg.latency=[.1 .45]-.5;
    elseif ll==3
      cfg.latency=[.1 .45]-.07;
    elseif ll==4
      cfg.latency=[.1 .45]-.02;
    elseif ll==5
      cfg.latency=[.1 .45];
    elseif ll==6
      cfg.latency=[.12 .47];
    elseif ll==7
      cfg.latency=[.17 .52];
    elseif ll==9
      cfg.latency=[.6 .95];
    end
    
    statt_aud_early_N1N2{ll}=ft_timelockstatistics(cfg,grind_audtlock_stageN1N2_save{ll,tt,11},grind_audtlock_stageN1N2_save{ll,tt,12});
    
    cfg.design=zeros(2,2*nsuball);
    cfg.design(1,:)=[ones(1,nsuball) 2*ones(1,nsuball)];
    cfg.design(2,:)=[1:nsuball 1:nsuball];
    
    statt_aud_early_WN1{ll}=ft_timelockstatistics(cfg,grind_audtlock_stageAll_save{ll,tt,11},grind_audtlock_stageAll_save{ll,tt,10});
    
    cfg.latency=[.1 .45];
    
    statt_tac_early_WN1{ll}=ft_timelockstatistics(cfg,grind_tactlock_stageAll_save{ll,tt,11},grind_tactlock_stageAll_save{ll,tt,10});
    
    cfg.design=zeros(2,2*nsubN1N2);
    cfg.design(1,:)=[ones(1,nsubN1N2) 2*ones(1,nsubN1N2)];
    cfg.design(2,:)=[1:nsubN1N2 1:nsubN1N2];
    
    statt_tac_early_N1N2{ll}=ft_timelockstatistics(cfg,grind_tactlock_stageN1N2_save{ll,tt,11},grind_tactlock_stageN1N2_save{ll,tt,12});
    
    save([edir 'tlock_statmc_STAGECOMP_sleep' num2str(sleep) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'stat*');
    save([edir 'tlock_grind_STAGECOMP_sleep' num2str(sleep) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'grind*');
    
    cfg=[];
    cfg.neighbours=neighbours;
    cfg.parameter='individual';
    cfg.method='montecarlo';
    cfg.numrandomization=1000;
    cfg.correctm='cluster';
    cfg.clusteralpha = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan = 2;
    cfg.statistic='depsamplesFunivariate';
    cfg.ivar=1;
    cfg.uvar=2;
    cfg.tail=1;
    cfg.design=zeros(2,3*nsuball);
    cfg.design(1,:)=[ones(1,nsuball) 2*ones(1,nsuball) 3*ones(1,nsuball)];
    cfg.design(2,:)=[1:nsuball 1:nsuball 1:nsuball];
    
    if ll==1 || ll==3 || ll==4 || ll==5
      cfg.latency=[.1 .45];
    elseif ll==6
      cfg.latency=[.12 .47];
    elseif ll==7
      cfg.latency=[.17 .52];
    elseif ll==9
      cfg.latency=[.6 .95];
    end
    
    statt_mc_All{ll}=ft_timelockstatistics(cfg,grind_TPA_MSPN_stageAll{ll,tt,10},grind_TPA_MSPN_stageAll{ll,tt,11},grind_TPA_MSPN_stageAll{ll,tt,12});
    statt_ms_early_All{ll}=ft_timelockstatistics(cfg,grind_MStlock_stageAll_save{ll,tt,10},grind_MStlock_stageAll_save{ll,tt,11},grind_MStlock_stageAll_save{ll,tt,12});
    
    if ll==1
      cfg.latency=[.1 .45]-.5;
    elseif ll==3
      cfg.latency=[.1 .45]-.07;
    elseif ll==4
      cfg.latency=[.1 .45]-.02;
    elseif ll==5
      cfg.latency=[.1 .45];
    elseif ll==6
      cfg.latency=[.12 .47];
    elseif ll==7
      cfg.latency=[.17 .52];
    elseif ll==9
      cfg.latency=[.6 .95];
    end
    
    statt_aud_early_All{ll}=ft_timelockstatistics(cfg,grind_audtlock_stageAll_save{ll,tt,10},grind_audtlock_stageAll_save{ll,tt,11},grind_audtlock_stageAll_save{ll,tt,12});
    
    cfg.latency=[.1 .45];
    
    statt_tac_early_All{ll}=ft_timelockstatistics(cfg,grind_tactlock_stageAll_save{ll,tt,10},grind_tactlock_stageAll_save{ll,tt,11},grind_tactlock_stageAll_save{ll,tt,12});
    
    save([edir 'tlock_statmc_STAGECOMP_sleep' num2str(sleep) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'stat*');
  end
  
end  %ll

%%  Something else

%%  Comparing awake to asleep

printflag=0;
plotflag=0;
dostats=0;

for tt=2
  for ll=[1 3 4 5 6 7 9]
    clearvars -except ll tt sub edir ddir sdir iiuse plotflag printflag dostats fstats*
    
    for sleep=[0 1]
      
      subuse=iiuse;
      
      submin=subuse(1)-1;
      subuseind=0;
      for ii=subuse
        %       for ii=setdiff(subuse,[8 9 10 12 14 15 16 17 18])
        cd([edir sub{ii} ])
        load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '.mat'])
        
        tk0=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(0) '.mat']);
        tk1=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '.mat']);
        
        % THis is preliminary...how best to included all stages later on?
        if sleep==0
          ss=10; % awake
        elseif sleep==1
          ss=23; % this is concatenation of N2 and N3
        end
        
        %         for ss=ssuse
        if sleep==0 %load for both at sleep0 then will still be in memory for sleep1
          numtrt(ll,tt,10,ii-submin)=numt_trials(ll,tt,10); % does this need to be per 'sleep01' as well?
          if audtacflag
            numtra(ll,tt,10,ii-submin)=numa_trials(ll,tt,10);
          end
          load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(1) '.mat'],'num*trials')
          numtrt(ll,tt,23,ii-submin)=numt_trials(ll,tt,23); % does this need to be per 'sleep01' as well?
          if audtacflag
            numtra(ll,tt,23,ii-submin)=numa_trials(ll,tt,23);
          end
        end
        %         end
        
        % discard from both if either sleep/wake doesn't have 'enough' (what is enough? 20?)
        if numtrt(ll,tt,10,ii-submin)<20 || numtrt(ll,tt,23,ii-submin)<20
          subuse=setdiff(subuse,ii);
        else
          subuseind=subuseind+1;
          tlock_tacPaud_each{subuseind}=tlock_tacPaud{ll,tt,ss};
          tlock_audPtac_each{subuseind}=tlock_audPtac{ll,tt,ss};
          tlock_tacMSpN_each{subuseind}=tlock_tacMSpN{ll,tt,ss};
          tlock_audMSpN_each{subuseind}=tlock_audMSpN{ll,tt,ss};
          if sleep==1
            fstats_tpa(:,ll,tt,:,subuseind)=featurestats_tacPaud(:,ll,tt,[10 23],ii);
            fstats_apt(:,ll,tt,:,subuseind)=featurestats_audPtac(:,ll,tt,[10 23],ii);
            fstats_tmspn(:,ll,tt,:,subuseind)=featurestats_tacMSpN(:,ll,tt,[10 23],ii);
            fstats_amspn(:,ll,tt,:,subuseind)=featurestats_audMSpN(:,ll,tt,[10 23],ii);
          end
        end
        %       clear *_s0
        clear tlock*N tlock*tac tlock*aud
      end
      subuseindfinal=subuseind;
      
      for ii=1:subuseindfinal
        cfg=[];
        cfg.operation='subtract';
        cfg.parameter='avg';
        tlock_TPA_MSPN{ii}=ft_math(cfg,tlock_tacPaud_each{ii},tlock_tacMSpN_each{ii});
        tlock_APT_MSPN{ii}=ft_math(cfg,tlock_audPtac_each{ii},tlock_audMSpN_each{ii});
      end
      cfg=[];
      cfg.keepindividual='yes';
      grind_TPA_MSPN{ss}=ft_timelockgrandaverage(cfg,tlock_TPA_MSPN{:});
      grind_APT_MSPN{ss}=ft_timelockgrandaverage(cfg,tlock_APT_MSPN{:});
      cfg=[];
      grave_TPA_MSPN{ss}=ft_timelockgrandaverage(cfg,tlock_TPA_MSPN{:});
      grave_APT_MSPN{ss}=ft_timelockgrandaverage(cfg,tlock_APT_MSPN{:});
      
      if plotflag
        figure(20);
        for ii=1:subuseindfinal
          if ss==23
            subplot(3,subuseindfinal,ii);               imagesc(tlock_TPA_MSPN{ii}.time,1:63,tlock_TPA_MSPN{ii}.avg);caxis([-6 6])
          elseif ss==10
            subplot(3,subuseindfinal,subuseindfinal+ii);imagesc(tlock_TPA_MSPN{ii}.time,1:63,tlock_TPA_MSPN{ii}.avg);caxis([-6 6])
          end
        end
        figure(21);
        for ii=1:subuseindfinal
          if ss==23
            subplot(3,subuseindfinal,ii);               imagesc(tlock_APT_MSPN{ii}.time,1:63,tlock_APT_MSPN{ii}.avg);caxis([-6 6])
          elseif ss==10
            subplot(3,subuseindfinal,subuseindfinal+ii);imagesc(tlock_APT_MSPN{ii}.time,1:63,tlock_APT_MSPN{ii}.avg);caxis([-6 6])
          end
        end
      end
      
      
    end % end ss
    
    cfg=[];
    cfg.operation='subtract';
    cfg.parameter='individual';
    grind_TPA_MSPN_sleepdiff=ft_math(cfg,grind_TPA_MSPN{23},grind_TPA_MSPN{10})
    grind_APT_MSPN_sleepdiff=ft_math(cfg,grind_APT_MSPN{23},grind_APT_MSPN{10})
    
    if plotflag
      figure(20);
      for ii=1:subuseindfinal
        subplot(3,subuseindfinal,2*subuseindfinal+ii);imagesc(grind_TPA_MSPN_sleepdiff.time,1:63,squeeze(grind_TPA_MSPN_sleepdiff.individual(ii,:,:)));caxis([-6 6])
      end
      figure(21);
      for ii=1:subuseindfinal
        subplot(3,subuseindfinal,2*subuseindfinal+ii);imagesc(grind_APT_MSPN_sleepdiff.time,1:63,squeeze(grind_APT_MSPN_sleepdiff.individual(ii,:,:)));caxis([-6 6])
      end
      if printflag
        print(20,['D:\audtac\figs\indiv_tpa_mspn_N23_W_diff_cond' num2str(ll) num2str(tt) '.png'],'-dpng')
        print(21,['D:\audtac\figs\indiv_apt_mspn_N23_W_diff_cond' num2str(ll) num2str(tt) '.png'],'-dpng')
      end
    end
    
    if plotflag
      topoplot_highlight(11,grave_TPA_MSPN{23},[-0.5 0.6],[]);
      topoplot_highlight(12,grave_TPA_MSPN{10},[-0.5 0.6],[]);
      topoplot_highlight(13,grave_APT_MSPN{23},[-0.1 1.1],[]);
      topoplot_highlight(14,grave_APT_MSPN{10},[-0.1 1.1],[]);
      
      if printflag
        print(11,['D:\audtac\figs\grave_tpamspnN23_topoOverTime_cond' num2str(ll) num2str(tt) '.png'],'-dpng')
        print(12,['D:\audtac\figs\grave_tpamspnW_topoOverTime_cond' num2str(ll) num2str(tt) '.png'],'-dpng')
        print(13,['D:\audtac\figs\grave_aptmspnN23_topoOverTime_cond' num2str(ll) num2str(tt) '.png'],'-dpng')
        print(14,['D:\audtac\figs\grave_aptmspnW_topoOverTime_cond' num2str(ll) num2str(tt) '.png'],'-dpng')
      end
    end
    
    if dostats
      load eeg1010_neighb
      
      nsub=subuseindfinal;
      
      cfg=[];
      cfg.latency=[-.1 .5];
      cfg.neighbours=neighbours;
      % cfg.parameter='avg';
      cfg.parameter='individual';
      cfg.method='montecarlo';
      % cfg.method='analytic';
      cfg.numrandomization=2000;
      % cfg.correctm='holm';
      cfg.correctm='cluster';
      cfg.clusteralpha = 0.05;
      cfg.clusterstatistic = 'maxsum';
      cfg.minnbchan = 2;
      cfg.statistic='depsamplesT';
      % cfg.statistic='indepsamplesregrT';
      % cfg.statistic='indepsamplesT';
      cfg.design=zeros(2,2*nsub);
      cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
      cfg.design(2,:)=[1:nsub 1:nsub];
      cfg.ivar=1;
      cfg.uvar=2;
      %     statt_mc=ft_timelockstatistics(cfg, grind_tacPaud, grind_tacMSpN);
      %     stata_mc=ft_timelockstatistics(cfg, grind_audPtac, grind_audMSpN);
      statt_mc=ft_timelockstatistics(cfg, grind_TPA_MSPN{23}, grind_TPA_MSPN{10});
      stata_mc=ft_timelockstatistics(cfg, grind_APT_MSPN{23}, grind_APT_MSPN{10});
      
      grave_TPA_MSPN_sleepdiff.avg=squeeze(mean(grind_TPA_MSPN_sleepdiff.individual,1));
      grave_TPA_MSPN_sleepdiff.time=grind_TPA_MSPN_sleepdiff.time;
      grave_TPA_MSPN_sleepdiff.label=grind_TPA_MSPN_sleepdiff.label;
      grave_TPA_MSPN_sleepdiff.dimord='chan_time';
      grave_APT_MSPN_sleepdiff.avg=squeeze(mean(grind_APT_MSPN_sleepdiff.individual,1));
      grave_APT_MSPN_sleepdiff.time=grind_APT_MSPN_sleepdiff.time;
      grave_APT_MSPN_sleepdiff.label=grind_APT_MSPN_sleepdiff.label;
      grave_APT_MSPN_sleepdiff.dimord='chan_time';
      
      if plotflag
        topoplot_highlight(22,grave_TPA_MSPN_sleepdiff,[statt_mc.time(1) statt_mc.time(end)],statt_mc);
        topoplot_highlight(23,grave_APT_MSPN_sleepdiff,[stata_mc.time(1) stata_mc.time(end)],stata_mc);
        if printflag
          print(22,['D:\audtac\figs\grSleepdiff_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) '.png'],'-dpng')
          print(23,['D:\audtac\figs\grSleepdiff_topoOverTime_mc_at_cond' num2str(ll) num2str(tt) '.png'],'-dpng')
        end
      end
    end
    
  end
  % do some processing on fstats* here, over all 'll' but within a 'tt'
  
  fstats_tpa(3,:,tt,:,:)
  
end
save([edir 'tlockSLEEP01_numtrlltt.mat'],'numtr*','grind*','fstats*');

cd(sdir)
spss_exporter(squeeze([fstats_tmspn(3,[1 3 4 5 6 7 9],2,2,:); fstats_tpa(3,[1 3 4 5 6 7 9],2,2,:)]),{'tKc','MultSens','SOA'},0);  % no
spss_exporter(squeeze([fstats_tmspn(4,[1 3 4 5 6 7 9],2,2,:); fstats_tpa(4,[1 3 4 5 6 7 9],2,2,:)]),{'tDelta','MultSens','SOA'},0) % no
spss_exporter(squeeze([fstats_tmspn(5,[1 3 4 5 6 7 9],2,2,:); fstats_tpa(5,[1 3 4 5 6 7 9],2,2,:)]),{'tSW','MultSens','SOA'},0) % no
spss_exporter(squeeze([fstats_tmspn(6,[1 3 4 5 6 7 9],2,2,:); fstats_tpa(6,[1 3 4 5 6 7 9],2,2,:)]),{'tSpindleF','MultSens','SOA'},0) % possible soa and/or interaction
spss_exporter(squeeze([fstats_tmspn(7,[1 3 4 5 6 7 9],2,2,:); fstats_tpa(7,[1 3 4 5 6 7 9],2,2,:)]),{'tSpindleS','MultSens','SOA'},0) % error

spss1factor_exporter(squeeze(fstats_tmspn(3,[1 3 4 5 6 7 9],2,2,:)-fstats_tpa(3,[1 3 4 5 6 7 9],2,2,:)),{'tKcDiff','SOA','MultSens'},0) % no
spss1factor_exporter(squeeze(fstats_tmspn(4,[1 3 4 5 6 7 9],2,2,:)-fstats_tpa(4,[1 3 4 5 6 7 9],2,2,:)),{'tDeltaDiff','SOA','MultSens'},0) % no
spss1factor_exporter(squeeze(fstats_tmspn(5,[1 3 4 5 6 7 9],2,2,:)-fstats_tpa(5,[1 3 4 5 6 7 9],2,2,:)),{'tSWDiff','SOA','MultSens'},0) % no
spss1factor_exporter(squeeze(fstats_tmspn(6,[1 3 4 5 6 7 9],2,2,:)-fstats_tpa(6,[1 3 4 5 6 7 9],2,2,:)),{'tSpindleFDiff','SOA','MultSens'},0) % yes
spss1factor_exporter(squeeze(fstats_tmspn(7,[1 3 4 5 6 7 9],2,2,:)-fstats_tpa(7,[1 3 4 5 6 7 9],2,2,:)),{'tSpindleSDiff','SOA','MultSens'},0) % no

spss_exporter(squeeze([fstats_amspn(3,[1 3 4 5 6 7 9],2,2,:); fstats_apt(3,[1 3 4 5 6 7 9],2,2,:)]),{'aKc','MultSens','SOA'},0)
spss_exporter(squeeze([fstats_amspn(4,[1 3 4 5 6 7 9],2,2,:); fstats_apt(4,[1 3 4 5 6 7 9],2,2,:)]),{'aDelta','MultSens','SOA'},0)
spss_exporter(squeeze([fstats_amspn(5,[1 3 4 5 6 7 9],2,2,:); fstats_apt(5,[1 3 4 5 6 7 9],2,2,:)]),{'aSW','MultSens','SOA'},0)
spss_exporter(squeeze([fstats_amspn(6,[1 3 4 5 6 7 9],2,2,:); fstats_apt(6,[1 3 4 5 6 7 9],2,2,:)]),{'aSpindleF','MultSens','SOA'},0)
spss_exporter(squeeze([fstats_amspn(7,[1 3 4 5 6 7 9],2,2,:); fstats_apt(7,[1 3 4 5 6 7 9],2,2,:)]),{'aSpindleS','MultSens','SOA'},0)

spss1factor_exporter(squeeze(fstats_amspn(3,[1 3 4 5 6 7 9],2,2,:)-fstats_apt(3,[1 3 4 5 6 7 9],2,2,:)),{'aKcDiff','SOA','MultSens'},0)
spss1factor_exporter(squeeze(fstats_amspn(4,[1 3 4 5 6 7 9],2,2,:)-fstats_apt(4,[1 3 4 5 6 7 9],2,2,:)),{'aDeltaDiff','SOA','MultSens'},0)
spss1factor_exporter(squeeze(fstats_amspn(5,[1 3 4 5 6 7 9],2,2,:)-fstats_apt(5,[1 3 4 5 6 7 9],2,2,:)),{'aSWDiff','SOA','MultSens'},0)
spss1factor_exporter(squeeze(fstats_amspn(6,[1 3 4 5 6 7 9],2,2,:)-fstats_apt(6,[1 3 4 5 6 7 9],2,2,:)),{'aSpindleFDiff','SOA','MultSens'},0) % no
spss1factor_exporter(squeeze(fstats_amspn(7,[1 3 4 5 6 7 9],2,2,:)-fstats_apt(7,[1 3 4 5 6 7 9],2,2,:)),{'aSpindleSDiff','SOA','MultSens'},0)

spss_exporter(squeeze([mean(fstats_tmspn(3:5,[1 3 4 5 6 7 9],2,2,:),1); mean(fstats_tpa(3:5,[1 3 4 5 6 7 9],2,2,:),1)]),{'tBigWaves','MultSens','SOA'},0) % no
spss_exporter(squeeze([mean(fstats_tmspn(6:7,[1 3 4 5 6 7 9],2,2,:),1); mean(fstats_tpa(6:7,[1 3 4 5 6 7 9],2,2,:),1)]),{'tSpindleAll','MultSens','SOA'},0) % possible soa

spss1factor_exporter(squeeze(mean(fstats_tmspn(3:5,[1 3 4 5 6 7 9],2,2,:)-fstats_tpa(3:5,[1 3 4 5 6 7 9],2,2,:),1)),{'tBigWavesDiff','SOA','MultSens'},0) % no
spss1factor_exporter(squeeze(mean(fstats_tmspn(6:7,[1 3 4 5 6 7 9],2,2,:)-fstats_tpa(6:7,[1 3 4 5 6 7 9],2,2,:),1)),{'tSpindleAllDiff','SOA','MultSens'},0) % no

spss_exporter(squeeze([mean(fstats_amspn(3:5,[1 3 4 5 6 7 9],2,2,:),1); mean(fstats_apt(3:5,[1 3 4 5 6 7 9],2,2,:),1)]),{'aBigWaves','MultSens','SOA'},0)
spss_exporter(squeeze([mean(fstats_amspn(6:7,[1 3 4 5 6 7 9],2,2,:),1); mean(fstats_apt(6:7,[1 3 4 5 6 7 9],2,2,:),1)]),{'aSpindleAll','MultSens','SOA'},0)

spss1factor_exporter(squeeze(mean(fstats_amspn(3:5,[1 3 4 5 6 7 9],2,2,:)-fstats_apt(3:5,[1 3 4 5 6 7 9],2,2,:),1)),{'aBigWavesDiff','SOA','MultSens'},0)
spss1factor_exporter(squeeze(mean(fstats_amspn(6:7,[1 3 4 5 6 7 9],2,2,:)-fstats_apt(6:7,[1 3 4 5 6 7 9],2,2,:),1)),{'aSpindleAllDiff','SOA','MultSens'},0)

% test correlations of BigWaves and Spindles to behavioural outcomes
load([bdir 'rtgroup_pcb_diffms.mat'])
pcb=pcb(~isnan(pcb));
diffms=diffms(~isnan(diffms));

%collapsing over Bigwaves, diff of SOA456 to SOA19, split by medican of PCB
[hh,pp]=ttest2(squeeze([mean(mean(fstats_tmspn(3:5,[4 5 6],2,2,pcb>median(pcb)),2),1) - mean(mean(fstats_tmspn(3:5,[1 9],2,2,pcb>median(pcb)),2),1)]),squeeze([mean(mean(fstats_tmspn(3:5,[4 5 6],2,2,pcb<median(pcb)),2),1) - mean(mean(fstats_tmspn(3:5,[1 9],2,2,pcb<median(pcb)),2),1)]))
[hh,pp]=ttest2(squeeze([mean(mean(fstats_tmspn(6:7,[4 5 6],2,2,pcb>median(pcb)),2),1) - mean(mean(fstats_tmspn(6:7,[1 9],2,2,pcb>median(pcb)),2),1)]),squeeze([mean(mean(fstats_tmspn(6:7,[4 5 6],2,2,pcb<median(pcb)),2),1) - mean(mean(fstats_tmspn(6:7,[1 9],2,2,pcb<median(pcb)),2),1)]))

%collapsing over Bigwaves, diff of SOA456 to SOA19, split by medican of DiffMs
[hh,pp]=ttest2(squeeze([mean(mean(fstats_tmspn(3:5,[4 5 6],2,2,diffms>median(diffms)),2),1) - mean(mean(fstats_tmspn(3:5,[1 9],2,2,diffms>median(diffms)),2),1)]),squeeze([mean(mean(fstats_tmspn(3:5,[4 5 6],2,2,diffms<median(diffms)),2),1) - mean(mean(fstats_tmspn(3:5,[1 9],2,2,diffms<median(diffms)),2),1)]))
[hh,pp]=ttest2(squeeze([mean(mean(fstats_tmspn(6:7,[4 5 6],2,2,diffms>median(diffms)),2),1) - mean(mean(fstats_tmspn(6:7,[1 9],2,2,diffms>median(diffms)),2),1)]),squeeze([mean(mean(fstats_tmspn(6:7,[4 5 6],2,2,diffms<median(diffms)),2),1) - mean(mean(fstats_tmspn(6:7,[1 9],2,2,diffms<median(diffms)),2),1)]))

