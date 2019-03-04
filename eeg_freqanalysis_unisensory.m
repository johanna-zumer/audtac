eeg_legomagic_preamble

%% Checking if Null, Tac-alone, and Aud-alone from each of 7 SOA conditions not signficantly diff from each other.
load elec1010_neighb.mat

% allcond_sameN=1;
soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];
tacbasemax=min(-.15,soades-.15);
audbasemax=min(-.15,fliplr(soades)-.15);
msstart=max(0,soades);

ylimlo=[4 7; 8 12; 14 30];
ylimhi=[30 45; 55 80];

soadesc={'Aud first by 500ms' '' 'Aud first by 70ms' 'Aud first by 20ms' 'Simultaneous' 'Tac first by 20ms' 'Tac first by 70ms' '' 'Tac first by 500ms'};
plotflag=0;
printflag=0;
audtacflag=0;
fftaddflag=0;
synchasynch=0;
statspowflag=1;
statsplvflag=1;
itcdepsampflag=1;
savegrindflag=0;

soalist=[1 3 4 5 6 7 9];

chanuse_sleep0={'all' '-F4'};
chanuse_sleep1={'all' '-AF7' '-AF3' '-Fp1'};
tt=3;

sleep=0;

clearvars -except ll tt sub *dir ii*use sleep *flag figind soadesc soalist chanuse* ylim* neigh* *basemax
if sleep
  chanuse=chanuse_sleep1;
  subuseall=iiBuse;
  iteruse=11;
  trialkc=-1;   % CHANGE ME
  ssuse=10:12; % all stages
  sleepcond='Sleep N2';
  usetr=1;
  medsplitflag=1;
else
  chanuse=chanuse_sleep0;
  subuseall=setdiff(iiSuse,[]);
  iteruse=31;
  trialkc=-1;
  ssuse=10; % awake
  sleepcond='Awake W';
  usetr=3;
  medsplitflag=0;
end
%     ss=ssuse;

submin=subuseall(1)-1;
subuseind=0;
iteruse

for ii=subuseall
  cd([edir sub{ii} ])
  %       load(['freq_diffs_averef_' sub{ii} '.mat']);
  %       try
  %         load(['freq_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_tt' num2str(tt) '.mat'])
  %         load(['freq_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(0) '_tt' num2str(tt) '.mat'])
  %       catch
  %         if tt==2
  %           load(['freq_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '.mat'])
  %           load(['freq_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(0) '.mat'])
  %         end
  %       end
  
  % NOTE: this below is with trialkc=0 implicit!!
  %         load(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_tt' num2str(tt) '.mat'])
  %       try
  try
    load(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iteruse) '_trialkc' num2str(trialkc) '.mat'])
  catch
    load(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_tt' num2str(tt) '_trialkc' num2str(trialkc) '.mat'])
  end
  %       catch
  %         if trialkc==0
  %           load(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_tt' num2str(tt) '.mat'])
  %         end
  %       end
  
  
  %       try
  %         tkt=load(['trialkeptTFR_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iteruse) '_usetr' num2str(usetr) '_trialkc' num2str(trialkc) '.mat']);
  %       catch
  %         tkt=load(['trialkeptTFR_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '.mat']);
  %       end
  %         for ss=ssuse
  %       subuse=subuseall; % reset to all for each sleep stage
  subuseind=subuseind+1;
  
  for ss=ssuse
    for ll=soalist
      
      numtrt(ll,tt,ss,subuseind)=numt_trials(ll,tt,ss);
      if numtrt(ll,tt,ss,subuseind)<20
        subuse(ll,ss,subuseind)=nan;
      else
        subuse(ll,ss,subuseind)=1;
      end
      
      if ~isnan(subuse(ll,ss,subuseind))
        freqhiall_tNulAlone_comb{ll,ss}{subuseind}=freqhi_tNulAlone_comb{ll,tt,ss};
        freqloall_tNulAlone_comb{ll,ss}{subuseind}=freqlo_tNulAlone_comb{ll,tt,ss};
        freqhiall_tTacAlone_comb{ll,ss}{subuseind}=freqhi_tTacAlone_comb{ll,tt,ss};
        freqloall_tTacAlone_comb{ll,ss}{subuseind}=freqlo_tTacAlone_comb{ll,tt,ss};
        freqhiall_tAudAlone_comb{ll,ss}{subuseind}=freqhi_tAudAlone_comb{ll,tt,ss};
        freqloall_tAudAlone_comb{ll,ss}{subuseind}=freqlo_tAudAlone_comb{ll,tt,ss};
        freqhiall_tMSAlone_comb{ll,ss}{subuseind} =freqhi_tMSAlone_comb{ll,tt,ss};
        freqloall_tMSAlone_comb{ll,ss}{subuseind} =freqlo_tMSAlone_comb{ll,tt,ss};
        freqloall_tNulAlone_comb{ll,ss}{subuseind}.dimord='chan_freq_time';
        freqloall_tTacAlone_comb{ll,ss}{subuseind}.dimord='chan_freq_time';
        freqloall_tAudAlone_comb{ll,ss}{subuseind}.dimord='chan_freq_time';
        freqloall_tMSAlone_comb{ll,ss}{subuseind}.dimord='chan_freq_time';
        freqhiall_tNulAlone_comb{ll,ss}{subuseind}.dimord='chan_freq_time';
        freqhiall_tTacAlone_comb{ll,ss}{subuseind}.dimord='chan_freq_time';
        freqhiall_tAudAlone_comb{ll,ss}{subuseind}.dimord='chan_freq_time';
        freqhiall_tMSAlone_comb{ll,ss}{subuseind}.dimord='chan_freq_time';
        freqhiall_tNulAlone_comb{ll,ss}{subuseind}.plvabs=abs(freqhi_tNulAlone_comb{ll,tt,ss}.plvspctrm);
        freqloall_tNulAlone_comb{ll,ss}{subuseind}.plvabs=abs(freqlo_tNulAlone_comb{ll,tt,ss}.plvspctrm);
        freqhiall_tTacAlone_comb{ll,ss}{subuseind}.plvabs=abs(freqhi_tTacAlone_comb{ll,tt,ss}.plvspctrm);
        freqloall_tTacAlone_comb{ll,ss}{subuseind}.plvabs=abs(freqlo_tTacAlone_comb{ll,tt,ss}.plvspctrm);
        freqhiall_tAudAlone_comb{ll,ss}{subuseind}.plvabs=abs(freqhi_tAudAlone_comb{ll,tt,ss}.plvspctrm);
        freqloall_tAudAlone_comb{ll,ss}{subuseind}.plvabs=abs(freqlo_tAudAlone_comb{ll,tt,ss}.plvspctrm);
        freqhiall_tMSAlone_comb{ll,ss}{subuseind}.plvabs=abs(freqhi_tMSAlone_comb{ll,tt,ss}.plvspctrm);
        freqloall_tMSAlone_comb{ll,ss}{subuseind}.plvabs=abs(freqlo_tMSAlone_comb{ll,tt,ss}.plvspctrm);
      end
    end %ll
  end % ss
end %ii
subuseindfinal=subuseind;


for ss=ssuse
  
  cfg=[];
  cfg.keepindividual='yes';
  cfg.parameter={'powspctrm' 'plvspctrm' 'plvabs'};
  
  for ll=soalist
    usesub=~isnan(subuse(ll,ss,:));
    grindlo_tNulAlone{ll,ss}=ft_freqgrandaverage(cfg,freqloall_tNulAlone_comb{ll,ss}{usesub});
    grindhi_tNulAlone{ll,ss}=ft_freqgrandaverage(cfg,freqhiall_tNulAlone_comb{ll,ss}{usesub});
    grindlo_tTacAlone{ll,ss}=ft_freqgrandaverage(cfg,freqloall_tTacAlone_comb{ll,ss}{usesub});
    grindhi_tTacAlone{ll,ss}=ft_freqgrandaverage(cfg,freqhiall_tTacAlone_comb{ll,ss}{usesub});
    grindlo_tAudAlone{ll,ss}=ft_freqgrandaverage(cfg,freqloall_tAudAlone_comb{ll,ss}{usesub});
    grindhi_tAudAlone{ll,ss}=ft_freqgrandaverage(cfg,freqhiall_tAudAlone_comb{ll,ss}{usesub});
    grindlo_tMSAlone{ll,ss} =ft_freqgrandaverage(cfg,freqloall_tMSAlone_comb{ll,ss}{usesub});
    grindhi_tMSAlone{ll,ss} =ft_freqgrandaverage(cfg,freqhiall_tMSAlone_comb{ll,ss}{usesub});
  end % ll
  
  if savegrindflag
    save(['grindTFR_UniMsNul_sleep' num2str(sleep) '_iter' num2str(iteruse) '_trialkc' num2str(trialkc) '.mat'],'grind*');
  end
  
  cfg=[];
  cfg.neighbours=neighbours;
  cfg.method='montecarlo';
  cfg.numrandomization=2000;
  cfg.correctm='cluster';
  cfg.clusteralpha = 0.05;
  cfg.clusterstatistic = 'maxsum';
  cfg.minnbchan = 2;
  cfg.ivar=1;
  
  if statspowflag
    cfg.statistic='depsamplesT';
    cfg.uvar=2;
    cfg.parameter='powspctrm';
    for ll=soalist
      nsub=size(grindlo_tNulAlone{ll,ss}.powspctrm,1);
      cfg.design=zeros(2,2*nsub);
      cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
      cfg.design(2,:)=[1:nsub 1:nsub];
      cfg.latency=[-.15-audbasemax(5); .35-audbasemax(5)];
      stattl_mc_TacVsNul_short{ll,ss}=ft_freqstatistics(cfg, grindlo_tTacAlone{ll,ss}, grindlo_tNulAlone{ll,ss});
      cfg.latency=[-.15-audbasemax(5); 1.05-audbasemax(5)];
      stattl_mc_TacVsNul_long{ll,ss}=ft_freqstatistics(cfg, grindlo_tTacAlone{ll,ss}, grindlo_tNulAlone{ll,ss});
      if ll>=5
        cfg.latency=[-.15-audbasemax(ll); .35-audbasemax(ll)];
      else
        cfg.latency=[.15+tacbasemax(ll); .65+tacbasemax(ll)];
      end
      stattl_mc_AudVsNul_short{ll,ss}=ft_freqstatistics(cfg, grindlo_tAudAlone{ll,ss}, grindlo_tNulAlone{ll,ss});
      if ll>=5
        cfg.latency=[-.15-audbasemax(ll); 1.05-audbasemax(ll)];
      else
        cfg.latency=[.15+tacbasemax(ll); 1.35+tacbasemax(ll)];
      end
      stattl_mc_AudVsNul_long{ll,ss}=ft_freqstatistics(cfg, grindlo_tAudAlone{ll,ss}, grindlo_tNulAlone{ll,ss});
      cfg.latency=[-.15-audbasemax(ll); .35-audbasemax(ll)];
      stattl_mc_MSVsNul_short{ll,ss} =ft_freqstatistics(cfg, grindlo_tMSAlone{ll,ss}, grindlo_tNulAlone{ll,ss});
      cfg.latency=[-.15-audbasemax(ll); 1.05-audbasemax(ll)];
      stattl_mc_MSVsNul_long{ll,ss} =ft_freqstatistics(cfg, grindlo_tMSAlone{ll,ss}, grindlo_tNulAlone{ll,ss});
    end  % ll
    
  end
  
  if statsplvflag
    if itcdepsampflag
      cfg.statistic='depsamplesT';
      cfg.uvar=2;
      cfg.parameter='plvabs';
    else
      cfg.statistic='diff_itc';
      cfg.uvar=[];
      cfg.parameter='plvspctrm';
      cfg.complex='diffabs'; % default; not sensitive to phase differences between conditions
      cfg.clusterthreshold = 'nonparametric_individual'; % or 'nonparametric_individual' see clusterstat.m
    end
    for ll=soalist
      nsub=size(grindlo_tNulAlone{ll,ss}.powspctrm,1);
      cfg.design=zeros(2,2*nsub);
      cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
      cfg.design(2,:)=[1:nsub 1:nsub];
      if itcdepsampflag
        cfg.latency=[-.15-audbasemax(5); .35-audbasemax(5)];
        stattl_mcplv_TacVsNul_short_itcdepT{ll,ss}=ft_freqstatistics(cfg, grindlo_tTacAlone{ll,ss}, grindlo_tNulAlone{ll,ss});
        cfg.latency=[-.15-audbasemax(5); 1.05-audbasemax(5)];
        stattl_mcplv_TacVsNul_long_itcdepT{ll,ss}=ft_freqstatistics(cfg, grindlo_tTacAlone{ll,ss}, grindlo_tNulAlone{ll,ss});
        if ll>=5
          cfg.latency=[-.15-audbasemax(ll); .35-audbasemax(ll)];
        else
          cfg.latency=[.15+tacbasemax(ll); .65+tacbasemax(ll)];
        end
        stattl_mcplv_AudVsNul_short_itcdepT{ll,ss}=ft_freqstatistics(cfg, grindlo_tAudAlone{ll,ss}, grindlo_tNulAlone{ll,ss});
        if ll>=5
          cfg.latency=[-.15-audbasemax(ll); 1.05-audbasemax(ll)];
        else
          cfg.latency=[.15+tacbasemax(ll); 1.35+tacbasemax(ll)];
        end
        stattl_mcplv_AudVsNul_long_itcdepT{ll,ss}=ft_freqstatistics(cfg, grindlo_tAudAlone{ll,ss}, grindlo_tNulAlone{ll,ss});
        cfg.latency=[-.15-audbasemax(ll); .35-audbasemax(ll)];
        stattl_mcplv_MSVsNul_short_itcdepT{ll,ss} =ft_freqstatistics(cfg, grindlo_tMSAlone{ll,ss}, grindlo_tNulAlone{ll,ss});
        cfg.latency=[-.15-audbasemax(ll); 1.05-audbasemax(ll)];
        stattl_mcplv_MSVsNul_long_itcdepT{ll,ss} =ft_freqstatistics(cfg, grindlo_tMSAlone{ll,ss}, grindlo_tNulAlone{ll,ss});
      else
        error('helpitc')
        %         stattl_mcplv_TacVsNul{ll,ss}=ft_freqstatistics(cfg, grindlo_tTacAlone{ll,ss}, grindlo_tNulAlone{ll,ss});
        %         if ll>=5
        %           cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
        %         else
        %           cfg.latency=[.15+tacbasemax(ll); .85+tacbasemax(ll)];
        %         end
        %         stattl_mcplv_AudVsNul{ll,ss}=ft_freqstatistics(cfg, grindlo_tAudAlone{ll,ss}, grindlo_tNulAlone{ll,ss});
        %         cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
        %         stattl_mcplv_MSVsNul{ll,ss} =ft_freqstatistics(cfg, grindlo_tMSAlone{ll,ss}, grindlo_tNulAlone{ll,ss});
      end
      
    end %ll
  end
  
  
  % NOTE: without trialkc label, it refers to trialkc=0;
  %     if ~exist([edir 'stats_TFR_UniMSNul_' num2str(tt) num2str(ss) num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'file')
  %       %     save([edir 'statsgrave_TFR_UniMSNul_' num2str(tt) num2str(ss) num2str(sleep) '.mat'],'stat*','grave*');
  %       save([edir 'stats_TFR_UniMSNul_' num2str(tt) num2str(ss) num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'stat*');
  %     else
  %       save([edir 'stats_TFR_UniMSNul_' num2str(tt) num2str(ss) num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'stat*','-append');
  %     end
  if statspowflag || statsplvflag
    if itcdepsampflag
      save([edir 'stats_TFR_UniMSNul_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '_itc.mat'],'stat*','subuse');
    else
      if ~exist([edir 'stats_TFR_UniMSNul_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'file')
        %     save([edir 'statsgrave_TFR_UniMSNul_' num2str(tt) num2str(ss) num2str(sleep) '.mat'],'stat*','grave*');
        save([edir 'stats_TFR_UniMSNul_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'stat*','subuse');
      else
        save([edir 'stats_TFR_UniMSNul_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'stat*','subuse','-append');
      end
    end
  end
  
  if ~medsplitflag
    continue
  else
    
    load([edir 'sortTacWP200.mat'],'iiuse_*halfWwideP2N1');
    
    for ll=soalist
      subind_top=dsearchn(iiBuse',iiuse_tophalfWwideP2N1');
      subind_bot=dsearchn(iiBuse',iiuse_bothalfWwideP2N1');
      usesub=find(squeeze(~isnan(subuse(ll,ss,:))));
      
      grindlo_tTacAlone_top{ll,ss}=grindlo_tTacAlone{ll,ss};
      grindlo_tTacAlone_top{ll,ss}.powspctrm=grindlo_tTacAlone{ll,ss}.powspctrm(dsearchn(usesub,intersect(usesub,subind_top)),:,:,:);
      grindlo_tTacAlone_top{ll,ss}.plvspctrm=grindlo_tTacAlone{ll,ss}.plvspctrm(dsearchn(usesub,intersect(usesub,subind_top)),:,:,:);
      grindlo_tTacAlone_top{ll,ss}.plvabs=grindlo_tTacAlone{ll,ss}.plvabs(dsearchn(usesub,intersect(usesub,subind_top)),:,:,:);
      grindlo_tTacAlone_bot{ll,ss}=grindlo_tTacAlone{ll,ss};
      grindlo_tTacAlone_bot{ll,ss}.powspctrm=grindlo_tTacAlone{ll,ss}.powspctrm(dsearchn(usesub,intersect(usesub,subind_bot)),:,:,:);
      grindlo_tTacAlone_bot{ll,ss}.plvspctrm=grindlo_tTacAlone{ll,ss}.plvspctrm(dsearchn(usesub,intersect(usesub,subind_bot)),:,:,:);
      grindlo_tTacAlone_bot{ll,ss}.plvabs=grindlo_tTacAlone{ll,ss}.plvabs(dsearchn(usesub,intersect(usesub,subind_bot)),:,:,:);
      grindlo_tNulAlone_top{ll,ss}=grindlo_tNulAlone{ll,ss};
      grindlo_tNulAlone_top{ll,ss}.powspctrm=grindlo_tNulAlone{ll,ss}.powspctrm(dsearchn(usesub,intersect(usesub,subind_top)),:,:,:);
      grindlo_tNulAlone_top{ll,ss}.plvspctrm=grindlo_tNulAlone{ll,ss}.plvspctrm(dsearchn(usesub,intersect(usesub,subind_top)),:,:,:);
      grindlo_tNulAlone_top{ll,ss}.plvabs=grindlo_tNulAlone{ll,ss}.plvabs(dsearchn(usesub,intersect(usesub,subind_top)),:,:,:);
      grindlo_tNulAlone_bot{ll,ss}=grindlo_tNulAlone{ll,ss};
      grindlo_tNulAlone_bot{ll,ss}.powspctrm=grindlo_tNulAlone{ll,ss}.powspctrm(dsearchn(usesub,intersect(usesub,subind_bot)),:,:,:);
      grindlo_tNulAlone_bot{ll,ss}.plvspctrm=grindlo_tNulAlone{ll,ss}.plvspctrm(dsearchn(usesub,intersect(usesub,subind_bot)),:,:,:);
      grindlo_tNulAlone_bot{ll,ss}.plvabs=grindlo_tNulAlone{ll,ss}.plvabs(dsearchn(usesub,intersect(usesub,subind_bot)),:,:,:);
      grindlo_tAudAlone_top{ll,ss}=grindlo_tAudAlone{ll,ss};
      grindlo_tAudAlone_top{ll,ss}.powspctrm=grindlo_tAudAlone{ll,ss}.powspctrm(dsearchn(usesub,intersect(usesub,subind_top)),:,:,:);
      grindlo_tAudAlone_top{ll,ss}.plvspctrm=grindlo_tAudAlone{ll,ss}.plvspctrm(dsearchn(usesub,intersect(usesub,subind_top)),:,:,:);
      grindlo_tAudAlone_top{ll,ss}.plvabs=grindlo_tAudAlone{ll,ss}.plvabs(dsearchn(usesub,intersect(usesub,subind_top)),:,:,:);
      grindlo_tAudAlone_bot{ll,ss}=grindlo_tAudAlone{ll,ss};
      grindlo_tAudAlone_bot{ll,ss}.powspctrm=grindlo_tAudAlone{ll,ss}.powspctrm(dsearchn(usesub,intersect(usesub,subind_bot)),:,:,:);
      grindlo_tAudAlone_bot{ll,ss}.plvspctrm=grindlo_tAudAlone{ll,ss}.plvspctrm(dsearchn(usesub,intersect(usesub,subind_bot)),:,:,:);
      grindlo_tAudAlone_bot{ll,ss}.plvabs=grindlo_tAudAlone{ll,ss}.plvabs(dsearchn(usesub,intersect(usesub,subind_bot)),:,:,:);
      grindlo_tMSAlone_top{ll,ss}=grindlo_tMSAlone{ll,ss};
      grindlo_tMSAlone_top{ll,ss}.powspctrm=grindlo_tMSAlone{ll,ss}.powspctrm(dsearchn(usesub,intersect(usesub,subind_top)),:,:,:);
      grindlo_tMSAlone_top{ll,ss}.plvspctrm=grindlo_tMSAlone{ll,ss}.plvspctrm(dsearchn(usesub,intersect(usesub,subind_top)),:,:,:);
      grindlo_tMSAlone_top{ll,ss}.plvabs=grindlo_tMSAlone{ll,ss}.plvabs(dsearchn(usesub,intersect(usesub,subind_top)),:,:,:);
      grindlo_tMSAlone_bot{ll,ss}=grindlo_tMSAlone{ll,ss};
      grindlo_tMSAlone_bot{ll,ss}.powspctrm=grindlo_tMSAlone{ll,ss}.powspctrm(dsearchn(usesub,intersect(usesub,subind_bot)),:,:,:);
      grindlo_tMSAlone_bot{ll,ss}.plvspctrm=grindlo_tMSAlone{ll,ss}.plvspctrm(dsearchn(usesub,intersect(usesub,subind_bot)),:,:,:);
      grindlo_tMSAlone_bot{ll,ss}.plvabs=grindlo_tMSAlone{ll,ss}.plvabs(dsearchn(usesub,intersect(usesub,subind_bot)),:,:,:);
    end
    
    cfg=[];
    cfg.neighbours=neighbours;
    cfg.method='montecarlo';
    cfg.numrandomization=2000;
    % cfg.correctm='holm';
    cfg.correctm='cluster';
    cfg.clusteralpha = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan = 2;
    cfg.ivar=1;
    
    
    if statspowflag
      cfg.parameter='powspctrm';
      for ll=soalist
        % top
        cfg.statistic='depsamplesT';
        cfg.uvar=2;
        nsub=size(grindlo_tNulAlone_top{ll,ss}.powspctrm,1);
        if nsub>1
          cfg.design=zeros(2,2*nsub);
          cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
          cfg.design(2,:)=[1:nsub 1:nsub];
          error('changed nonmedsplit to _short_ and _long_')
          cfg.latency=[-.15-audbasemax(5); .55-audbasemax(5)];
          stattl_mc_TacVsNul_top{ll,ss}=ft_freqstatistics(cfg, grindlo_tTacAlone_top{ll,ss}, grindlo_tNulAlone_top{ll,ss});
          if ll>=5
            cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
          else
            cfg.latency=[.15+tacbasemax(ll); .85+tacbasemax(ll)];
          end
          stattl_mc_AudVsNul_top{ll,ss}=ft_freqstatistics(cfg, grindlo_tAudAlone_top{ll,ss}, grindlo_tNulAlone_top{ll,ss});
          cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
          stattl_mc_MSVsNul_top{ll,ss} =ft_freqstatistics(cfg, grindlo_tMSAlone_top{ll,ss}, grindlo_tNulAlone_top{ll,ss});
        else
          stattl_mc_TacVsNul_top{ll,ss}=[];
          stattl_mc_AudVsNul_top{ll,ss}=[];
          stattl_mc_MSVsNul_top{ll,ss} =[];
        end
        
        
        % bottom
        cfg.statistic='depsamplesT';
        cfg.uvar=2;
        nsub=size(grindlo_tNulAlone_bot{ll,ss}.powspctrm,1);
        if nsub>1
          cfg.design=zeros(2,2*nsub);
          cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
          cfg.design(2,:)=[1:nsub 1:nsub];
          
          cfg.latency=[-.15-audbasemax(5); .55-audbasemax(5)];
          stattl_mc_TacVsNul_bot{ll,ss}=ft_freqstatistics(cfg, grindlo_tTacAlone_bot{ll,ss}, grindlo_tNulAlone_bot{ll,ss});
          if ll>=5
            cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
          else
            cfg.latency=[.15+tacbasemax(ll); .85+tacbasemax(ll)];
          end
          stattl_mc_AudVsNul_bot{ll,ss}=ft_freqstatistics(cfg, grindlo_tAudAlone_bot{ll,ss}, grindlo_tNulAlone_bot{ll,ss});
          cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
          stattl_mc_MSVsNul_bot{ll,ss} =ft_freqstatistics(cfg, grindlo_tMSAlone_bot{ll,ss}, grindlo_tNulAlone_bot{ll,ss});
        else
          stattl_mc_TacVsNul_bot{ll,ss}=[];
          stattl_mc_AudVsNul_bot{ll,ss}=[];
          stattl_mc_MSVsNul_bot{ll,ss} =[];
        end
        
        % top vs bottom
        cfg.statistic='indepsamplesT';
        cfg.uvar=[];
        nsubt=size(grindlo_tNulAlone_top{ll,ss}.powspctrm,1);
        nsubb=size(grindlo_tNulAlone_bot{ll,ss}.powspctrm,1);
        if nsubb>1 && nsubt>1
          %             cfg.design=zeros(2,nsubb + nsubt);
          %             cfg.design(1,:)=[ones(1,nsubt) 2*ones(1,nsubb)];
          %             cfg.design(2,:)=[1:nsubt 1:nsubb];
          cfg.design=zeros(1,nsubb + nsubt);
          cfg.design(1,:)=[ones(1,nsubt) 2*ones(1,nsubb)];
          
          cfg.latency=[-.15-audbasemax(5); .55-audbasemax(5)];
          stattl_mc_TacTVsTacB{ll,ss}=ft_freqstatistics(cfg, grindlo_tTacAlone_top{ll,ss}, grindlo_tTacAlone_bot{ll,ss});
          if ll>=5
            cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
          else
            cfg.latency=[.15+tacbasemax(ll); .85+tacbasemax(ll)];
          end
          stattl_mc_AudTVsAudB{ll,ss}=ft_freqstatistics(cfg, grindlo_tAudAlone_top{ll,ss}, grindlo_tAudAlone_bot{ll,ss});
          cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
          stattl_mc_MSTVsMSB{ll,ss}=ft_freqstatistics(cfg, grindlo_tMSAlone_top{ll,ss}, grindlo_tMSAlone_bot{ll,ss});
        else
          stattl_mc_TacTVsTacB{ll,ss}=[];
          stattl_mc_AudTVsAudB{ll,ss}=[];
          stattl_mc_MSTVsMSB{ll,ss}=[];
        end
      end %ll
    end
    
    if statsplvflag
      for ll=soalist
        if itcdepsampflag
          cfg.statistic='depsamplesT';
          cfg.parameter='plvabs';
          cfg.uvar=2;
        else
          
          cfg.statistic='diff_itc';  % use this for all PLV, whether dep or indep
          cfg.uvar=[];
          cfg.parameter='plvspctrm';
          cfg.complex='diffabs'; % default; not sensitive to phase differences between conditions
          cfg.clusterthreshold = 'nonparametric_individual'; % or 'nonparametric_individual' see clusterstat.m
        end
        
        
        
        % top
        nsub=size(grindlo_tNulAlone_top{ll,ss}.powspctrm,1);
        if nsub>1
          cfg.design=zeros(2,2*nsub);
          cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
          cfg.design(2,:)=[1:nsub 1:nsub];
          cfg.latency=[-.15-audbasemax(5); .55-audbasemax(5)];
          if itcdepsampflag
            stattl_mcplv_TacVsNul_top_itcdepT{ll,ss}=ft_freqstatistics(cfg, grindlo_tTacAlone_top{ll,ss}, grindlo_tNulAlone_top{ll,ss});
            if ll>=5
              cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            else
              cfg.latency=[.15+tacbasemax(ll); .85+tacbasemax(ll)];
            end
            stattl_mcplv_AudVsNul_top_itcdepT{ll,ss}=ft_freqstatistics(cfg, grindlo_tAudAlone_top{ll,ss}, grindlo_tNulAlone_top{ll,ss});
            cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            stattl_mcplv_MSVsNul_top_itcdepT{ll,ss} =ft_freqstatistics(cfg, grindlo_tMSAlone_top{ll,ss}, grindlo_tNulAlone_top{ll,ss});
          else
            stattl_mcplv_TacVsNul_top{ll,ss}=ft_freqstatistics(cfg, grindlo_tTacAlone_top{ll,ss}, grindlo_tNulAlone_top{ll,ss});
            if ll>=5
              cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            else
              cfg.latency=[.15+tacbasemax(ll); .85+tacbasemax(ll)];
            end
            stattl_mcplv_AudVsNul_top{ll,ss}=ft_freqstatistics(cfg, grindlo_tAudAlone_top{ll,ss}, grindlo_tNulAlone_top{ll,ss});
            cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            stattl_mcplv_MSVsNul_top{ll,ss} =ft_freqstatistics(cfg, grindlo_tMSAlone_top{ll,ss}, grindlo_tNulAlone_top{ll,ss});
          end
        else
          if itcdepsampflag
            stattl_mcplv_TacVsNul_top_itcdepT{ll,ss}=[];
            stattl_mcplv_AudVsNul_top_itcdepT{ll,ss}=[];
            stattl_mcplv_MSVsNul_top_itcdepT{ll,ss} =[];
          else
            stattl_mcplv_TacVsNul_top{ll,ss}=[];
            stattl_mcplv_AudVsNul_top{ll,ss}=[];
            stattl_mcplv_MSVsNul_top{ll,ss} =[];
          end
        end
        
        % bottom
        nsub=size(grindlo_tNulAlone_bot{ll,ss}.powspctrm,1);
        if nsub>1
          cfg.design=zeros(2,2*nsub);
          cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
          cfg.design(2,:)=[1:nsub 1:nsub];
          cfg.latency=[-.15-audbasemax(5); .55-audbasemax(5)];
          if itcdepsampflag
            stattl_mcplv_TacVsNul_bot_itcdepT{ll,ss}=ft_freqstatistics(cfg, grindlo_tTacAlone_bot{ll,ss}, grindlo_tNulAlone_bot{ll,ss});
            if ll>=5
              cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            else
              cfg.latency=[.15+tacbasemax(ll); .85+tacbasemax(ll)];
            end
            stattl_mcplv_AudVsNul_bot_itcdepT{ll,ss}=ft_freqstatistics(cfg, grindlo_tAudAlone_bot{ll,ss}, grindlo_tNulAlone_bot{ll,ss});
            cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            stattl_mcplv_MSVsNul_bot_itcdepT{ll,ss} =ft_freqstatistics(cfg, grindlo_tMSAlone_bot{ll,ss}, grindlo_tNulAlone_bot{ll,ss});
          else
            stattl_mcplv_TacVsNul_bot{ll,ss}=ft_freqstatistics(cfg, grindlo_tTacAlone_bot{ll,ss}, grindlo_tNulAlone_bot{ll,ss});
            if ll>=5
              cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            else
              cfg.latency=[.15+tacbasemax(ll); .85+tacbasemax(ll)];
            end
            stattl_mcplv_AudVsNul_bot{ll,ss}=ft_freqstatistics(cfg, grindlo_tAudAlone_bot{ll,ss}, grindlo_tNulAlone_bot{ll,ss});
            cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            stattl_mcplv_MSVsNul_bot{ll,ss} =ft_freqstatistics(cfg, grindlo_tMSAlone_bot{ll,ss}, grindlo_tNulAlone_bot{ll,ss});
          end
        else
          if itcdepsampflag
            stattl_mcplv_TacVsNul_bot_itcdepT{ll,ss}=[];
            stattl_mcplv_AudVsNul_bot_itcdepT{ll,ss}=[];
            stattl_mcplv_MSVsNul_bot_itcdepT{ll,ss} =[];
          else
            stattl_mcplv_TacVsNul_bot{ll,ss}=[];
            stattl_mcplv_AudVsNul_bot{ll,ss}=[];
            stattl_mcplv_MSVsNul_bot{ll,ss} =[];
          end
        end
        
        if itcdepsampflag
          cfg.statistic='indepsamplesT';
          cfg.parameter='plvabs';
          cfg.uvar=[];
        else
          cfg.statistic='diff_itc';  % use this for all PLV, whether dep or indep
          cfg.uvar=[];
          cfg.parameter='plvspctrm';
          cfg.complex='diffabs'; % default; not sensitive to phase differences between conditions
          cfg.clusterthreshold = 'nonparametric_individual'; % or 'nonparametric_individual' see clusterstat.m
        end          % top vs bottom
        nsubt=size(grindlo_tNulAlone_top{ll,ss}.powspctrm,1);
        nsubb=size(grindlo_tNulAlone_bot{ll,ss}.powspctrm,1);
        if nsubt>1 && nsubb>1
          %             cfg.design=zeros(2,nsubb + nsubt);
          %             cfg.design(1,:)=[ones(1,nsubt) 2*ones(1,nsubb)];
          %             cfg.design(2,:)=[1:nsubt 1:nsubb];
          cfg.design=zeros(1,nsubb + nsubt);
          cfg.design(1,:)=[ones(1,nsubt) 2*ones(1,nsubb)];
          cfg.latency=[-.15-audbasemax(5); .55-audbasemax(5)];
          if itcdepsampflag
            
            stattl_mcplv_TacTVsTacB_itcdepT{ll,ss}=ft_freqstatistics(cfg, grindlo_tTacAlone_top{ll,ss}, grindlo_tTacAlone_bot{ll,ss});
            if ll>=5
              cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            else
              cfg.latency=[.15+tacbasemax(ll); .85+tacbasemax(ll)];
            end
            stattl_mcplv_AudTVsAudB_itcdepT{ll,ss}=ft_freqstatistics(cfg, grindlo_tAudAlone_top{ll,ss}, grindlo_tAudAlone_bot{ll,ss});
            cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            stattl_mcplv_MSTVsMSB_itcdepT{ll,ss}=ft_freqstatistics(cfg, grindlo_tMSAlone_top{ll,ss}, grindlo_tMSAlone_bot{ll,ss});
          else
            stattl_mcplv_TacTVsTacB{ll,ss}=ft_freqstatistics(cfg, grindlo_tTacAlone_top{ll,ss}, grindlo_tTacAlone_bot{ll,ss});
            if ll>=5
              cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            else
              cfg.latency=[.15+tacbasemax(ll); .85+tacbasemax(ll)];
            end
            stattl_mcplv_AudTVsAudB{ll,ss}=ft_freqstatistics(cfg, grindlo_tAudAlone_top{ll,ss}, grindlo_tAudAlone_bot{ll,ss});
            cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            stattl_mcplv_MSTVsMSB{ll,ss}=ft_freqstatistics(cfg, grindlo_tMSAlone_top{ll,ss}, grindlo_tMSAlone_bot{ll,ss});
          end
        else
        end
      end %ll
      
    end
    
    % NOTE: without trialkc label, it refers to trialkc=0;
    %     save([edir 'stats_TFR_UniMSNul_' num2str(tt) num2str(ss) num2str(sleep) '.mat'],'stat*');
    %     save([edir 'stats_TFR_UniMSNul_' num2str(tt) num2str(ss) num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'stat*','-append');
    if statspowflag || statsplvflag
      if itcdepsampflag
        save([edir 'stats_TFR_UniMSNul_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '_itc.mat'],'stat*','-append');
      else
        save([edir 'stats_TFR_UniMSNul_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'stat*','-append');
      end
    end
  end % medsplitflag
end % ss

%   end % tt
% end % sleep

%%
% figures for paper:
trialkc=-1;
load(['stats_TFR_UniMSNul_sleep1_trialkc' num2str(trialkc) '.mat'])

% Auditory vs Null
figure(1);imagesc(squeeze(sum(stattl_mc_AudVsNul{5,12}.mask,1)));axis xy
ax=gca;
ax.XTick=ax.XTick+1;
ax.XTickLabel={'0.1' '0.2' '0.3' '0.4' '0.5' '0.6' '0.7'};
ax.YTick=[4 8 12];
ax.YTickLabel={'10' '18' '26'}
title('# channels significant: Aud. vs Null')
colorbar;   caxis([0 63])
print(1,[fdir 'TFP_AudVsNul_numchan_trialkc' num2str(trialkc) '.png'],'-dpng');
print(1,[fdir 'TFP_AudVsNul_numchan_trialkc' num2str(trialkc) '.eps'],'-depsc');

figure(2);imagesc(squeeze(mean(stattl_mc_AudVsNul{5,12}.stat,1)));axis xy
ax=gca;
ax.XTick=ax.XTick+1;
ax.XTickLabel={'0.1' '0.2' '0.3' '0.4' '0.5' '0.6' '0.7'};
ax.YTick=[4 8 12];
ax.YTickLabel={'10' '18' '26'}
title('T-value contrast: Aud. vs Null')
colorbar;  caxis([-4 4])
print(2,[fdir 'TFP_AudVsNul_meanstat_trialkc' num2str(trialkc) '.png'],'-dpng');
print(2,[fdir 'TFP_AudVsNul_meanstat_trialkc' num2str(trialkc) '.eps'],'-depsc');

figure(3);imagesc(squeeze(sum(stattl_mc_TacVsNul{5,12}.mask,1)));axis xy
ax=gca;
ax.XTick=ax.XTick+1;
ax.XTickLabel={'0.1' '0.2' '0.3' '0.4' '0.5' '0.6' '0.7'};
ax.YTick=[4 8 12];
ax.YTickLabel={'10' '18' '26'}
title('# channels significant: Aud. vs Null')
colorbar;   caxis([0 63])
print(3,[fdir 'TFP_TacVsNul_numchan_trialkc' num2str(trialkc) '.png'],'-dpng');
print(3,[fdir 'TFP_TacVsNul_numchan_trialkc' num2str(trialkc) '.eps'],'-depsc');


figure(4);imagesc(squeeze(mean(stattl_mc_TacVsNul{5,12}.stat,1)));axis xy
ax=gca;
ax.XTick=ax.XTick+1;
ax.XTickLabel={'0.1' '0.2' '0.3' '0.4' '0.5' '0.6' '0.7'};
ax.YTick=[4 8 12];
ax.YTickLabel={'10' '18' '26'}
title('T-value contrast: Tac. vs Null')
colorbar;   caxis([-4 4])
print(4,[fdir 'TFP_TacVsNul_meanstat_trialkc' num2str(trialkc) '.png'],'-dpng');
print(4,[fdir 'TFP_TacVsNul_meanstat_trialkc' num2str(trialkc) '.eps'],'-depsc');

% Auditory vs Null
% figure(5);imagesc(squeeze(sum(stattl_mcplv_AudVsNul{5,12}.mask,1)));axis xy
figure(5);imagesc(squeeze(sum(stattl_mcplv_AudVsNul_itcdepT{5,12}.mask,1)));axis xy
ax=gca;
ax.XTick=ax.XTick+1;
ax.XTickLabel={'0.1' '0.2' '0.3' '0.4' '0.5' '0.6' '0.7'};
ax.YTick=[4 8 12];
ax.YTickLabel={'10' '18' '26'}
title('# channels significant: Aud. vs Null')
colorbar;  caxis([0 63])
print(5,[fdir 'PLF_AudVsNul_numchan_trialkc' num2str(trialkc) '.png'],'-dpng');
print(5,[fdir 'PLF_AudVsNul_numchan_trialkc' num2str(trialkc) '.eps'],'-depsc');

% figure(6);imagesc(squeeze(mean(stattl_mcplv_AudVsNul{5,12}.stat,1)));axis xy
figure(6);imagesc(squeeze(mean(stattl_mcplv_AudVsNul_itcdepT{5,12}.stat,1)));axis xy
ax=gca;
ax.XTick=ax.XTick+1;
ax.XTickLabel={'0.1' '0.2' '0.3' '0.4' '0.5' '0.6' '0.7'};
ax.YTick=[4 8 12];
ax.YTickLabel={'10' '18' '26'}
title('T-value contrast: Aud. vs Null')
colorbar;  caxis([-4 4])
print(6,[fdir 'PLF_AudVsNul_meanstat_trialkc' num2str(trialkc) '.png'],'-dpng');
print(6,[fdir 'PLF_AudVsNul_meanstat_trialkc' num2str(trialkc) '.eps'],'-depsc');

% Tactile vs Null
% figure(7);imagesc(squeeze(sum(stattl_mcplv_TacVsNul{5,12}.mask,1)));axis xy
figure(7);imagesc(squeeze(sum(stattl_mcplv_TacVsNul_itcdepT{5,12}.mask,1)));axis xy
ax=gca;
ax.XTick=ax.XTick+1;
ax.XTickLabel={'0.1' '0.2' '0.3' '0.4' '0.5' '0.6' '0.7'};
ax.YTick=[4 8 12];
ax.YTickLabel={'10' '18' '26'}
title('# channels significant: Tac. vs Null')
colorbar;  caxis([0 63])
print(7,[fdir 'PLF_TacVsNul_numchan_trialkc' num2str(trialkc) '.png'],'-dpng');
print(7,[fdir 'PLF_TacVsNul_numchan_trialkc' num2str(trialkc) '.eps'],'-depsc');

figure(8);imagesc(squeeze(mean(stattl_mcplv_TacVsNul_itcdepT{5,12}.stat,1)));axis xy
ax=gca;
ax.XTick=ax.XTick+1;
ax.XTickLabel={'0.1' '0.2' '0.3' '0.4' '0.5' '0.6' '0.7'};
ax.YTick=[4 8 12];
ax.YTickLabel={'10' '18' '26'}
title('T-value contrast: Tac. vs Null')
colorbar;  caxis([-4 4])
print(8,[fdir 'PLF_TacVsNul_meanstat_trialkc' num2str(trialkc) '.png'],'-dpng');
print(8,[fdir 'PLF_TacVsNul_meanstat_trialkc' num2str(trialkc) '.eps'],'-depsc');

% line plots
