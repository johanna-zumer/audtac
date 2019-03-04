eeg_legomagic_preamble

%% MultSens contrast with Median-split

load elec1010_neighb.mat

soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];
tacbasemax=min(-.15,soades-.15);
audbasemax=min(-.15,fliplr(soades)-.15);

ylimlo=[4 7; 8 12; 14 30];
ylimhi=[30 45; 55 80];

soadesc={'Aud first by 500ms' '' 'Aud first by 70ms' 'Aud first by 20ms' 'Simultaneous' 'Tac first by 20ms' 'Tac first by 70ms' '' 'Tac first by 500ms'};
comb2flag=1;
mcseed=13;  % montecarlo cfg.randomseed
resetusetr=0;
soalist=[1 3 4 5 6 7 9];
chanuse_sleep0={'all' '-F4'};
chanuse_sleep1={'all' '-AF7' '-AF3' '-Fp1'};

statspowflag=1;
statsplvflag=1;
savegrindflag=1;
itcdepsampflag=1;


tt=3;
sleep=1;
clearvars -except ll tt sub *dir ii*use sleep *flag figind soades* soalist chanuse* ylim* neigh* *basemax mcseed
if sleep
  usetr=1; % 1 for sleep
  chanuse=chanuse_sleep1;
  subuseall=iiBuse;
  iteruse=11;
  trialkc=-1;
  %       ssuse=12;
  sleepcond='Sleep N2';
else
  usetr=3; % 1 with 27, or 2 or 3 with 31 or 32
  chanuse=chanuse_sleep0;
  subuseall=setdiff(iiSuse,[]);
  iteruse=27;
  trialkc=-1;
  %       ssuse=10; % awake
  sleepcond='Awake W';
end
%     ss=ssuse;

submin=subuseall(1)-1;
subuseind=0;
iteruse
for ii=subuseall
  cd([edir sub{ii} ])
%   try
    load(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iteruse) '_trialkc' num2str(trialkc) '.mat'])
%   catch
%     load(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_tt' num2str(tt) '_trialkc' num2str(trialkc) '.mat'])
%   end
  subuseind=subuseind+1;
  
  for ss=10:12
    for ll=soalist
      
      numtrt(ll,tt,ss,subuseind)=numt_trials(ll,tt,ss);
      if numtrt(ll,tt,ss,subuseind)<20
        subuse(ll,ss,subuseind)=nan;
      else
        subuse(ll,ss,subuseind)=1;
      end
      
      %       if sleep==0
      %         ssuse=10; % awake
      %         sleepcond='Awake W';
      %       elseif sleep==1
      %         ssuse=12; % N2
      %         sleepcond='Sleep N2';
      %       end
      try
        if usetr==2
          % load usetr=0 here; then later down load usetr=2
          tkt=load(['trialkeptTFR_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iteruse) '_usetr' num2str(3) '_trialkc' num2str(trialkc) '.mat']);
        else
          tkt=load(['trialkeptTFR_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iteruse) '_usetr' num2str(usetr) '.mat']);
        end
      catch
        tkt=load(['trialkeptTFR_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iteruse) '_usetr' num2str(usetr) '_trialkc' num2str(trialkc) '.mat']);
      end
      
      
      %       ss=ssuse;
      %         for ss=ssuse
      subuse=subuseall; % reset to all for each sleep stage
      numtrt(ll,tt,ss,ii-submin)=tkt.numcondtfinal(ll,tt,ss);
      
      if min(numtrt(ll,tt,ss,ii-submin),[],1)<20 % what is best number to use here?
        subuse=setdiff(subuse,ii);
      else
        %         subuseind=subuseind+1;
        freqloall_tacPaud_comb1{ll,ss}{subuseind,1}=freqlo_tacPaud_comb{ll,tt,ss,1};
        freqloall_tacMSpN_comb1{ll,ss}{subuseind,1}=freqlo_tacMSpN_comb{ll,tt,ss,1};
        freqhiall_tacPaud_comb1{ll,ss}{subuseind,1}=freqhi_tacPaud_comb{ll,tt,ss,1};
        freqhiall_tacMSpN_comb1{ll,ss}{subuseind,1}=freqhi_tacMSpN_comb{ll,tt,ss,1};
        freqloall_tacPaud_comb1{ll,ss}{subuseind,1}.dimord='chan_freq_time';
        freqloall_tacMSpN_comb1{ll,ss}{subuseind,1}.dimord='chan_freq_time';
        freqhiall_tacPaud_comb1{ll,ss}{subuseind,1}.dimord='chan_freq_time';
        freqhiall_tacMSpN_comb1{ll,ss}{subuseind,1}.dimord='chan_freq_time';
        freqloall_tacPaud_comb1{ll,ss}{subuseind,1}.plvabs=abs(freqlo_tacPaud_comb{ll,tt,ss,1}.plvspctrm);
        freqloall_tacMSpN_comb1{ll,ss}{subuseind,1}.plvabs=abs(freqlo_tacMSpN_comb{ll,tt,ss,1}.plvspctrm);
        freqhiall_tacPaud_comb1{ll,ss}{subuseind,1}.plvabs=abs(freqhi_tacPaud_comb{ll,tt,ss,1}.plvspctrm);
        freqhiall_tacMSpN_comb1{ll,ss}{subuseind,1}.plvabs=abs(freqhi_tacMSpN_comb{ll,tt,ss,1}.plvspctrm);
        if comb2flag
          freqloall_tacPaud_comb2{ll,ss}{subuseind,1}=freqlo_tacPaud_comb{ll,tt,ss,2};
          freqloall_tacMSpN_comb2{ll,ss}{subuseind,1}=freqlo_tacMSpN_comb{ll,tt,ss,2};
          freqhiall_tacPaud_comb2{ll,ss}{subuseind,1}=freqhi_tacPaud_comb{ll,tt,ss,2};
          freqhiall_tacMSpN_comb2{ll,ss}{subuseind,1}=freqhi_tacMSpN_comb{ll,tt,ss,2};
          freqloall_tacPaud_comb2{ll,ss}{subuseind,1}.dimord='chan_freq_time';
          freqloall_tacMSpN_comb2{ll,ss}{subuseind,1}.dimord='chan_freq_time';
          freqhiall_tacPaud_comb2{ll,ss}{subuseind,1}.dimord='chan_freq_time';
          freqhiall_tacMSpN_comb2{ll,ss}{subuseind,1}.dimord='chan_freq_time';
          freqloall_tacPaud_comb2{ll,ss}{subuseind,1}.plvabs=abs(freqlo_tacPaud_comb{ll,tt,ss,2}.plvspctrm);
          freqloall_tacMSpN_comb2{ll,ss}{subuseind,1}.plvabs=abs(freqlo_tacMSpN_comb{ll,tt,ss,2}.plvspctrm);
          freqhiall_tacPaud_comb2{ll,ss}{subuseind,1}.plvabs=abs(freqhi_tacPaud_comb{ll,tt,ss,2}.plvspctrm);
          freqhiall_tacMSpN_comb2{ll,ss}{subuseind,1}.plvabs=abs(freqhi_tacMSpN_comb{ll,tt,ss,2}.plvspctrm);
        end
      end
      
    end % ll
  end % ss
  clear freqlo_* freqhi_*
end % ii
subuseindfinal=subuseind;

for ii=1:subuseind
  for ss=10:12
    for ll=soalist
      cfg=[];
      cfg.operation='subtract';
      cfg.parameter={'powspctrm' 'plvspctrm' 'plvabs'};
      for iterind=1
        if ii<=size(freqloall_tacPaud_comb1{ll,ss},1) && ~isempty(freqloall_tacPaud_comb1{ll,ss}{ii,iterind})
          freqloall_TPA_MSPN_comb1{ll,ss}{ii,iterind}=ft_math(cfg,freqloall_tacPaud_comb1{ll,ss}{ii,iterind},freqloall_tacMSpN_comb1{ll,ss}{ii,iterind});
          freqhiall_TPA_MSPN_comb1{ll,ss}{ii,iterind}=ft_math(cfg,freqhiall_tacPaud_comb1{ll,ss}{ii,iterind},freqhiall_tacMSpN_comb1{ll,ss}{ii,iterind});
        end
      end
      if comb2flag
        for iterind=1
          if ii<=size(freqloall_tacPaud_comb2{ll,ss},1) && ~isempty(freqloall_tacPaud_comb2{ll,ss}{ii,iterind})
            freqloall_TPA_MSPN_comb2{ll,ss}{ii,iterind}=ft_math(cfg,freqloall_tacPaud_comb2{ll,ss}{ii,iterind},freqloall_tacMSpN_comb2{ll,ss}{ii,iterind});
            freqhiall_TPA_MSPN_comb2{ll,ss}{ii,iterind}=ft_math(cfg,freqhiall_tacPaud_comb2{ll,ss}{ii,iterind},freqhiall_tacMSpN_comb2{ll,ss}{ii,iterind});
          end
        end
      end
      
    end % ll
  end % ss
  
end % ii


load([edir 'sortTacWP200.mat'],'iiuse_*halfWwideP2N1');
subind_top=dsearchn(iiBuse',iiuse_tophalfWwideP2N1');
subind_bot=dsearchn(iiBuse',iiuse_bothalfWwideP2N1');
sublab=nan(19,1);
sublab(subind_top)=1;
sublab(subind_bot)=-1;

cfg=[];
cfg.parameter={'powspctrm' 'plvspctrm' 'plvabs'};


for ss=10:12
  for ll=soalist
    clear cellfull
    for ii=1:length(freqloall_TPA_MSPN_comb1{ll,ss}),cellfull(ii)=~isempty(freqloall_TPA_MSPN_comb1{ll,ss}{ii});end
    sublabuse{ll,ss}=sublab(cellfull);
    
    cfg.keepindividual='yes';
    grindlo_TPA_MSPN_comb1{ll,ss}=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb1{ll,ss}{cellfull});
    grindhi_TPA_MSPN_comb1{ll,ss}=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb1{ll,ss}{cellfull});
    if comb2flag
      grindlo_TPA_MSPN_comb2{ll,ss}=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb2{ll,ss}{cellfull});
      grindhi_TPA_MSPN_comb2{ll,ss}=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb2{ll,ss}{cellfull});
    end    
    grindlo_tacPaud_comb1{ll,ss}=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb1{ll,ss}{cellfull});
    grindlo_tacMSpN_comb1{ll,ss}=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb1{ll,ss}{cellfull});
    grindhi_tacPaud_comb1{ll,ss}=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb1{ll,ss}{cellfull});
    grindhi_tacMSpN_comb1{ll,ss}=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb1{ll,ss}{cellfull});
    if comb2flag
      grindlo_tacPaud_comb2{ll,ss}=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb2{ll,ss}{cellfull});
      grindlo_tacMSpN_comb2{ll,ss}=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb2{ll,ss}{cellfull});
      grindhi_tacPaud_comb2{ll,ss}=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb2{ll,ss}{cellfull});
      grindhi_tacMSpN_comb2{ll,ss}=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb2{ll,ss}{cellfull});
    end
    
    cfg.keepindividual='no';
    gravelo_TPA_MSPN_comb1{ll,ss}=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb1{ll,ss}{cellfull});
    gravehi_TPA_MSPN_comb1{ll,ss}=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb1{ll,ss}{cellfull});
    if comb2flag
      gravelo_TPA_MSPN_comb2{ll,ss}=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb2{ll,ss}{cellfull});
      gravehi_TPA_MSPN_comb2{ll,ss}=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb2{ll,ss}{cellfull});
    end    
    gravelo_tacPaud_comb1{ll,ss}=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb1{ll,ss}{cellfull});
    gravelo_tacMSpN_comb1{ll,ss}=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb1{ll,ss}{cellfull});
    gravehi_tacPaud_comb1{ll,ss}=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb1{ll,ss}{cellfull});
    gravehi_tacMSpN_comb1{ll,ss}=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb1{ll,ss}{cellfull});
    if comb2flag
      gravelo_tacPaud_comb2{ll,ss}=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb2{ll,ss}{cellfull});
      gravelo_tacMSpN_comb2{ll,ss}=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb2{ll,ss}{cellfull});
      gravehi_tacPaud_comb2{ll,ss}=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb2{ll,ss}{cellfull});
      gravehi_tacMSpN_comb2{ll,ss}=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb2{ll,ss}{cellfull});
    end
    
    freqloall_TPA_MSPN_comb1{ll,ss}=[];
    freqloall_TPA_MSPN_comb2{ll,ss}=[];
    freqhiall_TPA_MSPN_comb1{ll,ss}=[];
    freqhiall_TPA_MSPN_comb2{ll,ss}=[];
    freqloall_tacPaud_comb1{ll,ss}=[];
    freqloall_tacPaud_comb2{ll,ss}=[];
    freqhiall_tacPaud_comb1{ll,ss}=[];
    freqhiall_tacPaud_comb2{ll,ss}=[];
    freqloall_tacMSpN_comb1{ll,ss}=[];
    freqloall_tacMSpN_comb2{ll,ss}=[];
    freqhiall_tacMSpN_comb1{ll,ss}=[];
    freqhiall_tacMSpN_comb2{ll,ss}=[];
  end % ll
end % ss


clear freq*all*comb*


for ss=10:12
  for ll=soalist
    %     usesub=find(squeeze(~isnan(subuse(ll,ss,:))));
    
    % Make grind and grave for top/bottom subselections    
    %[grindout,graveout]=grindgravetopbot(grindorig,sublist,subval)
    
    [grindlo_tacPaud_top{ll,ss,1},gravelo_tacPaud_top{ll,ss,1}]=grindgravetopbot(grindlo_tacPaud_comb1{ll,ss},sublabuse{ll,ss},1);
    [grindlo_tacMSpN_top{ll,ss,1},gravelo_tacMSpN_top{ll,ss,1}]=grindgravetopbot(grindlo_tacMSpN_comb1{ll,ss},sublabuse{ll,ss},1);
    [grindhi_tacPaud_top{ll,ss,1},gravehi_tacPaud_top{ll,ss,1}]=grindgravetopbot(grindhi_tacPaud_comb1{ll,ss},sublabuse{ll,ss},1);
    [grindhi_tacMSpN_top{ll,ss,1},gravehi_tacMSpN_top{ll,ss,1}]=grindgravetopbot(grindhi_tacMSpN_comb1{ll,ss},sublabuse{ll,ss},1);
    
    [grindlo_tacPaud_bot{ll,ss,1},gravelo_tacPaud_bot{ll,ss,1}]=grindgravetopbot(grindlo_tacPaud_comb1{ll,ss},sublabuse{ll,ss},-1);
    [grindlo_tacMSpN_bot{ll,ss,1},gravelo_tacMSpN_bot{ll,ss,1}]=grindgravetopbot(grindlo_tacMSpN_comb1{ll,ss},sublabuse{ll,ss},-1);
    [grindhi_tacPaud_bot{ll,ss,1},gravehi_tacPaud_bot{ll,ss,1}]=grindgravetopbot(grindhi_tacPaud_comb1{ll,ss},sublabuse{ll,ss},-1);
    [grindhi_tacMSpN_bot{ll,ss,1},gravehi_tacMSpN_bot{ll,ss,1}]=grindgravetopbot(grindhi_tacMSpN_comb1{ll,ss},sublabuse{ll,ss},-1);
    
    if comb2flag
      [grindlo_tacPaud_top{ll,ss,2},gravelo_tacPaud_top{ll,ss,2}]=grindgravetopbot(grindlo_tacPaud_comb2{ll,ss},sublabuse{ll,ss},1);
      [grindlo_tacMSpN_top{ll,ss,2},gravelo_tacMSpN_top{ll,ss,2}]=grindgravetopbot(grindlo_tacMSpN_comb2{ll,ss},sublabuse{ll,ss},1);
      [grindhi_tacPaud_top{ll,ss,2},gravehi_tacPaud_top{ll,ss,2}]=grindgravetopbot(grindhi_tacPaud_comb2{ll,ss},sublabuse{ll,ss},1);
      [grindhi_tacMSpN_top{ll,ss,2},gravehi_tacMSpN_top{ll,ss,2}]=grindgravetopbot(grindhi_tacMSpN_comb2{ll,ss},sublabuse{ll,ss},1);
      
      [grindlo_tacPaud_bot{ll,ss,2},gravelo_tacPaud_bot{ll,ss,2}]=grindgravetopbot(grindlo_tacPaud_comb2{ll,ss},sublabuse{ll,ss},-1);
      [grindlo_tacMSpN_bot{ll,ss,2},gravelo_tacMSpN_bot{ll,ss,2}]=grindgravetopbot(grindlo_tacMSpN_comb2{ll,ss},sublabuse{ll,ss},-1);
      [grindhi_tacPaud_bot{ll,ss,2},gravehi_tacPaud_bot{ll,ss,2}]=grindgravetopbot(grindhi_tacPaud_comb2{ll,ss},sublabuse{ll,ss},-1);
      [grindhi_tacMSpN_bot{ll,ss,2},gravehi_tacMSpN_bot{ll,ss,2}]=grindgravetopbot(grindhi_tacMSpN_comb2{ll,ss},sublabuse{ll,ss},-1);
      
    end
    
    % begin stats
    cfg=[];
    cfg.neighbours=neighbours;
    cfg.method='montecarlo';
    cfg.numrandomization=2000;
    cfg.correctm='cluster';
    cfg.clusteralpha = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan = 2;
    cfg.ivar=1;
    cfg.latency=[-.15-audbasemax(ll); 1.05-audbasemax(ll)];
    
    
    if comb2flag
      combmax=2;
    else
      combmax=1;
    end
    
    for comb=1:combmax
      nsubt=size(grindlo_tacPaud_top{ll,ss,comb}.powspctrm,1);
      nsubb=size(grindlo_tacPaud_bot{ll,ss,comb}.powspctrm,1);
      
      if statspowflag          % power
        
        cfg.parameter='powspctrm';
        cfg.statistic='depsamplesT';
        cfg.uvar=2;
        
        if nsubt>1
          cfg.design=zeros(2,2*nsubt);
          cfg.design(1,:)=[ones(1,nsubt) 2*ones(1,nsubt)];
          cfg.design(2,:)=[1:nsubt 1:nsubt];
          stattl_mc_TpAMSpN_top{ll,ss,comb}=ft_freqstatistics(cfg, grindlo_tacPaud_top{ll,ss,comb}, grindlo_tacMSpN_top{ll,ss,comb});
          statth_mc_TpAMSpN_top{ll,ss,comb}=ft_freqstatistics(cfg, grindhi_tacPaud_top{ll,ss,comb}, grindhi_tacMSpN_top{ll,ss,comb});
        else
          stattl_mc_TpAMSpN_top{ll,ss,comb}=[];
          statth_mc_TpAMSpN_top{ll,ss,comb}=[];
        end
        
        if nsubb>1
          cfg.design=zeros(2,2*nsubb);
          cfg.design(1,:)=[ones(1,nsubb) 2*ones(1,nsubb)];
          cfg.design(2,:)=[1:nsubb 1:nsubb];
          stattl_mc_TpAMSpN_bot{ll,ss,comb}=ft_freqstatistics(cfg, grindlo_tacPaud_bot{ll,ss,comb}, grindlo_tacMSpN_bot{ll,ss,comb});
          statth_mc_TpAMSpN_bot{ll,ss,comb}=ft_freqstatistics(cfg, grindhi_tacPaud_bot{ll,ss,comb}, grindhi_tacMSpN_bot{ll,ss,comb});
        else
          stattl_mc_TpAMSpN_bot{ll,ss,comb}=[];
          statth_mc_TpAMSpN_bot{ll,ss,comb}=[];
        end
        
        cfg.statistic='indepsamplesT';
        cfg.uvar=[];
        if nsubt>1 && nsubb>1
          cfg.design=zeros(1,nsubb + nsubt);
          cfg.design(1,:)=[ones(1,nsubt) 2*ones(1,nsubb)];
          stattl_mc_TpAMSpN_TvB{ll,ss,comb}=ft_freqstatistics(cfg, grindlo_tacMSpN_top{ll,ss,comb}, grindlo_tacMSpN_bot{ll,ss,comb});
          statth_mc_TpAMSpN_TvB{ll,ss,comb}=ft_freqstatistics(cfg, grindhi_tacMSpN_top{ll,ss,comb}, grindhi_tacMSpN_bot{ll,ss,comb});
        else
          stattl_mc_TpAMSpN_TvB{ll,ss,comb}=[];
          statth_mc_TpAMSpN_TvB{ll,ss,comb}=[];
        end
      end % statspowflag
      
      
      if statsplvflag  % plf / plv / itc
        
        if itcdepsampflag
          cfg.parameter='plvabs';
          cfg.statistic='depsamplesT';
          cfg.uvar=2;
        else
          cfg.parameter='plvspctrm';
          cfg.statistic='diff_itc';
          cfg.uvar=[];
          cfg.complex='diffabs'; % default; not sensitive to phase differences between conditions
          cfg.clusterthreshold = 'nonparametric_individual'; % or 'nonparametric_individual' see clusterstat.m
        end
        
        
        if nsubt>1
          cfg.design=zeros(2,2*nsubt);
          cfg.design(1,:)=[ones(1,nsubt) 2*ones(1,nsubt)];
          cfg.design(2,:)=[1:nsubt 1:nsubt];
          if itcdepsampflag
            stattl_mcplv_TpAMSpN_top_itcdepT{ll,ss,comb}=ft_freqstatistics(cfg, grindlo_tacPaud_top{ll,ss,comb}, grindlo_tacMSpN_top{ll,ss,comb});
            statth_mcplv_TpAMSpN_top_itcdepT{ll,ss,comb}=ft_freqstatistics(cfg, grindhi_tacPaud_top{ll,ss,comb}, grindhi_tacMSpN_top{ll,ss,comb});
          else
            stattl_mcplv_TpAMSpN_top{ll,ss,comb}=ft_freqstatistics(cfg, grindlo_tacPaud_top{ll,ss,comb}, grindlo_tacMSpN_top{ll,ss,comb});
            statth_mcplv_TpAMSpN_top{ll,ss,comb}=ft_freqstatistics(cfg, grindhi_tacPaud_top{ll,ss,comb}, grindhi_tacMSpN_top{ll,ss,comb});
          end
        end
        
        if nsubb>1
          cfg.design=zeros(2,2*nsubb);
          cfg.design(1,:)=[ones(1,nsubb) 2*ones(1,nsubb)];
          cfg.design(2,:)=[1:nsubb 1:nsubb];
          if itcdepsampflag
            stattl_mcplv_TpAMSpN_bot_itcdepT{ll,ss,comb}=ft_freqstatistics(cfg, grindlo_tacPaud_bot{ll,ss,comb}, grindlo_tacMSpN_bot{ll,ss,comb});
            statth_mcplv_TpAMSpN_bot_itcdepT{ll,ss,comb}=ft_freqstatistics(cfg, grindhi_tacPaud_bot{ll,ss,comb}, grindhi_tacMSpN_bot{ll,ss,comb});
          else
            stattl_mcplv_TpAMSpN_bot{ll,ss,comb}=ft_freqstatistics(cfg, grindlo_tacPaud_bot{ll,ss,comb}, grindlo_tacMSpN_bot{ll,ss,comb});
            statth_mcplv_TpAMSpN_bot{ll,ss,comb}=ft_freqstatistics(cfg, grindhi_tacPaud_bot{ll,ss,comb}, grindhi_tacMSpN_bot{ll,ss,comb});
          end
        end
        
        if itcdepsampflag
          cfg.parameter='plvabs';
          cfg.statistic='indepsamplesT';
          cfg.uvar=[];
        else
          cfg.parameter='plvspctrm';
          cfg.statistic='diff_itc';
          cfg.uvar=[];
          cfg.complex='diffabs'; % default; not sensitive to phase differences between conditions
          cfg.clusterthreshold = 'nonparametric_individual'; % or 'nonparametric_individual' see clusterstat.m
        end
        
        if nsubt>1 && nsubb>1
          cfg.design=zeros(1,nsubb + nsubt);
          cfg.design(1,:)=[ones(1,nsubt) 2*ones(1,nsubb)];
          if itcdepsampflag
            stattl_mcplv_TpAMSpN_TvB_itcdepT{ll,ss,comb}=ft_freqstatistics(cfg, grindlo_tacMSpN_top{ll,ss,comb}, grindlo_tacMSpN_bot{ll,ss,comb});
            statth_mcplv_TpAMSpN_TvB_itcdepT{ll,ss,comb}=ft_freqstatistics(cfg, grindhi_tacMSpN_top{ll,ss,comb}, grindhi_tacMSpN_bot{ll,ss,comb});
          else
            stattl_mcplv_TpAMSpN_TvB{ll,ss,comb}=ft_freqstatistics(cfg, grindlo_tacMSpN_top{ll,ss,comb}, grindlo_tacMSpN_bot{ll,ss,comb});
            statth_mcplv_TpAMSpN_TvB{ll,ss,comb}=ft_freqstatistics(cfg, grindhi_tacMSpN_top{ll,ss,comb}, grindhi_tacMSpN_bot{ll,ss,comb});
          end
        end
      end % statsplvflag
      
    end % comb
    
  end % ll
end % ss


if itcdepsampflag
  try
    save([edir 'stats_TFR_MScon_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '_itc.mat'],'stat*','grave*','-append');
  catch
    save([edir 'stats_TFR_MScon_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '_itc.mat'],'stat*','grave*');
  end
else
  try
    save([edir 'stats_TFR_MScon_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'stat*','grave*','-append');
  catch
    save([edir 'stats_TFR_MScon_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'stat*','grave*');
  end
end

%% 
% figures to assess above
load stats_TFR_MScon_sleep1_trialkc-1_itc % takes a long time!

soalist=[1 3 4 5 6 7 9];
chanplot{1}={'Fz' 'Cz' 'F1' 'F2' 'FC1' 'FC2' 'C1' 'C2'}; % frontocentral
chpl{1}=match_str(stattl_mcplv_TpAMSpN_bot_itcdepT{7,12,1}.label,chanplot{1});
fb{1}=1:2;fb{2}=3:5;fb{3}=6:10;fb{4}=1:10;
fname{1}='theta';fname{2}='alpha';fname{3}='beta';fname{4}='gamma';

for ll=soalist
  for ff=1:3
    figure(ll+10*ff);plot(squeeze(mean(mean(stattl_mcplv_TpAMSpN_bot_itcdepT{ll,12,1}.stat(chpl{1},fb{ff},:),1),2)),'b');ylim([-1.5 1.5])
    hold on;
    plot(squeeze(mean(mean(stattl_mcplv_TpAMSpN_top_itcdepT{ll,12,1}.stat(chpl{1},fb{ff},:),1),2)),'k');ylim([-1.5 1.5])
    plot(squeeze(mean(mean(stattl_mcplv_TpAMSpN_TvB_itcdepT{ll,12,1}.stat(chpl{1},fb{ff},:),1),2)),'r');ylim([-1.5 1.5])
    sigpoints=find(squeeze(min(min(stattl_mcplv_TpAMSpN_bot_itcdepT{ll,12,1}.prob(:,fb{ff},:),[],1),[],2))<.05);
    if ~isempty(sigpoints)
      plot(sigpoints,squeeze(mean(mean(stattl_mcplv_TpAMSpN_bot_itcdepT{ll,12,1}.stat(chpl{1},fb{ff},sigpoints),1),2)),'b*');
    end
    sigpoints=find(squeeze(min(min(stattl_mcplv_TpAMSpN_top_itcdepT{ll,12,1}.prob(:,fb{ff},:),[],1),[],2))<.05);
    if ~isempty(sigpoints)
      plot(sigpoints,squeeze(mean(mean(stattl_mcplv_TpAMSpN_top_itcdepT{ll,12,1}.stat(chpl{1},fb{ff},sigpoints),1),2)),'k*');
    end
    sigpoints=find(squeeze(min(min(stattl_mcplv_TpAMSpN_TvB_itcdepT{ll,12,1}.prob(:,fb{ff},:),[],1),[],2))<.05);
    if ~isempty(sigpoints)
      plot(sigpoints,squeeze(mean(mean(stattl_mcplv_TpAMSpN_TvB_itcdepT{ll,12,1}.stat(chpl{1},fb{ff},sigpoints),1),2)),'r*');
    end
    legend({'Bottom' 'Top' 'TvB'})
    title(fname{ff})
    print(ll+10*ff,[fdir 'statl_plvitc_' fname{ff} '_' num2str(ll) '.eps'],'-painters','-depsc');
  end
end

for ll=soalist
  for ff=1:3
    figure(ll+10*ff+30);plot(squeeze(mean(mean(stattl_mc_TpAMSpN_bot{ll,12,1}.stat(chpl{1},fb{ff},:),1),2)),'b');
    hold on;
    plot(squeeze(mean(mean(stattl_mc_TpAMSpN_top{ll,12,1}.stat(chpl{1},fb{ff},:),1),2)),'k');
    plot(squeeze(mean(mean(stattl_mc_TpAMSpN_TvB{ll,12,1}.stat(chpl{1},fb{ff},:),1),2)),'r');
    sigpoints=find(squeeze(min(min(stattl_mc_TpAMSpN_bot{ll,12,1}.prob(:,fb{ff},:),[],1),[],2))<.05);
    if ~isempty(sigpoints)
      plot(sigpoints,squeeze(mean(mean(stattl_mc_TpAMSpN_bot{ll,12,1}.stat(chpl{1},fb{ff},sigpoints),1),2)),'b*');
    end
    sigpoints=find(squeeze(min(min(stattl_mc_TpAMSpN_top{ll,12,1}.prob(:,fb{ff},:),[],1),[],2))<.05);
    if ~isempty(sigpoints)
      plot(sigpoints,squeeze(mean(mean(stattl_mc_TpAMSpN_top{ll,12,1}.stat(chpl{1},fb{ff},sigpoints),1),2)),'k*');
    end
    sigpoints=find(squeeze(min(min(stattl_mc_TpAMSpN_TvB{ll,12,1}.prob(:,fb{ff},:),[],1),[],2))<.05);
    if ~isempty(sigpoints)
      plot(sigpoints,squeeze(mean(mean(stattl_mc_TpAMSpN_TvB{ll,12,1}.stat(chpl{1},fb{ff},sigpoints),1),2)),'r*');
    end
    ylim([-2 2])
    legend({'Bottom' 'Top' 'TvB'})
    title(fname{ff})
    print(ll+10*ff+30,[fdir 'statl_tfp_' fname{ff} '_' num2str(ll) '.eps'],'-painters','-depsc');
  end
end

ff=4; % gamma
for ll=soalist
    figure(ll+10*ff+60);plot(squeeze(mean(mean(statth_mcplv_TpAMSpN_bot_itcdepT{ll,12,1}.stat(chpl{1},fb{ff},:),1),2)),'b');
    hold on;
    plot(squeeze(mean(mean(statth_mcplv_TpAMSpN_top_itcdepT{ll,12,1}.stat(chpl{1},fb{ff},:),1),2)),'k');
    plot(squeeze(mean(mean(statth_mcplv_TpAMSpN_TvB_itcdepT{ll,12,1}.stat(chpl{1},fb{ff},:),1),2)),'r');
    sigpoints=find(squeeze(min(min(statth_mcplv_TpAMSpN_bot_itcdepT{ll,12,1}.prob(:,fb{ff},:),[],1),[],2))<.05);
    if ~isempty(sigpoints)
      plot(sigpoints,squeeze(mean(mean(statth_mcplv_TpAMSpN_bot_itcdepT{ll,12,1}.stat(chpl{1},fb{ff},sigpoints),1),2)),'b*');
    end
    sigpoints=find(squeeze(min(min(statth_mcplv_TpAMSpN_top_itcdepT{ll,12,1}.prob(:,fb{ff},:),[],1),[],2))<.05);
    if ~isempty(sigpoints)
      plot(sigpoints,squeeze(mean(mean(statth_mcplv_TpAMSpN_top_itcdepT{ll,12,1}.stat(chpl{1},fb{ff},sigpoints),1),2)),'k*');
    end
    sigpoints=find(squeeze(min(min(statth_mcplv_TpAMSpN_TvB_itcdepT{ll,12,1}.prob(:,fb{ff},:),[],1),[],2))<.05);
    if ~isempty(sigpoints)
      plot(sigpoints,squeeze(mean(mean(statth_mcplv_TpAMSpN_TvB_itcdepT{ll,12,1}.stat(chpl{1},fb{ff},sigpoints),1),2)),'r*');
    end
    ylim([-1.5 1.5])
    legend({'Bottom' 'Top' 'TvB'})
    title(fname{ff})
    print(ll+10*ff+60,[fdir 'stath_plvitc_' fname{ff} '_' num2str(ll) '.eps'],'-painters','-depsc');
end
for ll=soalist
    figure(ll+10*ff+70);plot(squeeze(mean(mean(statth_mc_TpAMSpN_bot{ll,12,1}.stat(chpl{1},fb{ff},:),1),2)),'b');
    hold on;
    plot(squeeze(mean(mean(statth_mc_TpAMSpN_top{ll,12,1}.stat(chpl{1},fb{ff},:),1),2)),'k');
    plot(squeeze(mean(mean(statth_mc_TpAMSpN_TvB{ll,12,1}.stat(chpl{1},fb{ff},:),1),2)),'r');
    sigpoints=find(squeeze(min(min(statth_mc_TpAMSpN_bot{ll,12,1}.prob(:,fb{ff},:),[],1),[],2))<.05);
    if ~isempty(sigpoints)
      plot(sigpoints,squeeze(mean(mean(statth_mc_TpAMSpN_bot{ll,12,1}.stat(chpl{1},fb{ff},sigpoints),1),2)),'b*');
    end
    sigpoints=find(squeeze(min(min(statth_mc_TpAMSpN_top{ll,12,1}.prob(:,fb{ff},:),[],1),[],2))<.05);
    if ~isempty(sigpoints)
      plot(sigpoints,squeeze(mean(mean(statth_mc_TpAMSpN_top{ll,12,1}.stat(chpl{1},fb{ff},sigpoints),1),2)),'k*');
    end
    sigpoints=find(squeeze(min(min(statth_mc_TpAMSpN_TvB{ll,12,1}.prob(:,fb{ff},:),[],1),[],2))<.05);
    if ~isempty(sigpoints)
      plot(sigpoints,squeeze(mean(mean(statth_mc_TpAMSpN_TvB{ll,12,1}.stat(chpl{1},fb{ff},sigpoints),1),2)),'r*');
    end
    ylim([-2 2])
    legend({'Bottom' 'Top' 'TvB'})
    title(fname{ff})
    print(ll+10*ff+70,[fdir 'stath_tfp_' fname{ff} '_' num2str(ll) '.eps'],'-painters','-depsc');
end


