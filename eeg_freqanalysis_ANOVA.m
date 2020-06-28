eeg_legomagic_preamble

%% PCA/ICA and ANOVA 

grsaveloadflag=1;  % =0 compute and save grlo*;   =1 load previously computed grlo*
statflag=1;  % =0 compute and save freqstat;   =1 load previously computed freqstat
sleep=0;
if sleep
  iter=11;
  usetr=1;
else
  iter=31;
  usetr=3;
end
soalist=[1 3 4 5 6 7 9];

if grsaveloadflag==0
  llind=1;
  for ll=soalist
    load(['grindTFR_cond' num2str(ll) '_sleep' num2str(sleep) '_iter' num2str(iter) '_trialkc-1_usetr' num2str(usetr) '_mcseed13_itc.mat']);
    clear *comb2    

    cfg=[];
    cfg.operation='subtract';
    cfg.parameter={'powspctrm' 'plvabs'};
    grindlo_TPA_MSPN{llind}=ft_math(cfg,grindlo_tacPaud_comb1,grindlo_tacMSpN_comb1);
    grindhi_TPA_MSPN{llind}=ft_math(cfg,grindhi_tacPaud_comb1,grindhi_tacMSpN_comb1);
    if sleep
      grindlo_TPA_MSPN_nKD{llind}=ft_math(cfg,grindlo_tacPaud_nKD,grindlo_tacMSpN_nKD);
      grindlo_TPA_MSPN_nSD{llind}=ft_math(cfg,grindlo_tacPaud_nSD,grindlo_tacMSpN_nSD);
      grindlo_TPA_MSPN_nKD_nSD{llind}=ft_math(cfg,grindlo_tacPaud_nKD_nSD,grindlo_tacMSpN_nKD_nSD);
      grindhi_TPA_MSPN_nKD{llind}=ft_math(cfg,grindhi_tacPaud_nKD,grindhi_tacMSpN_nKD);
      grindhi_TPA_MSPN_nSD{llind}=ft_math(cfg,grindhi_tacPaud_nSD,grindhi_tacMSpN_nSD);
      grindhi_TPA_MSPN_nKD_nSD{llind}=ft_math(cfg,grindhi_tacPaud_nKD_nSD,grindhi_tacMSpN_nKD_nSD);
    end
    if ll>5
      grindlo_TPA_MSPN{llind}.time=grindlo_TPA_MSPN{llind}.time-soades(ll);
      grindhi_TPA_MSPN{llind}.time=grindhi_TPA_MSPN{llind}.time-soades(ll);
      if sleep
        grindlo_TPA_MSPN_nKD{llind}.time=grindlo_TPA_MSPN_nKD{llind}.time-soades(ll);
        grindlo_TPA_MSPN_nSD{llind}.time=grindlo_TPA_MSPN_nSD{llind}.time-soades(ll);
        grindlo_TPA_MSPN_nKD_nSD{llind}.time=grindlo_TPA_MSPN_nKD_nSD{llind}.time-soades(ll);
        grindhi_TPA_MSPN_nKD{llind}.time=grindhi_TPA_MSPN_nKD{llind}.time-soades(ll);
        grindhi_TPA_MSPN_nSD{llind}.time=grindhi_TPA_MSPN_nSD{llind}.time-soades(ll);
        grindhi_TPA_MSPN_nKD_nSD{llind}.time=grindhi_TPA_MSPN_nKD_nSD{llind}.time-soades(ll);
      end
    end
    cfg=[];
    cfg.latency=[0 1.2];  % okay for both sleep and awake
    grindlo_TPA_MSPN_sel{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN{llind});
    grindhi_TPA_MSPN_sel{llind}=ft_selectdata(cfg,grindhi_TPA_MSPN{llind});
    if sleep
      grindlo_TPA_MSPN_nKD_sel{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_nKD{llind});
      grindlo_TPA_MSPN_nSD_sel{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_nSD{llind});
      grindlo_TPA_MSPN_nKD_nSD_sel{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_nKD_nSD{llind});
      grindhi_TPA_MSPN_nKD_sel{llind}=ft_selectdata(cfg,grindhi_TPA_MSPN_nKD{llind});
      grindhi_TPA_MSPN_nSD_sel{llind}=ft_selectdata(cfg,grindhi_TPA_MSPN_nSD{llind});
      grindhi_TPA_MSPN_nKD_nSD_sel{llind}=ft_selectdata(cfg,grindhi_TPA_MSPN_nKD_nSD{llind});
    end
    llind=llind+1;
  end
  cfg=[];
  cfg.avgoverrpt='yes';
  cfg.latency=[0 1.05];
  for llind=1:7
    gravelo_TPA_MSPN_avg{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_sel{llind});
    gravehi_TPA_MSPN_avg{llind}=ft_selectdata(cfg,grindhi_TPA_MSPN_sel{llind});
    if sleep
      gravelo_TPA_MSPN_nKD_avg{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_nKD_sel{llind});
      gravelo_TPA_MSPN_nSD_avg{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_nSD_sel{llind});
      gravelo_TPA_MSPN_nKD_nSD_avg{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_nKD_nSD_sel{llind});
      gravehi_TPA_MSPN_nKD_avg{llind}=ft_selectdata(cfg,grindhi_TPA_MSPN_nKD_sel{llind});
      gravehi_TPA_MSPN_nSD_avg{llind}=ft_selectdata(cfg,grindhi_TPA_MSPN_nSD_sel{llind});
      gravehi_TPA_MSPN_nKD_nSD_avg{llind}=ft_selectdata(cfg,grindhi_TPA_MSPN_nKD_nSD_sel{llind});
    end
  end
  save(['grlo_TPA_MSPN_sleep' num2str(sleep) '.mat'],'gravelo_TPA_MSPN_*avg','grindlo_TPA_MSPN_*sel')
  save(['grhi_TPA_MSPN_sleep' num2str(sleep) '.mat'],'gravehi_TPA_MSPN_*avg','grindhi_TPA_MSPN_*sel')
elseif grsaveloadflag==1
  load(['grlo_TPA_MSPN_sleep' num2str(sleep) '.mat']);
  if sleep
    load(['grhi_TPA_MSPN_sleep' num2str(sleep) '.mat']);  % note this not calculated yet for sleep=0;
  end
end


% Average over band / band-specific
cfg=[];
cfg.avgoverfreq='yes';
for llind=1:7
  cfg.frequency=[4 6.5];
  gravelo_TPA_MSPN_avg_theta{llind}=ft_selectdata(cfg,gravelo_TPA_MSPN_avg{llind});
  grindlo_TPA_MSPN_sel_theta{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_sel{llind});
  if sleep
    gravelo_TPA_MSPN_nKD_avg_theta{llind}=ft_selectdata(cfg,gravelo_TPA_MSPN_nKD_avg{llind});
    grindlo_TPA_MSPN_nKD_sel_theta{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_nKD_sel{llind});
    gravelo_TPA_MSPN_nSD_avg_theta{llind}=ft_selectdata(cfg,gravelo_TPA_MSPN_nSD_avg{llind});
    grindlo_TPA_MSPN_nSD_sel_theta{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_nSD_sel{llind});
    gravelo_TPA_MSPN_nKD_nSD_avg_theta{llind}=ft_selectdata(cfg,gravelo_TPA_MSPN_nKD_nSD_avg{llind});
    grindlo_TPA_MSPN_nKD_nSD_sel_theta{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_nKD_nSD_sel{llind});
  end
  cfg.frequency=[8 12];
  gravelo_TPA_MSPN_avg_alpha{llind}=ft_selectdata(cfg,gravelo_TPA_MSPN_avg{llind});
  grindlo_TPA_MSPN_sel_alpha{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_sel{llind});
  if sleep
    gravelo_TPA_MSPN_nKD_avg_alpha{llind}=ft_selectdata(cfg,gravelo_TPA_MSPN_nKD_avg{llind});
    grindlo_TPA_MSPN_nKD_sel_alpha{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_nKD_sel{llind});
    gravelo_TPA_MSPN_nSD_avg_alpha{llind}=ft_selectdata(cfg,gravelo_TPA_MSPN_nSD_avg{llind});
    grindlo_TPA_MSPN_nSD_sel_alpha{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_nSD_sel{llind});
    gravelo_TPA_MSPN_nKD_nSD_avg_alpha{llind}=ft_selectdata(cfg,gravelo_TPA_MSPN_nKD_nSD_avg{llind});
    grindlo_TPA_MSPN_nKD_nSD_sel_alpha{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_nKD_nSD_sel{llind});
  end
  cfg.frequency=[14 30];
  gravelo_TPA_MSPN_avg_beta{llind}=ft_selectdata(cfg,gravelo_TPA_MSPN_avg{llind});
  grindlo_TPA_MSPN_sel_beta{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_sel{llind});
  if sleep
    gravelo_TPA_MSPN_nKD_avg_beta{llind}=ft_selectdata(cfg,gravelo_TPA_MSPN_nKD_avg{llind});
    grindlo_TPA_MSPN_nKD_sel_beta{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_nKD_sel{llind});
    gravelo_TPA_MSPN_nSD_avg_beta{llind}=ft_selectdata(cfg,gravelo_TPA_MSPN_nSD_avg{llind});
    grindlo_TPA_MSPN_nSD_sel_beta{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_nSD_sel{llind});
    gravelo_TPA_MSPN_nKD_nSD_avg_beta{llind}=ft_selectdata(cfg,gravelo_TPA_MSPN_nKD_nSD_avg{llind});
    grindlo_TPA_MSPN_nKD_nSD_sel_beta{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_nKD_nSD_sel{llind});
  end
  if sleep
    cfg.frequency=[40 80];
    gravehi_TPA_MSPN_avg_gamma{llind}=ft_selectdata(cfg,gravehi_TPA_MSPN_avg{llind});
    grindhi_TPA_MSPN_sel_gamma{llind}=ft_selectdata(cfg,grindhi_TPA_MSPN_sel{llind});
    if sleep
      gravehi_TPA_MSPN_nKD_avg_gamma{llind}=ft_selectdata(cfg,gravehi_TPA_MSPN_nKD_avg{llind});
      grindhi_TPA_MSPN_nKD_sel_gamma{llind}=ft_selectdata(cfg,grindhi_TPA_MSPN_nKD_sel{llind});
      gravehi_TPA_MSPN_nSD_avg_gamma{llind}=ft_selectdata(cfg,gravehi_TPA_MSPN_nSD_avg{llind});
      grindhi_TPA_MSPN_nSD_sel_gamma{llind}=ft_selectdata(cfg,grindhi_TPA_MSPN_nSD_sel{llind});
      gravehi_TPA_MSPN_nKD_nSD_avg_gamma{llind}=ft_selectdata(cfg,gravehi_TPA_MSPN_nKD_nSD_avg{llind});
      grindhi_TPA_MSPN_nKD_nSD_sel_gamma{llind}=ft_selectdata(cfg,grindhi_TPA_MSPN_nKD_nSD_sel{llind});
    end
  end
  cfg.frequency=[8 20];
  gravelo_TPA_MSPN_avg_alfbet{llind}=ft_selectdata(cfg,gravelo_TPA_MSPN_avg{llind});
  grindlo_TPA_MSPN_sel_alfbet{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_sel{llind});
  cfg.frequency=[14 20];
  gravelo_TPA_MSPN_avg_beta1{llind}=ft_selectdata(cfg,gravelo_TPA_MSPN_avg{llind});
  grindlo_TPA_MSPN_sel_beta1{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_sel{llind});
  if sleep  % spindle range
    gravelo_TPA_MSPN_nKD_avg_beta1{llind}=ft_selectdata(cfg,gravelo_TPA_MSPN_nKD_avg{llind});
    grindlo_TPA_MSPN_nKD_sel_beta1{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_nKD_sel{llind});
    gravelo_TPA_MSPN_nSD_avg_beta1{llind}=ft_selectdata(cfg,gravelo_TPA_MSPN_nSD_avg{llind});
    grindlo_TPA_MSPN_nSD_sel_beta1{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_nSD_sel{llind});
    gravelo_TPA_MSPN_nKD_nSD_avg_beta1{llind}=ft_selectdata(cfg,gravelo_TPA_MSPN_nKD_nSD_avg{llind});
    grindlo_TPA_MSPN_nKD_nSD_sel_beta1{llind}=ft_selectdata(cfg,grindlo_TPA_MSPN_nKD_nSD_sel{llind});
  end
end

nsub=size(grindlo_TPA_MSPN_sel{1}.powspctrm,1);
if sleep
  nsub_nKD=size(grindlo_TPA_MSPN_nKD_sel{1}.powspctrm,1);
  nsub_nSD=size(grindlo_TPA_MSPN_nSD_sel{1}.powspctrm,1);
  nsub_nKD_nSD=size(grindlo_TPA_MSPN_nKD_nSD_sel{1}.powspctrm,1);
end
load eeg1010_neighb

if statflag==0  % need to run stats the first time, after that just load them.
  % Do band-averaged only.
  cfg=[];
  cfg.latency=[0 1.05];
  cfg.frequency='all';
  cfg.method='montecarlo';
  cfg.parameter='powspctrm';
  cfg.neighbours=neighbours;
  cfg.correctm='cluster';
  cfg.numrandomization=2000;
  cfg.tail=1; % only 1 makes sense for F-test
  cfg.ivar=1;
  cfg.uvar=2;
  cfg.design=set_cfg_design_depF(nsub);
  cfg.randomseed=13;
  cfg.clusterstatistic='maxsum';
  cfg.statistic='depsamplesFunivariate';
  cfg.computecritval='yes';
  cfg.comptueprob='yes';
  freqstat_TPA_MSPN_1wayANOVA_lo = ft_freqstatistics(cfg,grindlo_TPA_MSPN_sel{:}); % 0.0305   0.0995  *
  freqstat_TPA_MSPN_1wayANOVA_hi = ft_freqstatistics(cfg,grindhi_TPA_MSPN_sel{:}); % 
  freqstat_TPA_MSPN_1wayANOVA_theta = ft_freqstatistics(cfg,grindlo_TPA_MSPN_sel_theta{:}); %  0.0235  0.0995   *
  freqstat_TPA_MSPN_1wayANOVA_alpha = ft_freqstatistics(cfg,grindlo_TPA_MSPN_sel_alpha{:}); % 0.2574
  freqstat_TPA_MSPN_1wayANOVA_beta = ft_freqstatistics(cfg,grindlo_TPA_MSPN_sel_beta{:}); % 0.001 [but 0.1504 if 12-30 Hz]
  freqstat_TPA_MSPN_1wayANOVA_gamma = ft_freqstatistics(cfg,grindhi_TPA_MSPN_sel_gamma{:}); 
  freqstat_TPA_MSPN_1wayANOVA_alfbet = ft_freqstatistics(cfg,grindlo_TPA_MSPN_sel_alfbet{:}); % 0.0790
  freqstat_TPA_MSPN_1wayANOVA_beta1 = ft_freqstatistics(cfg,grindlo_TPA_MSPN_sel_beta1{:}); % 0.0865
  % cfg.frequency=[4 7];
  % freqstat_TPA_MSPN_1wayANOVA_thetaall = ft_freqstatistics(cfg,grindlo_TPA_MSPN_sel{:}); % 0.0015   *
  % cfg.frequency=[8 12];
  % freqstat_TPA_MSPN_1wayANOVA_alphaall = ft_freqstatistics(cfg,grindlo_TPA_MSPN_sel{:}); % 0.6222
  % cfg.frequency=[8 20];
  % freqstat_TPA_MSPN_1wayANOVA_alfbetall = ft_freqstatistics(cfg,grindlo_TPA_MSPN_sel{:}); % 0.0785
  % cfg.frequency=[14 30];
  % freqstat_TPA_MSPN_1wayANOVA_betaall = ft_freqstatistics(cfg,grindlo_TPA_MSPN_sel{:}); % 0.0010   *
  if sleep
    cfg.design=set_cfg_design_depF(nsub_nKD);
    freqstat_TPA_MSPN_1wayANOVA_nKD_lo = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nKD_sel{:}); 
    freqstat_TPA_MSPN_1wayANOVA_nKD_hi = ft_freqstatistics(cfg,grindhi_TPA_MSPN_nKD_sel{:}); 
    freqstat_TPA_MSPN_1wayANOVA_nKD_theta = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nKD_sel_theta{:}); 
    freqstat_TPA_MSPN_1wayANOVA_nKD_alpha = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nKD_sel_alpha{:}); 
    freqstat_TPA_MSPN_1wayANOVA_nKD_beta = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nKD_sel_beta{:}); 
    freqstat_TPA_MSPN_1wayANOVA_nKD_beta1 = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nKD_sel_beta1{:}); 
    freqstat_TPA_MSPN_1wayANOVA_nKD_gamma = ft_freqstatistics(cfg,grindhi_TPA_MSPN_nKD_sel_gamma{:}); 
%     cfg.design=set_cfg_design_depF(nsub_nSD); % commented out b/c ll=3 has only nsub_nSD=18 (rest are =19).
%     freqstat_TPA_MSPN_1wayANOVA_nSD_theta = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nSD_sel_theta{:}); 
%     freqstat_TPA_MSPN_1wayANOVA_nSD_alpha = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nSD_sel_theta{:}{:}); 
%     freqstat_TPA_MSPN_1wayANOVA_nSD_beta = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nSD_sel_beta{:}); 
    cfg.design=set_cfg_design_depF(nsub_nKD_nSD);
    freqstat_TPA_MSPN_1wayANOVA_nKD_nSD_lo = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nKD_nSD_sel{:}); 
    freqstat_TPA_MSPN_1wayANOVA_nKD_nSD_hi = ft_freqstatistics(cfg,grindhi_TPA_MSPN_nKD_nSD_sel{:}); 
    freqstat_TPA_MSPN_1wayANOVA_nKD_nSD_theta = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nKD_nSD_sel_theta{:}); 
    freqstat_TPA_MSPN_1wayANOVA_nKD_nSD_alpha = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nKD_nSD_sel_alpha{:}); 
    freqstat_TPA_MSPN_1wayANOVA_nKD_nSD_beta = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nKD_nSD_sel_beta{:}); 
    freqstat_TPA_MSPN_1wayANOVA_nKD_nSD_beta1 = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nKD_nSD_sel_beta1{:}); 
    freqstat_TPA_MSPN_1wayANOVA_nKD_nSD_gamma = ft_freqstatistics(cfg,grindhi_TPA_MSPN_nKD_nSD_sel_gamma{:}); 
  end
  
  cfg.design=set_cfg_design_depF(nsub);
  cfg.parameter='plvabs';
  cfg.frequency='all';
  plvstat_TPA_MSPN_1wayANOVA_lo = ft_freqstatistics(cfg,grindlo_TPA_MSPN_sel{:}); 
  plvstat_TPA_MSPN_1wayANOVA_hi = ft_freqstatistics(cfg,grindhi_TPA_MSPN_sel{:}); 
  plvstat_TPA_MSPN_1wayANOVA_theta = ft_freqstatistics(cfg,grindlo_TPA_MSPN_sel_theta{:}); % 0.0005**
  plvstat_TPA_MSPN_1wayANOVA_alpha = ft_freqstatistics(cfg,grindlo_TPA_MSPN_sel_alpha{:}); % 0.3318
  plvstat_TPA_MSPN_1wayANOVA_beta = ft_freqstatistics(cfg,grindlo_TPA_MSPN_sel_beta{:}); % 0.3433
  plvstat_TPA_MSPN_1wayANOVA_gamma = ft_freqstatistics(cfg,grindhi_TPA_MSPN_sel_gamma{:}); % 0.3433
  plvstat_TPA_MSPN_1wayANOVA_alfbet = ft_freqstatistics(cfg,grindlo_TPA_MSPN_sel_alfbet{:}); % 0.1624
  plvstat_TPA_MSPN_1wayANOVA_beta1 = ft_freqstatistics(cfg,grindlo_TPA_MSPN_sel_beta1{:}); % 0.1274
  % cfg.frequency=[4 7];
  % plvstat_TPA_MSPN_1wayANOVA_thetaall = ft_freqstatistics(cfg,grindlo_TPA_MSPN_sel{:}); % 0.0015   *
  % cfg.frequency=[8 12];
  % plvstat_TPA_MSPN_1wayANOVA_alphaall = ft_freqstatistics(cfg,grindlo_TPA_MSPN_sel{:}); % 0.6222
  % cfg.frequency=[8 20];
  % plvstat_TPA_MSPN_1wayANOVA_alfbetall = ft_freqstatistics(cfg,grindlo_TPA_MSPN_sel{:}); % 0.0785
  % cfg.frequency=[14 30];
  % plvstat_TPA_MSPN_1wayANOVA_betaall = ft_freqstatistics(cfg,grindlo_TPA_MSPN_sel{:}); % 0.0010   *
  if sleep
    cfg.design=set_cfg_design_depF(nsub_nKD);
    plvstat_TPA_MSPN_1wayANOVA_nKD_lo = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nKD_sel{:});
    plvstat_TPA_MSPN_1wayANOVA_nKD_hi = ft_freqstatistics(cfg,grindhi_TPA_MSPN_nKD_sel{:});
    plvstat_TPA_MSPN_1wayANOVA_nKD_theta = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nKD_sel_theta{:}); % 0.0005**
    plvstat_TPA_MSPN_1wayANOVA_nKD_alpha = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nKD_sel_alpha{:}); % 0.3318
    plvstat_TPA_MSPN_1wayANOVA_nKD_beta = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nKD_sel_beta{:}); % 0.3433
    plvstat_TPA_MSPN_1wayANOVA_nKD_beta1 = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nKD_sel_beta1{:}); 
    plvstat_TPA_MSPN_1wayANOVA_nKD_gamma = ft_freqstatistics(cfg,grindhi_TPA_MSPN_nKD_sel_gamma{:}); 
%     cfg.design=set_cfg_design_depF(nsub_nSD);
%     plvstat_TPA_MSPN_1wayANOVA_nSD_theta = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nSD_sel_theta{:}); % 0.0005**
%     plvstat_TPA_MSPN_1wayANOVA_nSD_alpha = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nSD_sel_alpha{:}); % 0.3318
%     plvstat_TPA_MSPN_1wayANOVA_nSD_beta = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nSD_sel_beta{:}); % 0.3433
    cfg.design=set_cfg_design_depF(nsub_nKD_nSD);
    plvstat_TPA_MSPN_1wayANOVA_nKD_nSD_lo = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nKD_nSD_sel{:});
    plvstat_TPA_MSPN_1wayANOVA_nKD_nSD_hi = ft_freqstatistics(cfg,grindhi_TPA_MSPN_nKD_nSD_sel{:});
    plvstat_TPA_MSPN_1wayANOVA_nKD_nSD_theta = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nKD_nSD_sel_theta{:}); % 0.0005**
    plvstat_TPA_MSPN_1wayANOVA_nKD_nSD_alpha = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nKD_nSD_sel_alpha{:}); % 0.3318
    plvstat_TPA_MSPN_1wayANOVA_nKD_nSD_beta = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nKD_nSD_sel_beta{:}); % 0.3433
    plvstat_TPA_MSPN_1wayANOVA_nKD_nSD_beta1 = ft_freqstatistics(cfg,grindlo_TPA_MSPN_nKD_nSD_sel_beta1{:}); 
    plvstat_TPA_MSPN_1wayANOVA_nKD_nSD_gamma = ft_freqstatistics(cfg,grindhi_TPA_MSPN_nKD_nSD_sel_gamma{:}); 
  end
  
  save(['freqstat_sleep' num2str(sleep) '.mat'],'freqstat*','plvstat*')
elseif statflag==1
  load(['freqstat_sleep' num2str(sleep) '.mat'])
end

if sleep==1
  any(unique(freqstat_TPA_MSPN_1wayANOVA_lo.prob)<.05)  %*
  any(unique(freqstat_TPA_MSPN_1wayANOVA_hi.prob)<.05) 
  any(unique(freqstat_TPA_MSPN_1wayANOVA_alfbet.prob)<.05)  %*
  any(unique(freqstat_TPA_MSPN_1wayANOVA_alpha.prob)<.05)  %*
  any(unique(freqstat_TPA_MSPN_1wayANOVA_beta.prob)<.05)  %*
  any(unique(freqstat_TPA_MSPN_1wayANOVA_beta1.prob)<.05)  %*
  any(unique(freqstat_TPA_MSPN_1wayANOVA_theta.prob)<.05)  %*
  any(unique(freqstat_TPA_MSPN_1wayANOVA_gamma.prob)<.05) 
  any(unique(plvstat_TPA_MSPN_1wayANOVA_lo.prob)<.05)
  any(unique(plvstat_TPA_MSPN_1wayANOVA_hi.prob)<.05)
  any(unique(plvstat_TPA_MSPN_1wayANOVA_alfbet.prob)<.05)
  any(unique(plvstat_TPA_MSPN_1wayANOVA_alpha.prob)<.05)
  any(unique(plvstat_TPA_MSPN_1wayANOVA_beta.prob)<.05)
  any(unique(plvstat_TPA_MSPN_1wayANOVA_beta1.prob)<.05)  %*
  any(unique(plvstat_TPA_MSPN_1wayANOVA_theta.prob)<.05)
  any(unique(plvstat_TPA_MSPN_1wayANOVA_gamma.prob)<.05)  %*
% %   any(unique(freqstat_TPA_MSPN_1wayANOVA_nKD_alpha.prob)<.05)
% %   any(unique(freqstat_TPA_MSPN_1wayANOVA_nKD_beta.prob)<.05)
% %   any(unique(freqstat_TPA_MSPN_1wayANOVA_nKD_theta.prob)<.05)
end

% % figure;imagesc(squeeze(sum(freqstat_TPA_MSPN_1wayANOVA.mask,1)));axis xy
% figure;imagesc(squeeze(freqstat_TPA_MSPN_1wayANOVA_theta.mask))
% figure;imagesc(squeeze(freqstat_TPA_MSPN_1wayANOVA_beta.mask))
% % figure;imagesc(squeeze(sum(freqstat_TPA_MSPN_1wayANOVA_thetaall.mask,2)))
% % figure;imagesc(squeeze(sum(plvstat_TPA_MSPN_1wayANOVA.mask,1)));axis xy
% figure;imagesc(squeeze(plvstat_TPA_MSPN_1wayANOVA_theta.mask))
% % figure;imagesc(squeeze(sum(plvstat_TPA_MSPN_1wayANOVA_thetaall.mask,2)))

%%
if sleep==0
  % re-order channels to be sensible front to back (final for paper)
  induse=[];
  cfg=[];cfg.channel='A*';tmp=ft_selectdata(cfg,freqstat_TPA_MSPN_1wayANOVA_theta);induse=[induse; match_str(freqstat_TPA_MSPN_1wayANOVA_theta.label,tmp.label)];
  cfg=[];cfg.channel='F*';tmp=ft_selectdata(cfg,freqstat_TPA_MSPN_1wayANOVA_theta);induse=[induse; match_str(freqstat_TPA_MSPN_1wayANOVA_theta.label,tmp.label)]
  cfg=[];cfg.channel='T*';tmp=ft_selectdata(cfg,freqstat_TPA_MSPN_1wayANOVA_theta);induse=[induse; match_str(freqstat_TPA_MSPN_1wayANOVA_theta.label,tmp.label)]
  cfg=[];cfg.channel='C*';tmp=ft_selectdata(cfg,freqstat_TPA_MSPN_1wayANOVA_theta);induse=[induse; match_str(freqstat_TPA_MSPN_1wayANOVA_theta.label,tmp.label)]
  cfg=[];cfg.channel='P*';tmp=ft_selectdata(cfg,freqstat_TPA_MSPN_1wayANOVA_theta);induse=[induse; match_str(freqstat_TPA_MSPN_1wayANOVA_theta.label,tmp.label)]
  cfg=[];cfg.channel='O*';tmp=ft_selectdata(cfg,freqstat_TPA_MSPN_1wayANOVA_theta);induse=[induse; match_str(freqstat_TPA_MSPN_1wayANOVA_theta.label,tmp.label)]
  % Figure for paper:
  figure(43);imagesc(freqstat_TPA_MSPN_1wayANOVA_theta.time,1:63,squeeze(freqstat_TPA_MSPN_1wayANOVA_theta.mask(induse,:)));
  xlim([-.01 1.051]);
  set(gca,'FontSize',30);
  set(gca,'XTick',[-.6:.1:1.8])
  set(gca,'XTickLabel',{ '' '-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'  '' ' ' ''  ' ' '1.5' '' ' ' '' })
  set(gca,'YTick',4:11:63)
  set(gca,'YTickLabel',freqstat_TPA_MSPN_1wayANOVA_theta.label(induse(4:11:end)))
%   set(gca,'YTickMode','auto')
  print(43,[fdir 'PowThetastatANOVA_sleep' num2str(sleep) '.eps'],'-painters','-depsc')
  figure(44);imagesc(freqstat_TPA_MSPN_1wayANOVA_theta.time,1:63,squeeze(freqstat_TPA_MSPN_1wayANOVA_beta.mask(induse,:)));
  xlim([-.01 1.051]);
  set(gca,'FontSize',30);
  set(gca,'XTick',[-.6:.1:1.8])
  set(gca,'XTickLabel',{ '' '-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'  '' ' ' ''  ' ' '1.5' '' ' ' '' })
  set(gca,'YTick',4:11:63)
  set(gca,'YTickLabel',freqstat_TPA_MSPN_1wayANOVA_theta.label(induse(4:11:end)))
%   set(gca,'YTickMode','auto')
  print(44,[fdir 'PowBetastatANOVA_sleep' num2str(sleep) '.eps'],'-painters','-depsc')
  figure(45);imagesc(freqstat_TPA_MSPN_1wayANOVA_theta.time,1:63,squeeze(plvstat_TPA_MSPN_1wayANOVA_theta.mask(induse,:)));
  xlim([-.01 1.051]);
  set(gca,'FontSize',30);
  set(gca,'XTick',[-.6:.1:1.8])
  set(gca,'XTickLabel',{ '' '-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'  '' ' ' ''  ' ' '1.5' '' ' ' '' })
  set(gca,'YTick',4:11:63)
  set(gca,'YTickLabel',freqstat_TPA_MSPN_1wayANOVA_theta.label(induse(4:11:end)))
%   set(gca,'YTickMode','auto')
  print(45,[fdir 'PlvThetastatANOVA_sleep' num2str(sleep) '.eps'],'-painters','-depsc')
elseif sleep==1
  % Figure for paper:
  figure(43);imagesc(freqstat_TPA_MSPN_1wayANOVA_theta.time,1:63,squeeze(freqstat_TPA_MSPN_1wayANOVA_theta.mask));
  xlim([-.01 1.051]);
  set(gca,'FontSize',30);
  set(gca,'XTick',[-.6:.1:1.8])
  set(gca,'XTickLabel',{ '' '-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'  '' ' ' ''  ' ' '1.5' '' ' ' '' })
  set(gca,'YTick',1:63)
  set(gca,'YTickLabel',freqstat_TPA_MSPN_1wayANOVA_theta.label)
  set(gca,'YTickMode','auto')
  print(43,[fdir 'PowThetastatANOVA_sleep' num2str(sleep) '.eps'],'-painters','-depsc')
  figure(44);imagesc(freqstat_TPA_MSPN_1wayANOVA_theta.time,1:63,squeeze(freqstat_TPA_MSPN_1wayANOVA_beta.mask));
  xlim([-.01 1.051]);
  set(gca,'FontSize',30);
  set(gca,'XTick',[-.6:.1:1.8])
  set(gca,'XTickLabel',{ '' '-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'  '' ' ' ''  ' ' '1.5' '' ' ' '' })
  set(gca,'YTick',1:63)
  set(gca,'YTickLabel',freqstat_TPA_MSPN_1wayANOVA_theta.label)
  set(gca,'YTickMode','auto')
  print(44,[fdir 'PowBetastatANOVA_sleep' num2str(sleep) '.eps'],'-painters','-depsc')
  figure(45);imagesc(freqstat_TPA_MSPN_1wayANOVA_theta.time,1:63,squeeze(plvstat_TPA_MSPN_1wayANOVA_theta.mask));
  xlim([-.01 1.051]);
  set(gca,'FontSize',30);
  set(gca,'XTick',[-.6:.1:1.8])
  set(gca,'XTickLabel',{ '' '-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'  '' ' ' ''  ' ' '1.5' '' ' ' '' })
  set(gca,'YTick',1:63)
  set(gca,'YTickLabel',freqstat_TPA_MSPN_1wayANOVA_theta.label)
  set(gca,'YTickMode','auto')
  print(45,[fdir 'PlvThetastatANOVA_sleep' num2str(sleep) '.eps'],'-painters','-depsc')
end

%%
% time window of each to use:
reltimepoints_thetaTFP=find(mean(squeeze(freqstat_TPA_MSPN_1wayANOVA_theta.mask))>.1); % note this one requires smaller threshold 
diffpoints_thetaTFP=find(diff(reltimepoints_thetaTFP)>1)
thetaTFP_anova_twall=[reltimepoints_thetaTFP(1):reltimepoints_thetaTFP(end)];

reltimepoints_beta=find(mean(squeeze(freqstat_TPA_MSPN_1wayANOVA_beta.mask))>.2);
diffpoints_beta=find(diff(reltimepoints_beta)>1)
betaTFP_anova_twall=[reltimepoints_beta(1):reltimepoints_beta(end)];
betaTFP_anova_tw1=[reltimepoints_beta(1):reltimepoints_beta(6)];
betaTFP_anova_tw2=[reltimepoints_beta(7):reltimepoints_beta(end)];

reltimepoints_thetaITC=find(mean(squeeze(plvstat_TPA_MSPN_1wayANOVA_theta.mask))>.2);
diffpoints_thetaITC=find(diff(reltimepoints_thetaITC)>1)
thetaITC_anova_twall=[reltimepoints_thetaITC(1):reltimepoints_thetaITC(end)];



% Plotting topo of stat for figures
freqstat_TPA_MSPN_1wayANOVA_theta.stat1=freqstat_TPA_MSPN_1wayANOVA_theta.stat.*freqstat_TPA_MSPN_1wayANOVA_theta.mask;
freqstat_TPA_MSPN_1wayANOVA_beta.stat1=freqstat_TPA_MSPN_1wayANOVA_beta.stat.*freqstat_TPA_MSPN_1wayANOVA_beta.mask;
plvstat_TPA_MSPN_1wayANOVA_theta.stat1=plvstat_TPA_MSPN_1wayANOVA_theta.stat.*plvstat_TPA_MSPN_1wayANOVA_theta.mask;

figure;imagesc(squeeze(freqstat_TPA_MSPN_1wayANOVA_theta.stat.*freqstat_TPA_MSPN_1wayANOVA_theta.mask))
figure;imagesc(squeeze(freqstat_TPA_MSPN_1wayANOVA_beta.stat.*freqstat_TPA_MSPN_1wayANOVA_beta.mask))
figure;imagesc(squeeze(plvstat_TPA_MSPN_1wayANOVA_theta.stat.*plvstat_TPA_MSPN_1wayANOVA_theta.mask))

% since beta might have 2 sub-clusters
mask_use=freqstat_TPA_MSPN_1wayANOVA_beta.mask;
mask_use(:,23:end)=0;
freqstat_TPA_MSPN_1wayANOVA_beta.stat_tw1=freqstat_TPA_MSPN_1wayANOVA_beta.stat.*mask_use;

mask_use=freqstat_TPA_MSPN_1wayANOVA_beta.mask;
mask_use(:,1:22)=0;
freqstat_TPA_MSPN_1wayANOVA_beta.stat_tw2=freqstat_TPA_MSPN_1wayANOVA_beta.stat.*mask_use;



cfg=[];
cfg.zlim='zeromax';
cfg.layout='eeg1010';
cfg.highlight          = 'on';
cfg.highlightsize = 16;
cfg.markersymbol = 'o';
cfg.parameter='stat';
cfg.colorbar='yes';

cfg.xlim=[.1 .37]; % reltimepoints_thetaTFP
cfg.highlightchannel   =  freqstat_TPA_MSPN_1wayANOVA_theta.label(find(nanmean(freqstat_TPA_MSPN_1wayANOVA_theta.stat1,3)));
figure(1);ft_topoplotTFR(cfg,freqstat_TPA_MSPN_1wayANOVA_theta)
print(1,[fdir 'ANOVA_PowTheta_topo1_allF_sleep' num2str(sleep) '.eps'],'-painters','-depsc')

cfg.xlim=[.13 .28]; % time
cfg.highlightchannel   =  freqstat_TPA_MSPN_1wayANOVA_beta.label(find(nanmean(freqstat_TPA_MSPN_1wayANOVA_beta.stat1,3)));
figure(2);ft_topoplotTFR(cfg,freqstat_TPA_MSPN_1wayANOVA_beta)
print(2,[fdir 'ANOVA_PowBeta_topoAll_allF_sleep' num2str(sleep) '.eps'],'-painters','-depsc')
cfg.xlim=[.13 .18]; % time
cfg.highlightchannel   =  freqstat_TPA_MSPN_1wayANOVA_beta.label(find(nanmean(freqstat_TPA_MSPN_1wayANOVA_beta.stat_tw1,3)));
figure(22);ft_topoplotTFR(cfg,freqstat_TPA_MSPN_1wayANOVA_beta)
print(22,[fdir 'ANOVA_PowBeta_topo1_allF_sleep' num2str(sleep) '.eps'],'-painters','-depsc')
cfg.xlim=[.25 .28]; % time
cfg.highlightchannel   =  freqstat_TPA_MSPN_1wayANOVA_beta.label(find(nanmean(freqstat_TPA_MSPN_1wayANOVA_beta.stat_tw2,3)));
figure(23);ft_topoplotTFR(cfg,freqstat_TPA_MSPN_1wayANOVA_beta)
print(23,[fdir 'ANOVA_PowBeta_topo2_allF_sleep' num2str(sleep) '.eps'],'-painters','-depsc')

cfg.xlim=[0 .44]; % time
cfg.highlightchannel   =  plvstat_TPA_MSPN_1wayANOVA_theta.label(find(nanmean(plvstat_TPA_MSPN_1wayANOVA_theta.stat1,3)));
figure(3);ft_topoplotTFR(cfg,plvstat_TPA_MSPN_1wayANOVA_theta)
print(3,[fdir 'ANOVA_PlvTheta_topo1_allF_sleep' num2str(sleep) '.eps'],'-painters','-depsc')


clear fullpow fullplv
for llind=1:7
%   fullpow(llind,:,:)=reshape(gravelo_TPA_MSPN_avg{llind}.powspctrm,[63 14*size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
%   fullplv(llind,:,:)=reshape(gravelo_TPA_MSPN_avg{llind}.plvabs,[63 14*size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
  thetapow(llind,:,:)=squeeze(gravelo_TPA_MSPN_avg_theta{llind}.powspctrm);
  thetaplv(llind,:,:)=squeeze(gravelo_TPA_MSPN_avg_theta{llind}.plvabs);
%   alphapow(llind,:,:)=squeeze(gravelo_TPA_MSPN_avg_alpha{llind}.powspctrm);
%   alphaplv(llind,:,:)=squeeze(gravelo_TPA_MSPN_avg_alpha{llind}.plvabs);
  betapow(llind,:,:)=squeeze(gravelo_TPA_MSPN_avg_beta{llind}.powspctrm);
%   betaplv(llind,:,:)=squeeze(gravelo_TPA_MSPN_avg_beta{llind}.plvabs);
end
% fullpow_reshape=reshape(fullpow,[7 63*14*size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
% fullplv_reshape=reshape(fullplv,[7 63*14*size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
thetapow_reshape=reshape(thetapow,[7 63*size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
thetaplv_reshape=reshape(thetaplv,[7 63*size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
% alphapow_reshape=reshape(alphapow,[7 63*size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
% alphaplv_reshape=reshape(alphaplv,[7 63*size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
betapow_reshape=reshape(betapow,[7 63*size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
% betaplv_reshape=reshape(betaplv,[7 63*size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);



% Masked-PCA 
mask_use=squeeze(freqstat_TPA_MSPN_1wayANOVA_theta.mask);
name=['PowTheta_sleep' num2str(sleep) ];
[aa1,aaa1_powtheta,vvv1_powtheta,bb1]=pca_masked(mask_use,thetapow_reshape,gravelo_TPA_MSPN_avg{llind}.label,1,gravelo_TPA_MSPN_avg{4}.time,[-.1 1.15],fdir,name,[1 1]);
% percent variance explained:
bb1(1,1)/sum(bb1(:))

mask_use=squeeze(freqstat_TPA_MSPN_1wayANOVA_beta.mask);
name=['PowBeta_all_sleep' num2str(sleep)];
[aa1,aaa1_powbeta,vvv1_powbeta,bb1]=pca_masked(mask_use,betapow_reshape,gravelo_TPA_MSPN_avg{llind}.label,1,gravelo_TPA_MSPN_avg{4}.time,[-.1 1.15],fdir,name,[1 0]);
bb1(1,1)/sum(bb1(:))

% mask_use=squeeze(freqstat_TPA_MSPN_1wayANOVA_beta.mask);
% mask_use(:,23:end)=0;
% name='PowBeta_tw1';
% [aa1,aaa1,vvv1]=pca_masked(mask_use,betapow_reshape,gravelo_TPA_MSPN_avg{llind}.label,1,gravelo_TPA_MSPN_avg{4}.time,[-.1 1.15],fdir,name);
% 
% mask_use=squeeze(freqstat_TPA_MSPN_1wayANOVA_beta.mask);
% mask_use(:,1:22)=0;
% name='PowBeta_tw2';
% [aa1,aaa1,vvv1]=pca_masked(mask_use,betapow_reshape,gravelo_TPA_MSPN_avg{llind}.label,1,gravelo_TPA_MSPN_avg{4}.time,[-.1 1.15],fdir,name);

mask_use=squeeze(plvstat_TPA_MSPN_1wayANOVA_theta.mask);
name=['PlvTheta_sleep' num2str(sleep)];
[aa1,aaa1_plvtheta,vvv1_plvtheta,bb1]=pca_masked(mask_use,thetaplv_reshape,gravelo_TPA_MSPN_avg{llind}.label,1,gravelo_TPA_MSPN_avg{4}.time,[-.1 1.15],fdir,name,[1 1]);
bb1(1,1)/sum(bb1(:))

%% Correlation of subfigure C to subfigure D (Figure 4)
load tfrsave.mat

aaa1_Freqlabel=gravelo_TPA_MSPN_avg{llind}.label;
tfrsave_freqlabel=aaa1_Freqlabel;

[topocorr_powtheta,topopval_powtheta]=corr(aaa1_powtheta(:,1),tfrsaveOT{7}.topo{1});
[topocorr_powbeta,topopval_powbeta]  =corr(aaa1_powbeta(:,1), tfrsaveOB{7}.topo{1});
[topocorr_plvtheta,topopval_plvtheta]=corr(aaa1_plvtheta(:,1),tfrsaveLT{7}.topo{1});

%  figure;plot([aaa1_powtheta(:,1),tfrsaveOT{7}.topo{1}])
 figure;plot([aaa1_plvtheta(:,1),tfrsaveLT{7}.topo{1}])

[tc_corr_powtheta,tc_pval_powtheta]=corr(vvv1_powtheta(:,1),tfrsaveOT{7}.courseFC(23:128)');
[tc_corr_powbeta,tc_pval_powbeta] =corr(vvv1_powbeta(:,1), tfrsaveOB{7}.courseFC(23:128)');
[tc_corr_plvtheta,tc_pval_plvtheta]=corr(vvv1_plvtheta(:,1),tfrsaveLT{7}.courseFC(23:128)');

% Now - across the rows within a figure (Figure 4)
load tmpd5_forcorr.mat; % see eeg_legomagic_erp_plot.m
load erp_2ndpca.mat; % see eeg_legomagic_erp_stats_anova.m

chanskip=match_str(aaa1_Freqlabel,setdiff(aaa1_Freqlabel,erppca_label));

chanskip=match_str(plvstat_TPA_MSPN_1wayANOVA_theta.label,setdiff(plvstat_TPA_MSPN_1wayANOVA_theta.label,erppca_label));

[topocorr_4Ci_4Cii,topopval_4Ci_4Cii]   = corr(aaa1_tw2(:,1),aaa1_plvtheta(setdiff(1:63,chanskip),1));
[topocorr_4Ci_4Ciii,topopval_4Ci_4Ciii]  = corr(aaa1_tw2(:,1),aaa1_powtheta(setdiff(1:63,chanskip),1));
[topocorr_4Cii_4Ciii,topopval_4Cii_4Ciii] = corr(aaa1_plvtheta(setdiff(1:63,chanskip),1),aaa1_powtheta(setdiff(1:63,chanskip),1));
[topocorr_4Cii_4Ciii,topopval_4Cii_4Ciii] = corr(aaa1_plvtheta(:,1),aaa1_powtheta(:,1));

figure;plot([aaa1_tw2(:,1),aaa1_plvtheta(setdiff(1:63,chanskip),1)])

chanskipD=match_str(tfrsave_freqlabel,setdiff(tfrsave_freqlabel,tmpd5_forcorr_tw2{3}.label));

[topocorr_4Di_4Dii,topopval_4Di_4Dii]     = corr(tmpd5_forcorr_tw2{3}.avg,tfrsaveLT{7}.topo{1}(setdiff(1:63,chanskipD)))
[topocorr_4Di_4Diii,topopval_4Di_4Diii]   = corr(tmpd5_forcorr_tw2{3}.avg,tfrsaveOT{7}.topo{1}(setdiff(1:63,chanskipD)))
[topocorr_4Dii_4Diii,topopval_4Dii_4Diii] = corr(tfrsaveLT{7}.topo{1}(setdiff(1:63,chanskipD)),tfrsaveOT{7}.topo{1}(setdiff(1:63,chanskipD)))

figure;plot([tmpd5_forcorr_tw2{3}.avg,10* tfrsaveLT{7}.topo{1}(setdiff(1:63,chanskipD))])
figure;plot(tmpd5_forcorr_tw2{3}.avg,10* tfrsaveLT{7}.topo{1}(setdiff(1:63,chanskipD)),'o')

% temporal correlation
[tc_corr_4Ci_4Cii,tc_pval_4Ci_4Cii]     = corr(vvv1_tw2(1:10:end,1),vvv1_plvtheta(1:51,1))
[tc_corr_4Ci_4Ciii,tc_pval_4Ci_4Ciii]   = corr(vvv1_tw2(1:10:end,1),vvv1_powtheta(1:51,1))
[tc_corr_4Cii_4Ciii,tc_pval_4Cii_4Ciii] = corr(vvv1_plvtheta(1:51,1),vvv1_powtheta(1:51,1))

figure;plot([-vvv1_tw2(1:10:end,1),vvv1_plvtheta(1:51,1)])

[tc_corr_4Di_4Dii,tc_pval_4Di_4Dii]     = corr(tmpd5_tc{3,1}.avg(51:10:551)',tfrsaveLT{7}.courseFC(23:73)')
[tc_corr_4Di_4Diii,tc_pval_4Di_4Diii]   = corr(tmpd5_tc{3,1}.avg(51:10:551)',tfrsaveOT{7}.courseFC(23:73)')
[tc_corr_4Dii_4Diii,tc_pval_4Dii_4Diii] = corr(tfrsaveLT{7}.courseFC(23:73)',tfrsaveOT{7}.courseFC(23:73)')


%%
% [aao,bbo,vvo]=svd(fullpow_reshape,'econ');
% [aal,bbl,vvl]=svd(fullplv_reshape,'econ');
[aaoT,bboT,vvoT]=svd(thetapow_reshape,'econ');
[aalT,bblT,vvlT]=svd(thetaplv_reshape,'econ');
% [aaoA,bboA,vvoA]=svd(alphapow_reshape,'econ');
% [aalA,bblA,vvlA]=svd(alphaplv_reshape,'econ');
[aaoB,bboB,vvoB]=svd(betapow_reshape,'econ');
% [aalB,bblB,vvlB]=svd(betaplv_reshape,'econ');

clear data4ica
data4ica.dimord='chan_time';
data4ica.time{1}=.001:.001:93.492;
data4ica.label={'1' '2' '3' '4' '5' '6' '7'};
cfg=[];
cfg.randomseed=13;
cfg.method='runica';
% data4ica.trial{1}=fullpow_reshape;
% runicaPow=ft_componentanalysis(cfg, data4ica);
% data4ica.trial{1}=fullplv_reshape;
% runicaPlv=ft_componentanalysis(cfg, data4ica);
data4ica.time{1}=.001:.001:6.678;
data4ica.trial{1}=thetapow_reshape;
runicaPowT=ft_componentanalysis(cfg, data4ica);
data4ica.trial{1}=thetaplv_reshape;
runicaPlvT=ft_componentanalysis(cfg, data4ica);
% data4ica.trial{1}=alphapow_reshape;
% runicaPowA=ft_componentanalysis(cfg, data4ica);
% data4ica.trial{1}=alphaplv_reshape;
% runicaPlvA=ft_componentanalysis(cfg, data4ica);
data4ica.trial{1}=betapow_reshape;
runicaPowB=ft_componentanalysis(cfg, data4ica);
% data4ica.trial{1}=betaplv_reshape;
% runicaPlvB=ft_componentanalysis(cfg, data4ica);

% runicaPlvT.topo(:,1)=-runicaPlvT.topo(:,1);
% runicaPlvT.trial{1}(1,:)=-runicaPlvT.trial{1}(1,:);



for llind=1:7
  %   rOica(:,:,llind)=reshape(runicaPow.trial{1}(llind,:),[63 14*size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
  %   rLica(:,:,llind)=reshape(runicaPlv.trial{1}(llind,:),[63 14*size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
  % %   pOica(:,:,llind)=reshape(vvo(:,llind)',[63 14*size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
  % %   pLica(:,:,llind)=reshape(vvl(:,llind)',[63 14*size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
  rOTica(:,:,llind)=reshape(runicaPowT.trial{1}(llind,:),[63 size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
  rLTica(:,:,llind)=reshape(runicaPlvT.trial{1}(llind,:),[63 size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
  %   pOTica(:,:,llind)=reshape(vvoT(:,llind)',[63 size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
  %   pLTica(:,:,llind)=reshape(vvlT(:,llind)',[63 size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
  %   rOAica(:,:,llind)=reshape(runicaPowA.trial{1}(llind,:),[63 size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
  %   rLAica(:,:,llind)=reshape(runicaPlvA.trial{1}(llind,:),[63 size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
  % %   pOAica(:,:,llind)=reshape(vvoA(:,llind)',[63 size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
  % %   pLAica(:,:,llind)=reshape(vvlA(:,llind)',[63 size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
  rOBica(:,:,llind)=reshape(runicaPowB.trial{1}(llind,:),[63 size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
  %   rLBica(:,:,llind)=reshape(runicaPlvB.trial{1}(llind,:),[63 size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
  %   pOBica(:,:,llind)=reshape(vvoB(:,llind)',[63 size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
  %   pLBica(:,:,llind)=reshape(vvlB(:,llind)',[63 size(gravelo_TPA_MSPN_avg{llind}.powspctrm,3)]);
  
  statplot.time=1;
  statplot.dimord='chan_time';
  statplot.label=gravelo_TPA_MSPN_avg{llind}.label;
  
  [af,bf,cf]=svd(rOTica(:,:,llind),'econ');
  if llind==1 || llind==2
    af=-af;cf=-cf;
  end
  statplot.avg=af(:,1);
  if 0 % for ICA alone figure
    figure(220);subplot(7,3,llind*3-2);bar(runicaPowT.topo(:,llind)); ylim([-2 2])
    figure(220);subplot(7,3,llind*3-1);cfg=[];cfg.latency=1;cfg.zlim='maxabs';cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
    figure(220);subplot(7,3,llind*3);plot(cf(:,1));ylim([-max(abs(cf(:,1)))-.04 max(abs(cf(:,1)))+.04])
  else % For ICA-ERP figure in paper
    close all
    if llind==3 % TA70
      timeadd=.07;
    elseif llind==2 % TA30
      timeadd=.02;
    else
      timeadd=0;
    end
    plot_ica(cf,runicaPowT,statplot,[-.1 1.15],timeadd,freqstat_TPA_MSPN_1wayANOVA_thetaall.time,llind,'PowThICA',fdir)
    icasaveOT{llind}.topo=statplot.avg;
    icasaveOT{llind}.bar=runicaPowT.topo(:,llind);
    icasaveOT{llind}.time=freqstat_TPA_MSPN_1wayANOVA_thetaall.time+timeadd;
    icasaveOT{llind}.course=cf(:,1);    
  end  
  
  [af,bf,cf]=svd(rLTica(:,:,llind),'econ');
%   if llind==1
%     af=-af;cf=-cf;
%   end
% % % P=polyfit(gravelo_TPA_MSPN_comb1_theta.plvabs(chanplot{1},1,[time 0 to 1.05]),cf(:,1),1);
  statplot.avg=af(:,1);
  if 0 % for ICA alone figure
    figure(230);subplot(7,3,llind*3-2);bar(runicaPlvT.topo(:,llind)); ylim([-2.2 2.2])
    figure(230);subplot(7,3,llind*3-1);cfg=[];cfg.latency=1;cfg.zlim='maxabs';cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
    figure(230);subplot(7,3,llind*3);plot(cf(:,1));ylim([-max(abs(cf(:,1)))-.04 max(abs(cf(:,1)))+.04])
  else % For ICA-ERP figure in paper
    close all
    timeadd=0; % primary one is llind=1 with AT70
    plot_ica(cf,runicaPlvT,statplot,[-.1 1.15],timeadd,freqstat_TPA_MSPN_1wayANOVA_thetaall.time,llind,'ItcThICA',fdir)
    icasaveLT{llind}.topo=statplot.avg;
    icasaveLT{llind}.bar=runicaPowT.topo(:,llind);
    icasaveLT{llind}.time=freqstat_TPA_MSPN_1wayANOVA_thetaall.time+timeadd;
    icasaveLT{llind}.course=cf(:,1);    
  end
  
  [af,bf,cf]=svd(rOBica(:,:,llind),'econ');
  statplot.avg=af(:,1);
  if 0 % for ICA alone figure
    figure(420);subplot(7,3,llind*3-2);bar(runicaPowB.topo(:,llind)); ylim([-2 2])
    figure(420);subplot(7,3,llind*3-1);cfg=[];cfg.latency=1;cfg.zlim='maxabs';cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
    figure(420);subplot(7,3,llind*3);plot(cf(:,1));ylim([-max(abs(cf(:,1)))-.04 max(abs(cf(:,1)))+.04])
  else % For ICA-ERP figure in paper
    close all
    timeadd=0; % primary one is llind=2 with AT70
    plot_ica(cf,runicaPowB,statplot,[-.1 1.15],timeadd,freqstat_TPA_MSPN_1wayANOVA_thetaall.time,llind,'PowBeICA',fdir)
    icasaveOB{llind}.topo=statplot.avg;
    icasaveOB{llind}.bar=runicaPowT.topo(:,llind);
    icasaveOB{llind}.time=freqstat_TPA_MSPN_1wayANOVA_thetaall.time+timeadd;
    icasaveOB{llind}.course=cf(:,1);    
  end
  
end

save('icasave.mat','icasave*','-append')

% Correlate (non-masked) ICA with data
for llind=1:7
  [af,bf,cf]=svd(rOTica(:,:,llind),'econ');
  pcarOTica=af(:,1)*cf(:,1)';
  [af,bf,cf]=svd(rOBica(:,:,llind),'econ');
  pcarOBica=af(:,1)*cf(:,1)';
  [af,bf,cf]=svd(rLTica(:,:,llind),'econ');
  pcarLTica=af(:,1)*cf(:,1)';
  for llcond=1:7
    rOTica_grave_corr(llind,llcond)=corr(reshape(rOTica(:,:,llind),[63*106 1]),reshape(squeeze(thetapow_reshape(llcond,:,:)),[63*106 1]));
    pcarOTica_grave_corr(llind,llcond)=corr(reshape(pcarOTica,[63*106 1]),reshape(squeeze(thetapow_reshape(llcond,:,:)),[63*106 1]));
    rOBica_grave_corr(llind,llcond)=corr(reshape(rOBica(:,:,llind),[63*106 1]),reshape(squeeze(betapow_reshape(llcond,:,:)),[63*106 1]));
    pcarOBica_grave_corr(llind,llcond)=corr(reshape(pcarOBica,[63*106 1]),reshape(squeeze(betapow_reshape(llcond,:,:)),[63*106 1]));
    rLTica_grave_corr(llind,llcond)=corr(reshape(rLTica(:,:,llind),[63*106 1]),reshape(squeeze(thetaplv_reshape(llcond,:,:)),[63*106 1]));
    pcarLTica_grave_corr(llind,llcond)=corr(reshape(pcarLTica,[63*106 1]),reshape(squeeze(thetaplv_reshape(llcond,:,:)),[63*106 1]));
  end
end



% 2) Power Theta averaged
mask_use=squeeze(freqstat_TPA_MSPN_1wayANOVA_theta.mask);
[rho,rho_nw,rhocondOT]=mask_corr(mask_use,runicaPowT,thetapow_reshape,0);
for llind=1:7
  figure;imagesc(rOTica(:,:,llind).*mask_use); colorbar; caxis([-0.2 0.2])
  figure;imagesc(mask_use.*squeeze(gravelo_TPA_MSPN_avg_theta{llind}.powspctrm));colorbar; caxis([-.4 .4])
end
% rho=    0.3420    0.2915    0.5854    0.3324    0.1777    0.0955    0.2303
    
% 2) Power Beta averaged
mask_use=squeeze(freqstat_TPA_MSPN_1wayANOVA_beta.mask);
[rho,rho_nw,rhocondOB]=mask_corr(mask_use,runicaPowB,betapow_reshape,0);
for llind=1:7
  figure;imagesc(rOBica(:,:,llind).*mask_use); colorbar; caxis([-0.2 0.2])
  figure;imagesc(mask_use.*squeeze(gravelo_TPA_MSPN_avg_beta{llind}.powspctrm));colorbar; caxis([-.4 .4])
end
% rho=     0.3602    0.6250    0.2434    0.2492    0.4613    0.3667    0.1864

% 3) ITC Theta averaged
mask_use=squeeze(plvstat_TPA_MSPN_1wayANOVA_theta.mask);
[rho,rho_nw,rhocondLT]=mask_corr(mask_use,runicaPlvT,thetaplv_reshape,0);
for llind=1:7
  figure;imagesc(rLTica(:,:,llind).*squeeze(plvstat_TPA_MSPN_1wayANOVA_theta.mask)); colorbar; caxis([-0.04 0.04])
  figure;imagesc(squeeze(plvstat_TPA_MSPN_1wayANOVA_theta.mask.*gravelo_TPA_MSPN_avg_theta{llind}.plvabs));colorbar; caxis([-.08 .08])
end
% rho = *0.7132    0.3244    0.4329    0.3616    0.1857    0.2565    0.1867

%%
% For Methods figure:
figure;imagesc(squeeze(plvstat_TPA_MSPN_1wayANOVA_theta.mask));colormap gray
figure;imagesc(runicaPlvT.topo(:,1));caxis([-1.5 1.5])
figure;imagesc(runicaPlvT.trial{1}(1,:));caxis([-.05 .05])
figure; imagesc(rLTica(:,:,1)); caxis([-.05 .05]);
figure;for llind=1:7,subplot(7,1,llind),imagesc(squeeze(plvstat_TPA_MSPN_1wayANOVA_theta.mask).*squeeze(thetaplv(llind,:,:)));caxis([-.15 .15]);end


load tfrsave.mat
load icasave.mat

topocorrOT=nan(7,9);
courseFCcorrOT=nan(7,9);
courseOPcorrOT=nan(7,9);
topocorrOB=nan(7,9);
courseFCcorrOB=nan(7,9);
courseOPcorrOB=nan(7,9);
topocorrLT=nan(7,9);
courseFCcorrLT=nan(7,9);
courseOPcorrLT=nan(7,9);
for ll=soalist
  for llind=1:7
    if isfield(tfrsaveOT{ll},'topo')
      topocorrOT(llind,ll)=corr(icasaveOT{llind}.topo,tfrsaveOT{ll}.topo{1});
    end
    time1=nearest(tfrsaveOT{ll}.time,icasaveOT{llind}.time(1));
    time2=nearest(tfrsaveOT{ll}.time,icasaveOT{llind}.time(end));
    courseFCcorrOT(llind,ll)=corr(icasaveOT{llind}.course,tfrsaveOT{ll}.courseFC(time1:time2)');
    courseOPcorrOT(llind,ll)=corr(icasaveOT{llind}.course,tfrsaveOT{ll}.courseOP(time1:time2)');
    if isfield(tfrsaveOB{ll},'topo')
      topocorrOB(llind,ll)=corr(icasaveOB{llind}.topo,tfrsaveOB{ll}.topo{1});
    end
    time1=nearest(tfrsaveOB{ll}.time,icasaveOB{llind}.time(1));
    time2=nearest(tfrsaveOB{ll}.time,icasaveOB{llind}.time(end));
    courseFCcorrOB(llind,ll)=corr(icasaveOB{llind}.course,tfrsaveOB{ll}.courseFC(time1:time2)');
    courseOPcorrOB(llind,ll)=corr(icasaveOB{llind}.course,tfrsaveOB{ll}.courseOP(time1:time2)');
    if isfield(tfrsaveLT{ll},'topo')
      topocorrLT(llind,ll)=corr(icasaveLT{llind}.topo,tfrsaveLT{ll}.topo{1});
    end
    time1=nearest(tfrsaveLT{ll}.time,icasaveLT{llind}.time(1));
    time2=nearest(tfrsaveLT{ll}.time,icasaveLT{llind}.time(end));
    courseFCcorrLT(llind,ll)=corr(icasaveLT{llind}.course,tfrsaveLT{ll}.courseFC(time1:time2)');
    courseOPcorrLT(llind,ll)=corr(icasaveLT{llind}.course,tfrsaveLT{ll}.courseOP(time1:time2)');
  end
end
figure;imagescc(topocorrOT)
figure;imagescc(courseFCcorrOT)
figure;imagescc(courseOPcorrOT)
figure;imagescc(topocorrOB)
figure;imagescc(courseFCcorrOB)
figure;imagescc(courseOPcorrOB)
figure;imagescc(topocorrLT)
figure;imagescc(courseFCcorrLT)
figure;imagescc(courseOPcorrLT)

abs(topocorrOT)>.45 & (abs(courseFCcorrOT)>.7 | abs(courseOPcorrOT)>.7)

max(abs(courseFCcorrOT),abs(courseOPcorrOT))


figure;imagesc((abs(topocorrOT)+max(abs(courseFCcorrOT),abs(courseOPcorrOT)))/2);caxis([0 1]);colormap('gray')
figure;imagesc((abs(topocorrOB)+max(abs(courseFCcorrOB),abs(courseOPcorrOB)))/2);caxis([0 1]);colormap('gray')
figure;imagesc((abs(topocorrLT)+max(abs(courseFCcorrLT),abs(courseOPcorrLT)))/2);caxis([0 1]);colormap('gray')

figure;imagesc([(abs(topocorrOT)+max(abs(courseFCcorrOT),abs(courseOPcorrOT)))/2]>.6)
figure;imagesc([(abs(topocorrOB)+max(abs(courseFCcorrOB),abs(courseOPcorrOB)))/2]>.6)
figure;imagesc([(abs(topocorrLT)+max(abs(courseFCcorrLT),abs(courseOPcorrLT)))/2]>.6)


% for llind=1:7
%   [af,bf,cf]=svd(rOica(:,:,llind),'econ');
%   statplot.avg=af(:,1);
%   statplot.time=1;
%   statplot.dimord='chan_time';
%   statplot.label=gravelo_TPA_MSPN_avg{llind}.label;
%   figure(120);subplot(7,3,llind*3-2);bar(runicaPow.topo(:,llind));%ylim([-3 3])
%   figure(120);subplot(7,3,llind*3-1);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
%   figure(120);subplot(7,3,llind*3);imagesc(reshape(cf(:,1),[14 106])); axis xy
% 
%   [af,bf,cf]=svd(rLica(:,:,llind),'econ');
%   statplot.avg=af(:,1);
%   figure(130);subplot(7,3,llind*3-2);bar(runicaPlv.topo(:,llind));%ylim([-3 3])
%   figure(130);subplot(7,3,llind*3-1);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
%   figure(130);subplot(7,3,llind*3);imagesc(reshape(cf(:,1),[14 106])); axis xy

%   [af,bf,cf]=svd(pOica(:,:,llind),'econ');
%   statplot.avg=af(:,1);
%   figure(140);subplot(7,3,llind*3-2);bar(aao(:,llind));
%   figure(140);subplot(7,3,llind*3-1);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
%   figure(140);subplot(7,3,llind*3);imagesc(reshape(cf(:,1),[14 106])); axis xy
% 
%   [af,bf,cf]=svd(pLica(:,:,llind),'econ');
%   statplot.avg=af(:,1);
%   figure(150);subplot(7,3,llind*3-2);bar(aal(:,llind));
%   figure(150);subplot(7,3,llind*3-1);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
%   figure(150);subplot(7,3,llind*3);imagesc(reshape(cf(:,1),[14 106])); axis xy

  
%   [af,bf,cf]=svd(pOTica(:,:,llind),'econ');
%   statplot.avg=af(:,1);
%   figure(240);subplot(7,3,llind*3-2);bar(aaoT(:,llind));
%   figure(240);subplot(7,3,llind*3-1);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
%   figure(240);subplot(7,3,llind*3);plot(cf(:,1));
% 
%   [af,bf,cf]=svd(pLTica(:,:,llind),'econ');
%   statplot.avg=af(:,1);
%   figure(250);subplot(7,3,llind*3-2);bar(aalT(:,llind));
%   figure(250);subplot(7,3,llind*3-1);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
%   figure(250);subplot(7,3,llind*3);plot(cf(:,1));

  
%   [af,bf,cf]=svd(rOAica(:,:,llind),'econ');
%   statplot.avg=af(:,1);
%   figure(320);subplot(7,3,llind*3-2);bar(runicaPowA.topo(:,llind));%ylim([-3 3])
%   figure(320);subplot(7,3,llind*3-1);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
%   figure(320);subplot(7,3,llind*3);plot(cf(:,1));
% 
%   [af,bf,cf]=svd(rLAica(:,:,llind),'econ');
%   statplot.avg=af(:,1);
%   figure(330);subplot(7,3,llind*3-2);bar(runicaPlvA.topo(:,llind));%ylim([-3 3])
%   figure(330);subplot(7,3,llind*3-1);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
%   figure(330);subplot(7,3,llind*3);plot(cf(:,1));

%   [af,bf,cf]=svd(pOAica(:,:,llind),'econ');
%   statplot.avg=af(:,1);
%   figure(340);subplot(7,3,llind*3-2);bar(aaoA(:,llind));
%   figure(340);subplot(7,3,llind*3-1);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
%   figure(340);subplot(7,3,llind*3);plot(cf(:,1));
% 
%   [af,bf,cf]=svd(pLAica(:,:,llind),'econ');
%   statplot.avg=af(:,1);
%   figure(350);subplot(7,3,llind*3-2);bar(aalA(:,llind));
%   figure(350);subplot(7,3,llind*3-1);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
%   figure(350);subplot(7,3,llind*3);plot(cf(:,1));
  
%   [af,bf,cf]=svd(rLBica(:,:,llind),'econ');
%   statplot.avg=af(:,1);
%   figure(430);subplot(7,3,llind*3-2);bar(runicaPlvB.topo(:,llind));%ylim([-3 3])
%   figure(430);subplot(7,3,llind*3-1);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
%   figure(430);subplot(7,3,llind*3);plot(cf(:,1));

%   [af,bf,cf]=svd(pOBica(:,:,llind),'econ');
%   statplot.avg=af(:,1);
%   figure(440);subplot(7,3,llind*3-2);bar(aaoB(:,llind));
%   figure(440);subplot(7,3,llind*3-1);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
%   figure(440);subplot(7,3,llind*3);plot(cf(:,1));
% 
%   [af,bf,cf]=svd(pLBica(:,:,llind),'econ');
%   statplot.avg=af(:,1);
%   figure(450);subplot(7,3,llind*3-2);bar(aalB(:,llind));
%   figure(450);subplot(7,3,llind*3-1);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
%   figure(450);subplot(7,3,llind*3);plot(cf(:,1));
% end



% %  ITC all
% mask_vector=reshape(plvstat_TPA_MSPN_1wayANOVA.mask,[1 63*14*106]);
% mask_vector_rep=repmat(mask_vector,[1 7]);
% mask_induse=find(mask_vector);  % 2.7%
% mask_induse_rep=find(mask_vector_rep); 
% clear rica_masked*
% for llind=1:7
%   rica_masked(llind,:)=runicaPlv.trial{1}(llind,mask_induse);
%   rica_masked_w(llind,:)=[rica_masked(llind,:)*runicaPlv.topo(1,llind) rica_masked(llind,:)*runicaPlv.topo(2,llind) rica_masked(llind,:)*runicaPlv.topo(3,llind) rica_masked(llind,:)*runicaPlv.topo(4,llind) rica_masked(llind,:)*runicaPlv.topo(5,llind) rica_masked(llind,:)*runicaPlv.topo(6,llind) rica_masked(llind,:)*runicaPlv.topo(7,llind)];
%   rica_masked_nw(llind,:)=repmat(rica_masked(llind,:),[1 7]);
% end
% maskdatavec=reshape(fullplv_reshape(:,mask_induse)',[1 7*2550]);
% for llind=1:7
%   [rho(llind),pval(llind)]=corr(maskdatavec',rica_masked_w(llind,:)');
%   [rho_nw(llind),pval_nw(llind)]=corr(maskdatavec',rica_masked_nw(llind,:)');
% end
% for llind=1:7
%   maskedICA=reshape(rLica(:,:,llind),[63 14 106]).*plvstat_TPA_MSPN_1wayANOVA.mask;
%   maskedData=plvstat_TPA_MSPN_1wayANOVA.mask.*gravelo_TPA_MSPN_avg{llind}.plvabs;
%   figure;imagesc(squeeze(mean(maskedICA,1)));axis xy; colorbar;caxis([-.01 .01])
%   figure;imagesc(squeeze(mean(maskedData,1)));axis xy; colorbar;caxis([-.01 .01])
% end


