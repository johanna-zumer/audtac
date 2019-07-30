
%%  Plotting band-specific line plots results
stim2nd= [  0 nan    0    0 0 .02 .07 nan .5];
clear chanplot
chanplot{1}={'Fz' 'Cz' 'F1' 'F2' 'FC1' 'FC2' 'C1' 'C2'}; % frontocentral
chanplot{2}={'CP5' 'POz' 'Pz' 'P3' 'P4' 'C4' 'O1' 'O2' 'P7' 'PO7'}; % occipital; % same as for ERP
% chanplot{3}={'C1' 'Cz' 'C2' 'CP1' 'CPz' 'CP2' 'P1' 'Pz' 'P2'}; % Centred on CPz
% chanplot{4}={'P3' 'P1' 'Pz' 'P2' 'P4' 'PO3' 'POz' 'PO4' 'O1' 'O2' 'Oz'};
chanlabel{1}='Frontocentral electrodes';
chanlabel{2}='Occipital-parietal electrodes';
soalist=[1 3 4 5 6 7 9];
tt=3;
fftaddflag=0;
synchasynch=0;
tacbasemax=min(-.15,soades-.15);
baseline2=[tacbasemax(1) tacbasemax(1)+.08];
itcdepsampflag=1;
ftver=0;
mcseed=13;
sleep=0;
if sleep
  ss=12;
  trialkc=-1; % change me
  iter=11;
  usetr=1;
else
  ss=10;
  iter=31;
  usetr=3;
  trialkc=-1;
end
plotplvflag=0;
plottfrflag=0;
peakfindflag=0;
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
plotallmask=1; % =0 means each subplot will use significane mask relevant for those sensors & time plotted only
% =1 means each subplot will use mask relevant for all sensors even those not plotted.
colorblindD  =[204 101 0]/256;
colorblindApT=[5 165 255]/256;
colorblindMpN=[0 59 179]/256;
colorblindT  =[255 51 166]/256;
colorblindA  =[0 158 115]/256;
colorblindM  =[0 0 0]/256;
colorblindN  =[128 128 128]/256;


% Do unisensory plotting & peak finding first.

load(['grindTFR_UniMsNul_sleep' num2str(sleep) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'grind*');
load(['stats_TFR_UniMSNul_sleep0_trialkc-1_itc.mat'])

ll=5;
cfg=[];
cfg.avgoverrpt='yes';
tmp2=ft_selectdata(cfg,grindlo_tNulAlone{5,10});

tmp4=ft_selectdata(cfg,grindlo_tTacAlone{5,10});
tmp4.tfrmaskS=logical(zeros(size(tmp4.powspctrm)));
tmp4.plvmaskS=logical(zeros(size(tmp4.powspctrm)));
tmp4.tfrmaskL=logical(zeros(size(tmp4.powspctrm)));
tmp4.plvmaskL=logical(zeros(size(tmp4.powspctrm)));
tmp4.tfrmaskS(:,:,dsearchn(tmp4.time',stattl_mc_TacVsNul_short{ll,ss}.time(1)):dsearchn(tmp4.time',stattl_mc_TacVsNul_short{ll,ss}.time(end)))=stattl_mc_TacVsNul_short{ll,ss}.mask;
tmp4.tfrmaskL(:,:,dsearchn(tmp4.time',stattl_mc_TacVsNul_long{ll,ss}.time(1)): dsearchn(tmp4.time',stattl_mc_TacVsNul_long{ll,ss}.time(end))) =stattl_mc_TacVsNul_long{ll,ss}.mask;
tmp4.plvmaskS(:,:,dsearchn(tmp4.time',stattl_mcplv_TacVsNul_short_itcdepT{ll,ss}.time(1)):dsearchn(tmp4.time',stattl_mcplv_TacVsNul_short_itcdepT{ll,ss}.time(end)))=stattl_mcplv_TacVsNul_short_itcdepT{ll,ss}.mask;
tmp4.plvmaskL(:,:,dsearchn(tmp4.time',stattl_mcplv_TacVsNul_long_itcdepT{ll,ss}.time(1)): dsearchn(tmp4.time',stattl_mcplv_TacVsNul_long_itcdepT{ll,ss}.time(end))) =stattl_mcplv_TacVsNul_long_itcdepT{ll,ss}.mask;

tmp9=ft_selectdata(cfg,grindlo_tAudAlone{5,10});
tmp9.tfrmaskS=logical(zeros(size(tmp9.powspctrm)));
tmp9.plvmaskS=logical(zeros(size(tmp9.powspctrm)));
tmp9.tfrmaskL=logical(zeros(size(tmp9.powspctrm)));
tmp9.plvmaskL=logical(zeros(size(tmp9.powspctrm)));
tmp9.tfrmaskS(:,:,dsearchn(tmp9.time',stattl_mc_AudVsNul_short{ll,ss}.time(1)):dsearchn(tmp9.time',stattl_mc_AudVsNul_short{ll,ss}.time(end)))=stattl_mc_AudVsNul_short{ll,ss}.mask;
tmp9.tfrmaskL(:,:,dsearchn(tmp9.time',stattl_mc_AudVsNul_long{ll,ss}.time(1)): dsearchn(tmp9.time',stattl_mc_AudVsNul_long{ll,ss}.time(end))) =stattl_mc_AudVsNul_long{ll,ss}.mask;
tmp9.plvmaskS(:,:,dsearchn(tmp9.time',stattl_mcplv_AudVsNul_short_itcdepT{ll,ss}.time(1)):dsearchn(tmp9.time',stattl_mcplv_AudVsNul_short_itcdepT{ll,ss}.time(end)))=stattl_mcplv_AudVsNul_short_itcdepT{ll,ss}.mask;
tmp9.plvmaskL(:,:,dsearchn(tmp9.time',stattl_mcplv_AudVsNul_long_itcdepT{ll,ss}.time(1)): dsearchn(tmp9.time',stattl_mcplv_AudVsNul_long_itcdepT{ll,ss}.time(end))) =stattl_mcplv_AudVsNul_long_itcdepT{ll,ss}.mask;

clear pow*
pow{1}=tmp2;
pow{2}=tmp4;
pow{3}=tmp9;

stattimwinS=[stattl_mc_AudVsNul_short{ll,ss}.time(1) stattl_mc_AudVsNul_short{ll,ss}.time(end)];
stattimwinL=[stattl_mc_AudVsNul_long{ll,ss}.time(1) stattl_mc_AudVsNul_long{ll,ss}.time(end)];

% average over bands first
for pp=1:3
  for ff=1:2 % two diff freq bands
    cfg=[];
    if ff==1
      cfg.frequency   = [4 6.5];
    elseif ff==2
      cfg.frequency   = [14 30];
    elseif ff==3
      cfg.frequency   = [8 12]; % not used
    end
    cfg.avgoverfreq = 'yes';
    cfg.nanmean     = 'yes';
    powf{pp,ff}=ft_selectdata(cfg,pow{pp});
  end
end
clear pow
% baseline correct first, then remove baseline area for plotting

for ff=1:2
  base=powf{2,ff};
  base.powspctrm=repmat(nanmean(cat(3,nanmean(powf{1,ff}.powspctrm(:,:,49:53),3),nanmean(powf{2,ff}.powspctrm(:,:,49:53),3),nanmean(powf{3,ff}.powspctrm(:,:,49:53),3)),3),[1 1 length(powf{1,ff}.time)]);
  base.plvabs=repmat(nanmean(cat(3,nanmean(powf{1,ff}.plvabs(:,:,49:53),3),nanmean(powf{2,ff}.plvabs(:,:,49:53),3),nanmean(powf{3,ff}.plvabs(:,:,49:53),3)),3),[1 1 length(powf{1,ff}.time)]);
  base.tfrmaskL=repmat(nanmean(powf{2,ff}.tfrmaskL(:,:,49:53),3),[1 1 length(powf{1,ff}.time)]);  % use 'nul' as baseline for all three conditions for plotting
  base.plvmaskL=repmat(nanmean(powf{2,ff}.plvmaskL(:,:,49:53),3),[1 1 length(powf{1,ff}.time)]);  % use 'nul' as baseline for all three conditions for plotting
  for pp=1:3
    cfg=[];
    if pp>1
      cfg.parameter={'powspctrm' 'plvabs' 'tfrmaskL' 'plvmaskL'};
    else
      cfg.parameter={'powspctrm' 'plvabs'};
    end
    cfg.operation='subtract';
    powfb{pp,ff}=ft_math(cfg,powf{pp,ff},base)
    
    
    cfg=[];
    cfg.latency=[tacbasemax(ll) 1.3];
    powfb{pp,ff}=ft_selectdata(cfg,powfb{pp,ff});
    
    powfb{pp,ff}.tfravg=squeeze(powfb{pp,ff}.powspctrm);
    powfb{pp,ff}.plvavg=squeeze(powfb{pp,ff}.plvabs);
    if pp>1
      powfb{pp,ff}.tfrmaskL=squeeze(powfb{pp,ff}.tfrmaskL);
      powfb{pp,ff}.plvmaskL=squeeze(powfb{pp,ff}.plvmaskL);
    end
    powfb{pp,ff}.dimord='chan_time';
    powfb{pp,ff}=rmfield(powfb{pp,ff},'freq');
    powfb{pp,ff}=rmfield(powfb{pp,ff},'powspctrm');
    powfb{pp,ff}=rmfield(powfb{pp,ff},'plvabs');
    if plotallmask && pp>1
      powfb{pp,ff}.tfrmaskL=repmat(ceil(mean(powfb{pp,ff}.tfrmaskL,1)),[size(powfb{pp,ff}.tfrmaskL,1) 1]);
      powfb{pp,ff}.plvmaskL=repmat(ceil(mean(powfb{pp,ff}.plvmaskL,1)),[size(powfb{pp,ff}.plvmaskL,1) 1]);
    end
  end
end

% peak finding
chpl{1}=match_str(powfb{2,1}.label,chanplot{1});
chpl{2}=match_str(powfb{2,1}.label,chanplot{2});
[mx,mind]=max(nanmean(powfb{2,1}.tfravg(chpl{1},:),1));   % Tactile theta FC
timeTacThetaTfr=powfb{2,1}.time(mind); % 0.21
[mx,mind]=max(nanmean(powfb{3,1}.tfravg(chpl{1},:),1));   % Tactile theta FC
timeAudThetaTfr=powfb{3,1}.time(mind); % 0.16
[mx,mind]=max(nanmean(powfb{2,2}.tfravg(chpl{2},:),1));   % Tactile alpha beta OP
timeTacAlpbetTfr=powfb{2,2}.time(mind); % 0.94
%plv
[mx,mind]=max(nanmean(powfb{2,1}.plvavg(chpl{1},:),1));   % Tactile theta FC
timeTacThetaPlv=powfb{2,1}.time(mind); % 0.14
[mx,mind]=max(nanmean(powfb{3,1}.plvavg(chpl{1},:),1));   % Tactile theta FC
timeAudThetaPlv=powfb{3,1}.time(mind); % 0.13


for cg=1:length(chanplot)
  for ff=1:2
    cfg=[];
    cfg.layout='elec1010.lay';
    if ff==1
      cfg.ylim=[-2.5 2.5];  % or [1 2.5]?
    elseif ff==2 % beta
      cfg.ylim=[-1 1];
    elseif ff==3 % alpha
      cfg.ylim=[-2.5 2.5];
    end
    cfg.linewidth=3;
    cfg.xlim=[-.6 1.8];
    cfg.channel=chanplot{cg};
    cfg.graphcolor=[colorblindN; colorblindT; colorblindA];
    cfg.interactive='no';
    cfg.parameter='tfravg';
    figure(10*(cg+1)+(ff-1)*100)
    ft_singleplotER(cfg,powfb{1,ff},powfb{2,ff},powfb{3,ff});
    hold on;plot(powfb{1,ff}.time,0,'k');
    set(gca,'XTick',[-.6:.1:1.8])
    set(gca,'XTickLabel',{ '' '-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'  '' ' ' ''  ' ' '1.5' '' ' ' '' })
    set(gca,'FontSize',30)
    title([])
    plot([stattimwinL(1) stattimwinL(1)],cfg.ylim,'k--','LineWidth',6)
    plot([stattimwinL(2) stattimwinL(2)],cfg.ylim,'k--','LineWidth',6)
    plot([0 0],cfg.ylim,'Color',colorblindT,'LineWidth',6)
    plot([soades(ll) soades(ll)],cfg.ylim,'Color',colorblindA,'LineWidth',6)
    plot(cfg.xlim,[0 0],'Color','k')
    for pp=2:3
      highlight=logical(powfb{pp,ff}.tfrmaskL(17,:));
      begsample = find(diff([0 highlight 0])== 1);
      endsample = find(diff([0 highlight 0])==-1)-1;
      for ii=1:length(begsample)
        begx=powfb{pp,ff}.time(begsample(ii));
        begy=powfb{pp,ff}.time(endsample(ii));
        if pp==2
          plot([begx begy], [cfg.ylim(2)-.5 cfg.ylim(2)-.5],'Color',colorblindT,'LineWidth',6)
        elseif pp==3
          plot([begx begy], [cfg.ylim(2)-.3 cfg.ylim(2)-.3],'Color',colorblindA,'LineWidth',6)
        end
      end
    end % pp
    axis([cfg.xlim(1) cfg.xlim(2) cfg.ylim(1) cfg.ylim(2)])
    %     legend('Null','Tactile','Auditory')
  end % ff
end % cg
print(20,[fdir 'tfr_UniNul_FC_' num2str(ll) '_theta' '.eps'],'-depsc2')
print(20+100,[fdir 'tfr_UniNul_FC_' num2str(ll) '_beta' '.eps'],'-depsc2')
print(30,[fdir 'tfr_UniNul_OP_' num2str(ll) '_theta' '.eps'],'-depsc2')
print(30+100,[fdir 'tfr_UniNul_OP_' num2str(ll) '_beta' '.eps'],'-depsc2')


for cg=1:length(chanplot)
  for ff=1:2
    cfg=[];
    cfg.layout='elec1010.lay';
    if ff==1
      cfg.ylim=[-.1 .4];  % or [1 2.5]?
    elseif ff==2
      cfg.ylim=[-.1 .4];
    end
    cfg.linewidth=3;
    cfg.xlim=[-.6 1.8];
    cfg.channel=chanplot{cg};
    cfg.graphcolor=[colorblindN; colorblindT; colorblindA];
    cfg.interactive='no';
    cfg.parameter='plvavg';
    figure(10*(cg+3)+(ff-1)*100)
    ft_singleplotER(cfg,powfb{1,ff},powfb{2,ff},powfb{3,ff});
    hold on;plot(powfb{1,ff}.time,0,'k');
    set(gca,'XTick',[-.6:.1:1.8])
    set(gca,'XTickLabel',{ '' '-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'  '' ' ' ''  ' ' '1.5' '' ' ' '' })
    set(gca,'FontSize',30)
    title([])
    plot([stattimwinL(1) stattimwinL(1)],cfg.ylim,'k--','LineWidth',6)
    plot([stattimwinL(2) stattimwinL(2)],cfg.ylim,'k--','LineWidth',6)
    plot([0 0],cfg.ylim,'Color',colorblindT,'LineWidth',6)
    plot([soades(ll) soades(ll)],cfg.ylim,'Color',colorblindA,'LineWidth',6)
    plot(cfg.xlim,[0 0],'Color','k')
    for pp=2:3
      highlight=logical(powfb{pp,ff}.plvmaskL(17,:));
      begsample = find(diff([0 highlight 0])== 1);
      endsample = find(diff([0 highlight 0])==-1)-1;
      for ii=1:length(begsample)
        begx=powfb{pp,ff}.time(begsample(ii));
        begy=powfb{pp,ff}.time(endsample(ii));
        if pp==2
          plot([begx begy], [cfg.ylim(2)-.04 cfg.ylim(2)-.04],'Color',colorblindT,'LineWidth',6)
        elseif pp==3
          plot([begx begy], [cfg.ylim(2)-.02 cfg.ylim(2)-.02],'Color',colorblindA,'LineWidth',6)
        end
      end
    end % pp
    axis([cfg.xlim(1) cfg.xlim(2) cfg.ylim(1) cfg.ylim(2)])
    %     legend('Null','Tactile','Auditory')
  end % ff
end % cg
print(40,[fdir 'plv_UniNul_FC_' num2str(ll) '_theta' '.eps'],'-depsc2')
print(40+100,[fdir 'plv_UniNul_FC_' num2str(ll) '_alpbet' '.eps'],'-depsc2')
print(50,[fdir 'plv_UniNul_OP_' num2str(ll) '_theta' '.eps'],'-depsc2')
print(50+100,[fdir 'plv_UniNul_OP_' num2str(ll) '_alpbet' '.eps'],'-depsc2')


% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multisensory interaction

load(['freqstat_sleep' num2str(sleep) '.mat']);  % see eeg_freqanalysis_ANOVA.m

mask_use=freqstat_TPA_MSPN_1wayANOVA_theta.mask;
mask_use(:,1:10)=0;
mask_use(:,39:end)=0;  % see thetaTFP_anova_twall in eeg_freqanalysis_ANOVA.m
freqstat_TPA_MSPN_1wayANOVA_theta.stat_tw1=freqstat_TPA_MSPN_1wayANOVA_theta.stat.*mask_use;
freqstat_TPA_MSPN_1wayANOVA_theta.mask_tw1=mask_use;

mask_use=freqstat_TPA_MSPN_1wayANOVA_beta.mask;
mask_use(:,23:end)=0;
freqstat_TPA_MSPN_1wayANOVA_beta.stat_tw1=freqstat_TPA_MSPN_1wayANOVA_beta.stat.*mask_use;
freqstat_TPA_MSPN_1wayANOVA_beta.mask_tw1=mask_use;
mask_use=freqstat_TPA_MSPN_1wayANOVA_beta.mask;
mask_use(:,1:22)=0;
freqstat_TPA_MSPN_1wayANOVA_beta.stat_tw2=freqstat_TPA_MSPN_1wayANOVA_beta.stat.*mask_use;
freqstat_TPA_MSPN_1wayANOVA_beta.mask_tw2=mask_use;

mask_use=plvstat_TPA_MSPN_1wayANOVA_theta.mask;
mask_use(:,46:end)=0;  % see thetaITC_anova_twall in eeg_freqanalysis_ANOVA.m
plvstat_TPA_MSPN_1wayANOVA_theta.stat_tw1=plvstat_TPA_MSPN_1wayANOVA_theta.stat.*mask_use;
plvstat_TPA_MSPN_1wayANOVA_theta.mask_tw1=mask_use;

chanplot_theta{1}=freqstat_TPA_MSPN_1wayANOVA_theta.label(find(mean(freqstat_TPA_MSPN_1wayANOVA_theta.mask,3)));

chanplot_beta{1}=freqstat_TPA_MSPN_1wayANOVA_beta.label(find(mean(freqstat_TPA_MSPN_1wayANOVA_beta.mask,3)));
chanplot_beta{2}=freqstat_TPA_MSPN_1wayANOVA_beta.label(find(mean(freqstat_TPA_MSPN_1wayANOVA_beta.mask_tw1,3)));
chanplot_beta{3}=freqstat_TPA_MSPN_1wayANOVA_beta.label(find(mean(freqstat_TPA_MSPN_1wayANOVA_beta.mask_tw2,3)));

chanplot_theta_itc{1}=plvstat_TPA_MSPN_1wayANOVA_theta.label(find(mean(plvstat_TPA_MSPN_1wayANOVA_theta.mask,3)));

for ll=soalist
% for ll=[1 3 4 5 6]
% for ll=4
  close all
%   clear grave* stat* tmp*
  clear tmp*
  load([edir 'statsgrave_TFR_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '_usetr' num2str(usetr) '_mcseed' num2str(mcseed) '_itc' '.mat']);
  %   baseline1=[tacbasemax(ll) tacbasemax(ll)+.08];
  stattimwin=[stattl_mc_comb1.time(1) stattl_mc_comb1.time(end)];
  
  
  cfg=[];
  cfg.avgoverfreq='yes';
  cfg.frequency=[4 6.5];
  gravelo_tacPaud_comb1_theta=ft_selectdata(cfg,gravelo_tacPaud_comb1);
  gravelo_tacMSpN_comb1_theta=ft_selectdata(cfg,gravelo_tacMSpN_comb1);
  gravelo_TPA_MSPN_comb1_theta=ft_selectdata(cfg,gravelo_TPA_MSPN_comb1);
  cfg.frequency=[8 12.5];
  gravelo_tacPaud_comb1_alpha=ft_selectdata(cfg,gravelo_tacPaud_comb1);
  gravelo_tacMSpN_comb1_alpha=ft_selectdata(cfg,gravelo_tacMSpN_comb1);
  gravelo_TPA_MSPN_comb1_alpha=ft_selectdata(cfg,gravelo_TPA_MSPN_comb1);
  cfg.frequency=[14 30];
  gravelo_tacPaud_comb1_beta=ft_selectdata(cfg,gravelo_tacPaud_comb1);
  gravelo_tacMSpN_comb1_beta=ft_selectdata(cfg,gravelo_tacMSpN_comb1);
  gravelo_TPA_MSPN_comb1_beta=ft_selectdata(cfg,gravelo_TPA_MSPN_comb1);

  
  if 0
    for xlimind=1:2
      if xlimind==1
        xlimxlim=[-.1 1.15];
        purpose='ica'; % main figure alongside ICA
      elseif xlimind==2
        xlimxlim=[-.6 1.8];
        purpose='within';  % supplemntary figure
      end
      % Plotting Power
      tfrsaveOT{ll} = plottfr(gravelo_TPA_MSPN_comb1_theta,gravelo_tacPaud_comb1_theta,gravelo_tacMSpN_comb1_theta,xlimxlim,purpose,stattl_mc_comb1_theta,'theta','powspctrm',ll,tacbasemax,plotallmask,chanplot,stattimwin,soades,fdir,'mask');
      plottfr(gravelo_TPA_MSPN_comb1_alpha,gravelo_tacPaud_comb1_alpha,gravelo_tacMSpN_comb1_alpha,xlimxlim,purpose,stattl_mc_comb1_alpha,'alpha','powspctrm',ll,tacbasemax,plotallmask,chanplot,stattimwin,soades,fdir,'mask')
      tfrsaveOB{ll} = plottfr(gravelo_TPA_MSPN_comb1_beta, gravelo_tacPaud_comb1_beta, gravelo_tacMSpN_comb1_beta, xlimxlim,purpose,stattl_mc_comb1_beta, 'beta', 'powspctrm',ll,tacbasemax,plotallmask,chanplot,stattimwin,soades,fdir,'mask');
      
      % Plotting PLV
      tfrsaveLT{ll} = plottfr(gravelo_TPA_MSPN_comb1_theta,gravelo_tacPaud_comb1_theta,gravelo_tacMSpN_comb1_theta,xlimxlim,purpose,stattl_mc_comb1plvabsDepT_theta,'theta','plvabs',ll,tacbasemax,plotallmask,chanplot,stattimwin,soades,fdir,'mask');
      plottfr(gravelo_TPA_MSPN_comb1_alpha,gravelo_tacPaud_comb1_alpha,gravelo_tacMSpN_comb1_alpha,xlimxlim,purpose,stattl_mc_comb1plvabsDepT_alpha,'alpha','plvabs',ll,tacbasemax,plotallmask,chanplot,stattimwin,soades,fdir,'mask')
      plottfr(gravelo_TPA_MSPN_comb1_beta, gravelo_tacPaud_comb1_beta, gravelo_tacMSpN_comb1_beta, xlimxlim,purpose,stattl_mc_comb1plvabsDepT_beta, 'beta', 'plvabs',ll,tacbasemax,plotallmask,chanplot,stattimwin,soades,fdir,'mask')
    end
  end
  
  purpose='anova';
  xlimxlim=[-.1 1.15];
  statstructuse=freqstat_TPA_MSPN_1wayANOVA_theta;
  statstructuse.time=statstructuse.time+stim2nd(ll);
  plottfr(gravelo_TPA_MSPN_comb1_theta,gravelo_tacPaud_comb1_theta,gravelo_tacMSpN_comb1_theta,xlimxlim,purpose,statstructuse,'theta','powspctrm',ll,tacbasemax,plotallmask,chanplot_theta,stattimwin,soades,fdir,'mask');
  plottfr(gravelo_TPA_MSPN_comb1_theta,gravelo_tacPaud_comb1_theta,gravelo_tacMSpN_comb1_theta,xlimxlim,purpose,statstructuse,'theta','powspctrm',ll,tacbasemax,plotallmask,chanplot_theta,stattimwin,soades,fdir,'mask_tw1');
  
  statstructuse=freqstat_TPA_MSPN_1wayANOVA_beta;
  statstructuse.time=statstructuse.time+stim2nd(ll);
  plottfr(gravelo_TPA_MSPN_comb1_beta, gravelo_tacPaud_comb1_beta, gravelo_tacMSpN_comb1_beta, xlimxlim,purpose,statstructuse, 'beta', 'powspctrm',ll,tacbasemax,plotallmask,chanplot_beta,stattimwin,soades,fdir,'mask');
  plottfr(gravelo_TPA_MSPN_comb1_beta, gravelo_tacPaud_comb1_beta, gravelo_tacMSpN_comb1_beta, xlimxlim,purpose,statstructuse, 'beta', 'powspctrm',ll,tacbasemax,plotallmask,chanplot_beta,stattimwin,soades,fdir,'mask_tw1');
  plottfr(gravelo_TPA_MSPN_comb1_beta, gravelo_tacPaud_comb1_beta, gravelo_tacMSpN_comb1_beta, xlimxlim,purpose,statstructuse, 'beta', 'powspctrm',ll,tacbasemax,plotallmask,chanplot_beta,stattimwin,soades,fdir,'mask_tw2');
  
  statstructuse=plvstat_TPA_MSPN_1wayANOVA_theta;
  statstructuse.time=statstructuse.time+stim2nd(ll);
  plottfr(gravelo_TPA_MSPN_comb1_theta,gravelo_tacPaud_comb1_theta,gravelo_tacMSpN_comb1_theta,xlimxlim,purpose,statstructuse,'theta','plvabs',ll,tacbasemax,plotallmask,chanplot_theta_itc,stattimwin,soades,fdir,'mask');
  plottfr(gravelo_TPA_MSPN_comb1_theta,gravelo_tacPaud_comb1_theta,gravelo_tacMSpN_comb1_theta,xlimxlim,purpose,statstructuse,'theta','plvabs',ll,tacbasemax,plotallmask,chanplot_theta_itc,stattimwin,soades,fdir,'mask_tw1');
  
    % power
  if plottfrflag
    clear pow*
    
%     plottfr(gravelo_TPA_MSPN_comb1_theta,gravelo_tacPaud_comb1_theta,gravelo_tacMSpN_comb1_theta,[-.6 1.8],stattl_mc_comb2_theta,'theta','powspctrm',ll,tacbasemax,plotallmask,chanplot,stattimwin,soades,fdir)
    
    cfg=[];
    tmpA=ft_selectdata(cfg,gravelo_TPA_MSPN_comb1);
    tmpA.mask=logical(zeros(size(tmpA.powspctrm)));
    tmpA.mask(:,:,dsearchn(tmpA.time',stattl_mc_comb1.time(1)):dsearchn(tmpA.time',stattl_mc_comb1.time(end)))=stattl_mc_comb1.mask;
    tmpsA=ft_selectdata(cfg,gravelo_tacMSpN_comb1);
    tmpsA.mask=logical(zeros(size(tmpsA.powspctrm)));
    tmpsA.mask(:,:,dsearchn(tmpsA.time',stattl_mc_comb1.time(1)):dsearchn(tmpsA.time',stattl_mc_comb1.time(end)))=stattl_mc_comb1.mask;
    tmpuA=ft_selectdata(cfg,gravelo_tacPaud_comb1);
    tmpuA.mask=logical(zeros(size(tmpuA.powspctrm)));
    tmpuA.mask(:,:,dsearchn(tmpuA.time',stattl_mc_comb1.time(1)):dsearchn(tmpuA.time',stattl_mc_comb1.time(end)))=stattl_mc_comb1.mask;
    
    pow{1}=tmpuA;
    pow{2}=tmpsA;
    pow{3}=tmpA;
    
    % average over bands first
    for pp=1:3
      for ff=1:2 % two diff freq bands
        cfg=[];
        if ff==1 % 6 , 9 doesn't matter
          %           if ll==1 || ll==3
          %             cfg.frequency   = [4 8.5];
          %           elseif ll==5 || ll==7 || ll==4 || ll==6  || ll==9
          cfg.frequency   = [4 6.5];
          %           end
        elseif ff==2
          %           if ll==1 || ll==3
          %             cfg.frequency   = [10 20];
          %           elseif ll==5 || ll==7 || ll==4 || ll==6  || ll==9
          cfg.frequency   = [8 20];
          %           end
        else
          disp('help');keyboard
        end
        cfg.avgoverfreq = 'yes';
        cfg.nanmean     = 'yes';
        powf{pp,ff}=ft_selectdata(cfg,pow{pp});
      end
    end
    clear pow
    
    % baseline correct first, then remove baseline area for plotting
    for ff=1:2
      %       base=nanmean(cat(3,nanmean(powf{1,ff}.powspctrm(:,:,49:53),3),nanmean(powf{2,ff}.powspctrm(:,:,49:53),3)),3);
      baset=powf{2,ff};
      baset.powspctrm=repmat(nanmean(cat(3,nanmean(powf{1,ff}.powspctrm(:,:,49:53),3),nanmean(powf{2,ff}.powspctrm(:,:,49:53),3)),3),[1 1 length(powf{1,ff}.time)]);
      %       base.plvabs=repmat(nanmean(cat(3,nanmean(powf{1,ff}.plvabs(:,:,49:53),3),nanmean(powf{2,ff}.plvabs(:,:,49:53),3),nanmean(powf{3,ff}.plvabs(:,:,49:53),3)),3),[1 1 length(powf{1,ff}.time)]);
      baset.mask=repmat(nanmean(powf{2,ff}.mask(:,:,49:53),3),[1 1 length(powf{1,ff}.time)]);  % use 'nul' as baseline for all three conditions for plotting
      %         cfg=[];
      %         cfg.baseline=[gravelo_TPA_MSPN_comb1.time(49) gravelo_TPA_MSPN_comb1.time(53)];
      %         cfg.baselinetype='absolute';
      %         cfg.parameter={'powspctrm' 'mask'};
      %         powfb{pp,ff}=ft_freqbaseline(cfg, powf{pp,ff});
      cfg=[];
      cfg.parameter={'powspctrm' 'mask'};
      cfg.operation='subtract';
      powfb{1,ff}=ft_math(cfg,powf{1,ff},baset)
      powfb{2,ff}=ft_math(cfg,powf{2,ff},baset)
      powfb{3,ff}=powf{3,ff};
      
      for pp=1:3
        cfg=[];
        if ll==1
          if ff==1
            cfg.latency=[tacbasemax(ll) 1.2];
          elseif ff==2
            cfg.latency=[tacbasemax(ll) 1.2];
          end
        elseif ll==3
          cfg.latency=[tacbasemax(ll) 1.3];
        elseif ll==4
          cfg.latency=[tacbasemax(ll) 1.3];
        elseif ll==5
          cfg.latency=[tacbasemax(ll) 1.3];
        elseif ll==6
          cfg.latency=[tacbasemax(ll) 1.3+.02];
        elseif ll==7
          cfg.latency=[tacbasemax(ll) 1.3+.07];
        elseif ll==9
          if ff==1
            cfg.latency=[tacbasemax(ll) 1.3+.2];
          elseif ff==2
            cfg.latency=[tacbasemax(ll) 1.3+.50];
          end
        end
        powfb{pp,ff}=ft_selectdata(cfg,powfb{pp,ff});
        
        powfb{pp,ff}.avg=squeeze(powfb{pp,ff}.powspctrm);
        powfb{pp,ff}.mask=squeeze(powfb{pp,ff}.mask);
        powfb{pp,ff}.dimord='chan_time';
        powfb{pp,ff}=rmfield(powfb{pp,ff},'freq');
        powfb{pp,ff}=rmfield(powfb{pp,ff},'powspctrm');
        if plotallmask
          powfb{pp,ff}.mask=repmat(ceil(mean(powfb{pp,ff}.mask,1)),[size(powfb{pp,ff}.mask,1) 1]);
        end
        if ff==2 && ll==1 % remove stat window from alp/bet when it's only for 8 Hz (from theta)
          powfb{pp,ff}.mask=zeros(size(powfb{pp,ff}.mask));
        end
      end % pp
    end
    
    
    for cg=1:length(chanplot)
      for ff=1:2
        cfg=[];
        cfg.parameter='avg';
        cfg.layout='elec1010.lay';
        if ff==1
          cfg.ylim=[-2.5 2.5];  % or [1 2.5]?
        elseif ff==2
          cfg.ylim=[-2.5 2.5];
        end
        cfg.linewidth=3;
        cfg.xlim=[-.6 1.8];
        if cg>length(chanplot)
          cfg.channel=tmpd5.label(any(tmpd5.mask,2));
        else
          cfg.channel=chanplot{cg};
        end
        cfg.graphcolor=[colorblindApT; colorblindMpN; colorblindD];
        cfg.interactive='no';
        cfg.maskparameter='mask';
        cfg.maskstyle='box'; % default
        %     subplot(4,1,cg)
        %     subplot(1,2,cg)
        figure(ll+10*(cg+1)+(ff-1)*100)
        ft_singleplotER(cfg,powfb{1,ff},powfb{2,ff},powfb{3,ff});
        hold on;plot(powfb{1,ff}.time,0,'k');
        set(gca,'XTick',[-.6:.1:1.8])
        set(gca,'XTickLabel',{ '' '-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'  '' ' ' ''  ' ' '1.5' '' ' ' '' })
        set(gca,'FontSize',30)
        title([])
        plot([stattimwin(1) stattimwin(1)],cfg.ylim,'k--','LineWidth',6)
        plot([stattimwin(2) stattimwin(2)],cfg.ylim,'k--','LineWidth',6)
        plot([0 0],cfg.ylim,'Color',colorblindT,'LineWidth',6)
        plot([soades(ll) soades(ll)],cfg.ylim,'Color',colorblindA,'LineWidth',6)
        %         plot(cfg.xlim,[0 0],'Color',coloruse(10,:),'LineWidth',1)
        plot(cfg.xlim,[0 0],'Color','k')
        axis([cfg.xlim(1) cfg.xlim(2) cfg.ylim(1) cfg.ylim(2)])
        if cg==3
          legend('Sum Unisensory','MultSens + Null','SumUnisens - MultsensNull')
        end
      end % ff
    end % cg
    print(ll+20,[fdir 'tfr_tacPaud_MSpN_diff_FC_' num2str(ll) '_theta' '.png'],'-dpng')
    print(ll+20+100,[fdir 'tfr_tacPaud_MSpN_diff_FC_' num2str(ll) '_alpbet' '.png'],'-dpng')
    print(ll+30,[fdir 'tfr_tacPaud_MSpN_diff_OP_' num2str(ll) '_theta' '.png'],'-dpng')
    print(ll+30+100,[fdir 'tfr_tacPaud_MSpN_diff_OP_' num2str(ll) '_alpbet' '.png'],'-dpng')
    print(ll+20,[fdir 'tfr_tacPaud_MSpN_diff_FC_' num2str(ll) '_theta' '.eps'],'-depsc')
    print(ll+20+100,[fdir 'tfr_tacPaud_MSpN_diff_FC_' num2str(ll) '_alpbet' '.eps'],'-depsc')
    print(ll+30,[fdir 'tfr_tacPaud_MSpN_diff_OP_' num2str(ll) '_theta' '.eps'],'-depsc')
    print(ll+30+100,[fdir 'tfr_tacPaud_MSpN_diff_OP_' num2str(ll) '_alpbet' '.eps'],'-depsc')
    
    if peakfindflag
      if ll==3
        [mx,mindD3]=max(nanmean(powfb{3,1}.avg(chpl{1},:),1));
        timeTacThetaTFR_AT70=powfb{3,1}.time(mindD3); % 300 ms
      elseif ll==7
        [mx,mindD3]=max(nanmean(powfb{3,1}.avg(chpl{1},:),1));
        timeTacThetaTFR_TA70=powfb{3,1}.time(mindD3);
      end
      
      if exist('timeTacThetaTFR_AT70') && exist('timeTacThetaTFR_TA70')
        timeTacThetaTFR_Diff=mean([timeTacThetaTFR_AT70 timeTacThetaTFR_TA70]);
%         timeTacThetaTFR_Uni=mean([timeTacThetaTfr timeAudThetaTfr]);
        for pp=1:3,
          peakthetaTFR_Diff(pp,ll)=nanmean(powfb{pp,1}.avg(chpl{1},dsearchn(powfb{pp,1}.time',timeTacThetaTFR_Diff)),1);
          peakthetaTFR_Diff_range(pp,ll)=nanmean(nanmean(powfb{pp,1}.avg(chpl{1},[dsearchn(powfb{pp,1}.time',timeTacThetaTFR_Diff-.02):dsearchn(powfb{pp,1}.time',timeTacThetaTFR_Diff+.02)]),1),2);
%           peakthetaTFR_Uni(pp,ll) =nanmean(powfb{pp,1}.avg(chpl{1},dsearchn(powfb{pp,1}.time',timeTacThetaTFR_Uni)),1);
        end
      end
      
      % posterior alpha/beta ERS rebound
      for pp=1:3,
        peakAlphaTFR_Diff(pp,ll)=nanmean(nanmean(powfb{pp,2}.avg(chpl{2},end-40:end-10),1),2); % 900-1200ms
        peakAlphaTFR_Diff_short(pp,ll)=nanmean(nanmean(powfb{pp,2}.avg(chpl{2},end-27:end-23),1),2); % 1030-1070ms
%         peakalphaTFR_Tac(pp,ll) =nanmean(powfb{pp,2}.avg(chpl{2},dsearchn(powfb{pp,2}.time',timeTacAlpbet)),1);
      end
      
    end %peak
    
    
  end % plottfrflag
  
  if plotplvflag
    clear pow*
    cfg=[];
    tmpA=ft_selectdata(cfg,gravelo_TPA_MSPN_comb1);
    tmpA.mask=logical(zeros(size(tmpA.powspctrm)));
    tmpA.mask(:,:,dsearchn(tmpA.time',stattl_mc_comb1plvabsDepT.time(1)):dsearchn(tmpA.time',stattl_mc_comb1plvabsDepT.time(end)))=stattl_mc_comb1plvabsDepT.mask;
    tmpsA=ft_selectdata(cfg,gravelo_tacMSpN_comb1);
    tmpsA.mask=logical(zeros(size(tmpsA.powspctrm)));
    tmpsA.mask(:,:,dsearchn(tmpsA.time',stattl_mc_comb1plvabsDepT.time(1)):dsearchn(tmpsA.time',stattl_mc_comb1plvabsDepT.time(end)))=stattl_mc_comb1plvabsDepT.mask;
    tmpuA=ft_selectdata(cfg,gravelo_tacPaud_comb1);
    tmpuA.mask=logical(zeros(size(tmpuA.powspctrm)));
    tmpuA.mask(:,:,dsearchn(tmpuA.time',stattl_mc_comb1plvabsDepT.time(1)):dsearchn(tmpuA.time',stattl_mc_comb1plvabsDepT.time(end)))=stattl_mc_comb1plvabsDepT.mask;
    
    pow{1}=tmpuA;
    pow{2}=tmpsA;
    pow{3}=tmpA;
    
    % average over bands first
    for pp=1:3
      for ff=1:2 % two diff freq bands
        cfg=[];
        if ff==1
          cfg.frequency   = [4 8.5];
        elseif ff==2
          cfg.frequency   = [10 20];
        end
        cfg.avgoverfreq = 'yes';
        cfg.nanmean     = 'yes';
        powf{pp,ff}=ft_selectdata(cfg,pow{pp});
      end
    end
    clear pow
    
    % baseline correct first, then remove baseline area for plotting
    for ff=1:2
      %         cfg=[];
      %         cfg.baseline=[gravelo_TPA_MSPN_comb1.time(49) gravelo_TPA_MSPN_comb1.time(53)];
      %         cfg.baselinetype='absolute';
      %         cfg.parameter={'plvabs' 'mask'};
      %         powfb{pp,ff}=ft_freqbaseline(cfg, powf{pp,ff});
      basep=powf{2,ff};
      basep.plvabs=repmat(nanmean(cat(3,nanmean(powf{1,ff}.plvabs(:,:,49:53),3),nanmean(powf{2,ff}.plvabs(:,:,49:53),3)),3),[1 1 length(powf{1,ff}.time)]);
      basep.mask=repmat(nanmean(powf{2,ff}.mask(:,:,49:53),3),[1 1 length(powf{1,ff}.time)]);  % use 'nul' as baseline for all three conditions for plotting
      cfg=[];
      cfg.parameter={'plvabs' 'mask'};
      cfg.operation='subtract';
      powfb{1,ff}=ft_math(cfg,powf{1,ff},basep)
      powfb{2,ff}=ft_math(cfg,powf{2,ff},basep)
      powfb{3,ff}=powf{3,ff};
      
      for pp=1:3
        cfg=[];
        if ll==1
          cfg.latency=[tacbasemax(ll) 1.3];
        elseif ll==3
          cfg.latency=[tacbasemax(ll) 1.3];
        elseif ll==4
          cfg.latency=[tacbasemax(ll) 1.3];
        elseif ll==5
          cfg.latency=[tacbasemax(ll) 1.3];
        elseif ll==6
          cfg.latency=[tacbasemax(ll) 1.3+.02];
        elseif ll==7
          cfg.latency=[tacbasemax(ll) 1.3+.07];
        elseif ll==9
          cfg.latency=[tacbasemax(ll) 1.3+.50];
        end
        powfb{pp,ff}=ft_selectdata(cfg,powfb{pp,ff});
        
        powfb{pp,ff}.avg=squeeze(powfb{pp,ff}.plvabs);
        powfb{pp,ff}.mask=squeeze(powfb{pp,ff}.mask);
        powfb{pp,ff}.dimord='chan_time';
        powfb{pp,ff}=rmfield(powfb{pp,ff},'freq');
        powfb{pp,ff}=rmfield(powfb{pp,ff},'plvabs');
        if plotallmask
          powfb{pp,ff}.mask=repmat(ceil(mean(powfb{pp,ff}.mask,1)),[size(powfb{pp,ff}.mask,1) 1]);
        end
      end
    end
    
    for cg=1:length(chanplot)
      for ff=1:2
        
        cfg=[];
        cfg.parameter='avg';
        cfg.layout='elec1010.lay';
        if ff==1
          cfg.ylim=[-.1 .4];  % or [1 2.5]?
        elseif ff==2
          cfg.ylim=[-.1 .4];
        end
        cfg.linewidth=3;
        cfg.xlim=[-.6 1.8];
        if cg>length(chanplot)
          cfg.channel=tmpd5.label(any(tmpd5.mask,2));
        else
          cfg.channel=chanplot{cg};
        end
        cfg.graphcolor=[colorblindApT; colorblindMpN; colorblindD];
        cfg.interactive='no';
        cfg.maskparameter='mask';
        cfg.maskstyle='box'; % default
        %     subplot(4,1,cg)
        %     subplot(1,2,cg)
        figure(ll+10*(cg+1)+(ff-1)*100)
        ft_singleplotER(cfg,powfb{1,ff},powfb{2,ff},powfb{3,ff});
        hold on;plot(powfb{1,ff}.time,0,'k');
        set(gca,'XTick',[-.6:.1:1.8])
        set(gca,'XTickLabel',{ '' '-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'  '' ' ' ''  ' ' '1.5' '' ' ' '' })
        set(gca,'FontSize',30)
        title([])
        plot([stattimwin(1) stattimwin(1)],cfg.ylim,'k--','LineWidth',6)
        plot([stattimwin(2) stattimwin(2)],cfg.ylim,'k--','LineWidth',6)
        plot([0 0],cfg.ylim,'Color',colorblindT,'LineWidth',6)
        plot([soades(ll) soades(ll)],cfg.ylim,'Color',colorblindA,'LineWidth',6)
        %         plot(cfg.xlim,[0 0],'Color',coloruse(10,:),'LineWidth',1)
        plot(cfg.xlim,[0 0],'Color','k')
        axis([cfg.xlim(1) cfg.xlim(2) cfg.ylim(1) cfg.ylim(2)])
        if cg==3
          legend('A+T','MS+N','Difference')
%           legend('Sum Unisensory','MultSens + Null','SumUnisens - MultsensNull')
        end
      end % ff
    end % cg
    print(ll+20,[fdir 'plv_tacPaud_MSpN_diff_FC_' num2str(ll) '_theta' '.png'],'-dpng')
    print(ll+20+100,[fdir 'plv_tacPaud_MSpN_diff_FC_' num2str(ll) '_alpbet' '.png'],'-dpng')
    print(ll+30,[fdir 'plv_tacPaud_MSpN_diff_OP_' num2str(ll) '_theta' '.png'],'-dpng')
    print(ll+30+100,[fdir 'plv_tacPaud_MSpN_diff_OP_' num2str(ll) '_alpbet' '.png'],'-dpng')
    print(ll+20,[fdir 'plv_tacPaud_MSpN_diff_FC_' num2str(ll) '_theta' '.eps'],'-depsc')
    print(ll+20+100,[fdir 'plv_tacPaud_MSpN_diff_FC_' num2str(ll) '_alpbet' '.eps'],'-depsc')
    print(ll+30,[fdir 'plv_tacPaud_MSpN_diff_OP_' num2str(ll) '_theta' '.eps'],'-depsc')
    print(ll+30+100,[fdir 'plv_tacPaud_MSpN_diff_OP_' num2str(ll) '_alpbet' '.eps'],'-depsc')
    
    
    if peakfindflag
      timeThetaPLV_Uni=mean([timeTacThetaPlv timeAudThetaPlv]); %135 ms
      for pp=1:3
        peakthetaPLV_Uni(pp,ll) =nanmean(powfb{pp,1}.avg(chpl{1},dsearchn(powfb{pp,1}.time',timeThetaPLV_Uni)),1);
        peakthetaPLV_Uni_range(pp,ll) =nanmean(nanmean(powfb{pp,1}.avg(chpl{1},[dsearchn(powfb{pp,1}.time',timeThetaPLV_Uni-.02):dsearchn(powfb{pp,1}.time',timeThetaPLV_Uni+.02)]),1),2);
      end
    end % peak
  
    
  end %plvplotflag
  
  
end
save tfrsave.mat tfrsave*
save('peakTFvalues.mat','peak*','time*');



%% Old



%% plot peak bars

load('peakTFvalues.mat');

data=peakthetaTFR_Diff(:,soalist)';
starind=[1 2 3 6];
yminmax=[-.2 1.3];
figind=56;
figname=[fdir 'tfr_condDiff_theta.eps'];
plotlegsummarybars(data,figind,figname,starind,yminmax)

data=peakthetaTFR_Diff_range(:,soalist)';
figind=59;
figname=[fdir 'tfr_condDiff_thetarange.eps'];
plotlegsummarybars(data,figind,figname,starind,yminmax)

data=peakAlphaTFR_Diff(:,soalist)';
starind=[2 3 4 6];
yminmax=[-.7 2.3];
figind=57;
figname=[fdir 'tfr_condDiff_alpha.eps'];
plotlegsummarybars(data,figind,figname,starind,yminmax)

data=peakAlphaTFR_Diff_short(:,soalist)';
figind=60;
yminmax=[-1.2 2.4];
figname=[fdir 'tfr_condDiff_alphashort.eps'];
plotlegsummarybars(data,figind,figname,starind,yminmax)

data=peakthetaPLV_Uni(:,soalist)';
starind=[2 6];
yminmax=[-.07 .37];
figind=58;
figname=[fdir 'plv_condDiff_theta.eps'];
plotlegsummarybars(data,figind,figname,starind,yminmax)

data=peakthetaPLV_Uni_range(:,soalist)';
figind=61;
figname=[fdir 'plv_condDiff_thetarange.eps'];
plotlegsummarybars(data,figind,figname,starind,yminmax)


% figure(56); % preferred
% subplot(2,1,1);ph=bar(peakthetaTFR_Diff(1:2,soalist)');ylim([-0.2 1.3])
% set(ph(1),'FaceColor',coloruse(1,:))
% set(ph(2),'FaceColor',coloruse(10,:))
% set(gca,'XTickLabel',{ '' '' '' '' '' '' ''   })
% set(gca,'XTickLabelRotation',0)
% set(gca,'FontSize',18)
% title('TFP Theta, FC channels')
% subplot(2,1,2);ph=bar(peakthetaTFR_Diff(3,soalist),0.3);ylim([-0.2 1.3])
% set(ph(1),'FaceColor',coloruse(5,:))
% set(gca,'XTickLabel',{ 'AT500' 'AT70' 'AT20' 'AT0' 'TA20' 'TA70' 'TA500'   })
% set(gca,'XTickLabelRotation',45)
% set(gca,'FontSize',18)
% print(56,[fdir 'tfr_theta_peak.eps'],'-depsc')
% % figure(57);
% % subplot(2,1,1);ph=bar(peakthetaTFR_Uni(1:2,soalist)');ylim([-0.2 2.1])
% % set(ph(1),'FaceColor',coloruse(1,:))
% % set(ph(2),'FaceColor',coloruse(10,:))
% % subplot(2,1,2);ph=bar(peakthetaTFR_Uni(3,soalist),0.3);ylim([-0.2 2.1])
% % set(ph(1),'FaceColor',coloruse(5,:))
% 
% 
% figure(58);  % preferred
% subplot(2,1,1);ph=bar(peakAlphaTFR_Diff(1:2,soalist)');ylim([-0.7 2.3])
% set(ph(1),'FaceColor',coloruse(1,:))
% set(ph(2),'FaceColor',coloruse(10,:))
% set(gca,'XTickLabel',{ '' '' '' '' '' '' ''   })
% set(gca,'XTickLabelRotation',0)
% set(gca,'FontSize',18)
% title('TFP Alpha/Beta, Par. channels')
% subplot(2,1,2);ph=bar(peakAlphaTFR_Diff(3,soalist),0.3);ylim([-0.7 2.3])
% set(ph(1),'FaceColor',coloruse(5,:))
% set(gca,'XTickLabel',{ 'AT500' 'AT70' 'AT20' 'AT0' 'TA20' 'TA70' 'TA500'   })
% set(gca,'XTickLabelRotation',45)
% set(gca,'FontSize',18)
% print(58,[fdir 'tfr_alpha_peak.eps'],'-depsc')
% % figure(59);
% % subplot(2,1,1);ph=bar(peakalphaTFR_Tac(1:2,soalist)');ylim([-0.7 2.1])
% % set(ph(1),'FaceColor',coloruse(1,:))
% % set(ph(2),'FaceColor',coloruse(10,:))
% % subplot(2,1,2);ph=bar(peakalphaTFR_Tac(3,soalist),0.3);ylim([-0.7 2.1])
% % set(ph(1),'FaceColor',coloruse(5,:))
% 
% figure(60);
% subplot(2,1,1);ph=bar(peakthetaPLV_Uni(1:2,soalist)');ylim([-0.07 0.37])
% set(ph(1),'FaceColor',coloruse(1,:))
% set(ph(2),'FaceColor',coloruse(10,:))
% set(gca,'XTickLabel',{ '' '' '' '' '' '' ''   })
% set(gca,'XTickLabelRotation',0)
% set(gca,'FontSize',18)
% title('PLF Theta, FC channels')
% subplot(2,1,2);ph=bar(peakthetaPLV_Uni(3,soalist),0.3);ylim([-0.07 0.37])
% set(ph(1),'FaceColor',coloruse(5,:))
% set(gca,'XTickLabel',{ 'AT500' 'AT70' 'AT20' 'AT0' 'TA20' 'TA70' 'TA500'   })
% set(gca,'XTickLabelRotation',45)
% set(gca,'FontSize',18)
% print(60,[fdir 'plv_theta_peak.eps'],'-depsc')
% 
% % line lots to match ERP now instead
% try close(57); end;
% figure(57);
% x=[-.6 -.35 -.1 0 .1 .35 .6];
% [ph,h1,h2]=plotyy(x,peakthetaTFR_Diff(1,soalist),x,peakthetaTFR_Diff(3,soalist));
% line(x,peakthetaTFR_Diff(2,soalist),'Parent',ph(1));
% line([-.8 .8], [0 0], 'Parent', ph(2));
% set(h2,'Color',coloruse([5 ],:),'Marker','*','MarkerSize',15,'LineWidth',3,'LineStyle',':');
% set(ph(2).Children(1),'Color',coloruse(5,:))
% set(ph(1).Children(1),'Color',coloruse(10,:),'Marker','^','MarkerSize',15,'LineWidth',3,'LineStyle',':')
% set(ph(1).Children(2),'Color',coloruse(1,:),'Marker','v','MarkerSize',15,'LineWidth',3,'LineStyle',':')
% set(ph(1),'YColor',coloruse([10 ],:));
% set(ph(2),'YColor',coloruse([5 ],:));
% set(ph,'FontSize',15)
% set(ph,'LineWidth',3)
% xlim(ph(2), [-.8 .8])
% xlim(ph(1), [-.8 .8])
% ylim(ph(1), [-.2 1.3])
% ylim(ph(2), [-.2 1.3])
% set(gca,'Xtick',x)
% set(gca,'XTickLabel',{ '-500' '-70' '-20' '0' '20' '70' '500'})
% xlabel('Multisensory Asynchrony (ms)')
% legend('A+T','MS+N','Difference')
% ylabel('Theta Power (uV^2)')
% print(57,[fdir 'TFP_condDiff_theta_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')
% 
% 
% try close(58); end;
% figure(58);
% x=[-.6 -.35 -.1 0 .1 .35 .6];
% [ph,h1,h2]=plotyy(x,peakAlphaTFR_Diff(1,soalist),x,peakAlphaTFR_Diff(3,soalist));
% line(x,peakAlphaTFR_Diff(2,soalist),'Parent',ph(1));
% line([-.8 .8], [0 0], 'Parent', ph(2));
% set(h2,'Color',coloruse([5 ],:),'Marker','*','MarkerSize',15,'LineWidth',3,'LineStyle',':');
% set(ph(2).Children(1),'Color',coloruse(5,:))
% set(ph(1).Children(1),'Color',coloruse(10,:),'Marker','^','MarkerSize',15,'LineWidth',3,'LineStyle',':')
% set(ph(1).Children(2),'Color',coloruse(1,:),'Marker','v','MarkerSize',15,'LineWidth',3,'LineStyle',':')
% set(ph(1),'YColor',coloruse([10 ],:));
% set(ph(2),'YColor',coloruse([5 ],:));
% set(ph,'FontSize',15)
% set(ph,'LineWidth',3)
% xlim(ph(2), [-.8 .8])
% xlim(ph(1), [-.8 .8])
% ylim(ph(1), [-.7 2.3])
% ylim(ph(2), [-.7 2.3])
% set(gca,'Xtick',x)
% set(gca,'XTickLabel',{ '-500' '-70' '-20' '0' '20' '70' '500'})
% xlabel('Multisensory Asynchrony (ms)')
% legend('A+T','MS+N','Difference')
% ylabel('Alpha/Beta Power (uV^2)')
% print(58,[fdir 'TFP_condDiff_alpha_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')
% 
% 
% try close(59); end;
% figure(59);
% x=[-.6 -.35 -.1 0 .1 .35 .6];
% [ph,h1,h2]=plotyy(x,peakthetaPLV_Uni(1,soalist),x,peakthetaPLV_Uni(3,soalist));
% line(x,peakthetaPLV_Uni(2,soalist),'Parent',ph(1));
% line([-.8 .8], [0 0], 'Parent', ph(2));
% set(h2,'Color',coloruse([5 ],:),'Marker','*','MarkerSize',15,'LineWidth',3,'LineStyle',':');
% set(ph(2).Children(1),'Color',coloruse(5,:))
% set(ph(1).Children(1),'Color',coloruse(10,:),'Marker','^','MarkerSize',15,'LineWidth',3,'LineStyle',':')
% set(ph(1).Children(2),'Color',coloruse(1,:),'Marker','v','MarkerSize',15,'LineWidth',3,'LineStyle',':')
% set(ph(1),'YColor',coloruse([10 ],:));
% set(ph(2),'YColor',coloruse([5 ],:));
% set(ph,'FontSize',15)
% set(ph,'LineWidth',3)
% xlim(ph(2), [-.8 .8])
% xlim(ph(1), [-.8 .8])
% ylim(ph(1), [-.07 .37])
% ylim(ph(2), [-.07 .37])
% set(gca,'Xtick',x)
% set(gca,'XTickLabel',{ '-500' '-70' '-20' '0' '20' '70' '500'})
% xlabel('Multisensory Asynchrony (ms)')
% legend('A+T','MS+N','Difference')
% ylabel('Phase locking factor')
% print(59,[fdir 'PLV_condDiff_theta_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')

for ll=3
  close all
  clear grave* stat* tmp*
  load([edir 'statsgrave_TFR_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '_usetr' num2str(usetr) '_mcseed' num2str(mcseed) '_itc' '.mat']);
  %   baseline1=[tacbasemax(ll) tacbasemax(ll)+.08];
  %   stattimwin=[stattl_mc_comb1.time(1) stattl_mc_comb1.time(end)];
  
  
  cfg=[];
  cfg.avgoverfreq='yes';
  cfg.frequency=[4 6.5];
  grdiff_theta=ft_selectdata(cfg,gravelo_TPA_MSPN_comb1);
  cfg.frequency=[8 20];
  grdiff_alpha=ft_selectdata(cfg,gravelo_TPA_MSPN_comb1);
  
  
  cfg=[];
  cfg.layout='elec1010.lay';
  cfg.highlight          = 'on';
  cfg.highlightsymbol    = '*';
  cfg.highlightsize      = 10;
  cfg.comment='no';
  
  cfg.parameter='powspctrm';
%   cfg.xlim=[timeTacThetaTFR_Diff-.1 timeTacThetaTFR_Diff+.1];
  cfg.xlim=[timeTacThetaTFR_Diff-.02 timeTacThetaTFR_Diff+.02];
  cfg.zlim=[-1.1 1.1];
  cfg.highlightchannel   =  chanplot{1};
  figure(120+ll);
  ft_topoplotTFR(cfg,grdiff_theta);
  print(120+ll,[fdir 'TFP_topoDiff_Theta_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')
  
  cfg.parameter='powspctrm';
%   cfg.xlim=gravelo_TPA_MSPN_comb1.time([end-40 end-10]);
  cfg.xlim=gravelo_TPA_MSPN_comb1.time([end-27 end-23]);
  cfg.zlim=[-2.5 2.5];
  cfg.highlightchannel   =  chanplot{2};
  figure(130+ll);
  ft_topoplotTFR(cfg,grdiff_alpha);
  print(130+ll,[fdir 'TFP_topoDiff_Alpha_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')

  cfg.parameter='plvabs';
%   cfg.xlim=[timeThetaPLV_Uni-.03 timeThetaPLV_Uni+.03];
  cfg.xlim=[timeThetaPLV_Uni-.02 timeThetaPLV_Uni+.02];
  cfg.zlim=[-.2 .2];
  cfg.highlightchannel   =  chanplot{1};
  figure(140+ll);
  ft_topoplotTFR(cfg,grdiff_theta);
  print(140+ll,[fdir 'PLV_topoDiff_Theta_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-painters','-depsc')
  
end



%% Plotting time courses of TFP / PLF

chanplot{1}={'Fz' 'Cz' 'F1' 'F2' 'FC1' 'FC2' 'C1' 'C2'}; % frontocentral
chanplot{2}={'CP5' 'POz' 'Pz' 'P3' 'P4' 'C4' 'O1' 'O2' 'P7' 'PO7'}; % occipital; % same as for ERP
% coloruse=varycolor(10);
soalist=[1 3 4 5 6 7 9];
colorblindD  =[204 101 0]/256;
colorblindApT=[5 165 255]/256;
colorblindMpN=[0 59 179]/256;
colorblindT  =[255 51 166]/256;
colorblindA  =[0 158 115]/256;
colorblindM  =[0 0 0]/256;
colorblindN  =[128 128 128]/256;

sleep=1;
trialkc=-1;

llind=0;
for ll=soalist
  llind=llind+1;
  if sleep==1
    load(['statsgrave_TFR_cond' num2str(ll) '3121_iter11_trialkc' num2str(trialkc) '_usetr1_mcseed13_itc.mat'])
  elseif sleep==0
    load(['statsgrave_TFR_cond' num2str(ll) '3100_iter31_trialkc-1_usetr3_mcseed13_itc.mat'])
  end
  
  chpl{1}=match_str(gravelo_TPA_MSPN_comb1.label,chanplot{1});
  chpl{2}=match_str(gravelo_TPA_MSPN_comb1.label,chanplot{2});
  timeuse=gravelo_TPA_MSPN_comb1.time(21:end-20);
  %   timwin=[gravelo_TPA_MSPN_comb1.time(21) gravelo_TPA_MSPN_comb1.time(end-20)];
  stattimese=dsearchn(stattl_mc_comb1.time',[timeuse(1) timeuse(end)]');
  stattimeuse=stattl_mc_comb1.time(stattimese(1):stattimese(2));
  
  for cp=1:2
    figure(1)
    frequse=1:2;
    base1=squeeze(nanmean(nanmean(gravelo_tacPaud_comb1.powspctrm(chpl{cp},frequse,51),2),1));
    base2=squeeze(nanmean(nanmean(gravelo_tacMSpN_comb1.powspctrm(chpl{cp},frequse,51),2),1));
    base=mean([base1 base2]);
    %     subplot(2,7,(cp-1)*7+llind);
    subplot(7,2,(llind-1)*2+cp);
    pk=plot(timeuse,squeeze(nanmean(nanmean(gravelo_TPA_MSPN_comb1.powspctrm(chpl{cp},frequse,21:end-20),2),1)),'m')
    pk.Color=colorblindD;
    hold on;pk=plot(timeuse,squeeze(nanmean(nanmean(gravelo_tacPaud_comb1.powspctrm(chpl{cp},frequse,21:end-20),2),1))-base,'k')
    pk.Color=colorblindApT;
    hold on;pk=plot(timeuse,squeeze(nanmean(nanmean(gravelo_tacMSpN_comb1.powspctrm(chpl{cp},frequse,21:end-20),2),1))-base,'g')
    pk.Color=colorblindMpN;
    hold on;plot(stattimeuse,squeeze(nanmean(nanmean(stattl_mc_comb1.stat(chpl{cp},frequse,stattimese(1):stattimese(2)),2),1)),'b')
    ylim([-3.5 3.5])
    xlim([-0.6 1.8]);
    
    figure(2)
    frequse=3:6;
    base1=squeeze(nanmean(nanmean(gravelo_tacPaud_comb1.powspctrm(chpl{cp},frequse,51),2),1));
    base2=squeeze(nanmean(nanmean(gravelo_tacMSpN_comb1.powspctrm(chpl{cp},frequse,51),2),1));
    base=mean([base1 base2]);
    %     subplot(2,7,(cp-1)*7+llind);
    subplot(7,2,(llind-1)*2+cp);
    pk=plot(timeuse,squeeze(nanmean(nanmean(gravelo_TPA_MSPN_comb1.powspctrm(chpl{cp},frequse,21:end-20),2),1)),'m')
    pk.Color=colorblindD;
    hold on;pk=plot(timeuse,squeeze(nanmean(nanmean(gravelo_tacPaud_comb1.powspctrm(chpl{cp},frequse,21:end-20),2),1))-base,'k')
    pk.Color=colorblindApT;
    hold on;pk=plot(timeuse,squeeze(nanmean(nanmean(gravelo_tacMSpN_comb1.powspctrm(chpl{cp},frequse,21:end-20),2),1))-base,'g')
    pk.Color=colorblindMpN;
    hold on;plot(stattimeuse,squeeze(nanmean(nanmean(stattl_mc_comb1.stat(chpl{cp},frequse,stattimese(1):stattimese(2)),2),1)),'b')
    ylim([-3 3])
    xlim([-0.6 1.8]);
    
    figure(3)
    frequse=6:14;
    base1=squeeze(nanmean(nanmean(gravelo_tacPaud_comb1.powspctrm(chpl{cp},frequse,51),2),1));
    base2=squeeze(nanmean(nanmean(gravelo_tacMSpN_comb1.powspctrm(chpl{cp},frequse,51),2),1));
    base=mean([base1 base2]);
    %     subplot(2,7,(cp-1)*7+llind);
    subplot(7,2,(llind-1)*2+cp);
    pk=plot(timeuse,squeeze(nanmean(nanmean(gravelo_TPA_MSPN_comb1.powspctrm(chpl{cp},frequse,21:end-20),2),1)),'m')
    pk.Color=colorblindD;
    hold on;pk=plot(timeuse,squeeze(nanmean(nanmean(gravelo_tacPaud_comb1.powspctrm(chpl{cp},frequse,21:end-20),2),1))-base,'k')
    pk.Color=colorblindApT;
    hold on;pk=plot(timeuse,squeeze(nanmean(nanmean(gravelo_tacMSpN_comb1.powspctrm(chpl{cp},frequse,21:end-20),2),1))-base,'g')
    pk.Color=colorblindMpN;
    hold on;plot(stattimeuse,squeeze(nanmean(nanmean(stattl_mc_comb1.stat(chpl{cp},frequse,stattimese(1):stattimese(2)),2),1)),'b')
    ylim([-0.7 0.7])
    xlim([-0.6 1.8]);
    
    figure(4)
    frequse=3:9;
    base1=squeeze(nanmean(nanmean(gravelo_tacPaud_comb1.powspctrm(chpl{cp},frequse,51),2),1));
    base2=squeeze(nanmean(nanmean(gravelo_tacMSpN_comb1.powspctrm(chpl{cp},frequse,51),2),1));
    base=mean([base1 base2]);
    %     subplot(2,7,(cp-1)*7+llind);
    subplot(7,2,(llind-1)*2+cp);
    pk=plot(timeuse,squeeze(nanmean(nanmean(gravelo_TPA_MSPN_comb1.powspctrm(chpl{cp},frequse,21:end-20),2),1)),'m')
    pk.Color=colorblindD;
    hold on;pk=plot(timeuse,squeeze(nanmean(nanmean(gravelo_tacPaud_comb1.powspctrm(chpl{cp},frequse,21:end-20),2),1))-base,'k')
    pk.Color=colorblindApT;
    hold on;pk=plot(timeuse,squeeze(nanmean(nanmean(gravelo_tacMSpN_comb1.powspctrm(chpl{cp},frequse,21:end-20),2),1))-base,'g')
    pk.Color=colorblindMpN;
    hold on;plot(stattimeuse,squeeze(nanmean(nanmean(stattl_mc_comb1.stat(chpl{cp},frequse,stattimese(1):stattimese(2)),2),1)),'b')
    ylim([-2.5 2.5])
    xlim([-0.6 1.8]);
    
  end
  
  %   print(ll,[fdir 'TFR_tacPaud_MSpN_diff_' num2str(ll) '.png'],'-dpng')
  print(1,[fdir 'TFR_tacPaud_MSpN_diff_theta.png'],'-dpng')
  print(2,[fdir 'TFR_tacPaud_MSpN_diff_alpha.png'],'-dpng')
  print(3,[fdir 'TFR_tacPaud_MSpN_diff_beta.png'],'-dpng')
  print(4,[fdir 'TFR_tacPaud_MSpN_diff_alfbet.png'],'-dpng')
  print(1,[fdir 'TFR_tacPaud_MSpN_diff_theta.eps'],'-depsc')
  print(2,[fdir 'TFR_tacPaud_MSpN_diff_alpha.eps'],'-depsc')
  print(3,[fdir 'TFR_tacPaud_MSpN_diff_beta.eps'],'-depsc')
  print(4,[fdir 'TFR_tacPaud_MSpN_diff_alfbet.eps'],'-depsc')
  
end
