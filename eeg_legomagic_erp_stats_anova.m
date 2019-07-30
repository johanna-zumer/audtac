eeg_legomagic_preamble

%% ANOVA F-test (initially added/begun in November 2018 after JoN comments)

sleep=0;

soalist=[1 3 4 5 6 7 9];
soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];
if sleep
  ss=12;
else
  ss=10;
end

runstatflag=0;

%%
if sleep
  load tlock_grind_sleep1_ss12_iter11_trialkc-1_statwinorig0_ftver0.mat;
  iicnt=0;
  for ii=iiBuse
    iicnt=iicnt+1;
    load([sub{ii} '/trlfeat_' sub{ii} '_iter11.mat']);
    subind_nKD(iicnt,:)=min_nKD>20;
    subind_nSD(iicnt,:)=min_nSD>20;
    subind_nKD_nSD(iicnt,:)=min_nKD_nSD>20;
  end
  clear min* trlkeep* featind*
  keyboard  % this easist with manual editing;  recheck next few lines if tlockdiff  has been re-run
  for ll=setdiff(soalist,3)
    grind_tacPaud_nSD_save{ll,3,ss}.individual=grind_tacPaud_nSD_save{ll,3,ss}.individual(subind_nSD(:,3),:,:);
    grind_tacMSpN_nSD_save{ll,3,ss}.individual=grind_tacMSpN_nSD_save{ll,3,ss}.individual(subind_nSD(:,3),:,:);
  end
else
  if ~exist('grind_tacPaud_save','var')
    load tlock_grind_sleep0_iter27_statwinorig0_ftver0.mat;
  end
end
cfg=[];
cfg.operation='subtract';
cfg.parameter='individual';
llind=1;
clear grind_TPA_MSPN
for ll=soalist
  [grind_TPA_MSPN{llind},grind_tacPaud{llind},grind_tacMSpN{llind}]=grind_anova_helper(cfg,grind_tacPaud_save{ll,3,ss},grind_tacMSpN_save{ll,3,ss},ll);
  if sleep
    [grind_TPA_MSPN_nKD{llind},grind_tacPaud_nKD{llind},grind_tacMSpN_nKD{llind}]=grind_anova_helper(cfg,grind_tacPaud_nKD_save{ll,3,ss},grind_tacMSpN_nKD_save{ll,3,ss},ll);
    [grind_TPA_MSPN_nSD{llind},grind_tacPaud_nSD{llind},grind_tacMSpN_nSD{llind}]=grind_anova_helper(cfg,grind_tacPaud_nSD_save{ll,3,ss},grind_tacMSpN_nSD_save{ll,3,ss},ll);
    [grind_TPA_MSPN_nKD_nSD{llind},grind_tacPaud_nKD_nSD{llind},grind_tacMSpN_nKD_nSD{llind}]=grind_anova_helper(cfg,grind_tacPaud_nKD_nSD_save{ll,3,ss},grind_tacMSpN_nKD_nSD_save{ll,3,ss},ll);
  end
  %   grind_TPA_MSPN{llind}=ft_math(cfg,grind_tacPaud_save{ll,3,ss},grind_tacMSpN_save{ll,3,ss});
  %   grind_tacPaud{llind}=grind_tacPaud_save{ll,3,ss};
  %   grind_tacMSpN{llind}=grind_tacMSpN_save{ll,3,ss};
  %   if ll>5
  %     grind_TPA_MSPN{llind}.time=grind_TPA_MSPN{llind}.time-soades(ll);
  %     grind_tacPaud{llind}.time=grind_tacPaud{llind}.time-soades(ll);
  %     grind_tacMSpN{llind}.time=grind_tacMSpN{llind}.time-soades(ll);
  %   end
  llind=llind+1;
end

%%
nsub=size(grind_TPA_MSPN{1}.individual,1);
if sleep
  nsub_nKD=size(grind_TPA_MSPN_nKD{1}.individual,1);
  nsub_nSD=size(grind_TPA_MSPN_nSD{1}.individual,1);
  nsub_nKD_nSD=size(grind_TPA_MSPN_nKD_nSD{1}.individual,1);
end
load eeg1010_neighb


if 0
  cfg=[];
  cfg.operation='add';
  cfg.parameter='individual';
  grind_TPA_MSPN_allsummed=ft_math(cfg,grind_TPA_MSPN{:});
  grind_tacPaud_allsummed=ft_math(cfg,grind_tacPaud{:});
  grind_tacMSpN_allsummed=ft_math(cfg,grind_tacMSpN{:});
  cfg=[];
  cfg.latency=[0 0.5];
  cfg.method='montecarlo';
  cfg.parameter='individual';
  cfg.neighbours=neighbours;
  cfg.correctm='cluster';
  cfg.numrandomization=2000;
  cfg.tail=0; % 0 for two-sided t-test
  cfg.ivar=1;
  cfg.uvar=2;
  cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
  cfg.design(2,:)=[1:nsub 1:nsub];
  cfg.randomseed=13;
  cfg.correcttail='alpha';
  cfg.clusterstatistic='maxsum';
  cfg.statistic='depsamplesT';
  cfg.computecritval='yes';
  cfg.computeprob='yes';
  stat_TPA_MSPN_allsummed = ft_timelockstatistics(cfg,grind_tacPaud_allsummed,grind_tacMSpN_allsummed);
  % Nothing significant  (p>0.62)
end

if runstatflag
  cfg=[];
  cfg.latency=[0 0.5];
  cfg.method='montecarlo';
  cfg.parameter='individual';
  cfg.neighbours=neighbours;
  cfg.correctm='cluster';
  cfg.numrandomization=2000;
  cfg.tail=1; % only 1 makes sense for F-test
  cfg.ivar=1;
  cfg.uvar=2;
  cfg.design = set_cfg_design_depF(nsub);
  % cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub) 3*ones(1,nsub) 4*ones(1,nsub) 5*ones(1,nsub) 6*ones(1,nsub) 7*ones(1,nsub)];
  % cfg.design(2,:)=[1:nsub 1:nsub 1:nsub 1:nsub 1:nsub 1:nsub 1:nsub];
  cfg.randomseed=13;
  cfg.clusterstatistic='maxsum';
  cfg.statistic='depsamplesFunivariate';
  cfg.computecritval='yes';
  cfg.computeprob='yes';
  stat_TPA_MSPN_1wayANOVA = ft_timelockstatistics(cfg,grind_TPA_MSPN{:});
  if sleep
    cfg.latency=[0 0.5];
    cfg.design = set_cfg_design_depF(nsub_nKD);
    stat_TPA_MSPN_nKD_1wayANOVA = ft_timelockstatistics(cfg,grind_TPA_MSPN_nKD{:});
    cfg.design = set_cfg_design_depF(nsub_nSD);
    stat_TPA_MSPN_nSD_1wayANOVA = ft_timelockstatistics(cfg,grind_TPA_MSPN_nSD{:});
    cfg.design = set_cfg_design_depF(nsub_nKD_nSD);
    stat_TPA_MSPN_nKD_nSD_1wayANOVA = ft_timelockstatistics(cfg,grind_TPA_MSPN_nKD_nSD{:});
    cfg.design = set_cfg_design_depF(nsub);
    cfg.latency=[0.1 0.8];
    stat_TPA_MSPN_1wayANOVA_long = ft_timelockstatistics(cfg,grind_TPA_MSPN{:});
    cfg.design = set_cfg_design_depF(nsub_nKD);
    stat_TPA_MSPN_nKD_1wayANOVA_long = ft_timelockstatistics(cfg,grind_TPA_MSPN_nKD{:});
    cfg.design = set_cfg_design_depF(nsub_nSD);
    stat_TPA_MSPN_nSD_1wayANOVA_long = ft_timelockstatistics(cfg,grind_TPA_MSPN_nSD{:});
    cfg.design = set_cfg_design_depF(nsub_nKD_nSD);
    stat_TPA_MSPN_nKD_nSD_1wayANOVA_long = ft_timelockstatistics(cfg,grind_TPA_MSPN_nKD_nSD{:});
    cfg.design = set_cfg_design_depF(nsub);
  end
  save([edir 'ANOVA_sleep' num2str(sleep) 'ERP.mat'],'stat*')
  
else
  load([edir 'ANOVA_sleep' num2str(sleep) '_ERP.mat'])
end

if sleep
  keyboard
else
  posthoctimwin=[min(stat_TPA_MSPN_1wayANOVA.time(find(ceil(mean(stat_TPA_MSPN_1wayANOVA.mask,1))))) max(stat_TPA_MSPN_1wayANOVA.time(find(ceil(mean(stat_TPA_MSPN_1wayANOVA.mask,1)))))];
  clustermat=stat_TPA_MSPN_1wayANOVA.posclusterslabelmat;
  clustermat(clustermat>=3)=0;
  
  reltimepoints=find(mean(clustermat)>.2);
  diffpoints=find(diff(reltimepoints)>1)
  erp_anova_tw1=reltimepoints(1:89);
  erp_anova_tw2=reltimepoints(91:185);
  erp_anova_tw3=reltimepoints(189:239);
  erp_anova_tw4=reltimepoints(240:end);
end

%%
% For awake paper:
figure(42);imagesc(0:.001:.5,1:62,clustermat);
xlim([-.01 .51]);
set(gca,'FontSize',30);
set(gca,'XTick',[-.6:.1:1.8])
set(gca,'XTickLabel',{ '' '-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'  '' ' ' ''  ' ' '1.5' '' ' ' '' })
set(gca,'YTick',1:62)
set(gca,'YTickLabel',stat_TPA_MSPN_1wayANOVA.label)
set(gca,'YTickMode','auto')
print(42,[fdir 'ERPstatANOVA.eps'],'-painters','-depsc')

cfg=[];
cfg.latency=[0 0.5];
for llind=1:7
  grind_TPA_MSPN_sel{llind}=ft_selectdata(cfg,grind_TPA_MSPN{llind});
end
cfg=[];
cfg.latency=[0 0.5];
cfg.avgoverrpt='yes';
for llind=1:7
  grave_TPA_MSPN_avg{llind}=ft_selectdata(cfg,grind_TPA_MSPN{llind});
end
% cfg.latency=posthoctimwin;
% for llind=1:7
%   grave_TPA_MSPN_timwin{llind}=ft_selectdata(cfg,grind_TPA_MSPN{llind});
% end

for llind=1:7
  fullconmat(llind,:,:)=grave_TPA_MSPN_avg{llind}.individual;
  %   fullconmat_tw(llind,:,:)=grave_TPA_MSPN_timwin{llind}.individual;
  fullconmat_grind(llind,:,:,:)=grind_TPA_MSPN_sel{llind}.individual;
end
fullconmat_reshape=reshape(fullconmat,[7 62*501]);
% fullconmat_tw_reshape=reshape(fullconmat_tw,[7 62*353]);

fullconmat_tmp=reshape(fullconmat_grind,[7 22 62*501]);
fullconmat_grind_reshape=reshape(fullconmat_tmp,[7 22*62*501]);

% cluster 1 (early)
mask_use=clustermat; % contains 2nd cluster at p=0.065
mask_use(mask_use==1)=0;
mask_use(mask_use==2)=1;
name='ERP1';
[aa1,aaa1,vvv1,bb1]=pca_masked(mask_use,fullconmat_reshape,grave_TPA_MSPN_avg{1}.label,1,[0:.001:.5],[-.1 .6],fdir,name,[1 0]);
bb1(1,1)/sum(bb1(:))

% Cluster 2
mask_use=stat_TPA_MSPN_1wayANOVA.mask;
mask_use(:,1:147)=0;
mask_use(:,301:end)=0;
name='ERP2';
[aa1,aaa1,vvv1,bb1]=pca_masked(mask_use,fullconmat_reshape,grave_TPA_MSPN_avg{1}.label,1,[0:.001:.5],[-.1 .6],fdir,name,[1 0]);
bb1(1,1)/sum(bb1(:))

% Cluster 3
mask_use=stat_TPA_MSPN_1wayANOVA.mask;
mask_use(:,1:350)=0;
mask_use(:,422:end)=0;
name='ERP3';
[aa1,aaa1,vvv1,bb1]=pca_masked(mask_use,fullconmat_reshape,grave_TPA_MSPN_avg{1}.label,1,[0:.001:.5],[-.1 .6],fdir,name,[0 1]);
bb1(1,1)/sum(bb1(:))

% Cluster 4
mask_use=stat_TPA_MSPN_1wayANOVA.mask;
mask_use(:,1:450)=0;
name='ERP4';
[aa1,aaa1,vvv1,bb1]=pca_masked(mask_use,fullconmat_reshape,grave_TPA_MSPN_avg{1}.label,1,[0:.001:.5],[-.1 .6],fdir,name,[0 1]);
bb1(1,1)/sum(bb1(:))

% Cluster 34
mask_use=stat_TPA_MSPN_1wayANOVA.mask;
mask_use(:,1:350)=0;
name='ERP34';
[aa1,aaa1,vvv1,bb1]=pca_masked(mask_use,fullconmat_reshape,grave_TPA_MSPN_avg{1}.label,1,[0:.001:.5],[-.1 .6],fdir,name);
bb1(1,1)/sum(bb1(:))


% % % %%%%%%%%%%%%%%%%
% PCA
% [aa1,bb1,vv1]=svd(fullconmat_reshape,'econ');
% [aa2,bb2,vv2]=svd(fullconmat_tw_reshape,'econ');
% [aa3,bb3,vv3]=svd(fullconmat_grind_reshape,'econ');
% [aa4,bb4,vv4]=svd(fullconmat_pcared_reshape,'econ');

% figure;for llind=1:7,subplot(4,2,llind);bar(aa4(:,llind));end
%
% [V,rotdata] = varimax(aa3);
%
%
% COEFF1 = pca(fullconmat_reshape);
% COEFF2 = pca(fullconmat_grind_reshape);

clear data4ica
data4ica.trial{1}=fullconmat_reshape;
% data4ica.trial{1}=fullconmat_pcared_ind_reshape;
data4ica.dimord='chan_time';
data4ica.time{1}=.001:.001:31.062;
% data4ica.time{1}=.001:.001:683.364;
data4ica.label={'1' '2' '3' '4' '5' '6' '7'};
cfg=[];
cfg.randomseed=13;
cfg.method='runica';
runicaout=ft_componentanalysis(cfg, data4ica);
cfg.method='fastica';
fasticaout_ind=ft_componentanalysis(cfg, data4ica);
% data4ica.trial{1}=fullconmat_pcared_avg_reshape;
% fasticaout_avg=ft_componentanalysis(cfg, data4ica);


runicaout.topo(:,1)=-runicaout.topo(:,1);
runicaout.trial{1}(1,:)=-runicaout.trial{1}(1,:);

for llind=1:7
  rica(:,:,llind)=reshape(runicaout.trial{1}(llind,:),[62 501]);
  %   pica(:,:,llind)=reshape(vv1(:,llind)',[62 501]);
  
  [af,bf,cf]=svd(rica(:,:,llind),'econ');
  if llind==4 || llind==1
    af=-af;cf=-cf;
  end
  statplot.avg=af(:,1);
  statplot.time=1;
  statplot.dimord='chan_time';
  statplot.label=grave_TPA_MSPN_avg{llind}.label;
  if 0 % for ICA alone figure
    figure(20);subplot(7,3,llind*3-2);bar(runicaout.topo(:,llind));ylim([-3 3])
    xticklabels({'AT500' 'AT70' 'AT20' 'AT0' 'TA20' 'TA70' 'TA500'})
    set(gca,'XTickLabelRotation',45)
    %   figure(20);subplot(7,3,llind*3-2);barh(runicaout.topo(:,llind));xlim([-3 3])
    figure(20);subplot(7,3,llind*3-1);cfg=[];cfg.xlim=1;cfg.zlim='maxabs';cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
    figure(20);subplot(7,3,llind*3);plot(cf(:,1));ylim([-max(abs(cf(:,1)))-.02 max(abs(cf(:,1)))+.02])
  else % For ICA-ERP figure in paper
    close all
    if llind==4 % TA70
      timeadd=.07;
    else
      timeadd=0;
    end
    plot_ica(cf,runicaout,statplot,[-.1 .6],timeadd,[0:.001:.5],llind,'erpica',fdir)
    icasaveERP{llind}.label=statplot.label;
    icasaveERP{llind}.topo=statplot.avg;
    icasaveERP{llind}.bar=runicaout.topo(:,llind);
    icasaveERP{llind}.time=[0:.001:.5];  % don't put 'timeadd' here; it will be added later when correlating to ERP
    icasaveERP{llind}.course=cf(:,1);
  end
  
  
  %   [af,bf,cf]=svd(pica(:,:,llind),'econ');
  %   statplot.avg=af(:,1);
  %   figure(40);subplot(7,3,llind*3-2);bar(aa1(:,llind));
  %   figure(40);subplot(7,3,llind*3-1);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
  %   figure(40);subplot(7,3,llind*3);plot(cf(:,1));
end
save('icasave.mat','icasave*','-append')

% Correlate (non-masked) ICA with data
for llind=1:7
  [af,bf,cf]=svd(rica(:,:,llind),'econ');
  pcarica=af(:,1)*cf(:,1)';
  for llcond=1:7
    rica_grave_corr(llind,llcond)=corr(reshape(rica(:,:,llind),[62*501 1]),reshape(squeeze(fullconmat(llcond,:,:)),[62*501 1]));
    pcarica_grave_corr(llind,llcond)=corr(reshape(pcarica,[62*501 1]),reshape(squeeze(fullconmat(llcond,:,:)),[62*501 1]));
  end
end
figure;imagescc(rica_grave_corr.*[abs(rica_grave_corr)>.7])

% Correlate sub-masked sections with sub-masked components
figure;plot(mean(stat_TPA_MSPN_1wayANOVA.mask,1))
% First 'peak' of cluster is 148-300ms
% Second 'peak' is 351-421
% Third 'peak' is 451-501

% cluster 2 (148-300ms)
mask_use=stat_TPA_MSPN_1wayANOVA.mask;
mask_use(:,1:147)=0;
mask_use(:,301:end)=0;
[rho,rho_nw,rhocond2]=mask_corr(mask_use,runicaout,fullconmat_reshape,0);
for llind=1:7
  figure;imagesc(rica(:,:,llind).*mask_use);caxis([-1.4 1.4])
  figure;imagesc(mask_use.*grave_TPA_MSPN_avg{llind}.individual);caxis([-4 4])
end
% rho = 0.3571    0.3532    *0.6089    *0.5210    0.1721    0.2454    0.0532

stat_TPA_MSPN_1wayANOVA.stat2=stat_TPA_MSPN_1wayANOVA.stat.*mask_use;

% cluster 3 (351-421ms)
mask_use=stat_TPA_MSPN_1wayANOVA.mask;
mask_use(:,1:350)=0;
mask_use(:,422:end)=0;
[rho,rho_nw,rhocond3]=mask_corr(mask_use,runicaout,fullconmat_reshape,0);
for llind=1:7
  figure;imagesc(rica(:,:,llind).*mask_use);caxis([-1.4 1.4])
  figure;imagesc(mask_use.*grave_TPA_MSPN_avg{llind}.individual);caxis([-4 4])
end
% rho = *0.7694    0.2988    0.4354    0.1108    *0.4713    0.1926    0.0820

stat_TPA_MSPN_1wayANOVA.stat3=stat_TPA_MSPN_1wayANOVA.stat.*mask_use;

% cluster 4 (451-501ms)
mask_use=stat_TPA_MSPN_1wayANOVA.mask;
mask_use(:,1:450)=0;
[rho,rho_nw,rhocond4]=mask_corr(mask_use,runicaout,fullconmat_reshape,0);
for llind=1:7
  figure;imagesc(rica(:,:,llind).*mask_use);caxis([-1.4 1.4])
  figure;imagesc(mask_use.*grave_TPA_MSPN_avg{llind}.individual);caxis([-4 4])
end
% rho = 0.3342   *0.6531    0.3425    0.3038    0.2538    0.3992    0.1928

stat_TPA_MSPN_1wayANOVA.stat4=stat_TPA_MSPN_1wayANOVA.stat.*mask_use;
stat_TPA_MSPN_1wayANOVA.stat4rm=repmat(nanmean(stat_TPA_MSPN_1wayANOVA.stat4,2),[62 501]);

% cluster 1 (early)
mask_use=clustermat; % contains 2nd cluster at p=0.065
mask_use(mask_use==1)=0;
mask_use(mask_use==2)=1;
[rho,rho_nw,rhocond1]=mask_corr(mask_use,runicaout,fullconmat_reshape,0);
for llind=1:7
  figure;imagesc(rica(:,:,llind).*mask_use);caxis([-1.4 1.4])
  figure;imagesc(mask_use.*grave_TPA_MSPN_avg{llind}.individual);caxis([-4 4])
end
% rho= *0.5890    0.2918    0.1013    0.2599    *0.7392    0.2761    0.1745

stat_TPA_MSPN_1wayANOVA.stat1=stat_TPA_MSPN_1wayANOVA.stat.*mask_use;
stat_TPA_MSPN_1wayANOVA.stat1rm=repmat(nanmean(stat_TPA_MSPN_1wayANOVA.stat1,2),[62 501]);


mask_use=zeros(size(clustermat));
mask_use(:,erp_anova_tw1(1):erp_anova_tw1(end))=1;
mask_use=mask_use.*clustermat/2;
stat_TPA_MSPN_1wayANOVA.stat1thresh=stat_TPA_MSPN_1wayANOVA.stat.*mask_use;

mask_use=zeros(size(clustermat));
mask_use(:,erp_anova_tw2(1):erp_anova_tw2(end))=1;
mask_use=mask_use.*clustermat;
stat_TPA_MSPN_1wayANOVA.stat2thresh=stat_TPA_MSPN_1wayANOVA.stat.*mask_use;

mask_use=zeros(size(clustermat));
mask_use(:,erp_anova_tw3(1):erp_anova_tw3(end))=1;
mask_use=mask_use.*clustermat;
stat_TPA_MSPN_1wayANOVA.stat3thresh=stat_TPA_MSPN_1wayANOVA.stat.*mask_use;

mask_use=zeros(size(clustermat));
mask_use(:,erp_anova_tw4(1):erp_anova_tw4(end))=1;
mask_use=mask_use.*clustermat;
stat_TPA_MSPN_1wayANOVA.stat4thresh=stat_TPA_MSPN_1wayANOVA.stat.*mask_use;

cfg=[];
cfg.layout='eeg1010';
cfg.highlight          = 'on';
cfg.highlightsize = 16;
cfg.markersymbol = 'o';
cfg.parameter='stat';

cfg.xlim=[erp_anova_tw1(1) erp_anova_tw1(end)]/1000;
cfg.highlightchannel   =  stat_TPA_MSPN_1wayANOVA.label(find(mean(stat_TPA_MSPN_1wayANOVA.stat1thresh,2)));
erp_anova_chan1 = cfg.highlightchannel;
cfg.zlim=[0 6.8]; % max of all three
figure(1);ft_topoplotER(cfg,stat_TPA_MSPN_1wayANOVA)
print(1,[fdir 'ANOVA_ERP_topo1.eps'],'-painters','-depsc')
cfg.zlim='zeromax';
figure(11);ft_topoplotER(cfg,stat_TPA_MSPN_1wayANOVA)
print(11,[fdir 'ANOVA_ERP_topo1.eps'],'-painters','-depsc')

cfg.xlim=[erp_anova_tw2(1) erp_anova_tw2(end)]/1000;
cfg.highlightchannel   =  stat_TPA_MSPN_1wayANOVA.label(find(mean(stat_TPA_MSPN_1wayANOVA.stat2thresh,2)));
erp_anova_chan2 = cfg.highlightchannel;
cfg.zlim=[0 6.8]; % max of all three
figure(2);ft_topoplotER(cfg,stat_TPA_MSPN_1wayANOVA)
print(2,[fdir 'ANOVA_ERP_topo2.eps'],'-painters','-depsc')
cfg.zlim='zeromax';
figure(22);ft_topoplotER(cfg,stat_TPA_MSPN_1wayANOVA)
print(22,[fdir 'ANOVA_ERP_topo2.eps'],'-painters','-depsc')

cfg.xlim=[erp_anova_tw3(1) erp_anova_tw3(end)]/1000;
cfg.highlightchannel   =  stat_TPA_MSPN_1wayANOVA.label(find(mean(stat_TPA_MSPN_1wayANOVA.stat3thresh,2)));
erp_anova_chan3 = cfg.highlightchannel;
cfg.zlim=[0 6.8]; % max of all three
figure(3);ft_topoplotER(cfg,stat_TPA_MSPN_1wayANOVA)
print(3,[fdir 'ANOVA_ERP_topo3.eps'],'-painters','-depsc')
cfg.zlim='zeromax';
figure(33);ft_topoplotER(cfg,stat_TPA_MSPN_1wayANOVA)
print(33,[fdir 'ANOVA_ERP_topo3.eps'],'-painters','-depsc')

cfg.xlim=[erp_anova_tw4(1) erp_anova_tw4(end)]/1000;
cfg.highlightchannel   =  stat_TPA_MSPN_1wayANOVA.label(find(mean(stat_TPA_MSPN_1wayANOVA.stat4thresh,2)));
erp_anova_chan4 = cfg.highlightchannel;
cfg.zlim=[0 6.8]; % max of all three
figure(4);ft_topoplotER(cfg,stat_TPA_MSPN_1wayANOVA)
print(4,[fdir 'ANOVA_ERP_topo4.eps'],'-painters','-depsc')
cfg.zlim='zeromax';
figure(44);ft_topoplotER(cfg,stat_TPA_MSPN_1wayANOVA)
print(44,[fdir 'ANOVA_ERP_topo4.eps'],'-painters','-depsc')

save([edir 'ANOVA_sleep0_ERP.mat'],'stat*','erp_anova*','-append')


% cfg.xlim=[.37 .39];
% cfg.highlightchannel   =  stat_TPA_MSPN_1wayANOVA.label(find(mean(stat_TPA_MSPN_1wayANOVA.stat4,2)));
% figure(4);ft_topoplotER(cfg,stat_TPA_MSPN_1wayANOVA)
% print(4,[fdir 'ANOVA_ERP_topo4.eps'],'-painters','-depsc')


% Older way: this makes false impression of F-value at zero when it is outside mask
% cfg=[];
% cfg.zlim='maxabs';
% cfg.layout='eeg1010';
% cfg.highlight          = 'on';
% cfg.highlightsize = 16;
% cfg.markersymbol = 'o';
% cfg.parameter='stat1rm';
% cfg.xlim=[.5 1.3];
% cfg.highlightchannel   =  stat_TPA_MSPN_1wayANOVA.label(find(mean(stat_TPA_MSPN_1wayANOVA.stat1,2)));
% figure(1);ft_topoplotER(cfg,stat_TPA_MSPN_1wayANOVA)
% print(1,[fdir 'ANOVA_ERP_topo1.eps'],'-painters','-depsc')
%
% cfg.parameter='stat2';
% cfg.xlim=[.21 .23];
% cfg.highlightchannel   =  stat_TPA_MSPN_1wayANOVA.label(find(mean(stat_TPA_MSPN_1wayANOVA.stat2,2)));
% figure(2);ft_topoplotER(cfg,stat_TPA_MSPN_1wayANOVA)
% print(2,[fdir 'ANOVA_ERP_topo2.eps'],'-painters','-depsc')
% cfg.parameter='stat3';
% cfg.xlim=[.37 .39];
% cfg.highlightchannel   =  stat_TPA_MSPN_1wayANOVA.label(find(mean(stat_TPA_MSPN_1wayANOVA.stat3,2)));
% figure(3);ft_topoplotER(cfg,stat_TPA_MSPN_1wayANOVA)
% print(3,[fdir 'ANOVA_ERP_topo3.eps'],'-painters','-depsc')
% cfg.parameter='stat4rm';
% cfg.xlim=[.37 .39];
% cfg.highlightchannel   =  stat_TPA_MSPN_1wayANOVA.label(find(mean(stat_TPA_MSPN_1wayANOVA.stat4,2)));
% figure(4);ft_topoplotER(cfg,stat_TPA_MSPN_1wayANOVA)
% print(4,[fdir 'ANOVA_ERP_topo4.eps'],'-painters','-depsc')




% Correlate masked data with masked components

[rho,rho_nw]=mask_corr(stat_TPA_MSPN_1wayANOVA.mask,runicaout,fullconmat_reshape);
for llind=1:7
  figure;imagesc(rica(:,:,llind).*stat_TPA_MSPN_1wayANOVA.mask);caxis([-1.4 1.4])
  figure;imagesc(stat_TPA_MSPN_1wayANOVA.mask.*grave_TPA_MSPN_avg{llind}.individual);caxis([-4 4])
end

% % PCA denoise each condition first
% for compkeep=1:9
%   for llind=1:7
%     [aa,bb,cc]=svd(grave_TPA_MSPN_avg{llind}.individual,'econ');
%     fullconmat_pcared_avg_reshape(llind,:)=reshape(aa(:,1:compkeep)*bb(1:compkeep,1:compkeep)*cc(:,1:compkeep)',[1 62*501]);
%   end
%   maskdatavec_pcared=reshape(fullconmat_pcared_avg_reshape(:,mask_induse)',[1 7*5647]); % yes transpose is critical here
%   for llind=1:7
%     [rhoc(llind,compkeep),pvalc(llind,compkeep)]=corr(maskdatavec_pcared',rica_masked_w(llind,:)');
%     [rhoc_nw(llind,compkeep),pvalc_nw(llind,compkeep)]=corr(maskdatavec_pcared',rica_masked_nw(llind,:)');
%   end
% end


% PCA denoise each condition first
for compkeep=1:9
  for llind=1:7
    [aa,bb,cc]=svd(reshape(permute(grind_TPA_MSPN_sel{llind}.individual,[2 1 3]),[62 22*501]),'econ');
    fullconmat_pcared_ind_reshape(llind,:)=reshape(squeeze(mean(reshape(aa(:,1:compkeep)*bb(1:compkeep,1:compkeep)*cc(:,1:compkeep)',[62 22 501]),2)),[1 62*501]);
    
    [aa,bb,cc]=svd(grave_TPA_MSPN_avg{llind}.individual,'econ');
    fullconmat_pcared_avg_reshape(llind,:)=reshape(aa(:,1:compkeep)*bb(1:compkeep,1:compkeep)*cc(:,1:compkeep)',[1 62*501]);
  end
  
  %   [aa4,bb4,vv4]=svd(fullconmat_pcared_reshape,'econ');
  
  clear data4ica
  % data4ica.trial{1}=fullconmat_reshape;
  data4ica.trial{1}=fullconmat_pcared_ind_reshape;
  data4ica.dimord='chan_time';
  data4ica.time{1}=.001:.001:31.062;
  % data4ica.time{1}=.001:.001:683.364;
  data4ica.label={'1' '2' '3' '4' '5' '6' '7'};
  cfg=[];
  cfg.randomseed=13;
  %   cfg.method='runica';
  %   runicaout=ft_componentanalysis(cfg, data4ica);
  cfg.method='fastica';
  fasticaout_ind=ft_componentanalysis(cfg, data4ica);
  data4ica.trial{1}=fullconmat_pcared_avg_reshape;
  fasticaout_avg=ft_componentanalysis(cfg, data4ica);
  
  % figure;for llind=1:7,subplot(4,2,llind);bar(runicaout.topo(:,llind));end
  % figure;for llind=1:7,subplot(4,2,llind);bar(fasticaout.topo(:,llind));end
  
  
  for llind=1:7
    fiica(:,:,llind)=reshape(fasticaout_ind.trial{1}(llind,:),[62 501]);
    faica(:,:,llind)=reshape(fasticaout_avg.trial{1}(llind,:),[62 501]);
    %     rica(:,:,llind)=reshape(runicaout.trial{1}(llind,:),[62 501]);
    %     pica(:,:,llind)=reshape(vv4(:,llind)',[62 501]);
    
    %     clear data4ica
    %     data4ica.dimord='chan_time';
    %     data4ica.time{1}=0:.001:.5;
    %     data4ica.label=grave_TPA_MSPN_avg{llind}.label;
    %     data4ica.trial{1}=fica(:,:,llind);
    %     cfg=[];
    %     cfg.randomseed=13;
    %     cfg.method='fastica';
    %       fAW{llind}=ft_componentanalysis(cfg, data4ica);
    %     data4ica.trial{1}=rica(:,:,llind);
    %     cfg.method='runica';
    %       rAW{llind}=ft_componentanalysis(cfg, data4ica);
    
    [af,bf,cf]=svd(fiica(:,:,llind),'econ');
    statplot.avg=af(:,1);
    statplot.time=1;
    statplot.dimord='chan_time';
    statplot.label=grave_TPA_MSPN_avg{llind}.label;
    figure(30+compkeep);subplot(7,3,llind*3-2);bar(fasticaout_ind.topo(:,llind));
    figure(30+compkeep);subplot(7,3,llind*3-1);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
    figure(30+compkeep);subplot(7,3,llind*3);plot(cf(:,1));
    
    [af,bf,cf]=svd(faica(:,:,llind),'econ');
    statplot.avg=af(:,1);
    figure(50+compkeep);subplot(7,3,llind*3-2);bar(fasticaout_avg.topo(:,llind));
    figure(50+compkeep);subplot(7,3,llind*3-1);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
    figure(50+compkeep);subplot(7,3,llind*3);plot(cf(:,1));
    
    
    %     [af,bf,cf]=svd(pica(:,:,llind),'econ');
    %     statplot.avg=af(:,1);
    %     figure(40+compkeep);subplot(7,3,llind*3-2);bar(aa4(:,llind));
    %     figure(40+compkeep);subplot(7,3,llind*3-1);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
    %     figure(40+compkeep);subplot(7,3,llind*3);plot(cf(:,1));
    
  end
end % compkeep

vv11=reshape(vv1(:,1),[62 501]); % 1st space-time map, reshaped
[u1,d1,v1]=svd(vv11,'econ'); % get best chan topography and time course
vv12=reshape(vv1(:,2),[62 501]); % 2nd space-time map, reshaped
[u2,d2,v2]=svd(vv12,'econ'); % get best chan topography and time course
vv13=reshape(vv1(:,3),[62 501]); % 3rd space-time map, reshaped
[u3,d3,v3]=svd(vv13,'econ'); % get best chan topography and time course
vv14=reshape(vv1(:,4),[62 501]); % 3rd space-time map, reshaped
[u4,d4,v4]=svd(vv14,'econ'); % get best chan topography and time course
vv15=reshape(vv1(:,5),[62 501]); % 3rd space-time map, reshaped
[u5,d5,v5]=svd(vv15,'econ'); % get best chan topography and time course
vv16=reshape(vv1(:,6),[62 501]); % 3rd space-time map, reshaped
[u6,d6,v6]=svd(vv16,'econ'); % get best chan topography and time course
vv17=reshape(vv1(:,7),[62 501]); % 3rd space-time map, reshaped
[u7,d7,v7]=svd(vv17,'econ'); % get best chan topography and time course

vv21=reshape(permute(reshape(vv3(:,1),[22 62 501]),[2 1 3]),[62 22*501]); % 1st space-time map, reshaped
[u1,d1,v1]=svd(vv21,'econ'); % get best chan topography and time course
vv22=reshape(permute(reshape(vv3(:,2),[22 62 501]),[2 1 3]),[62 22*501]); % 1st space-time map, reshaped
[u2,d2,v2]=svd(vv22,'econ'); % get best chan topography and time course
vv23=reshape(permute(reshape(vv3(:,3),[22 62 501]),[2 1 3]),[62 22*501]); % 1st space-time map, reshaped
[u3,d3,v3]=svd(vv23,'econ'); % get best chan topography and time course
vv24=reshape(permute(reshape(vv3(:,4),[22 62 501]),[2 1 3]),[62 22*501]); % 1st space-time map, reshaped
[u4,d4,v4]=svd(vv24,'econ'); % get best chan topography and time course
vv25=reshape(permute(reshape(vv3(:,5),[22 62 501]),[2 1 3]),[62 22*501]); % 1st space-time map, reshaped
[u5,d5,v5]=svd(vv25,'econ'); % get best chan topography and time course
vv26=reshape(permute(reshape(vv3(:,6),[22 62 501]),[2 1 3]),[62 22*501]); % 1st space-time map, reshaped
[u6,d6,v6]=svd(vv25,'econ'); % get best chan topography and time course
vv27=reshape(permute(reshape(vv3(:,7),[22 62 501]),[2 1 3]),[62 22*501]); % 1st space-time map, reshaped
[u7,d7,v7]=svd(vv25,'econ'); % get best chan topography and time course

% plot topography
statplot.avg=u1(:,1);
statplot.time=1;
statplot.dimord='chan_time';
statplot.label=grave_TPA_MSPN_avg{llind}.label;
figure(1);subplot(3,1,2);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
statplot.avg=u2(:,1);
figure(2);subplot(3,1,2);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
statplot.avg=u3(:,1);
figure(3);subplot(3,1,2);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
statplot.avg=u4(:,1);
figure(4);subplot(3,1,2);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
statplot.avg=u5(:,1);
figure(5);subplot(3,1,2);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
statplot.avg=u6(:,1);
figure(6);subplot(3,1,2);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
statplot.avg=u7(:,1);
figure(7);subplot(3,1,2);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)


for llind=1:7
  statplot.avg=fAW{llind}.topo(:,1);
  figure(10+llind);subplot(3,1,2);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
  statplot.avg=rAW{llind}.topo(:,1);
  figure(20+llind);subplot(3,1,2);cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
end

% % plot ICA topography
% statplot.avg=A1(:,1);
% statplot.time=1;
% statplot.dimord='chan_time';
% statplot.label=grave_TPA_MSPN_avg{llind}.label;
% figure;cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
% statplot.avg=A2(:,1);
% figure;cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
% statplot.avg=A3(:,1);
% figure;cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
% statplot.avg=A4(:,1);
% figure;cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
% statplot.avg=A5(:,1);
% figure;cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)

% plot time course
figure(1);subplot(3,1,3);plot(v1(:,1))
figure(2);subplot(3,1,3);plot(v2(:,1))
figure(3);subplot(3,1,3);plot(v3(:,1))
figure(4);subplot(3,1,3);plot(v4(:,1))
figure(5);subplot(3,1,3);plot(v5(:,1))
figure(6);subplot(3,1,3);plot(v6(:,1))
figure(7);subplot(3,1,3);plot(v7(:,1))
figure(1);subplot(3,1,3);plot(mean(reshape(v1(:,1),[22 501]),1))
figure(2);subplot(3,1,3);plot(mean(reshape(v2(:,1),[22 501]),1))
figure(3);subplot(3,1,3);plot(mean(reshape(v3(:,1),[22 501]),1))
figure(4);subplot(3,1,3);plot(mean(reshape(v4(:,1),[22 501]),1))
figure(5);subplot(3,1,3);plot(mean(reshape(v5(:,1),[22 501]),1))
figure(6);subplot(3,1,3);plot(mean(reshape(v6(:,1),[22 501]),1))
figure(7);subplot(3,1,3);plot(mean(reshape(v7(:,1),[22 501]),1))
for llind=1:7
  figure(10+llind);subplot(3,1,3);plot(fAW{llind}.trial{1}(1,:));
  figure(20+llind);subplot(3,1,3);plot(rAW{llind}.trial{1}(1,:));
end

% plot condition loading
ymax=.8;
figure(1);subplot(3,1,1);bar(aa1(:,1));ylim([-ymax ymax])
figure(2);subplot(3,1,1);bar(aa1(:,2));ylim([-ymax ymax])
figure(3);subplot(3,1,1);bar(aa1(:,3));ylim([-ymax ymax])
figure(4);subplot(3,1,1);bar(aa1(:,4));ylim([-ymax ymax])
figure(5);subplot(3,1,1);bar(aa1(:,5));ylim([-ymax ymax])
figure(6);subplot(3,1,1);bar(aa1(:,6));ylim([-ymax ymax])
figure(7);subplot(3,1,1);bar(aa1(:,7));ylim([-ymax ymax])

for llind=1:7
  ymax=0.8;
  figure(10+llind);subplot(3,1,1);bar(fasticaout.topo(:,llind));ylim([-ymax ymax])
  ymax=2.2;
  figure(20+llind);subplot(3,1,1);bar(runicaout.topo(:,llind));ylim([-ymax ymax])
end

%

% erpname={'erp_AT500' 'erp_AT70' 'erp_AT20' 'erp_AT0' 'erp_TA20' 'erp_TA70' 'erp_TA500'};
% for llind=1:7
%   spm_eeg_ft2spm(grind_TPA_MSPN{llind},erpname{llind})
% end

%% PCA of main results
% obsolete now

load tlock_statmc_sleep0_iter27_statwinorig0_ftver0mcseed13
soalist=[1 3 4 5 6 7 9];
for ind=1:length(soalist)
  statmat(ind,:)=reshape(statt_mc{soalist(ind),3,ss}.stat,[1 61*501]);
end
[uu,dd,vv]=svd(statmat,'econ');

vv1=reshape(vv(:,1),[61 501]); % 1st space-time map, reshaped
[u1,d1,v1]=svd(vv1,'econ'); % get best chan topography and time course
vv2=reshape(vv(:,2),[61 501]); % 2nd space-time map, reshaped
[u2,d2,v2]=svd(vv2,'econ'); % get best chan topography and time course
vv3=reshape(vv(:,3),[61 501]); % 3rd space-time map, reshaped
[u3,d3,v3]=svd(vv3,'econ'); % get best chan topography and time course
vv4=reshape(vv(:,4),[61 501]); % 3rd space-time map, reshaped
[u4,d4,v4]=svd(vv4,'econ'); % get best chan topography and time course
vv5=reshape(vv(:,5),[61 501]); % 3rd space-time map, reshaped
[u5,d5,v5]=svd(vv5,'econ'); % get best chan topography and time course

% plot topography
statplot.avg=u1(:,1);
statplot.time=1;
statplot.dimord='chan_time';
statplot.label=statt_mc{9,3,ss}.label;
cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
statplot.avg=u2(:,1);
cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
statplot.avg=u3(:,1);
cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
statplot.avg=u4(:,1);
cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)
statplot.avg=u5(:,1);
cfg=[];cfg.latency=1;cfg.layout='eeg1010';ft_topoplotER(cfg,statplot)

% plot time course
figure;plot(v1(:,1))
figure;plot(v2(:,1))
figure;plot(v3(:,1))
figure;plot(v4(:,1))
figure;plot(v5(:,1))

% plot condition loading
ymax=.8;
figure;bar(uu(:,1));ylim([-ymax ymax])
figure;bar(uu(:,2));ylim([-ymax ymax])
figure;bar(uu(:,3));ylim([-ymax ymax])
figure;bar(uu(:,4));ylim([-ymax ymax])
figure;bar(uu(:,5));ylim([-ymax ymax])

% % %
if sleep
  load tlock_grind_sleep1_ss12_iter11_trialkc-1_statwinorig0_ftver0.mat
else
  load tlock_grind_sleep0_iter27_statwinorig0_ftver0
end
grindDiff(:,:,:,1)=grind_tacPaud_save{1,3,ss}.individual(:,:,500:1000)-grind_tacMSpN_save{1,3,ss}.individual(:,:,500:1000);
grindDiff(:,:,:,2)=grind_tacPaud_save{3,3,ss}.individual(:,:,500:1000)-grind_tacMSpN_save{3,3,ss}.individual(:,:,500:1000);
grindDiff(:,:,:,3)=grind_tacPaud_save{4,3,ss}.individual(:,:,500:1000)-grind_tacMSpN_save{4,3,ss}.individual(:,:,500:1000);
grindDiff(:,:,:,4)=grind_tacPaud_save{5,3,ss}.individual(:,:,500:1000)-grind_tacMSpN_save{5,3,ss}.individual(:,:,500:1000);
grindDiff(:,:,:,5)=grind_tacPaud_save{6,3,ss}.individual(:,:,520:1020)-grind_tacMSpN_save{6,3,ss}.individual(:,:,520:1020);
grindDiff(:,:,:,6)=grind_tacPaud_save{7,3,ss}.individual(:,:,570:1070)-grind_tacMSpN_save{7,3,ss}.individual(:,:,570:1070);
grindDiff(:,:,:,7)=grind_tacPaud_save{9,3,ss}.individual(:,:,1000:1500)-grind_tacMSpN_save{9,3,ss}.individual(:,:,1000:1500);


