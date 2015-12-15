% preprocessing of EEG data from Hills, 64ch MR-cap
% clear all
close all
if ispc
  edir='D:\audtac\eeg_data\';
  ddir='D:\audtac\legomagic\diaries\';
else
  edir='/mnt/hgfs/D/audtac/eeg_data/';
  ddir='/mnt/hgfs/D/audtac/legomagic/diaries/';
end
cd(edir)

% sub{100}='p01'; 
sub{1}='e01'; 
sub{2}='e02'; 
sub{3}='e03'; 
sub{4}='e04'; 
%%%  above is pilot, below is real
sub{5}='e05'; 
sub{6}='e06'; 
sub{7}='e07'; 
sub{8}='e08'; 
sub{9}='e09';
sub{10}='e10';
sub{11}='e11';
sub{12}='e12';
sub{13}='e13';
sub{14}='e14';
sub{15}='e15';
sub{16}='e16';
sub{17}='e17';
sub{18}='e18';
sub{19}='e19';
sub{20}='e20';
sub{21}='e21';
sub{22}='e22';
sub{23}='e23';
sub{24}='e24';
sub{25}='e25';
sub{26}='e26';
sub{27}='e27';
sub{28}='e28';
sub{29}='e29';
sub{30}='e30';
sub{31}='e31';
sub{32}='e32';

ii=5;

cd([edir sub{ii} ])

if ispc
  rmpath(genpath('D:\matlab\spm8\external\fieldtrip\'))
  rmpath(genpath('D:\fieldtrip_svn\'))
  addpath('D:\fieldtrip_svn\')
else
  rmpath(genpath('/mnt/hgfs/D/matlab/spm8/external/fieldtrip/'))
  rmpath(genpath('/mnt/hgfs/D/fieldtrip_svn/'))
  addpath('/mnt/hgfs/D/fieldtrip_svn/')
end
which ft_defaults.m
ft_defaults;



%% freq filtering

load(['raw_each_rej_' sub{ii}],'raw_*');

cfg=[];
cfg.channel= {'all', '-ECG'};
cfg.reref='yes';
cfg.refchannel={'all', '-ECG'};
data_tac_ref=ft_preprocessing(cfg,raw_tac_rej);
data_aud_ref=ft_preprocessing(cfg,raw_aud_rej);
data_nul_ref=ft_preprocessing(cfg,raw_nul_rej);

clear raw_*

frq=[4 8; 8 13; 13 30; 30 48; 52 85];
if ii<8
  soalist=[3 4 5 6 7];
else
  soalist=[1 3 4 5 6 7 9];
end

for ff=5:size(frq,1)
    cfg=[];
    cfg.bpfilter='yes';
    cfg.bpfreq=frq(ff,:);
    cfg.padding=data_tac_ref.time{1}(end)-data_tac_ref.time{1}(1)+1;
    cfg.demean='yes';
    cfg.baselinewindow=[-.6 -.1];
    cfg.channel= {'all', '-ECG'}
    data_tac_filt{ff}=ft_preprocessing(cfg,data_tac_ref);
    data_aud_filt{ff}=ft_preprocessing(cfg,data_aud_ref);
    data_nul_filt{ff}=ft_preprocessing(cfg,data_nul_ref);
    cfg=[];
    cfg.hilbert='abs';
    data_tac_hilm{ff}=ft_preprocessing(cfg,data_tac_filt{ff});
    data_aud_hilm{ff}=ft_preprocessing(cfg,data_aud_filt{ff});
    data_nul_hilm{ff}=ft_preprocessing(cfg,data_nul_filt{ff});
    cfg=[];
    cfg.hilbert='angle';
    data_tac_hilp{ff}=ft_preprocessing(cfg,data_tac_filt{ff});
    data_aud_hilp{ff}=ft_preprocessing(cfg,data_aud_filt{ff});
    data_nul_hilp{ff}=ft_preprocessing(cfg,data_nul_filt{ff});

    for ll=1:length(soalist)
      for tt=3
        cfgs=[];
        cfgs.latency=[-.6 1.2];
        cfg=[];
        cfg.trials=[data_tac_filt{ff}.trialinfo(:,tt)==soalist(ll) & data_tac_filt{ff}.trialinfo(:,2)==0 & ~isnan(data_tac_filt{ff}.trialinfo(:,10))];
        tlock_tac{soalist(ll),ff}=ft_timelockanalysis(cfg,ft_selectdata_new(cfgs,data_tac_filt{ff}));
        tlock_tac_hm{soalist(ll),ff}=ft_timelockanalysis(cfg,ft_selectdata_new(cfgs,data_tac_hilm{ff}));
        cfg.keeptrials='yes';
        tlock_tac_hp{soalist(ll),ff}=ft_timelockanalysis(cfg,ft_selectdata_new(cfgs,data_tac_hilp{ff}));
        cfg=[];
%         cfg.trials=data_aud_filt{ff}.trialinfo(:,3)==soalist(ll);
        cfg.trials=[data_aud_filt{ff}.trialinfo(:,tt)==soalist(ll) & data_aud_filt{ff}.trialinfo(:,2)==0 & ~isnan(data_aud_filt{ff}.trialinfo(:,10))];
        tlock_aud{soalist(ll),ff}=ft_timelockanalysis(cfg,ft_selectdata_new(cfgs,data_aud_filt{ff}));
        tlock_aud_hm{soalist(ll),ff}=ft_timelockanalysis(cfg,ft_selectdata_new(cfgs,data_aud_hilm{ff}));
        cfg.keeptrials='yes';
        tlock_aud_hp{soalist(ll),ff}=ft_timelockanalysis(cfg,ft_selectdata_new(cfgs,data_aud_hilp{ff}));
      end
    end
    
    cfg=[];
%     cfg.trials=isnan(data_tac_filt{ff}.trialinfo(:,3));
    cfg.trials=[data_tac_filt{ff}.trialinfo(:,3)==-2 & data_tac_filt{ff}.trialinfo(:,2)==0];
    tlock_tac{10,ff}=ft_timelockanalysis(cfg,ft_selectdata_new(cfgs,data_tac_filt{ff}));
    tlock_tac_hm{10,ff}=ft_timelockanalysis(cfg,ft_selectdata_new(cfgs,data_tac_hilm{ff}));
    cfg.keeptrials='yes';
    tlock_tac_hp{10,ff}=ft_timelockanalysis(cfg,ft_selectdata_new(cfgs,data_tac_hilp{ff}));
    
    
    cfg=[];
%     cfg.trials=isnan(data_aud_filt{ff}.trialinfo(:,3));
    cfg.trials=[data_aud_filt{ff}.trialinfo(:,3)==-1 & data_aud_filt{ff}.trialinfo(:,2)==0];
    tlock_aud{10,ff}=ft_timelockanalysis(cfg,ft_selectdata_new(cfgs,data_aud_filt{ff}));
    tlock_aud_hm{10,ff}=ft_timelockanalysis(cfg,ft_selectdata_new(cfgs,data_aud_hilm{ff}));
    cfg.keeptrials='yes';
    tlock_aud_hp{10,ff}=ft_timelockanalysis(cfg,ft_selectdata_new(cfgs,data_aud_hilp{ff}));
    
    cfg=[];
    cfg.trials=data_nul_filt{ff}.trialinfo(:,2)==0;
    tlock_nul{10,ff}=ft_timelockanalysis(cfg,ft_selectdata_new(cfgs,data_nul_filt{ff}));
    tlock_nul_hm{10,ff}=ft_timelockanalysis(cfg,ft_selectdata_new(cfgs,data_nul_hilm{ff}));
    cfg.keeptrials='yes';
    tlock_nul_hp{10,ff}=ft_timelockanalysis(cfg,ft_selectdata_new(cfgs,data_nul_hilp{ff}));

    clear data_tac_filt data_tac_hilm data_tac_hilp
    clear data_aud_filt data_aud_hilm data_aud_hilp
    clear data_nul_filt data_nul_hilm data_nul_hilp

    save(['tlock_each_filthmhp_' num2str(frq(ff,1)) 'to' num2str(frq(ff,2)) '_' sub{ii}],'tlock_*');
end

clear data*ref

return

%%  timelock average

for ff=1:size(frq,1)
    tmp=load(['tlock_each_filthmhp_' num2str(frq(ff,1)) 'to' num2str(frq(ff,2)) '_' sub{ii}],'tlock_tac*');
    
    for ll=1:length(soalist)
        tlock_tac{soalist(ll),ff}=tmp.tlock_tac{soalist(ll),ff};
        tlock_tac_hm{soalist(ll),ff}=tmp.tlock_tac_hm{soalist(ll),ff};
        tlock_tac_hp{soalist(ll),ff}=tmp.tlock_tac_hp{soalist(ll),ff};
    end
    clear tmp
end


for ll=1:length(soalist)
    for ff=1:size(frq,1)
        hmavg(:,ff,:)=tlock_tac_hm{soalist(ll),ff}.avg;
    end
    hmfreq{soalist(ll)}.dimord='chan_freq_time';
    hmfreq{soalist(ll)}.powspctrm=hmavg;
    hmfreq{soalist(ll)}.label=tlock_tac_hm{soalist(ll),ff}.label;
    hmfreq{soalist(ll)}.freq=mean(frq,2);
    hmfreq{soalist(ll)}.time=tlock_tac_hm{soalist(ll),ff}.time;

    for ff=1:size(frq,1)
        hpavg(:,ff,:)=tlock_tac_hp{soalist(ll),ff}.avg;
    end
    hpfreq{soalist(ll)}.dimord='chan_freq_time';
    hpfreq{soalist(ll)}.powspctrm=hpavg;
    hpfreq{soalist(ll)}.label=tlock_tac_hp{soalist(ll),ff}.label;
    hpfreq{soalist(ll)}.freq=mean(frq,2);
    hpfreq{soalist(ll)}.time=tlock_tac_hp{soalist(ll),ff}.time;
end


cfg=[];
cfg.interactive='yes';
cfg.layout='EEG1010.lay';
cfg.baseline=[-.6 -.45];
cfg.baselinetype='relative';
ft_multiplotTFR(cfg,hmfreq{soalist(ll)});
cfg.baselinetype='absolute';
ft_multiplotTFR(cfg,hpfreq{soalist(ll)});




%%


% cfg=[];
% tlock_tac{11}=ft_timelockanalysis(cfg,data_tac_filt_ref);
% cfg=[];
% tlock_aud{11}=ft_timelockanalysis(cfg,data_aud_filt_ref);
% 
% save(['tlock_' sub{ii}],'tlock*');

% cfg=[];
% cfg.interactive='yes';
% cfg.layout='EEG1010.lay';
% ft_multiplotER(cfg,tlock_tac_hm{1});

% ft_multiplotER(cfg,tlock_tac{10});
% ft_multiplotER(cfg,tlock_aud{10});
% ft_multiplotER(cfg,tlock_nul{11});

figure;
for tt=1:floor([tlock_tac.time(end)-tlock_tac.time(1)]/.04)
    cfg=[];
    cfg.xlim=[[tlock_tac.time(1)+(tt-1)*.04]  [tlock_tac.time(1) + tt*.04] ];
    cfg.zlim=[-2 2];
    cfg.layout='EEG1010.lay';
    subplot(5,10,tt);
    ft_topoplotER(cfg,tlock_tac);
end

figure;
for tt=1:floor([tlock_aud.time(end)-tlock_aud.time(1)]/.04)
    cfg=[];
    cfg.xlim=[[tlock_aud.time(1)+(tt-1)*.04]  [tlock_aud.time(1) + tt*.04] ];
    cfg.zlim=[-2 2];
    cfg.layout='EEG1010.lay';
    subplot(5,10,tt);
    ft_topoplotER(cfg,tlock_aud);
end

figure;
for tt=1:floor([tlock_nul.time(end)-tlock_nul.time(1)]/.04)
    cfg=[];
    cfg.xlim=[[tlock_nul.time(1)+(tt-1)*.04]  [tlock_nul.time(1) + tt*.04] ];
    cfg.zlim=[-2 2];
    cfg.layout='EEG1010.lay';
    subplot(5,10,tt);
    ft_topoplotER(cfg,tlock_nul);
end


cfg=[];
cfg.operation='subtract';
cfg.parameter='avg';

tlock_tac_aud=ft_math(cfg,tlock_tac{10},tlock_aud{10});

tlock_tac_aud=ft_math(cfg,tlock_tac{5},tlock_tac{10});

tlock_tac_aud=ft_math(cfg,tlock_tac{5},tlock_tac{6});

tlock_tac_aud=ft_math(cfg,tlock_tac{5},tlock_tac{9});

tlock_tac_aud=ft_math(cfg,tlock_aud{5},tlock_aud{10});

tlock_tac_aud=ft_math(cfg,tlock_aud{5},tlock_aud{6});

tlock_tac_aud=ft_math(cfg,tlock_aud{5},tlock_aud{9});

cfg=[];
cfg.interactive='yes';
cfg.layout='EEG1010.lay';
ft_multiplotER(cfg,tlock_tac_aud);



cfg=[];
cfg.layout='EEG1010.lay';
cfg.xlim=[.05 .1];
cfg.zlim='maxabs';
ft_topoplotER(cfg, tlock_tac)

%%  compute each SOA timelock with keeptrials


soalist=[1 3 4 5 6 7 9];
for ll=1:7
    cfg=[];
    cfg.keeptrials='yes';
    cfg.trials=data_tac_filt_ref.trialinfo(:,3)==soalist(ll);
    tlock_tac{soalist(ll)}=ft_timelockanalysis(cfg,data_tac_filt_ref);
    cfg=[];
    cfg.keeptrials='yes';
    cfg.trials=data_aud_filt_ref.trialinfo(:,3)==soalist(ll);
    tlock_aud{soalist(ll)}=ft_timelockanalysis(cfg,data_aud_filt_ref);
end
cfg=[];
cfg.keeptrials='yes';
cfg.trials=isnan(data_tac_filt_ref.trialinfo(:,3));
tlock_tac{10}=ft_timelockanalysis(cfg,data_tac_filt_ref);
cfg=[];
cfg.keeptrials='yes';
cfg.trials=isnan(data_aud_filt_ref.trialinfo(:,3));
tlock_aud{10}=ft_timelockanalysis(cfg,data_aud_filt_ref);
cfg=[];
cfg.keeptrials='yes';
tlock_nul{10}=ft_timelockanalysis(cfg,data_nul_filt_ref);

%% sensor level stats

load elec1010_neighb.mat
cfg = [];
cfg.method = 'montecarlo';       % use the Monte Carlo Method to calculate the significance probability
cfg.statistic = 'indepsamplesT'; % use the independent samples T-statistic as a measure to 
                                 % evaluate the effect at the sample level
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;         % alpha level of the sample-specific test statistic that 
                                 % will be used for thresholding
cfg.clusterstatistic = 'maxsum'; % test statistic that will be evaluated under the 
                                 % permutation distribution. 
cfg.minnbchan = 2;               % minimum number of neighborhood channels that is 
                                 % required for a selected sample to be included 
                                 % in the clustering algorithm (default=0).
cfg.tail = 0;                    % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail = 0;
cfg.alpha = 0.025;               % alpha level of the permutation test
cfg.numrandomization = 100;      % number of draws from the permutation distribution
cfg.neighbours    = neighbours;  % the neighbours specify for each sensor with 
                                 % which other sensors it can form clusters
% cfg.channel       = {'MEG'};     % cell-array with selected channel labels
cfg.latency       = [-.5 1];       % time interval over which the experimental 
                                 % conditions must be compared (in seconds)
cfg.ivar  = 1;      

cfgm=[];
cfgm.operation='subtract';
cfgm.parameter='avg';

for ll=[1 3 4 5 6 7 9]
    data1=tlock_tac{ll}; % TA_ll
    data2=tlock_tac{10}; % T_0
    
    design = zeros(1,size(data1.trial,1) + size(data2.trial,1));
    design(1,1:size(data1.trial,1)) = 1;
    design(1,(size(data1.trial,1)+1):(size(data1.trial,1) + size(data2.trial,1)))= 2;
    
    cfg.design = design;             % design matrix
    
    stat{ll} = ft_timelockstatistics(cfg, data1, data2);
    
    aud_residual{ll}=ft_math(cfgm,ft_timelockanalysis([],data1),ft_timelockanalysis([],data2));

    data1=tlock_aud{ll}; % AT_ll
    data2=tlock_aud{10}; % T_0
    
    design = zeros(1,size(data1.trial,1) + size(data2.trial,1));
    design(1,1:size(data1.trial,1)) = 1;
    design(1,(size(data1.trial,1)+1):(size(data1.trial,1) + size(data2.trial,1)))= 2;
    
    cfg.design = design;             % design matrix
    
    stat{ll} = ft_timelockstatistics(cfg, data1, data2);
    
    tac_residual{ll}=ft_math(cfgm,ft_timelockanalysis([],data1),ft_timelockanalysis([],data2));
end

ll=9;
figure;imagesc(aud_residual{ll}.avg(:,dsearchn(aud_residual{ll}.time',stat{ll}.time(1)):dsearchn(aud_residual{ll}.time',stat{ll}.time(end))).*stat{ll}.mask);caxis([-6 6]);

figure; % auditory residual, with tactile at time 0 subtracted out
for ll=[1 3 4 5 6 7 9]
    subplot(9,1,ll)
    imagesc(stat{ll}.time,1:63,aud_residual{ll}.avg(:,dsearchn(aud_residual{ll}.time',stat{ll}.time(1)):dsearchn(aud_residual{ll}.time',stat{ll}.time(end))));caxis([-6 6]);
    ylabel(num2str(soatimes(ll)));
end
xlabel('Time 0 indicates tactile onset, subtracted out')

figure; % Tactile at time 0, with auditory at varying SOA
[sel1,sel2]=match_str(aud_residual{ll}.label,tlock_tac{ll}.label);
for ll=[1 3 4 5 6 7 9]
    subplot(9,1,ll)
    imagesc(stat{ll}.time,1:63,tlock_tac{ll}.avg(sel2,dsearchn(tlock_tac{ll}.time',stat{ll}.time(1)):dsearchn(tlock_tac{ll}.time',stat{ll}.time(end))));caxis([-6 6]);
    ylabel(num2str(soatimes(ll)));
end
xlabel('Time 0 indicates tactile onset')
subplot(9,1,8);
imagesc(stat{ll}.time,1:63,tlock_tac{10}.avg(sel2,dsearchn(tlock_tac{ll}.time',stat{ll}.time(1)):dsearchn(tlock_tac{ll}.time',stat{ll}.time(end))));caxis([-6 6]);
ylabel('Tac ALONE');

ll=9;
figure;imagesc(tac_residual{ll}.avg(:,dsearchn(tac_residual{ll}.time',stat{ll}.time(1)):dsearchn(tac_residual{ll}.time',stat{ll}.time(end))).*stat{ll}.mask);caxis([-6 6]);

figure; % tactile residual, with auditory at time 0 subtracted out
for ll=[9 7 6 5 4 3 1]
    subplot(9,1,(-(ll-5))+5)
    imagesc(stat{ll}.time,1:63,tac_residual{ll}.avg(:,dsearchn(tac_residual{ll}.time',stat{ll}.time(1)):dsearchn(tac_residual{ll}.time',stat{ll}.time(end))));caxis([-6 6]);
    ylabel(num2str(-soatimes(ll)));
end
xlabel('Time 0 indicates auditory onset, subtracted out')

figure; % Auditory at time 0, with tactile at varying SOA
[sel1,sel2]=match_str(tac_residual{ll}.label,tlock_aud{ll}.label);
for ll=[9 7 6 5 4 3 1]
    subplot(9,1,(-(ll-5))+5)
%     imagesc(stat{ll}.time,1:63,tlock_aud{ll}.avg(:,dsearchn(tlock_aud{ll}.time',stat{ll}.time(1)):dsearchn(tlock_aud{ll}.time',stat{ll}.time(end))));caxis([-6 6]);
    imagesc(stat{ll}.time,1:63,tlock_aud{ll}.avg(sel2,dsearchn(tlock_aud{ll}.time',stat{ll}.time(1)):dsearchn(tlock_aud{ll}.time',stat{ll}.time(end))));caxis([-6 6]);
    ylabel(num2str(-soatimes(ll)));
end
xlabel('Time 0 indicates auditory onset')
subplot(9,1,8);
imagesc(stat{ll}.time,1:63,tlock_aud{10}.avg(sel2,dsearchn(tlock_aud{ll}.time',stat{ll}.time(1)):dsearchn(tlock_aud{ll}.time',stat{ll}.time(end))));caxis([-6 6]);
ylabel('Aud ALONE');

figure;imagesc(data1_data2.avg(:,dsearchn(data1_data2.time',stat.time(1)):dsearchn(data1_data2.time',stat.time(end))).*stat.mask);caxis([-6 6]);

cfg=[];
cfg.xlim=[.08 .081];
cfg.layout='elec1010.lay';
cfg.zlim='maxabs';
ft_topoplotER(cfg,data1_data2);

% Make a vector of all p-values associated with the clusters from ft_timelockstatistics.
pos_cluster_pvals = [stat.posclusters(:).prob];
% Then, find which clusters are significant, outputting their indices as held in stat.posclusters
pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
% (stat.cfg.alpha is the alpha level we specified earlier for cluster comparisons; In this case, 0.025)
% make a boolean matrix of which (channel,time)-pairs are part of a significant cluster
pos = ismember(stat.posclusterslabelmat, pos_signif_clust);

% and now for the negative clusters...
neg_cluster_pvals = [stat.negclusters(:).prob];
neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
neg = ismember(stat.negclusterslabelmat, neg_signif_clust);

timestep = 0.03;		% timestep between time windows for each subplot (in seconds)
sampling_rate = data_tac_filt.fsample;	% Data has a temporal resolution of 300 Hz
sample_count = length(stat.time);
					% number of temporal samples in the statistics object
j = [0:timestep:1];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in MEEG samples

for k = 1:20;
     subplot(4,5,k);
     cfg = [];   
     cfg.xlim=[j(k) j(k+1)];   % time interval of the subplot
     cfg.zlim = [-2.5 2.5];
%      cfg.zlim = 'maxabs';
   % If a channel reaches this significance, then
   % the element of pos_int with an index equal to that channel
   % number will be set to 1 (otherwise 0).
   
   % Next, check which channels are significant over the
   % entire time interval of interest.
     pos_int = all(pos(:, m(k):m(k+1)), 2);
     neg_int = all(neg(:, m(k):m(k+1)), 2);

     cfg.highlight = 'on';
   % Get the index of each significant channel
     cfg.highlightchannel = find(pos_int | neg_int)
     cfg.comment = 'xlim';   
     cfg.commentpos = 'title';   
     cfg.layout = 'elec1010.lay';
     ft_topoplotER(cfg, data1_data2);   
end


