% statistics of EEG awake data

clear all
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


%%
% contrast of conditions in GLM: Aalone + Talone_shifted - ATmult - N

% ref='average';
for ii=5:18
  cd([edir sub{ii} ])
  clearvars -except ii sub edir ddir ref
  
%   while ~exist('ttrialkept','var')
%     load(['tlock_diffs_' sub{ii} '.mat']);
%     pause(10)
%   end
%   
%   load(['raw_each_rej_' sub{ii}],'raw_*');
%   
%   switch ref
%     case 'average'
%       cfg=[];
%       cfg.bpfilter='yes';
%       cfg.bpfreq=[1 40];
%       cfg.demean='yes';
%       cfg.baselinewindow=[-.6 -.1];
%       cfg.channel= {'all', '-ECG'};
%       cfg.reref='yes';
%       cfg.refchannel='all'; % necessary if doing source localisation
%       data_tac_filt_ref=ft_preprocessing(cfg,raw_tac_rej);
%       data_aud_filt_ref=ft_preprocessing(cfg,raw_aud_rej);
%       data_nul_filt_ref=ft_preprocessing(cfg,raw_nul_rej);
%     case 'ft'
%       cfg=[];
%       cfg.bpfilter='yes';
%       cfg.bpfreq=[1 40];
%       cfg.demean='yes';
%       cfg.baselinewindow=[-.6 -.1];
%       cfg.channel= {'all', '-ECG'}
%       data_tac_filt=ft_preprocessing(cfg,raw_tac_rej);
%       data_aud_filt=ft_preprocessing(cfg,raw_aud_rej);
%       data_nul_filt=ft_preprocessing(cfg,raw_nul_rej);
%       
%       % % Solve this: what should be reference channel for scalp ERP? fixme
%       cfg=[];
%       cfg.reref='yes';
%       %       cfg.refchannel='all'; % necessary if doing source localisation
%       cfg.refchannel={'FT9', 'FT10'}; % sort of like linked mastoids?
%       data_tac_filt_ref=ft_preprocessing(cfg,data_tac_filt);
%       data_aud_filt_ref=ft_preprocessing(cfg,data_aud_filt);
%       data_nul_filt_ref=ft_preprocessing(cfg,data_nul_filt);
%   end
%   
%   clear raw_*
%   
%   if ii>7
%     soalist=[1 3 4 5 6 7 9];
%   else
%     soalist=[3 4 5 6 7];
%   end
%   % include only trials during awake segments (in case participant went into N1 during sitting up portion)
%   for ll=1:length(soalist)
%     for tt=3:6 % different thresholds of inclusion of SOA variance around desired SOA (see concat_stim_files.m)
%       cfg=[];
%       cfg.vartrllength=2;
%       cfg.keeptrials='yes';
%       cfg.trials=[data_tac_filt_ref.trialinfo(:,tt)==soalist(ll) & data_tac_filt_ref.trialinfo(:,2)==0 & ~isnan(data_tac_filt_ref.trialinfo(:,10))];
%       tlock_tac_s0{soalist(ll),tt-2}=ft_timelockanalysis(cfg,data_tac_filt_ref); % only 'awake' state
%       %   cfg.trials=data_tac_filt_ref.trialinfo(:,3)==soalist(ll);
%       %   tlock_tac_sall{soalist(ll)}=ft_timelockanalysis(cfg,data_tac_filt_ref); % all trials, irrepective if fell asleep
%     end
%   end
%   for ll=1:length(soalist)
%     for tt=3:6 % different thresholds of inclusion of SOA variance around desired SOA (see concat_stim_files.m)
%       cfg=[];
%       cfg.vartrllength=2;
%       cfg.keeptrials='yes';
%       cfg.trials=[data_aud_filt_ref.trialinfo(:,tt)==soalist(ll) & data_aud_filt_ref.trialinfo(:,2)==0 & ~isnan(data_aud_filt_ref.trialinfo(:,10))];
%       tlock_aud_s0{soalist(ll),tt-2}=ft_timelockanalysis(cfg,data_aud_filt_ref);
%       %   cfg.trials=data_aud_filt_ref.trialinfo(:,3)==soalist(ll);
%       %   tlock_aud_sall{soalist(ll)}=ft_timelockanalysis(cfg,data_aud_filt_ref);
%     end
%   end
%   cfg=[]; % tac alone 38 (-2)
%   cfg.vartrllength=2;
%   cfg.keeptrials='yes';
%   cfg.trials=[data_tac_filt_ref.trialinfo(:,3)==-2 & data_tac_filt_ref.trialinfo(:,2)==0];
%   tlock_tac_s0{10,1}=ft_timelockanalysis(cfg,data_tac_filt_ref);
%   % cfg.trials=data_tac_filt_ref.trialinfo(:,3)==-2;
%   % tlock_tac_sall{10}=ft_timelockanalysis(cfg,data_tac_filt_ref);
%   cfg=[]; % aud alone 39 (-1)
%   cfg.vartrllength=2;
%   cfg.keeptrials='yes';
%   cfg.trials=[data_aud_filt_ref.trialinfo(:,3)==-1 & data_aud_filt_ref.trialinfo(:,2)==0];
%   tlock_aud_s0{10,1}=ft_timelockanalysis(cfg,data_aud_filt_ref);
%   % cfg.trials=data_aud_filt_ref.trialinfo(:,3)==-1;
%   % tlock_aud_sall{10}=ft_timelockanalysis(cfg,data_aud_filt_ref);
%   cfg=[];
%   cfg.vartrllength=2;
%   cfg.keeptrials='yes';
%   cfg.trials=data_nul_filt_ref.trialinfo(:,2)==0;
%   tlock_nul_s0{10,1}=ft_timelockanalysis(cfg,data_nul_filt_ref);
%   % cfg.trials='all';
%   % tlock_nul_sall{10}=ft_timelockanalysis(cfg,data_nul_filt_ref);
%   
%   clear data*
%   for ll=[soalist 10] % this loop is to remove trials with NaN in (e.g. if artifact rejection cut through middle of a trial)
%     for tt=1:4
%       if ll==10 && tt>1
%         continue
%       else
%         cfg=[];
%         cfg.keeptrials='yes';
%         cfg.trials=find(sum(squeeze(isnan(tlock_tac_s0{ll,tt}.trial(:,10,:)))')<.05*size(tlock_tac_s0{ll,tt}.trial,3));
%         tlock_tac_s0{ll,tt}=ft_timelockanalysis(cfg,tlock_tac_s0{ll,tt});
%         %   cfg.trials=find(sum(squeeze(isnan(tlock_tac_sall{ll}.trial(:,10,:)))')<.05*size(tlock_tac_sall{ll}.trial,3));
%         %   tlock_tac_sall{ll}=ft_timelockanalysis(cfg,tlock_tac_sall{ll});
%         cfg=[];
%         cfg.keeptrials='yes';
%         cfg.trials=find(sum(squeeze(isnan(tlock_aud_s0{ll,tt}.trial(:,10,:)))')<.05*size(tlock_aud_s0{ll,tt}.trial,3));
%         tlock_aud_s0{ll,tt}=ft_timelockanalysis(cfg,tlock_aud_s0{ll,tt});
%         %   cfg.trials=find(sum(squeeze(isnan(tlock_aud_sall{ll}.trial(:,10,:)))')<.05*size(tlock_aud_sall{ll}.trial,3));
%         %   tlock_aud_sall{ll}=ft_timelockanalysis(cfg,tlock_aud_sall{ll});
%       end
%       if ll==10 && tt==1
%         cfg=[];
%         cfg.keeptrials='yes';
%         cfg.trials=find(sum(squeeze(isnan(tlock_nul_s0{ll}.trial(:,10,:)))')<.05*size(tlock_nul_s0{ll}.trial,3));
%         tlock_nul_s0{ll}=ft_timelockanalysis(cfg,tlock_nul_s0{ll});
%         %     cfg.trials=find(sum(squeeze(isnan(tlock_nul_sall{ll}.trial(:,10,:)))')<.05*size(tlock_nul_sall{ll}.trial,3));
%         %     tlock_nul_sall{ll}=ft_timelockanalysis(cfg,tlock_nul_sall{ll});
%       end
%     end
%   end
%   
%   load([ddir sub{ii} '_audtac.mat']);
%   soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];
%   chanuse=match_str(tlock_tac_s0{3}.label,'Fz');
%   fsample=1/diff(tlock_tac_s0{3,1}.time(1:2)); % sampling rate
%   
  
  for ll=soalist
    for tt=1:4
      
%       % make sure the multisensory auditory-locked and ms tactile-locked have same trials; sometimes not exactly!
%       trcom=intersect(tlock_tac_s0{ll,tt}.trialinfo(:,11),tlock_aud_s0{ll,tt}.trialinfo(:,11));
%       [~,tbb,~]=intersect(tlock_tac_s0{ll,tt}.trialinfo(:,11),trcom);
%       [~,abb,~]=intersect(tlock_aud_s0{ll,tt}.trialinfo(:,11),trcom);
%       
%       cfg=[];
%       % I have thought it through, and no minus sign needed here.
%       %     cfg.offset=round(fsample*(tlock_tac_s0{ll,tt}.trialinfo(trialkept,10)-soades(ll))); % offset is in samples not seconds
%       cfg.offset=round(fsample*(tlock_tac_s0{ll,tt}.trialinfo(tbb(mtrialkept{ll,tt}),10)-soades(ll) -soades(ll))); % offset is in samples not seconds; -2*soades required for tactile not auditory
%       cfg.trials=ttrialkept{ll,tt};
%       tlock_tac_s0{ll+40,tt}=ft_redefinetrial(cfg,tlock_tac_s0{10}); % Talone with shift in time matching SOA jitter
%       cfg=[];
%       cfg.vartrllength=2;
%       tlock_tac_s0tlock{ll+40,tt}=ft_timelockanalysis(cfg,tlock_tac_s0{ll+40,tt}); % Talone with shift in time matching SOA jitter
%       
%       cfg=[];
%       % I have thought it through, and no minus sign needed here.
%       %     cfg.offset=round(fsample*(tlock_aud_s0{ll,tt}.trialinfo(trialkept,10)-soades(ll))); % tlock_aud_s0{ll,tt}.trialinfo(trialkept,10) is identical to tlock_tac_s0{ll,tt}.trialinfo(trialkept,10)
%       cfg.offset=round(fsample*(tlock_aud_s0{ll,tt}.trialinfo(abb(mtrialkept{ll,tt}),10)-soades(ll) +soades(ll))); % tlock_aud_s0{ll,tt}.trialinfo(trialkept,10) is identical to tlock_tac_s0{ll,tt}.trialinfo(trialkept,10)
%       cfg.trials=atrialkept{ll,tt};
%       tlock_aud_s0{ll+40,tt}=ft_redefinetrial(cfg,tlock_aud_s0{10}); % Aalone with shift in time matching SOA jitter
%       cfg=[];
%       cfg.vartrllength=2;
%       tlock_aud_s0tlock{ll+40,tt}=ft_timelockanalysis(cfg,tlock_aud_s0{ll+40,tt}); % Aalone with shift in time matching SOA jitter
%       
%       cfg=[];
%       cfg.trials=ntrialkept{ll,tt};
%       tlock_nul_s0{ll+40,tt}=ft_redefinetrial(cfg,tlock_nul_s0{10});


      % 1-40Hz bandpass filter and average

      % create sum of unisensory conditions
      cfg=[];
      cfg.operation='add';
      cfg.parameter='avg';
      % the call to ft_timelockanalysis is b/c .avg field not present after ft_redefinetrial
      tlock_tacPaud_s0{ll,tt}=ft_math(cfg,tlock_tac_s0{10},tlock_aud_s0tlock{ll+40,tt}); % keeping tac centred but jitter aud (to match tlock_tac_s0{X,tt})
      tlock_audPtac_s0{ll,tt}=ft_math(cfg,tlock_aud_s0{10},tlock_tac_s0tlock{ll+40,tt}); % keeping aud centred but jitter tac (to match tlock_aud_s0{X,tt})
      
      % create sum of multisensory plus nul conditions
      cfg=[];
      cfg.operation='add';
      cfg.parameter='avg';
      tlock_tacMSpN_s0{ll,tt}=ft_math(cfg,tlock_tac_s0{ll,tt},ft_timelockanalysis([],tlock_nul_s0{ll+40,tt})); % TA+N
      tlock_audMSpN_s0{ll,tt}=ft_math(cfg,tlock_aud_s0{ll,tt},ft_timelockanalysis([],tlock_nul_s0{ll+40,tt})); % AT+N
      
      %       cfg=[];
      %       cfg.operation='subtract';
      %       cfg.parameter='avg';
      %       tlock_tpa_mtamn_s0{ll,tt}=ft_math(cfg,tlock_tacPaud_s0{ll,tt},tlock_tacMSpN_s0{ll,tt}); % (T + As) - (TA + N)
      %       tlock_apt_matmn_s0{ll,tt}=ft_math(cfg,tlock_audPtac_s0{ll,tt},tlock_audMSpN_s0{ll,tt}); % (A + Ts) - (AT + N)
    end
  end
  switch ref
    case 'average'
      save(['tlock_diffs_averef_' sub{ii} '.mat'],'*trialkept','tlock_*P*','tlock_*N*')
    case 'ft'
      save(['tlock_diffs_' sub{ii} '.mat'],'*trialkept','tlock_*P*','tlock_*N*')
  end
  
  
  
end

return

%%  Individual t-score

time=-0.1:0.001:0.4;

ccc=nan(length(time),63,9,4,28);

for ii=5:18
  cd([edir sub{ii} ])
  load(['tlock_diffs_' sub{ii} '.mat']);
  if ii>7
    soalist=[1 3 4 5 6 7 9];
  else
    soalist=[3 4 5 6 7];
  end
  
  if ii==5
    labels=tlock_tac_s0{5,1}.label;
  end
  
  chanuse=match_str(labels,tlock_tac_s0{5,1}.label);
  
  for ll=soalist
    for tt=1:4
      
      
      for teatime=1:length(time)
        
        audtac=tlock_tac_s0{ll,tt}.trial(:,:,dsearchn(tlock_tac_s0{ll,tt}.time',time(teatime)));
        nul=tlock_nul_s0{ll+40,tt}.trial(:,:,dsearchn(tlock_nul_s0{ll+40,tt}.time',time(teatime)));
        tac=tlock_tac_s0{10,tt}.trial(ttrialkept{ll,tt},:,dsearchn(tlock_tac_s0{10,tt}.time',time(teatime)));
        aud=tlock_aud_s0{ll+40,tt}.trial(:,:,dsearchn(tlock_aud_s0{ll+40,tt}.time',time(teatime)));
        
        
        data=[tac; aud; audtac; nul];
        numtr=size(tac,1);
        design=zeros(size(data,1),5);
        design(:,1)=[zeros(0*numtr,1); ones(numtr,1); zeros(3*numtr,1)];
        design(:,2)=[zeros(1*numtr,1); ones(numtr,1); zeros(2*numtr,1)];
        design(:,3)=[zeros(2*numtr,1); ones(numtr,1); zeros(1*numtr,1)];
        design(:,4)=[zeros(3*numtr,1); ones(numtr,1); zeros(0*numtr,1)];
        design(:,5)=[zeros(0*numtr,1); ones(4*numtr,1); zeros(0*numtr,1)];
        cfg=[];
        cfg.glm.statistic='beta';
        cfg.glm.standardise=0;
        stat=ft_statfun_glm(cfg,data',design');
        beta=reshape(stat.stat,[size(design,2) size(data,2)]);
        con(teatime,chanuse,ll,tt,ii)=[1 1 -1 -1 0]*beta;
      end
    end
  end
  
  
  numtr=size(tlock_tac_s0{ll,tt}.trialinfo,1);
  mu=repmat(nanmean(tlock_tac_s0{ll,tt}.trial,1),[numtr 1 1]);
  sigma=repmat((nansum((tlock_tac_s0{ll,tt}.trial-mu).^2,1)/numtr).^0.5,[numtr 1 1]);
  zvalue=(tlock_tac_s0{ll,tt}.trial-mu)./sigma;
  zvalue=tlock_tac_s0{ll,tt}.trial./sigma;
  zvalue=squeeze(mu(1,:,:)./sigma(1,:,:));
  
  
  
  
  
  
  
  cfg.design=[];
  cfg.design(:,1)=[ ones(1,size(tlock_aud_s0{10}.trial,1)) ones(1,size(tlock_tac_s0{30+ll}.trial,1)) 2*ones(1,size(tlock_tac_s0{ll}.trial,1)) 2*ones(1,size(tlock_nul_s0{10}.trial,1))]';
  %   cfg.design(:,2)=[ ones(1,size(tlock_aud{10}.trial,1)) 2*ones(1,size(tlock_tac{30+ll}.trial,1)) 3*ones(1,size(tlock_tac{ll}.trial,1)) 4*ones(1,size(tlock_nul{10}.trial,1))]';
  stat1_s0{ll} = ft_timelockstatistics(cfg, tlock_aud_s0{10}, tlock_tac_s0{30+ll}, tlock_tac_s0{ll}, tlock_nul_s0{10})
  tlock_tpa_mtamn_s0{ll}.mask=stat1_s0{ll}.mask;
  
  cfg.design=[];
  cfg.design(:,1)=[ ones(1,size(tlock_aud_sall{10}.trial,1)) ones(1,size(tlock_tac_sall{30+ll}.trial,1)) 2*ones(1,size(tlock_tac_sall{ll}.trial,1)) 2*ones(1,size(tlock_nul_sall{10}.trial,1))]';
  %   cfg.design(:,2)=[ ones(1,size(tlock_aud{10}.trial,1)) 2*ones(1,size(tlock_tac{30+ll}.trial,1)) 3*ones(1,size(tlock_tac{ll}.trial,1)) 4*ones(1,size(tlock_nul{10}.trial,1))]';
  stat1_sall{ll} = ft_timelockstatistics(cfg, tlock_aud_sall{10}, tlock_tac_sall{30+ll}, tlock_tac_sall{ll}, tlock_nul_sall{10})
  tlock_tpa_mtamn_sall{ll}.mask=stat1_sall{ll}.mask;
end

% end
% thus, stat1_* is with Aud-alone at time zero, shifted tac, and AT with aud-first for 3, and tac-first for 7

for ll=soalist
  cfg.design=[];
  cfg.design(:,1)=[ ones(1,size(tlock_tac_s0{10}.trial,1)) ones(1,size(tlock_aud_s0{30+ll}.trial,1)) 2*ones(1,size(tlock_aud_s0{ll}.trial,1)) 2*ones(1,size(tlock_nul_s0{10}.trial,1))]';
  stat2_s0{ll} = ft_timelockstatistics(cfg, tlock_tac_s0{10}, tlock_aud_s0{30+ll}, tlock_aud_s0{ll}, tlock_nul_s0{10})
  
  cfg.design=[];
  cfg.design(:,1)=[ ones(1,size(tlock_tac_sall{10}.trial,1)) ones(1,size(tlock_aud_sall{30+ll}.trial,1)) 2*ones(1,size(tlock_aud_sall{ll}.trial,1)) 2*ones(1,size(tlock_nul_sall{10}.trial,1))]';
  stat2_sall{ll} = ft_timelockstatistics(cfg, tlock_tac_sall{10}, tlock_aud_sall{30+ll}, tlock_aud_sall{ll}, tlock_nul_sall{10})
end
% thus, stat2_* is with tac-alone at time zero, shifted aud, and AT with aud-first for 3, and tac-first for 7

save(['stat_erp_' sub{ii} '.mat'],'stat*')

figure;for ll=soalist,subplot(1,7,ll);imagesc(stat1_sall{ll}.time,1:62,stat1_sall{ll}.mask);end
figure;for ll=soalist,subplot(1,7,ll);imagesc(stat2_sall{ll}.time,1:62,stat2_sall{ll}.mask);end
figure;for ll=soalist,subplot(1,7,ll);imagesc(stat1_s0{ll}.time,1:62,stat1_s0{ll}.mask);end
figure;for ll=soalist,subplot(1,7,ll);imagesc(stat2_s0{ll}.time,1:62,stat2_s0{ll}.mask);end

cfg=[];
cfg.interactive='yes';
cfg.layout='EEG1010.lay';
ft_multiplotER(cfg,tlock_tpa_mtamn_sall{7});
ft_multiplotER(cfg,tlock_tpa_mtamn_sall{5});


clear *_s0
% end





%% Group level

for ll=[1 3 4 5 6 7 9]
  for tt=1:4
    
    if ll==1 | ll==9
      subuse=8:18;
      submin=7;
    else
      subuse=5:18;
      submin=4;
    end
    for ii=subuse
      cd([edir sub{ii} ])
      %   load(['tlock_diffs_' sub{ii} '.mat']);
      load(['tlock_diffs_averef_' sub{ii} '.mat']);
      tlock_tacPaud{ii-submin}=tlock_tacPaud_s0{ll,tt};
      tlock_audPtac{ii-submin}=tlock_audPtac_s0{ll,tt};
      tlock_tacMSpN{ii-submin}=tlock_tacMSpN_s0{ll,tt};
      tlock_audMSpN{ii-submin}=tlock_audMSpN_s0{ll,tt};
      clear *_s0
    end
    
    figure(20);
    for ii=1:length(subuse)
      cfg=[];
      cfg.parameter='avg';
      cfg.operation='subtract';
      diff{ii}=ft_math(cfg,tlock_tacPaud{ii},tlock_tacMSpN{ii});
      subplot(3,length(subuse),ii);imagesc(tlock_tacPaud{ii}.time,1:63,tlock_tacPaud{ii}.avg);caxis([-6 6])
      subplot(3,length(subuse),length(subuse)+ii);imagesc(tlock_tacMSpN{ii}.time,1:63,tlock_tacMSpN{ii}.avg);caxis([-6 6])
      subplot(3,length(subuse),2*length(subuse)+ii);imagesc(diff{ii}.time,1:63,diff{ii}.avg);caxis([-6 6])
    end
    print(20,['D:\audtac\figs\sumuni_mspn_diff_ta_cond' num2str(ll) num2str(tt) 'a.png'],'-dpng')
    figure(21);
    for ii=1:length(subuse)
      cfg=[];
      cfg.parameter='avg';
      cfg.operation='subtract';
      diff{ii}=ft_math(cfg,tlock_audPtac{ii},tlock_audMSpN{ii});
      subplot(3,length(subuse),ii);imagesc(tlock_audPtac{ii}.time,1:63,tlock_audPtac{ii}.avg);caxis([-6 6])
      subplot(3,length(subuse),length(subuse)+ii);imagesc(tlock_audMSpN{ii}.time,1:63,tlock_audMSpN{ii}.avg);caxis([-6 6])
      subplot(3,length(subuse),2*length(subuse)+ii);imagesc(diff{ii}.time,1:63,diff{ii}.avg);caxis([-6 6])
    end
    print(21,['D:\audtac\figs\sumuni_mspn_diff_at_cond' num2str(ll) num2str(tt) 'a.png'],'-dpng')
    
    % cfg=[];
    % grave_diff=ft_timelockgrandaverage(cfg,diff{:});
    % figure;
    % cfg=[];
    % cfg.layout='elec1010.lay';
    % ft_multiplotER(cfg,grave_diff);
    
    cfg=[];
    cfg.keepindividual='yes';
    grind_tacPaud=ft_timelockgrandaverage(cfg,tlock_tacPaud{:});
    grind_audPtac=ft_timelockgrandaverage(cfg,tlock_audPtac{:});
    grind_tacMSpN=ft_timelockgrandaverage(cfg,tlock_tacMSpN{:});
    grind_audMSpN=ft_timelockgrandaverage(cfg,tlock_audMSpN{:});
    
    % cfg=[];
    % grave_tacPaud=ft_timelockgrandaverage(cfg,tlock_tacPaud{:});
    % grave_audPtac=ft_timelockgrandaverage(cfg,tlock_audPtac{:});
    % grave_tacMSpN=ft_timelockgrandaverage(cfg,tlock_tacMSpN{:});
    % grave_audMSpN=ft_timelockgrandaverage(cfg,tlock_audMSpN{:});
    
    for ii=1:length(subuse)
      cfg=[];
      cfg.operation='subtract';
      cfg.parameter='avg';
      tlock_TPA_MSPN{ii}=ft_math(cfg,tlock_tacPaud{ii},tlock_tacMSpN{ii});
      tlock_APT_MSPN{ii}=ft_math(cfg,tlock_audPtac{ii},tlock_audMSpN{ii});
    end
    cfg=[];
    grave_TPA_MSPN=ft_timelockgrandaverage(cfg,tlock_TPA_MSPN{:});
    grave_APT_MSPN=ft_timelockgrandaverage(cfg,tlock_APT_MSPN{:});
    
    
    % figure;
    % % define parameters for plotting
    % timestep = 0.025;      %(in seconds)
    % fsample=1/(tlock_tacMSpN{1}.time(2)-tlock_tacMSpN{1}.time(1));
    % sample_count = length(grave_TPA_MSPN.time);
    % j = [-0.1:timestep:0.5];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
    % m = [1:timestep*fsample:sample_count];  % temporal endpoints in MEEG samples
    % % plot
    % for k = 1:24;
    %      subplot(4,6,k);
    %      cfg = [];
    %      cfg.xlim=[j(k) j(k+1)];
    %      cfg.zlim = [-1.5 1.5];
    % %      pos_int = all(pos(:, m(k):m(k+1)), 2);
    % %      cfg.highlight = 'on';
    % %      cfg.highlightchannel = find(pos_int);
    %      cfg.comment = 'xlim';
    %      cfg.commentpos = 'title';
    %      cfg.layout = 'elec1010.lay';
    %      ft_topoplotER(cfg, grave_TPA_MSPN);
    % end
    %
    % figure;
    % % define parameters for plotting
    % timestep = 0.025;      %(in seconds)
    % fsample=1/(tlock_tacMSpN{1}.time(2)-tlock_tacMSpN{1}.time(1));
    % sample_count = length(grave_TAmMS.time);
    % j = [-0.1:timestep:0.5];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
    % m = [1:timestep*fsample:sample_count];  % temporal endpoints in MEEG samples
    % % plot
    % for k = 1:24;
    %      subplot(4,6,k);
    %      cfg = [];
    %      cfg.xlim=[j(k) j(k+1)];
    %      cfg.zlim = [-1.5 1.5];
    % %      pos_int = all(pos(:, m(k):m(k+1)), 2);
    % %      cfg.highlight = 'on';
    % %      cfg.highlightchannel = find(pos_int);
    %      cfg.comment = 'xlim';
    %      cfg.commentpos = 'title';
    %      cfg.layout = 'elec1010.lay';
    %      ft_topoplotER(cfg, grave_APT_MSPN);
    % end
    %
    
    load eeg1010_neighb
    
    nsub=length(tlock_audMSpN);
    
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
    statt_mc=ft_timelockstatistics(cfg, grind_tacPaud, grind_tacMSpN);
    stata_mc=ft_timelockstatistics(cfg, grind_audPtac, grind_audMSpN);
    
    
    % get relevant (significant) values
    pos_cluster_pvals = [statt_mc.posclusters(:).prob];
    % pos_signif_clust = find(pos_cluster_pvals < statt.cfg.alpha);
    pos_signif_clust = find(pos_cluster_pvals < 0.05);
    pos = ismember(statt_mc.posclusterslabelmat, pos_signif_clust);
    % define parameters for plotting
    figure(22);
    timestep = 0.025;      %(in seconds)
    fsample=1/(tlock_tacMSpN{1}.time(2)-tlock_tacMSpN{1}.time(1));
    sample_count = length(grave_TPA_MSPN.time);
    j = [-0.1:timestep:0.5];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
    m = [1:timestep*fsample:sample_count];  % temporal endpoints in MEEG samples
    % plot
    for k = 1:24;
      subplot(4,6,k);
      cfg = [];
      cfg.xlim=[j(k) j(k+1)];
      cfg.zlim = [-1.5 1.5];
      %      pos_int = all(pos(:, m(k):m(k+1)), 2);
      pos_int = any(pos(:, m(k):m(k+1)), 2);
      cfg.highlight = 'on';
      cfg.highlightchannel = find(pos_int)
      %      keyboard
      cfg.comment = 'xlim';
      cfg.commentpos = 'title';
      cfg.layout = 'elec1010.lay';
      ft_topoplotER(cfg, grave_TPA_MSPN);
    end
    print(22,['D:\audtac\figs\grdiff_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) 'a.png'],'-dpng')
    
    % get relevant (significant) values
    pos_cluster_pvals = [stata_mc.posclusters(:).prob];
    % pos_signif_clust = find(pos_cluster_pvals < statt.cfg.alpha);
    pos_signif_clust = find(pos_cluster_pvals < 0.05);
    pos = ismember(stata_mc.posclusterslabelmat, pos_signif_clust);
    % define parameters for plotting
    figure(23);
    timestep = 0.025;      %(in seconds)
    fsample=1/(tlock_tacMSpN{1}.time(2)-tlock_tacMSpN{1}.time(1));
    sample_count = length(grave_APT_MSPN.time);
    j = [-0.1:timestep:0.5];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
    m = [1:timestep*fsample:sample_count];  % temporal endpoints in MEEG samples
    % plot
    for k = 1:24;
      subplot(4,6,k);
      cfg = [];
      cfg.xlim=[j(k) j(k+1)];
      cfg.zlim = [-1.5 1.5];
      %      pos_int = all(pos(:, m(k):m(k+1)), 2);
      pos_int = any(pos(:, m(k):m(k+1)), 2);
      cfg.highlight = 'on';
      cfg.highlightchannel = find(pos_int)
      %      keyboard
      cfg.comment = 'xlim';
      cfg.commentpos = 'title';
      cfg.layout = 'elec1010.lay';
      ft_topoplotER(cfg, grave_APT_MSPN);
    end
    print(23,['D:\audtac\figs\grdiff_topoOverTime_mc_at_cond' num2str(ll) num2str(tt) 'a.png'],'-dpng')
    
    
    cfg=[];
    cfg.latency=[-.1 .5];
    % cfg.neighbours=neighbours;
    % cfg.parameter='avg';
    cfg.parameter='individual';
    % cfg.method='montecarlo';
    cfg.method='analytic';
    cfg.alpha=0.05;
    % cfg.numrandomization=200;
    cfg.correctm='fdr';
    % cfg.correctm='cluster';
    % cfg.clusteralpha = 0.05;
    % cfg.clusterstatistic = 'maxsum';
    % cfg.minnbchan = 2;
    cfg.statistic='depsamplesT';
    % cfg.statistic='indepsamplesregrT';
    % cfg.statistic='indepsamplesT';
    cfg.design=zeros(2,2*nsub);
    cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
    cfg.design(2,:)=[1:nsub 1:nsub];
    cfg.ivar=1;
    cfg.uvar=2;
    statt_an=ft_timelockstatistics(cfg, grind_tacPaud, grind_tacMSpN);
    stata_an=ft_timelockstatistics(cfg, grind_audPtac, grind_audMSpN);
    
    % get relevant (significant) values
    pos = statt_an.mask;
    % define parameters for plotting
    figure(24);
    timestep = 0.025;      %(in seconds)
    fsample=1/(tlock_tacMSpN{1}.time(2)-tlock_tacMSpN{1}.time(1));
    sample_count = length(grave_TPA_MSPN.time);
    j = [-0.1:timestep:0.5];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
    m = [1:timestep*fsample:sample_count];  % temporal endpoints in MEEG samples
    % plot
    for k = 1:24;
      subplot(4,6,k);
      cfg = [];
      cfg.xlim=[j(k) j(k+1)];
      cfg.zlim = [-1.5 1.5];
      %      pos_int = all(pos(:, m(k):m(k+1)), 2);
      pos_int = any(pos(:, m(k):m(k+1)), 2);
      cfg.highlight = 'on';
      cfg.highlightchannel = find(pos_int);
      %      keyboard
      cfg.comment = 'xlim';
      cfg.commentpos = 'title';
      cfg.layout = 'elec1010.lay';
      ft_topoplotER(cfg, grave_TPA_MSPN);
    end
    print(24,['D:\audtac\figs\grdiff_topoOverTime_an_ta_cond' num2str(ll) num2str(tt) 'a.png'],'-dpng')
    
    % get relevant (significant) values
    pos = stata_an.mask;
    % define parameters for plotting
    figure(25);
    timestep = 0.025;      %(in seconds)
    fsample=1/(tlock_tacMSpN{1}.time(2)-tlock_tacMSpN{1}.time(1));
    sample_count = length(grave_APT_MSPN.time);
    j = [-0.1:timestep:0.5];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
    m = [1:timestep*fsample:sample_count];  % temporal endpoints in MEEG samples
    % plot
    for k = 1:24;
      subplot(4,6,k);
      cfg = [];
      cfg.xlim=[j(k) j(k+1)];
      cfg.zlim = [-1.5 1.5];
      %      pos_int = all(pos(:, m(k):m(k+1)), 2);
      pos_int = any(pos(:, m(k):m(k+1)), 2);
      cfg.highlight = 'on';
      cfg.highlightchannel = find(pos_int);
      %      keyboard
      cfg.comment = 'xlim';
      cfg.commentpos = 'title';
      cfg.layout = 'elec1010.lay';
      ft_topoplotER(cfg, grave_APT_MSPN);
    end
    print(25,['D:\audtac\figs\grdiff_topoOverTime_an_at_cond' num2str(ll) num2str(tt) 'a.png'],'-dpng')
    
  end
end

%%


%% Across subjects combining stats

for ii=2:4
  cd([edir sub{ii} ])
  stat{ii}=load(['stat_erp_' sub{ii} '.mat']);
  mmsu{ii}=load(['tlock_diffs_' sub{ii} '.mat']);
  if ii==2  % find channels in common to all subjects
    labelkeep=stat{ii}.stat1_s0{3}.label;
  else
    labelkeep=intersect(labelkeep,stat{ii}.stat1_s0{3}.label);
  end
end

chanuse=match_str(labelkeep,'Fz');


for ii=2:4
  for ll=3:7
    mask1_s0(:,:,ll,ii)=stat{ii}.stat1_s0{ll}.mask(match_str(stat{ii}.stat1_s0{ll}.label,labelkeep),:);
    mask2_s0(:,:,ll,ii)=stat{ii}.stat2_s0{ll}.mask(match_str(stat{ii}.stat2_s0{ll}.label,labelkeep),:);
    mask1_sall(:,:,ll,ii)=stat{ii}.stat1_sall{ll}.mask(match_str(stat{ii}.stat1_sall{ll}.label,labelkeep),:);
    mask2_sall(:,:,ll,ii)=stat{ii}.stat2_sall{ll}.mask(match_str(stat{ii}.stat2_sall{ll}.label,labelkeep),:);
    stat1_s0(:,:,ll,ii)=stat{ii}.stat1_s0{ll}.stat(match_str(stat{ii}.stat1_s0{ll}.label,labelkeep),:);
    stat2_s0(:,:,ll,ii)=stat{ii}.stat2_s0{ll}.stat(match_str(stat{ii}.stat2_s0{ll}.label,labelkeep),:);
    stat1_sall(:,:,ll,ii)=stat{ii}.stat1_sall{ll}.stat(match_str(stat{ii}.stat1_sall{ll}.label,labelkeep),:);
    stat2_sall(:,:,ll,ii)=stat{ii}.stat2_sall{ll}.stat(match_str(stat{ii}.stat2_sall{ll}.label,labelkeep),:);
    
    tpa_s0(:,:,ll,ii)=mmsu{ii}.tlock_tpa_mtamn_s0{ll}.avg(match_str(stat{ii}.stat1_s0{ll}.label,labelkeep),dsearchn(mmsu{ii}.tlock_tpa_mtamn_s0{ll}.time',-.2):dsearchn(mmsu{ii}.tlock_tpa_mtamn_s0{ll}.time',.5));
    tpa_sall(:,:,ll,ii)=mmsu{ii}.tlock_tpa_mtamn_sall{ll}.avg(match_str(stat{ii}.stat1_s0{ll}.label,labelkeep),dsearchn(mmsu{ii}.tlock_tpa_mtamn_sall{ll}.time',-.2):dsearchn(mmsu{ii}.tlock_tpa_mtamn_sall{ll}.time',.5));
    apt_s0(:,:,ll,ii)=mmsu{ii}.tlock_apt_matmn_s0{ll}.avg(match_str(stat{ii}.stat1_s0{ll}.label,labelkeep),dsearchn(mmsu{ii}.tlock_apt_matmn_s0{ll}.time',-.2):dsearchn(mmsu{ii}.tlock_apt_matmn_s0{ll}.time',.5));
    apt_sall(:,:,ll,ii)=mmsu{ii}.tlock_apt_matmn_sall{ll}.avg(match_str(stat{ii}.stat1_s0{ll}.label,labelkeep),dsearchn(mmsu{ii}.tlock_apt_matmn_sall{ll}.time',-.2):dsearchn(mmsu{ii}.tlock_apt_matmn_sall{ll}.time',.5));
    
  end
end

figure;
for ll=3:7
  subplot(1,5,ll-2);imagesc((nanmean(mask2_sall(:,:,ll,2:4),4)));caxis([0 1]);
end
figure;
for ll=3:7
  subplot(1,5,ll-2);imagesc((nanmean(mask1_sall(:,:,ll,2:4),4)));caxis([0 1]);
end
figure;
for ll=3:7
  subplot(1,5,ll-2);imagesc((nanmean(stat2_sall(:,:,ll,2:4),4)));caxis([-3 3]);
end
figure;
for ll=3:7
  subplot(1,5,ll-2);imagesc((nanmean(stat1_sall(:,:,ll,2:4),4)));caxis([-3 3]);
end

figure;for ll=3:7,subplot(1,5,ll-2);imagesc(-.2:.001:.5,1:length(labelkeep),mean(tpa_s0(:,:,ll,2:4),4));caxis([-4 4]);end
figure;for ll=3:7,subplot(1,5,ll-2);imagesc(-.2:.001:.5,1:length(labelkeep),mean(tpa_sall(:,:,ll,2:4),4));caxis([-4 4]);end
figure;for ll=3:7,subplot(1,5,ll-2);imagesc(-.2:.001:.5,1:length(labelkeep),mean(apt_s0(:,:,ll,2:4),4));caxis([-4 4]);end
figure;for ll=3:7,subplot(1,5,ll-2);imagesc(-.2:.001:.5,1:length(labelkeep),mean(apt_sall(:,:,ll,2:4),4));caxis([-4 4]);end


figure(111); % (Sum of Unisensory) minus (Multisensory plus Nul)
for ll=3:7
  subplot(2,5,ll-2);plot(-.2:.001:.5,mean(tpa_sall(chanuse,:,ll,2:4),4),'k');axis([-.2 .5 -7 7]);
  legend('(T+A)-(TA-N), time0 is tactile')
  subplot(2,5,ll-2+5);plot(-.2:.001:.5,mean(apt_sall(chanuse,:,ll,2:4),4),'k');axis([-.2 .5 -7 7])
  legend('(A+T)-(AT-N), time0 is auditory')
end
