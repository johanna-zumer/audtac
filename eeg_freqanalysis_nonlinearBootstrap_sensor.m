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
for ii=5:32
  
  cd([edir sub{ii} ])
  clearvars -except ii sub edir ddir
  
  % trial selection
  %   [tlock_tac_s0,tlock_aud_s0,tlock_nul_s0]=eeg_legomagic_trialSelection(ii);
  [tlock_tac_s0,tlock_aud_s0,tlock_nul_s0,tr]=eeg_legomagic_trialSelection1(ii);
  
  %  Frequency analysis
  if ii<8
    soalist=[3 4 5 6 7];
  else
    soalist=[1 3 4 5 6 7 9];
  end
  
  %   soanew=[soalist 10 soalist+40];
  
  for tt=1:4
    for ll=soalist
      cfg=[];
      cfg.trials=tr.t10trialkept{ll,tt};
      tlock_tac_s0{ll+20,tt}=ft_selectdata(cfg,tlock_tac_s0{10}); % create +20 here as 10 will be different for every ll,tt
      cfg=[];
      cfg.trials=tr.nllttrialkept{ll,tt};
      tlock_nul_s0{ll+50,tt}=ft_selectdata(cfg,tlock_nul_s0{10});
      
      cfg=[];
      cfg.trials=tr.a10trialkept{ll,tt};
      tlock_aud_s0{ll+20,tt}=ft_selectdata(cfg,tlock_aud_s0{10}); % create +20 here as 10 will be different for every ll,tt
      cfg=[];
      cfg.trials=tr.nllatrialkept{ll,tt};
      tlock_nul_s0{ll+60,tt}=ft_selectdata(cfg,tlock_nul_s0{10});
      numt_trials(ll,tt)=size(tlock_tac_s0{ll+20,tt}.trial,1); % this should be the same for tacll20, audll40, tacll, nulll50
      numa_trials(ll,tt)=size(tlock_aud_s0{ll+20,tt}.trial,1); % this should be the same for audll20, tacll40, audll, nulll60
    end
    
    % Inspired by Senkowski 2007 Exp Brain Res.
    for ll=soalist
      % tacPaud
      numcomb=numt_trials(ll,tt)^2;
      numtests=min(numcomb,1000);
      combindex=reshape(1:numcomb,numt_trials(ll,tt),numt_trials(ll,tt));
      for cc=1:numtests
        tmp=Shuffle(combindex(:));
        combuse=tmp(1:numt_trials(ll,tt));
        
        tlock_fake=tlock_tac_s0{ll+20,tt};
        for at=1:numt_trials(ll,tt)
          [tind,aind]=find(combindex==combuse(at));
          tlock_fake.trial(at,:,:)=tlock_tac_s0{ll+20,tt}.trial(tind,:,:)+tlock_aud_s0{ll+40,tt}.trial(aind,:,:);
        end
        
        cfg=[];
        cfg.method='mtmconvol';
        cfg.pad=4;
        cfg.foi=4:2:30;
        cfg.taper='hanning';
        cfg.toi=-1.2:0.1:1.3;
        cfg.t_ftimwin=4./cfg.foi;
        %         cfg.keeptrials='yes';  % not necessary at this stage for keeping trials (yes later for source analysis and phase resetting)
        freqlo_tacPaud_comb{ll,tt,cc}=ft_freqanalysis(cfg,tlock_fake);
        cfg=[];
        cfg.method='mtmconvol';
        cfg.pad=4;
        cfg.foi=[30:5:45 55:5:80];
        cfg.taper='dpss';
        cfg.tapsmofrq=7*ones(1,length(cfg.foi));
        cfg.toi=-1.2:0.1:1.3;
        cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
        freqhi_tacPaud_comb{ll,tt,cc}=ft_freqanalysis(cfg,tlock_fake);
      end
      
      % audPtac
      numcomb=numa_trials(ll,tt)^2;
      numtests=min(numcomb,1000);
      combindex=reshape(1:numcomb,numa_trials(ll,tt),numa_trials(ll,tt));
      for cc=1:numtests
        tmp=Shuffle(combindex(:));
        combuse=tmp(1:numa_trials(ll,tt));
        
        tlock_fake=tlock_aud_s0{ll+20,tt};
        for at=1:numa_trials(ll,tt)
          [aind,tind]=find(combindex==combuse(at));
          tlock_fake.trial(at,:,:)=tlock_aud_s0{ll+20,tt}.trial(aind,:,:)+tlock_tac_s0{ll+40,tt}.trial(tind,:,:);
        end
        
        cfg=[];
        cfg.method='mtmconvol';
        cfg.pad=4;
        cfg.foi=4:2:30;
        cfg.taper='hanning';
        cfg.toi=-1.2:0.1:1.3;
        cfg.t_ftimwin=4./cfg.foi;
        %         cfg.keeptrials='yes';  % not necessary at this stage for keeping trials (yes later for source analysis and phase resetting)
        freqlo_audPtad_comb{ll,tt,cc}=ft_freqanalysis(cfg,tlock_fake);
        cfg=[];
        cfg.method='mtmconvol';
        cfg.pad=4;
        cfg.foi=[30:5:45 55:5:80];
        cfg.taper='dpss';
        cfg.tapsmofrq=7*ones(1,length(cfg.foi));
        cfg.toi=-1.2:0.1:1.3;
        cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
        freqhi_audPtac_comb{ll,tt,cc}=ft_freqanalysis(cfg,tlock_fake);
      end
      
      
      % tacMSpN
      numcomb=numt_trials(ll,tt)^2;
      numtests=min(numcomb,1000);
      combindex=reshape(1:numcomb,numt_trials(ll,tt),numt_trials(ll,tt));
      for cc=1:numtests
        tmp=Shuffle(combindex(:));
        combuse=tmp(1:numt_trials(ll,tt));
        
        tlock_fake=tlock_tac_s0{ll,tt};
        for at=1:numt_trials(ll,tt)
          [tind,aind]=find(combindex==combuse(at));
          tlock_fake.trial(at,:,:)=tlock_tac_s0{ll,tt}.trial(tind,:,:)+tlock_nul_s0{ll+50,tt}.trial(aind,:,:);
        end
        cfg=[];
        cfg.method='mtmconvol';
        cfg.pad=4;
        cfg.foi=4:2:30;
        cfg.taper='hanning';
        cfg.toi=-1.2:0.1:1.3;
        cfg.t_ftimwin=4./cfg.foi;
        %         cfg.keeptrials='yes';  % not necessary at this stage for keeping trials (yes later for source analysis and phase resetting)
        freqlo_tacMSpN_comb{ll,tt,cc}=ft_freqanalysis(cfg,tlock_fake);
        cfg=[];
        cfg.method='mtmconvol';
        cfg.pad=4;
        cfg.foi=[30:5:45 55:5:80];
        cfg.taper='dpss';
        cfg.tapsmofrq=7*ones(1,length(cfg.foi));
        cfg.toi=-1.2:0.1:1.3;
        cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
        freqhi_tacMSpN_comb{ll,tt,cc}=ft_freqanalysis(cfg,tlock_fake);
      end
      
      
      % audMSpN
      numcomb=numa_trials(ll,tt)^2;
      numtests=min(numcomb,1000);
      combindex=reshape(1:numcomb,numa_trials(ll,tt),numa_trials(ll,tt));
      for cc=1:numtests
        tmp=Shuffle(combindex(:));
        combuse=tmp(1:numa_trials(ll,tt));
        
        tlock_fake=tlock_aud_s0{ll,tt};
        for at=1:numa_trials(ll,tt)
          [aind,tind]=find(combindex==combuse(at));
          tlock_fake.trial(at,:,:)=tlock_aud_s0{ll,tt}.trial(aind,:,:)+tlock_nul_s0{ll+60,tt}.trial(tind,:,:);
        end
        
        cfg=[];
        cfg.method='mtmconvol';
        cfg.pad=4;
        cfg.foi=4:2:30;
        cfg.taper='hanning';
        cfg.toi=-1.2:0.1:1.3;
        cfg.t_ftimwin=4./cfg.foi;
        %         cfg.keeptrials='yes';  % not necessary at this stage for keeping trials (yes later for source analysis and phase resetting)
        freqlo_audMSpN_comb{ll,tt,cc}=ft_freqanalysis(cfg,tlock_fake);
        cfg=[];
        cfg.method='mtmconvol';
        cfg.pad=4;
        cfg.foi=[30:5:45 55:5:80];
        cfg.taper='dpss';
        cfg.tapsmofrq=7*ones(1,length(cfg.foi));
        cfg.toi=-1.2:0.1:1.3;
        cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
        freqhi_audMSpN_comb{ll,tt,cc}=ft_freqanalysis(cfg,tlock_fake);
      end
    end
    
    
    
    
  end
  
  % Question: do baseline correction prior to enter in to stats?
  
  save(['freq_diffs_averef_comb_' sub{ii} '.mat'],'freq*P*','freq*N*','num*trials')
  
end


return



%% sensor level group stats

load elec1010_neighb.mat

allcond_sameN=1;


for ll=[1 3 4 5 6 7 9]
  for tt=1:4
    
    if ll==1 | ll==9
      subuse=8:32;
    else
      subuse=5:32;
    end
    if allcond_sameN
      subuse=8:32;
    end
    submin=subuse(1)-1;
    % Baseline correct each participant prior to entering to stats????
    for ii=subuse
      cd([edir sub{ii} ])
      load(['freq_diffs_averef_' sub{ii} '.mat']);
      freqloall_tacPaud{ii-submin}=freqlo_tacPaud{ll,tt};
      freqloall_audPtac{ii-submin}=freqlo_audPtac{ll,tt};
      freqloall_tacMSpN{ii-submin}=freqlo_tacMSpN{ll,tt};
      freqloall_audMSpN{ii-submin}=freqlo_audMSpN{ll,tt};
      freqhiall_tacPaud{ii-submin}=freqhi_tacPaud{ll,tt};
      freqhiall_audPtac{ii-submin}=freqhi_audPtac{ll,tt};
      freqhiall_tacMSpN{ii-submin}=freqhi_tacMSpN{ll,tt};
      freqhiall_audMSpN{ii-submin}=freqhi_audMSpN{ll,tt};
      clear freqlo_* freqhi_*
    end
    
    figure(20);
    for ii=1:length(subuse)
      cfg=[];
      cfg.parameter='powspctrm';
      cfg.operation='subtract';
      diff{ii}=ft_math(cfg,freqloall_tacPaud{ii},freqloall_tacMSpN{ii});
      subplot(3,length(subuse),ii);imagesc(freqloall_tacPaud{ii}.time,freqloall_tacPaud{ii}.freq,squeeze(mean(freqloall_tacPaud{ii}.powspctrm,1)));caxis([-6 6])
      subplot(3,length(subuse),length(subuse)+ii);imagesc(freqloall_tacMSpN{ii}.time,freqloall_tacPaud{ii}.freq,squeeze(mean(freqloall_tacMSpN{ii}.powspctrm,1)));caxis([-6 6])
      subplot(3,length(subuse),2*length(subuse)+ii);imagesc(diff{ii}.time,freqloall_tacPaud{ii}.freq,squeeze(mean(diff{ii}.powspctrm,1)));caxis([-6 6])
    end
    if allcond_sameN
      print(20,['D:\audtac\figs\sumuni_mspn_fl_diff_ta_cond' num2str(ll) num2str(tt) 'a25.png'],'-dpng')
    else
      print(20,['D:\audtac\figs\sumuni_mspn_fl_diff_ta_cond' num2str(ll) num2str(tt) 'a.png'],'-dpng')
    end
    figure(21);
    for ii=1:length(subuse)
      cfg=[];
      cfg.parameter='powspctrm';
      cfg.operation='subtract';
      diff{ii}=ft_math(cfg,freqloall_audPtac{ii},freqloall_audMSpN{ii});
      subplot(3,length(subuse),ii);imagesc(freqloall_audPtac{ii}.time,freqloall_tacPaud{ii}.freq,squeeze(mean(freqloall_audPtac{ii}.powspctrm,1)));caxis([-6 6])
      subplot(3,length(subuse),length(subuse)+ii);imagesc(freqloall_audMSpN{ii}.time,freqloall_tacPaud{ii}.freq,squeeze(mean(freqloall_audMSpN{ii}.powspctrm,1)));caxis([-6 6])
      subplot(3,length(subuse),2*length(subuse)+ii);imagesc(diff{ii}.time,freqloall_tacPaud{ii}.freq,squeeze(mean(diff{ii}.powspctrm,1)));caxis([-6 6])
    end
    if allcond_sameN
      print(21,['D:\audtac\figs\fsumuni_mspn_fl_diff_at_cond' num2str(ll) num2str(tt) 'a25.png'],'-dpng')
    else
      print(21,['D:\audtac\figs\fsumuni_mspn_fl_diff_at_cond' num2str(ll) num2str(tt) 'a.png'],'-dpng')
    end
    figure(22);
    for ii=1:length(subuse)
      cfg=[];
      cfg.parameter='powspctrm';
      cfg.operation='subtract';
      diff{ii}=ft_math(cfg,freqhiall_tacPaud{ii},freqhiall_tacMSpN{ii});
      subplot(3,length(subuse),ii);imagesc(freqhiall_tacPaud{ii}.time,freqloall_tacPaud{ii}.freq,squeeze(mean(freqhiall_tacPaud{ii}.powspctrm,1)));caxis([-1 1])
      subplot(3,length(subuse),length(subuse)+ii);imagesc(freqhiall_tacMSpN{ii}.time,freqloall_tacPaud{ii}.freq,squeeze(mean(freqhiall_tacMSpN{ii}.powspctrm,1)));caxis([-1 1])
      subplot(3,length(subuse),2*length(subuse)+ii);imagesc(diff{ii}.time,freqloall_tacPaud{ii}.freq,squeeze(mean(diff{ii}.powspctrm,1)));caxis([-.2 .2])
    end
    if allcond_sameN
      print(22,['D:\audtac\figs\fsumuni_mspn_fh_diff_ta_cond' num2str(ll) num2str(tt) 'a25.png'],'-dpng')
    else
      print(22,['D:\audtac\figs\fsumuni_mspn_fh_diff_ta_cond' num2str(ll) num2str(tt) 'a.png'],'-dpng')
    end
    figure(23);
    for ii=1:length(subuse)
      cfg=[];
      cfg.parameter='powspctrm';
      cfg.operation='subtract';
      diff{ii}=ft_math(cfg,freqhiall_audPtac{ii},freqhiall_audMSpN{ii});
      subplot(3,length(subuse),ii);imagesc(freqhiall_audPtac{ii}.time,freqloall_tacPaud{ii}.freq,squeeze(mean(freqhiall_audPtac{ii}.powspctrm,1)));caxis([-1 1])
      subplot(3,length(subuse),length(subuse)+ii);imagesc(freqhiall_audMSpN{ii}.time,freqloall_tacPaud{ii}.freq,squeeze(mean(freqhiall_audMSpN{ii}.powspctrm,1)));caxis([-1 1])
      subplot(3,length(subuse),2*length(subuse)+ii);imagesc(diff{ii}.time,freqloall_tacPaud{ii}.freq,squeeze(mean(diff{ii}.powspctrm,1)));caxis([-.2 .2])
    end
    if allcond_sameN
      print(23,['D:\audtac\figs\sumuni_mspn_fh_diff_at_cond' num2str(ll) num2str(tt) 'a25.png'],'-dpng')
    else
      print(23,['D:\audtac\figs\sumuni_mspn_fh_diff_at_cond' num2str(ll) num2str(tt) 'a.png'],'-dpng')
    end
    
    
    for ii=1:length(subuse)
      baselinetime=[dsearchn(freqloall_tacPaud_bn{ii}.time,-1.2) dsearchn(freqloall_tacPaud_bn{ii}.time,-0.7)];
      freqloall_tacPaud_bn{ii}=freqloall_tacPaud{ii};
      freqloall_tacPaud_bn{ii}.powspctrm=freqloall_tacPaud_bn{ii}.powspctrm./repmat(nanmean(freqloall_tacPaud_bn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud{ii}.time)]);
      freqloall_tacMSpN_bn{ii}=freqloall_tacMSpN{ii};
      freqloall_tacMSpN_bn{ii}.powspctrm=freqloall_tacMSpN_bn{ii}.powspctrm./repmat(nanmean(freqloall_tacMSpN_bn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud{ii}.time)]);
      freqloall_audPtac_bn{ii}=freqloall_audPtac{ii};
      freqloall_audPtac_bn{ii}.powspctrm=freqloall_audPtac_bn{ii}.powspctrm./repmat(nanmean(freqloall_audPtac_bn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud{ii}.time)]);
      freqloall_audMSpN_bn{ii}=freqloall_audMSpN{ii};
      freqloall_audMSpN_bn{ii}.powspctrm=freqloall_audMSpN_bn{ii}.powspctrm./repmat(nanmean(freqloall_audMSpN_bn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud{ii}.time)]);
      freqhiall_tacPaud_bn{ii}=freqhiall_tacPaud{ii};
      freqhiall_tacPaud_bn{ii}.powspctrm=freqhiall_tacPaud_bn{ii}.powspctrm./repmat(nanmean(freqhiall_tacPaud_bn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud{ii}.time)]);
      freqhiall_tacMSpN_bn{ii}=freqhiall_tacMSpN{ii};
      freqhiall_tacMSpN_bn{ii}.powspctrm=freqhiall_tacMSpN_bn{ii}.powspctrm./repmat(nanmean(freqhiall_tacMSpN_bn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud{ii}.time)]);
      freqhiall_audPtac_bn{ii}=freqhiall_audPtac{ii};
      freqhiall_audPtac_bn{ii}.powspctrm=freqhiall_audPtac_bn{ii}.powspctrm./repmat(nanmean(freqhiall_audPtac_bn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud{ii}.time)]);
      freqhiall_audMSpN_bn{ii}=freqhiall_audMSpN{ii};
      freqhiall_audMSpN_bn{ii}.powspctrm=freqhiall_audMSpN_bn{ii}.powspctrm./repmat(nanmean(freqhiall_audMSpN_bn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud{ii}.time)]);
    end
    
    cfg=[];
    cfg.keepindividual='yes';
    grindlo_tacPaud=ft_freqgrandaverage(cfg,freqloall_tacPaud{:});
    grindlo_tacMSpN=ft_freqgrandaverage(cfg,freqloall_tacMSpN{:});
    grindlo_audPtac=ft_freqgrandaverage(cfg,freqloall_audPtac{:});
    grindlo_audMSpN=ft_freqgrandaverage(cfg,freqloall_audMSpN{:});
    grindhi_tacPaud=ft_freqgrandaverage(cfg,freqhiall_tacPaud{:});
    grindhi_audPtac=ft_freqgrandaverage(cfg,freqhiall_audPtac{:});
    grindhi_tacMSpN=ft_freqgrandaverage(cfg,freqhiall_tacMSpN{:});
    grindhi_audMSpN=ft_freqgrandaverage(cfg,freqhiall_audMSpN{:});
    
    grindlo_tacPaud_bn=ft_freqgrandaverage(cfg,freqloall_tacPaud_bn{:});
    grindlo_tacMSpN_bn=ft_freqgrandaverage(cfg,freqloall_tacMSpN_bn{:});
    grindlo_audPtac_bn=ft_freqgrandaverage(cfg,freqloall_audPtac_bn{:});
    grindlo_audMSpN_bn=ft_freqgrandaverage(cfg,freqloall_audMSpN_bn{:});
    grindhi_tacPaud_bn=ft_freqgrandaverage(cfg,freqhiall_tacPaud_bn{:});
    grindhi_audPtac_bn=ft_freqgrandaverage(cfg,freqhiall_audPtac_bn{:});
    grindhi_tacMSpN_bn=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_bn{:});
    grindhi_audMSpN_bn=ft_freqgrandaverage(cfg,freqhiall_audMSpN_bn{:});
    
    for ii=1:length(subuse)
      cfg=[];
      cfg.operation='subtract';
      cfg.parameter='powspctrm';
      freqloall_TPA_MSPN{ii}=ft_math(cfg,freqloall_tacPaud{ii},freqloall_tacMSpN{ii});
      freqloall_APT_MSPN{ii}=ft_math(cfg,freqloall_audPtac{ii},freqloall_audMSpN{ii});
      freqhiall_TPA_MSPN{ii}=ft_math(cfg,freqhiall_tacPaud{ii},freqhiall_tacMSpN{ii});
      freqhiall_APT_MSPN{ii}=ft_math(cfg,freqhiall_audPtac{ii},freqhiall_audMSpN{ii});
      freqloall_TPA_MSPN_bn{ii}=ft_math(cfg,freqloall_tacPaud_bn{ii},freqloall_tacMSpN_bn{ii});
      freqloall_APT_MSPN_bn{ii}=ft_math(cfg,freqloall_audPtac_bn{ii},freqloall_audMSpN_bn{ii});
      freqhiall_TPA_MSPN_bn{ii}=ft_math(cfg,freqhiall_tacPaud_bn{ii},freqhiall_tacMSpN_bn{ii});
      freqhiall_APT_MSPN_bn{ii}=ft_math(cfg,freqhiall_audPtac_bn{ii},freqhiall_audMSpN_bn{ii});
    end
    cfg=[];
    gravelo_TPA_MSPN=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN{:});
    gravelo_APT_MSPN=ft_freqgrandaverage(cfg,freqloall_APT_MSPN{:});
    gravehi_TPA_MSPN=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN{:});
    gravehi_APT_MSPN=ft_freqgrandaverage(cfg,freqhiall_APT_MSPN{:});
    gravelo_TPA_MSPN_bn=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_bn{:});
    gravelo_APT_MSPN_bn=ft_freqgrandaverage(cfg,freqloall_APT_MSPN_bn{:});
    gravehi_TPA_MSPN_bn=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_bn{:});
    gravehi_APT_MSPN_bn=ft_freqgrandaverage(cfg,freqhiall_APT_MSPN_bn{:});
    
    
    nsub=length(freqloall_audMSpN);
    
    cfg=[];
    cfg.latency=[-.8 0.8];
    cfg.neighbours=neighbours;
    % cfg.parameter='avg';
    cfg.parameter='powspctrm';
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
    stattl_mc=ft_freqstatistics(cfg, grindlo_tacPaud, grindlo_tacMSpN);
    statal_mc=ft_freqstatistics(cfg, grindlo_audPtac, grindlo_audMSpN);
    statth_mc=ft_freqstatistics(cfg, grindhi_tacPaud, grindhi_tacMSpN);
    statah_mc=ft_freqstatistics(cfg, grindhi_audPtac, grindhi_audMSpN);
    stattl_mc_bn=ft_freqstatistics(cfg, grindlo_tacPaud_bn, grindlo_tacMSpN_bn);
    statal_mc_bn=ft_freqstatistics(cfg, grindlo_audPtac_bn, grindlo_audMSpN_bn);
    statth_mc_bn=ft_freqstatistics(cfg, grindhi_tacPaud_bn, grindhi_tacMSpN_bn);
    statah_mc_bn=ft_freqstatistics(cfg, grindhi_audPtac_bn, grindhi_audMSpN_bn);
    
    
    ylimlo=[4 7; 8 12; 14 30];
    for ss=1:size(ylimlo,1)
      ylim=ylimlo(ss,:);
      
      try
        fig=240+ss;
        topoplotTFR_highlight(fig,stattl_mc,gravelo_TPA_MSPN,ylim);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a.png'],'-dpng')
        end
      catch
      end
      try
        fig=340+ss;
        topoplotTFR_highlight(fig,stattl_mc_bn,gravelo_TPA_MSPN_bn,ylim);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25_bn.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a_bn.png'],'-dpng')
        end
      catch
      end
      try
        fig=260+ss;
        topoplotTFR_highlight(fig,statal_mc,gravelo_APT_MSPN,ylim);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_at_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_at_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a.png'],'-dpng')
        end
      catch
      end
      try
        fig=360+ss;
        topoplotTFR_highlight(fig,statal_mc_bn,gravelo_APT_MSPN_bn,ylim);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_at_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25_bn.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_at_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a_bn.png'],'-dpng')
        end
      catch
      end
    end
    
    ylimhi=[30 45; 55 80];
    for ss=1:size(ylimhi,1)
      ylim=ylimhi(ss,:);
      try
        fig=250+ss;
        topoplotTFR_highlight(fig,statth_mc,gravehi_TPA_MSPN,ylim);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a.png'],'-dpng')
        end
      catch
      end
      
      try
        fig=350+ss;
        topoplotTFR_highlight(fig,statth_mc_bn,gravehi_TPA_MSPN_bn,ylim);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25_bn.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a_bn.png'],'-dpng')
        end
      catch
      end
      try
        fig=270+ss;
        topoplotTFR_highlight(fig,statah_mc,gravehi_APT_MSPN,ylim);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_at_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_at_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a.png'],'-dpng')
        end
      catch
      end
      try
        fig=370+ss;
        topoplotTFR_highlight(fig,statah_mc_bn,gravehi_APT_MSPN_bn,ylim);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_at_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25_bn.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_at_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a_bn.png'],'-dpng')
        end
      catch
      end
      
    end
    
    
    
    
    
    
    
    
  end
end






