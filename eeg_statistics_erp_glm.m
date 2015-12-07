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

sub{100}='p01'; % ma.a. 03/04/14
sub{1}='e01'; % ab.m. 21/05/14
sub{2}='e02'; % ma.a. 04/06/14
sub{3}='e03'; % ag.m. 10/06/14
sub{4}='e04'; % re.g. 17/06/14
%%%  above is pilot, below is real
sub{5}='e05'; % ma.a. 25/06/14
sub{6}='e06'; % e.u.  01/07/14
sub{7}='e07'; % a.s.  07/07/14
sub{8}='e08'; % k.t.  09/07/14
sub{9}='e09';% d.a.  14/07/14
sub{10}='e10';% k.l.  15/07/14
sub{11}='e11';% ab.m.  16/07/14  % from here on, had EOG/EMG
sub{12}='e12';% b.s.  17/07/14
sub{13}='e13';% d.t.  21/07/14
sub{14}='e14';% f.g.  22/07/14
sub{15}='e15';% r.m.  23/07/14
sub{16}='e16';% t.p.  24/07/14 % from here on, had attempt Polhemus
sub{17}='e17';% t.t.  28/07/14
sub{18}='e18';% k.n.v.  29/07/14
sub{19}='e19';% j.b.  30/07/14
sub{20}='e20';% n.m.  31/07/14
sub{21}='e21';% l.c.  04/08/14
sub{22}='e22';% a.b.  05/08/14
sub{23}='e23';% r.c.
sub{24}='e24';% a.d.
sub{25}='e25';% j.c.
sub{26}='e26';% r.s.
sub{27}='e27';% a.p.
sub{28}='e28';% w.p.
sub{29}='e29';% i.r.
sub{30}='e30';% o.y.l.
sub{31}='e31';% r.b.
sub{32}='e32';% i.f.

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


%% contrast of conditions in GLM: Aalone + Talone_shifted - ATmult - N


time=-0.1:0.001:0.4;

contac=nan(length(time),63,9,4,28);
conaud=nan(length(time),63,9,4,28);

for ii=5:18
  clearvars -except ii sub edir ddir time conaud contac labels
  cd([edir sub{ii} ])
  load(['tlock_diffs_' sub{ii} '.mat']);
  
  
  
  load(['raw_each_rej_' sub{ii}],'raw_*');
  
  cfg=[];
  cfg.bpfilter='yes';
  cfg.bpfreq=[1 40];
  cfg.demean='yes';
  cfg.baselinewindow=[-.6 -.1];
  cfg.channel= {'all', '-ECG'};
  data_tac_filt=ft_preprocessing(cfg,raw_tac_rej);
  data_aud_filt=ft_preprocessing(cfg,raw_aud_rej);
  data_nul_filt=ft_preprocessing(cfg,raw_nul_rej);
  
  % % Solve this: what should be reference channel for scalp ERP? fixme
  cfg=[];
  cfg.reref='yes';
  cfg.refchannel='all'; % necessary if doing source localisation
  cfg.refchannel={'FT9', 'FT10'}; % sort of like linked mastoids?
  data_tac_filt_ref=ft_preprocessing(cfg,data_tac_filt);
  data_aud_filt_ref=ft_preprocessing(cfg,data_aud_filt);
  data_nul_filt_ref=ft_preprocessing(cfg,data_nul_filt);
  
  clear raw_*
  
  if ii>7
    soalist=[1 3 4 5 6 7 9];
  else
    soalist=[3 4 5 6 7];
  end
  % include only trials during awake segments (in case participant went into N1 during sitting up portion)
  for ll=1:length(soalist)
    for tt=3:6 % different thresholds of inclusion of SOA variance around desired SOA (see concat_stim_files.m)
      cfg=[];
      cfg.vartrllength=2;
      cfg.keeptrials='yes';
      cfg.trials=[data_tac_filt_ref.trialinfo(:,tt)==soalist(ll) & data_tac_filt_ref.trialinfo(:,2)==0 & ~isnan(data_tac_filt_ref.trialinfo(:,10))];
      tlock_tac_s0{soalist(ll),tt-2}=ft_timelockanalysis(cfg,data_tac_filt_ref); % only 'awake' state
      %   cfg.trials=data_tac_filt_ref.trialinfo(:,3)==soalist(ll);
      %   tlock_tac_sall{soalist(ll)}=ft_timelockanalysis(cfg,data_tac_filt_ref); % all trials, irrepective if fell asleep
    end
  end
  for ll=1:length(soalist)
    for tt=3:6 % different thresholds of inclusion of SOA variance around desired SOA (see concat_stim_files.m)
      cfg=[];
      cfg.vartrllength=2;
      cfg.keeptrials='yes';
      cfg.trials=[data_aud_filt_ref.trialinfo(:,tt)==soalist(ll) & data_aud_filt.trialinfo(:,2)==0 & ~isnan(data_aud_filt_ref.trialinfo(:,10))];
      tlock_aud_s0{soalist(ll),tt-2}=ft_timelockanalysis(cfg,data_aud_filt_ref);
      %   cfg.trials=data_aud_filt_ref.trialinfo(:,3)==soalist(ll);
      %   tlock_aud_sall{soalist(ll)}=ft_timelockanalysis(cfg,data_aud_filt_ref);
    end
  end
  cfg=[]; % tac alone 38 (-2)
  cfg.vartrllength=2;
  cfg.keeptrials='yes';
  cfg.trials=[data_tac_filt_ref.trialinfo(:,3)==-2 & data_tac_filt_ref.trialinfo(:,2)==0];
  tlock_tac_s0{10,1}=ft_timelockanalysis(cfg,data_tac_filt_ref);
  % cfg.trials=data_tac_filt_ref.trialinfo(:,3)==-2;
  % tlock_tac_sall{10}=ft_timelockanalysis(cfg,data_tac_filt_ref);
  cfg=[]; % aud alone 39 (-1)
  cfg.vartrllength=2;
  cfg.keeptrials='yes';
  cfg.trials=[data_aud_filt_ref.trialinfo(:,3)==-1 & data_aud_filt_ref.trialinfo(:,2)==0];
  tlock_aud_s0{10,1}=ft_timelockanalysis(cfg,data_aud_filt_ref);
  % cfg.trials=data_aud_filt_ref.trialinfo(:,3)==-1;
  % tlock_aud_sall{10}=ft_timelockanalysis(cfg,data_aud_filt_ref);
  cfg=[];
  cfg.vartrllength=2;
  cfg.keeptrials='yes';
  cfg.trials=data_nul_filt_ref.trialinfo(:,2)==0;
  tlock_nul_s0{10,1}=ft_timelockanalysis(cfg,data_nul_filt_ref);
  % cfg.trials='all';
  % tlock_nul_sall{10}=ft_timelockanalysis(cfg,data_nul_filt_ref);
  
  clear data*
  
  if ii==5
    labels=tlock_tac_s0{5,1}.label;
  end
  chanuse=match_str(labels,tlock_tac_s0{5,1}.label);

  for ll=[soalist 10] % this loop is to remove trials with NaN in (e.g. if artifact rejection cut through middle of a trial)
    for tt=1:4
      if ll==10 && tt>1
        continue
      else
        cfg=[];
        cfg.keeptrials='yes';
        cfg.trials=find(sum(squeeze(isnan(tlock_tac_s0{ll,tt}.trial(:,10,:)))')<.05*size(tlock_tac_s0{ll,tt}.trial,3));
        tlock_tac_s0{ll,tt}=ft_timelockanalysis(cfg,tlock_tac_s0{ll,tt});
        %   cfg.trials=find(sum(squeeze(isnan(tlock_tac_sall{ll}.trial(:,10,:)))')<.05*size(tlock_tac_sall{ll}.trial,3));
        %   tlock_tac_sall{ll}=ft_timelockanalysis(cfg,tlock_tac_sall{ll});
        cfg=[];
        cfg.keeptrials='yes';
        cfg.trials=find(sum(squeeze(isnan(tlock_aud_s0{ll,tt}.trial(:,10,:)))')<.05*size(tlock_aud_s0{ll,tt}.trial,3));
        tlock_aud_s0{ll,tt}=ft_timelockanalysis(cfg,tlock_aud_s0{ll,tt});
        %   cfg.trials=find(sum(squeeze(isnan(tlock_aud_sall{ll}.trial(:,10,:)))')<.05*size(tlock_aud_sall{ll}.trial,3));
        %   tlock_aud_sall{ll}=ft_timelockanalysis(cfg,tlock_aud_sall{ll});
      end
      if ll==10 && tt==1
        cfg=[];
        cfg.keeptrials='yes';
        cfg.trials=find(sum(squeeze(isnan(tlock_nul_s0{ll}.trial(:,10,:)))')<.05*size(tlock_nul_s0{ll}.trial,3));
        tlock_nul_s0{ll}=ft_timelockanalysis(cfg,tlock_nul_s0{ll});
        %     cfg.trials=find(sum(squeeze(isnan(tlock_nul_sall{ll}.trial(:,10,:)))')<.05*size(tlock_nul_sall{ll}.trial,3));
        %     tlock_nul_sall{ll}=ft_timelockanalysis(cfg,tlock_nul_sall{ll});
      end
    end
  end
  
  %  contrasting conditions
  
%   plotflag=1;
  load([ddir sub{ii} '_audtac.mat']);
  soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];
%   chanuse=match_str(tlock_tac_s0{3}.label,'Fz');
  
  % determine trial numbers per condition
  numuni=min([size(tlock_tac_s0{10,1}.trialinfo,1),size(tlock_aud_s0{10,1}.trialinfo,1),size(tlock_nul_s0{10,1}.trialinfo,1)]);
  for ll=soalist
    for tt=1:4 % limits of plus/minus:  [no-limit .01 .005 .002]
      trcom=intersect(tlock_tac_s0{ll,tt}.trialinfo(:,11),tlock_aud_s0{ll,tt}.trialinfo(:,11));
      nummul(ll,tt)=length(trcom);
      numcond(ll,tt)=min(nummul(ll,tt),numuni);
    end
  end
  numcondall=min(numcond(soalist,:)); % min over all conditions for a given SOA variance threshold
  
  fsample=1/diff(tlock_tac_s0{ll,tt}.time(1:2)); % sampling rate
  
  try close(10);end
  try close(11);end
  try close(19);end
  try close(20);end
  % create sum of unisensory, specific to the exact SOAs and actual jitters, for each multisensory condition
  for ll=soalist
    for tt=1:4
      
      % make sure the multisensory auditory-locked and ms tactile-locked have same trials; sometimes not exactly!
      trcom=intersect(tlock_tac_s0{ll,tt}.trialinfo(:,11),tlock_aud_s0{ll,tt}.trialinfo(:,11));
      [~,tbb,~]=intersect(tlock_tac_s0{ll,tt}.trialinfo(:,11),trcom);
      [~,abb,~]=intersect(tlock_aud_s0{ll,tt}.trialinfo(:,11),trcom);
      
      cfg=[];
      % I have thought it through, and no minus sign needed here.
      %     cfg.offset=round(fsample*(tlock_tac_s0{ll,tt}.trialinfo(trialkept,10)-soades(ll))); % offset is in samples not seconds
      cfg.offset=round(fsample*(tlock_tac_s0{ll,tt}.trialinfo(tbb(mtrialkept{ll,tt}),10)-soades(ll) -soades(ll))); % offset is in samples not seconds; -2*soades required for tactile not auditory
      cfg.trials=ttrialkept{ll,tt};
      tlock_tac_s0{ll+40,tt}=ft_redefinetrial(cfg,tlock_tac_s0{10}); % Talone with shift in time matching SOA jitter
      cfg=[];
      cfg.vartrllength=2;
      tlock_tac_s0tlock{ll+40,tt}=ft_timelockanalysis(cfg,tlock_tac_s0{ll+40,tt}); % Talone with shift in time matching SOA jitter
      
      cfg=[];
      % I have thought it through, and no minus sign needed here.
      %     cfg.offset=round(fsample*(tlock_aud_s0{ll,tt}.trialinfo(trialkept,10)-soades(ll))); % tlock_aud_s0{ll,tt}.trialinfo(trialkept,10) is identical to tlock_tac_s0{ll,tt}.trialinfo(trialkept,10)
      cfg.offset=round(fsample*(tlock_aud_s0{ll,tt}.trialinfo(abb(mtrialkept{ll,tt}),10)-soades(ll) +soades(ll))); % tlock_aud_s0{ll,tt}.trialinfo(trialkept,10) is identical to tlock_tac_s0{ll,tt}.trialinfo(trialkept,10)
      cfg.trials=atrialkept{ll,tt};
      tlock_aud_s0{ll+40,tt}=ft_redefinetrial(cfg,tlock_aud_s0{10}); % Aalone with shift in time matching SOA jitter
      cfg=[];
      cfg.vartrllength=2;
      tlock_aud_s0tlock{ll+40,tt}=ft_timelockanalysis(cfg,tlock_aud_s0{ll+40,tt}); % Aalone with shift in time matching SOA jitter
      
      cfg=[];
      cfg.trials=ntrialkept{ll,tt};
      tlock_nul_s0{ll+40,tt}=ft_redefinetrial(cfg,tlock_nul_s0{10});
      
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
      tlock_tac_s0{ll+20,tt}=ft_math(cfg,tlock_tac_s0{ll,tt},ft_timelockanalysis([],tlock_nul_s0{ll+40,tt})); % TA+N
      %   tlock_tac_sall{ll+20,tt}=ft_math(cfg,tlock_tac_sall{ll,tt},tlock_nul_sall{10}); % TA+N
      tlock_aud_s0{ll+20,tt}=ft_math(cfg,tlock_aud_s0{ll,tt},ft_timelockanalysis([],tlock_nul_s0{ll+40,tt})); % AT+N
      %   tlock_aud_sall{ll+20,tt}=ft_math(cfg,tlock_aud_sall{ll,tt},tlock_nul_sall{10}); % AT+N
      
      
      
      cfg=[];
      cfg.operation='subtract';
      cfg.parameter='avg';
      %   tlock_tpa_mtamn_s0{ll}=ft_math(cfg,ft_timelockanalysis([],tlock_tacPaud_s0{ll}),ft_timelockanalysis([],tlock_tac_s0{ll+20})); % (T + As) - (TA + N)
      %   tlock_tpa_mtamn_sall{ll}=ft_math(cfg,ft_timelockanalysis([],tlock_tacPaud_sall{ll}),ft_timelockanalysis([],tlock_tac_sall{ll+20})); % (T + As) - (TA + N)
      %   tlock_apt_matmn_s0{ll}=ft_math(cfg,ft_timelockanalysis([],tlock_audPtac_s0{ll}),ft_timelockanalysis([],tlock_aud_s0{ll+20})); % (A + Ts) - (AT + N)
      %   tlock_apt_matmn_sall{ll}=ft_math(cfg,ft_timelockanalysis([],tlock_audPtac_sall{ll}),ft_timelockanalysis([],tlock_aud_sall{ll+20})); % (A + Ts) - (AT + N)
      tlock_tpa_mtamn_s0{ll,tt}=ft_math(cfg,tlock_tacPaud_s0{ll,tt},tlock_tac_s0{ll+20,tt}); % (T + As) - (TA + N)
      %     tlock_tpa_mtamn_sall{ll}=ft_math(cfg,tlock_tacPaud_sall{ll},tlock_tac_sall{ll+20}); % (T + As) - (TA + N)
      tlock_apt_matmn_s0{ll,tt}=ft_math(cfg,tlock_audPtac_s0{ll,tt},tlock_aud_s0{ll+20,tt}); % (A + Ts) - (AT + N)
      %     tlock_apt_matmn_sall{ll}=ft_math(cfg,tlock_audPtac_sall{ll},tlock_aud_sall{ll+20}); % (A + Ts) - (AT + N)
      
    end
  end
  
  
  
  
  
  for ll=soalist
    for tt=1:4
      for teatime=1:length(time)
        audtac=tlock_tac_s0{ll,tt}.trial(mtrialkept{ll,tt},:,dsearchn(tlock_tac_s0{ll,tt}.time',time(teatime)));
        nul=tlock_nul_s0{ll+40,tt}.trial(:,:,dsearchn(tlock_nul_s0{ll+40,tt}.time',time(teatime)));
        tac=tlock_tac_s0{10,1}.trial(ttrialkept{ll,tt},:,dsearchn(tlock_tac_s0{10,1}.time',time(teatime)));
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
        cfg.glm.standardise=0; % this should be 1 (need function)
        stat=ft_statfun_glm(cfg,data',design');
        beta=reshape(stat.stat,[size(design,2) size(data,2)]);
        contac(teatime,chanuse,ll,tt,ii-4)=[1 1 -1 -1 0]*beta;
        
        
        audtac=tlock_aud_s0{ll,tt}.trial(mtrialkept{ll,tt},:,dsearchn(tlock_aud_s0{ll,tt}.time',time(teatime)));
        nul=tlock_nul_s0{ll+40,tt}.trial(:,:,dsearchn(tlock_nul_s0{ll+40,tt}.time',time(teatime)));
        aud=tlock_aud_s0{10,1}.trial(atrialkept{ll,tt},:,dsearchn(tlock_aud_s0{10,1}.time',time(teatime)));
        tac=tlock_tac_s0{ll+40,tt}.trial(:,:,dsearchn(tlock_tac_s0{ll+40,tt}.time',time(teatime)));
                
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
        cfg.glm.standardise=0; % this should be 1 (need function)
        stat=ft_statfun_glm(cfg,data',design');
        beta=reshape(stat.stat,[size(design,2) size(data,2)]);
        conaud(teatime,chanuse,ll,tt,ii-4)=[1 1 -1 -1 0]*beta;
      end
    end
  end
  save([edir 'conaudtac'],'con*','labels','time');
end
  
%%

soalist=[1 3 4 5 6 7 9];
for ll=soalist
  for tt=1:4
    for cc=1:63
      [hh,ppa(:,cc,ll,tt)]=ttest(squeeze(conaud(:,cc,ll,tt,:))');
      [hh,ppt(:,cc,ll,tt)]=ttest(squeeze(contac(:,cc,ll,tt,:))');
    end
  end
end

figure;
for ll=soalist
  for tt=1:4
    subplot(9,4,(ll-1)*4 +tt);
    imagesc(ppa(:,:,ll,tt)'<.001);
    imagesc(ppt(:,:,ll,tt)'<.001);
  end
end

  