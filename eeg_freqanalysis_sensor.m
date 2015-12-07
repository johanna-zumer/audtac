% preprocessing of EEG data from Hills, 64ch MR-cap
% clear all
close all
if ispc
  edir='D:\audtac\eeg_data\';
  ddir='D:\audtac\legomagic\diaries\';
  bdir='D:\audtac\behav_data\';
  sdir='D:\audtac\spss_stuff\';
else
  edir='/mnt/hgfs/D/audtac/eeg_data/';
  ddir='/mnt/hgfs/D/audtac/legomagic/diaries/';
  bdir='/mnt/hgfs/D/audtac/behav_data/';
  sdir='/mnt/hgfs/D/audtac/spss_stuff/';
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

%%
for ii=5:32
  
  cd([edir sub{ii} ])
  clearvars -except ii sub edir ddir
  pack
  
  
%   while ~exist(['tlock_trialSel_' sub{ii} '.mat'],'file')
%     pause(60);
%   end
  
%   try
%     load(['tlock_trialSel_' sub{ii} '.mat']);
%   catch % catch is currently redundant
%     [tlock_tac_s0,tlock_aud_s0,tlock_nul_s0,tr]=eeg_legomagic_trialSelection1(ii);
    [tlock_tac_s0,tlock_aud_s0,tlock_nul_s0,tr]=eeg_legomagic_trialSelection_freqwide(ii);
    save(['tlock_trialSel_' sub{ii} '.mat'],'tlock*s0','tr','-v7.3');
%   end
  
  %  Frequency analysis
  if ii<8
    soalist=[3 4 5 6 7];
  else
    soalist=[1 3 4 5 6 7 9];
  end
  
  
  for tt=1:4
    % Lower frequencies with Hanning taper
    for ll=[soalist soalist+20 soalist+40]
      
      if ll<10
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
      end
      if ll<10
        numt_trials(ll,tt)=size(tlock_tac_s0{ll+20,tt}.trial,1); % this should be the same for tacll20, audll40, tacll, nulll50
        numa_trials(ll,tt)=size(tlock_aud_s0{ll+20,tt}.trial,1); % this should be the same for audll20, tacll40, audll, nulll60
      end
      
      if ll>40
        tlock_tac_s0{ll,tt}=trim_nans(tlock_tac_s0{ll,tt});
        tlock_aud_s0{ll,tt}=trim_nans(tlock_aud_s0{ll,tt});
      end
      
      % Do it two ways:  FFT of each condition then add, or add then FFT
      
      % %%%%      % FFT then add    %%%%%%%%%%%%%
      if size(tlock_tac_s0{ll,tt}.trial,1)
        %         cfg=[];
        %         cfg.lpfilter='yes';
        %         cfg.lpfreq=40;
        %         cfg.demean='yes';
        %         cfg.baselinewindow=[-1.7 -0.6];
        cfg=[];
        cfg.method='mtmconvol';
        cfg.pad=4;
        cfg.foi=4:2:30;
        cfg.taper='hanning';
        cfg.toi=-1.2:0.1:1.3;
        cfg.t_ftimwin=4./cfg.foi;
        %         cfg.keeptrials='yes';  % not necessary at this stage for keeping trials (yes later for source analysis and phase resetting)
        freqlo_tac{ll,tt}=ft_freqanalysis(cfg,tlock_tac_s0{ll,tt});
        
        % Higher frequencies with dpss multitaper
        cfg=[];
        cfg.method='mtmconvol';
        cfg.pad=4;
        cfg.foi=[30:5:45 55:5:80];
        cfg.taper='dpss';
        cfg.tapsmofrq=7*ones(1,length(cfg.foi));
        cfg.toi=-1.2:0.1:1.3;
        cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
        freqhi_tac{ll,tt}=ft_freqanalysis(cfg,tlock_tac_s0{ll,tt});
      else
        freqlo_tac{ll,tt}=[];
        freqhi_tac{ll,tt}=[];
      end
      %       tlock_tac_s0{ll,tt}=[]; % clearing to help save memory as running
      
      if size(tlock_aud_s0{ll,tt}.trial,1)
        cfg=[];
        cfg.method='mtmconvol';
        cfg.pad=4;
        cfg.foi=4:2:30;
        cfg.taper='hanning';
        cfg.toi=-1.2:0.1:1.3;
        cfg.t_ftimwin=4./cfg.foi;
        freqlo_aud{ll,tt}=ft_freqanalysis(cfg,tlock_aud_s0{ll,tt});
        cfg=[];
        cfg.method='mtmconvol';
        cfg.pad=4;
        cfg.foi=[30:5:45 55:5:80];
        cfg.taper='dpss';
        cfg.tapsmofrq=7*ones(1,length(cfg.foi));
        cfg.toi=-1.2:0.1:1.3;
        cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
        freqhi_aud{ll,tt}=ft_freqanalysis(cfg,tlock_aud_s0{ll,tt});
      else
        freqlo_aud{ll,tt}=[];
        freqhi_aud{ll,tt}=[];
      end
      %       tlock_aud_s0{ll,tt}=[]; % clearing to help save memory as running
    end
    for ll=[soalist+50 soalist+60]
      if size(tlock_nul_s0{ll,tt}.trial,1)
        cfg=[];
        cfg.method='mtmconvol';
        cfg.pad=4;
        cfg.foi=4:2:30;
        cfg.taper='hanning';
        cfg.toi=-1.2:0.1:1.3;
        cfg.t_ftimwin=4./cfg.foi;
        freqlo_nul{ll,tt}=ft_freqanalysis(cfg,tlock_nul_s0{ll,tt});
        cfg=[];
        cfg.method='mtmconvol';
        cfg.pad=4;
        cfg.foi=[30:5:45 55:5:80];
        cfg.taper='dpss';
        cfg.tapsmofrq=7*ones(1,length(cfg.foi));
        cfg.toi=-1.2:0.1:1.3;
        cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
        freqhi_nul{ll,tt}=ft_freqanalysis(cfg,tlock_nul_s0{ll,tt});
      else
        freqlo_nul{ll,tt}=[];
        freqhi_nul{ll,tt}=[];
      end
      %       tlock_nul_s0{ll,tt}=[]; % clearing to help save memory as running
    end
    
    for ll=soalist
      
      % create sum of unisensory conditions
      cfg=[];
      cfg.operation='add';
      cfg.parameter='powspctrm';
      %       freqlo_tacPaud{ll,tt}=ft_math(cfg,ft_selectdata(cfgavg,freqlo_tac{10}),ft_selectdata(cfgavg,freqlo_aud{ll+40,tt})); % keeping tac centred but jitter aud (to match tlock_tac_s0{X,tt})
      %       freqlo_audPtac{ll,tt}=ft_math(cfg,ft_selectdata(cfgavg,freqlo_aud{10}),ft_selectdata(cfgavg,freqlo_tac{ll+40,tt})); % keeping aud centred but jitter tac (to match tlock_aud_s0{X,tt})
      %       freqhi_tacPaud{ll,tt}=ft_math(cfg,ft_selectdata(cfgavg,freqhi_tac{10}),ft_selectdata(cfgavg,freqhi_aud{ll+40,tt})); % keeping tac centred but jitter aud (to match tlock_tac_s0{X,tt})
      %       freqhi_audPtac{ll,tt}=ft_math(cfg,ft_selectdata(cfgavg,freqhi_aud{10}),ft_selectdata(cfgavg,freqhi_tac{ll+40,tt})); % keeping aud centred but jitter tac (to match tlock_aud_s0{X,tt})
      if numt_trials(ll,tt)
        freqlo_tacPaud_fftadd{ll,tt}=ft_math(cfg,freqlo_tac{ll+20,tt},freqlo_aud{ll+40,tt}); % keeping tac centred but jitter aud (to match tlock_tac_s0{X,tt})
        freqhi_tacPaud_fftadd{ll,tt}=ft_math(cfg,freqhi_tac{ll+20,tt},freqhi_aud{ll+40,tt}); % keeping tac centred but jitter aud (to match tlock_tac_s0{X,tt})
      else
        freqlo_tacPaud_fftadd{ll,tt}=[];
        freqhi_tacPaud_fftadd{ll,tt}=[];
      end
      if numa_trials(ll,tt)
        freqlo_audPtac_fftadd{ll,tt}=ft_math(cfg,freqlo_aud{ll+20,tt},freqlo_tac{ll+40,tt}); % keeping aud centred but jitter tac (to match tlock_aud_s0{X,tt})
        freqhi_audPtac_fftadd{ll,tt}=ft_math(cfg,freqhi_aud{ll+20,tt},freqhi_tac{ll+40,tt}); % keeping aud centred but jitter tac (to match tlock_aud_s0{X,tt})
      else
        freqlo_audPtac_fftadd{ll,tt}=[];
        freqhi_audPtac_fftadd{ll,tt}=[];
      end
      
      % create sum of mulitsensory and null conditions
      %       freqlo_tacMSpN{ll,tt}=ft_math(cfg,ft_selectdata(cfgavg,freqlo_tac{ll,tt}),ft_selectdata(cfgavg,freqlo_nul{ll+40,tt})); % keeping tac centred but jitter aud (to match tlock_tac_s0{X,tt})
      %       freqlo_audMSpN{ll,tt}=ft_math(cfg,ft_selectdata(cfgavg,freqlo_aud{ll,tt}),ft_selectdata(cfgavg,freqlo_nul{ll+40,tt})); % keeping aud centred but jitter tac (to match tlock_aud_s0{X,tt})
      %       freqhi_tacMSpN{ll,tt}=ft_math(cfg,ft_selectdata(cfgavg,freqhi_tac{ll,tt}),ft_selectdata(cfgavg,freqhi_nul{ll+40,tt})); % keeping tac centred but jitter aud (to match tlock_tac_s0{X,tt})
      %       freqhi_audMSpN{ll,tt}=ft_math(cfg,ft_selectdata(cfgavg,freqhi_aud{ll,tt}),ft_selectdata(cfgavg,freqhi_nul{ll+40,tt})); % keeping aud centred but jitter tac (to match tlock_aud_s0{X,tt})
      if numt_trials(ll,tt)
        freqlo_tacMSpN_fftadd{ll,tt}=ft_math(cfg,freqlo_tac{ll,tt},freqlo_nul{ll+50,tt}); % keeping tac centred but jitter aud (to match tlock_tac_s0{X,tt})
        freqhi_tacMSpN_fftadd{ll,tt}=ft_math(cfg,freqhi_tac{ll,tt},freqhi_nul{ll+50,tt}); % keeping tac centred but jitter aud (to match tlock_tac_s0{X,tt})
      else
        freqlo_tacMSpN_fftadd{ll,tt}=[];
        freqhi_tacMSpN_fftadd{ll,tt}=[];
      end
      if numa_trials(ll,tt)
        freqlo_audMSpN_fftadd{ll,tt}=ft_math(cfg,freqlo_aud{ll,tt},freqlo_nul{ll+60,tt}); % keeping aud centred but jitter tac (to match tlock_aud_s0{X,tt})
        freqhi_audMSpN_fftadd{ll,tt}=ft_math(cfg,freqhi_aud{ll,tt},freqhi_nul{ll+60,tt}); % keeping aud centred but jitter tac (to match tlock_aud_s0{X,tt})
      else
        freqlo_audMSpN_fftadd{ll,tt}=[];
        freqhi_audMSpN_fftadd{ll,tt}=[];
      end
    end
    
    
    %%%                     %%%%%
    %%% add first then FFT  %%%%%
    
    minnumcomb=1000; % don't let numcomb go higher than 1000.
    minnumcomb=1; % in effect just the normal combination of trials
    
    
    for ll=soalist
      % tacPaud
      numcomb=numt_trials(ll,tt)^2;
      numtests=min(numcomb,minnumcomb);
      combindex=reshape(1:numcomb,numt_trials(ll,tt),numt_trials(ll,tt));
      for cc=1:numtests
        tmp=Shuffle(combindex(:));
        combuse=tmp(1:numt_trials(ll,tt));
        for at=1:numt_trials(ll,tt)
          if cc==1  % do it as 'normal'
            tind=at;aind=at;
          else
            [tind,aind]=find(combindex==combuse(at));
          end
          cfg=[];
          cfg.trials=tind;
          tmpt=ft_selectdata(cfg,tlock_tac_s0{ll+20,tt});
          cfg.trials=aind;
          tmpa=ft_selectdata(cfg,tlock_aud_s0{ll+40,tt});
          cfg=[];
          cfg.operation='add';
          cfg.parameter='trial';
          tmpsum=ft_math(cfg,tmpt,tmpa);
          if at==1
            tlock_fake=tmpsum;
          end
          tlock_fake.trial(at,:,:)=tmpsum.trial(1,:,:);
        end
        
        cfg=[];
        cfg.method='mtmconvol';
        cfg.pad=4;
        cfg.foi=4:2:30;
        cfg.taper='hanning';
        cfg.toi=-1.2:0.1:1.3;
        cfg.t_ftimwin=4./cfg.foi;
        %         cfg.keeptrials='yes';  % not necessary at this stage for keeping trials (yes later for source analysis and phase resetting)
        freqlo_tacPaud_tmp{ll,tt,cc}=ft_freqanalysis(cfg,tlock_fake);
        cfg=[];
        cfg.method='mtmconvol';
        cfg.pad=4;
        cfg.foi=[30:5:45 55:5:80];
        cfg.taper='dpss';
        cfg.tapsmofrq=7*ones(1,length(cfg.foi));
        cfg.toi=-1.2:0.1:1.3;
        cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
        freqhi_tacPaud_tmp{ll,tt,cc}=ft_freqanalysis(cfg,tlock_fake);
      end
      clear tmp
      for cc=1:numtests,tmp(:,:,:,cc)=freqlo_tacPaud_tmp{ll,tt,cc}.powspctrm;end
      freqlo_tacPaud_comb{ll,tt,1}=freqlo_tacPaud_tmp{ll,tt,1};
      freqlo_tacPaud_comb{ll,tt,2}=freqlo_tacPaud_tmp{ll,tt,1};
      freqlo_tacPaud_comb{ll,tt,2}.powspctrm=mean(tmp,4);
      freqlo_tacPaud_comb{ll,tt,2}.stdspctrm=std(tmp,[],4);
      clear tmp
      for cc=1:numtests,tmp(:,:,:,cc)=freqhi_tacPaud_tmp{ll,tt,cc}.powspctrm;end
      freqhi_tacPaud_comb{ll,tt,1}=freqhi_tacPaud_tmp{ll,tt,1};
      freqhi_tacPaud_comb{ll,tt,2}=freqhi_tacPaud_tmp{ll,tt,1};
      freqhi_tacPaud_comb{ll,tt,2}.powspctrm=mean(tmp,4);
      freqhi_tacPaud_comb{ll,tt,2}.stdspctrm=std(tmp,[],4);
      clear freqlo_tacPaud_tmp freqhi_tacPaud_tmp
      
      % audPtac
      numcomb=numa_trials(ll,tt)^2;
      numtests=min(numcomb,minnumcomb);
      combindex=reshape(1:numcomb,numa_trials(ll,tt),numa_trials(ll,tt));
      for cc=1:numtests
        tmp=Shuffle(combindex(:));
        combuse=tmp(1:numa_trials(ll,tt));
        
        tlock_fake=tlock_aud_s0{ll+20,tt};
        for at=1:numa_trials(ll,tt)
          if cc==1  % do it as 'normal'
            tind=at;aind=at;
          else
            [aind,tind]=find(combindex==combuse(at));
          end
          cfg=[];
          cfg.trials=tind;
          tmpt=ft_selectdata(cfg,tlock_tac_s0{ll+40,tt});
          cfg.trials=aind;
          tmpa=ft_selectdata(cfg,tlock_aud_s0{ll+20,tt});
          cfg=[];
          cfg.operation='add';
          cfg.parameter='trial';
          tmpsum=ft_math(cfg,tmpt,tmpa);
          if at==1
            tlock_fake=tmpsum;
          end
          tlock_fake.trial(at,:,:)=tmpsum.trial(1,:,:);
%           tlock_fake.trial(at,:,:)=tlock_aud_s0{ll+20,tt}.trial(aind,:,:)+tlock_tac_s0{ll+40,tt}.trial(tind,:,:);
        end
        
        cfg=[];
        cfg.method='mtmconvol';
        cfg.pad=4;
        cfg.foi=4:2:30;
        cfg.taper='hanning';
        cfg.toi=-1.2:0.1:1.3;
        cfg.t_ftimwin=4./cfg.foi;
        %         cfg.keeptrials='yes';  % not necessary at this stage for keeping trials (yes later for source analysis and phase resetting)
        freqlo_audPtac_tmp{ll,tt,cc}=ft_freqanalysis(cfg,tlock_fake);
        cfg=[];
        cfg.method='mtmconvol';
        cfg.pad=4;
        cfg.foi=[30:5:45 55:5:80];
        cfg.taper='dpss';
        cfg.tapsmofrq=7*ones(1,length(cfg.foi));
        cfg.toi=-1.2:0.1:1.3;
        cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
        freqhi_audPtac_tmp{ll,tt,cc}=ft_freqanalysis(cfg,tlock_fake);
      end
      clear tmp
      for cc=1:numtests,tmp(:,:,:,cc)=freqlo_audPtac_tmp{ll,tt,cc}.powspctrm;end
      freqlo_audPtac_comb{ll,tt,1}=freqlo_audPtac_tmp{ll,tt,1};
      freqlo_audPtac_comb{ll,tt,2}=freqlo_audPtac_tmp{ll,tt,1};
      freqlo_audPtac_comb{ll,tt,2}.powspctrm=mean(tmp,4);
      freqlo_audPtac_comb{ll,tt,2}.stdspctrm=std(tmp,[],4);
      clear tmp
      for cc=1:numtests,tmp(:,:,:,cc)=freqhi_audPtac_tmp{ll,tt,cc}.powspctrm;end
      freqhi_audPtac_comb{ll,tt,1}=freqhi_audPtac_tmp{ll,tt,1};
      freqhi_audPtac_comb{ll,tt,2}=freqhi_audPtac_tmp{ll,tt,1};
      freqhi_audPtac_comb{ll,tt,2}.powspctrm=mean(tmp,4);
      freqhi_audPtac_comb{ll,tt,2}.stdspctrm=std(tmp,[],4);
      clear freqlo_audPtac_tmp freqhi_audPtac_tmp
      
      % tacMSpN
      numcomb=numt_trials(ll,tt)^2;
      numtests=min(numcomb,minnumcomb);
      combindex=reshape(1:numcomb,numt_trials(ll,tt),numt_trials(ll,tt));
      for cc=1:numtests
        tmp=Shuffle(combindex(:));
        combuse=tmp(1:numt_trials(ll,tt));
        
        tlock_fake=tlock_tac_s0{ll,tt};
        for at=1:numt_trials(ll,tt)
          if cc==1  % do it as 'normal'
            tind=at;aind=at;
          else
            [tind,aind]=find(combindex==combuse(at));
          end
          cfg=[];
          cfg.trials=tind;
          tmpt=ft_selectdata(cfg,tlock_tac_s0{ll,tt});
          cfg.trials=aind;
          tmpa=ft_selectdata(cfg,tlock_nul_s0{ll+50,tt});
          cfg=[];
          cfg.operation='add';
          cfg.parameter='trial';
          tmpsum=ft_math(cfg,tmpt,tmpa);
          if at==1
            tlock_fake=tmpsum;
          end
          tlock_fake.trial(at,:,:)=tmpsum.trial(1,:,:);
%           tlock_fake.trial(at,:,:)=tlock_tac_s0{ll,tt}.trial(tind,:,:)+tlock_nul_s0{ll+50,tt}.trial(aind,:,:);
        end
        cfg=[];
        cfg.method='mtmconvol';
        cfg.pad=4;
        cfg.foi=4:2:30;
        cfg.taper='hanning';
        cfg.toi=-1.2:0.1:1.3;
        cfg.t_ftimwin=4./cfg.foi;
        %         cfg.keeptrials='yes';  % not necessary at this stage for keeping trials (yes later for source analysis and phase resetting)
        freqlo_tacMSpN_tmp{ll,tt,cc}=ft_freqanalysis(cfg,tlock_fake);
        cfg=[];
        cfg.method='mtmconvol';
        cfg.pad=4;
        cfg.foi=[30:5:45 55:5:80];
        cfg.taper='dpss';
        cfg.tapsmofrq=7*ones(1,length(cfg.foi));
        cfg.toi=-1.2:0.1:1.3;
        cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
        freqhi_tacMSpN_tmp{ll,tt,cc}=ft_freqanalysis(cfg,tlock_fake);
      end
      clear tmp
      for cc=1:numtests,tmp(:,:,:,cc)=freqlo_tacMSpN_tmp{ll,tt,cc}.powspctrm;end
      freqlo_tacMSpN_comb{ll,tt,1}=freqlo_tacMSpN_tmp{ll,tt,1};
      freqlo_tacMSpN_comb{ll,tt,2}=freqlo_tacMSpN_tmp{ll,tt,1};
      freqlo_tacMSpN_comb{ll,tt,2}.powspctrm=mean(tmp,4);
      freqlo_tacMSpN_comb{ll,tt,2}.stdspctrm=std(tmp,[],4);
      clear tmp
      for cc=1:numtests,tmp(:,:,:,cc)=freqhi_tacMSpN_tmp{ll,tt,cc}.powspctrm;end
      freqhi_tacMSpN_comb{ll,tt,1}=freqhi_tacMSpN_tmp{ll,tt,1};
      freqhi_tacMSpN_comb{ll,tt,2}=freqhi_tacMSpN_tmp{ll,tt,1};
      freqhi_tacMSpN_comb{ll,tt,2}.powspctrm=mean(tmp,4);
      freqhi_tacMSpN_comb{ll,tt,2}.stdspctrm=std(tmp,[],4);
      clear freqlo_tacMSpN_tmp freqhi_tacMSpN_tmp
      
      % audMSpN
      numcomb=numa_trials(ll,tt)^2;
      numtests=min(numcomb,minnumcomb);
      combindex=reshape(1:numcomb,numa_trials(ll,tt),numa_trials(ll,tt));
      for cc=1:numtests
        tmp=Shuffle(combindex(:));
        combuse=tmp(1:numa_trials(ll,tt));
        
        tlock_fake=tlock_aud_s0{ll,tt};
        for at=1:numa_trials(ll,tt)
          if cc==1  % do it as 'normal'
            tind=at;aind=at;
          else
            [aind,tind]=find(combindex==combuse(at));
          end
          cfg=[];
          cfg.trials=tind;
          tmpt=ft_selectdata(cfg,tlock_nul_s0{ll+60,tt});
          cfg.trials=aind;
          tmpa=ft_selectdata(cfg,tlock_aud_s0{ll,tt});
          cfg=[];
          cfg.operation='add';
          cfg.parameter='trial';
          tmpsum=ft_math(cfg,tmpt,tmpa);
          if at==1
            tlock_fake=tmpsum;
          end
          tlock_fake.trial(at,:,:)=tmpsum.trial(1,:,:);
%           tlock_fake.trial(at,:,:)=tlock_aud_s0{ll,tt}.trial(aind,:,:)+tlock_nul_s0{ll+60,tt}.trial(tind,:,:);
        end
        cfg=[];
        cfg.method='mtmconvol';
        cfg.pad=4;
        cfg.foi=4:2:30;
        cfg.taper='hanning';
        cfg.toi=-1.2:0.1:1.3;
        cfg.t_ftimwin=4./cfg.foi;
        %         cfg.keeptrials='yes';  % not necessary at this stage for keeping trials (yes later for source analysis and phase resetting)
        freqlo_audMSpN_tmp{ll,tt,cc}=ft_freqanalysis(cfg,tlock_fake);
        cfg=[];
        cfg.method='mtmconvol';
        cfg.pad=4;
        cfg.foi=[30:5:45 55:5:80];
        cfg.taper='dpss';
        cfg.tapsmofrq=7*ones(1,length(cfg.foi));
        cfg.toi=-1.2:0.1:1.3;
        cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
        freqhi_audMSpN_tmp{ll,tt,cc}=ft_freqanalysis(cfg,tlock_fake);
      end
      clear tmp
      for cc=1:numtests,tmp(:,:,:,cc)=freqlo_audMSpN_tmp{ll,tt,cc}.powspctrm;end
      freqlo_audMSpN_comb{ll,tt,1}=freqlo_audMSpN_tmp{ll,tt,1};
      freqlo_audMSpN_comb{ll,tt,2}=freqlo_audMSpN_tmp{ll,tt,1};
      freqlo_audMSpN_comb{ll,tt,2}.powspctrm=mean(tmp,4);
      freqlo_audMSpN_comb{ll,tt,2}.stdspctrm=std(tmp,[],4);
      clear tmp
      for cc=1:numtests,tmp(:,:,:,cc)=freqhi_audMSpN_tmp{ll,tt,cc}.powspctrm;end
      freqhi_audMSpN_comb{ll,tt,1}=freqhi_audMSpN_tmp{ll,tt,1};
      freqhi_audMSpN_comb{ll,tt,2}=freqhi_audMSpN_tmp{ll,tt,1};
      freqhi_audMSpN_comb{ll,tt,2}.powspctrm=mean(tmp,4);
      freqhi_audMSpN_comb{ll,tt,2}.stdspctrm=std(tmp,[],4);
      clear freqlo_audMSpN_tmp freqhi_audMSpN_tmp
      
    end
    
    for ll=[soalist soalist+20 soalist+40]
      tlock_tac_s0{ll,tt}=[]; % clearing to help save memory as running
      tlock_aud_s0{ll,tt}=[]; % clearing to help save memory as running
      tlock_nul_s0{ll,tt}=[]; % clearing to help save memory as running
    end
  end
  
  % Question: do baseline correction prior to enter in to stats?
  
  save(['freq_diffs_averef_' sub{ii} '.mat'],'freq*P*','freq*N*','num*trials')
  
end


return



%% sensor level group stats

load elec1010_neighb.mat

allcond_sameN=1;
soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];
 tacbasemax=min(-.1,soades-.2);
 audbasemax=min(-.1,fliplr(soades)-.2);
 
 ylimlo=[4 7; 8 12; 14 30];
    ylimhi=[30 45; 55 80];


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
    subuseind=0;
    % Baseline correct each participant prior to entering to stats????
    for ii=subuse
      cd([edir sub{ii} ])
      load(['freq_diffs_averef_' sub{ii} '.mat']);
      numtrt(ll,tt,ii-submin)=numt_trials(ll,tt);
      numtra(ll,tt,ii-submin)=numa_trials(ll,tt);
      if numtrt(ll,tt,ii-submin)<20 % what is best number to use here?
        subuse=setdiff(subuse,ii);
      else
        subuseind=subuseind+1;
        freqloall_tacPaud_comb1{subuseind}=freqlo_tacPaud_comb{ll,tt,1};
        freqloall_audPtac_comb1{subuseind}=freqlo_audPtac_comb{ll,tt,1};
        freqloall_tacMSpN_comb1{subuseind}=freqlo_tacMSpN_comb{ll,tt,1};
        freqloall_audMSpN_comb1{subuseind}=freqlo_audMSpN_comb{ll,tt,1};
        freqhiall_tacPaud_comb1{subuseind}=freqhi_tacPaud_comb{ll,tt,1};
        freqhiall_audPtac_comb1{subuseind}=freqhi_audPtac_comb{ll,tt,1};
        freqhiall_tacMSpN_comb1{subuseind}=freqhi_tacMSpN_comb{ll,tt,1};
        freqhiall_audMSpN_comb1{subuseind}=freqhi_audMSpN_comb{ll,tt,1};
        freqloall_tacPaud_comb2{subuseind}=freqlo_tacPaud_comb{ll,tt,2};
        freqloall_audPtac_comb2{subuseind}=freqlo_audPtac_comb{ll,tt,2};
        freqloall_tacMSpN_comb2{subuseind}=freqlo_tacMSpN_comb{ll,tt,2};
        freqloall_audMSpN_comb2{subuseind}=freqlo_audMSpN_comb{ll,tt,2};
        freqhiall_tacPaud_comb2{subuseind}=freqhi_tacPaud_comb{ll,tt,2};
        freqhiall_audPtac_comb2{subuseind}=freqhi_audPtac_comb{ll,tt,2};
        freqhiall_tacMSpN_comb2{subuseind}=freqhi_tacMSpN_comb{ll,tt,2};
        freqhiall_audMSpN_comb2{subuseind}=freqhi_audMSpN_comb{ll,tt,2};
        freqloall_tacPaud_fftadd{subuseind}=freqlo_tacPaud_fftadd{ll,tt};
        freqloall_audPtac_fftadd{subuseind}=freqlo_audPtac_fftadd{ll,tt};
        freqloall_tacMSpN_fftadd{subuseind}=freqlo_tacMSpN_fftadd{ll,tt};
        freqloall_audMSpN_fftadd{subuseind}=freqlo_audMSpN_fftadd{ll,tt};
        freqhiall_tacPaud_fftadd{subuseind}=freqhi_tacPaud_fftadd{ll,tt};
        freqhiall_audPtac_fftadd{subuseind}=freqhi_audPtac_fftadd{ll,tt};
        freqhiall_tacMSpN_fftadd{subuseind}=freqhi_tacMSpN_fftadd{ll,tt};
        freqhiall_audMSpN_fftadd{subuseind}=freqhi_audMSpN_fftadd{ll,tt};
      end
      clear freqlo_* freqhi_*
    end
    
    if 0
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
    end
    
    
    for ii=1:length(subuse)
      
      cfg=[];
      cfg.baselinetype='relchange';
      cfg.baseline=[-1.3 tacbasemax(ll)];
      freqloall_tacPaud_fftaddbn{ii}=ft_freqbaseline(cfg,freqloall_tacPaud_fftadd{ii});
      freqloall_tacMSpN_fftaddbn{ii}=ft_freqbaseline(cfg,freqloall_tacMSpN_fftadd{ii});
      freqhiall_tacPaud_fftaddbn{ii}=ft_freqbaseline(cfg,freqhiall_tacPaud_fftadd{ii});
      freqhiall_tacMSpN_fftaddbn{ii}=ft_freqbaseline(cfg,freqhiall_tacMSpN_fftadd{ii});
      freqloall_tacPaud_comb1bn{ii}=ft_freqbaseline(cfg,freqloall_tacPaud_comb1{ii});
      freqloall_tacMSpN_comb1bn{ii}=ft_freqbaseline(cfg,freqloall_tacMSpN_comb1{ii});
      freqhiall_tacPaud_comb1bn{ii}=ft_freqbaseline(cfg,freqhiall_tacPaud_comb1{ii});
      freqhiall_tacMSpN_comb1bn{ii}=ft_freqbaseline(cfg,freqhiall_tacMSpN_comb1{ii});

      cfg.baseline=[-1.3 audbasemax(ll)];
      freqloall_audPtac_fftaddbn{ii}=ft_freqbaseline(cfg,freqloall_audPtac_fftadd{ii});
      freqloall_audMSpN_fftaddbn{ii}=ft_freqbaseline(cfg,freqloall_audMSpN_fftadd{ii});
      freqhiall_audPtac_fftaddbn{ii}=ft_freqbaseline(cfg,freqhiall_audPtac_fftadd{ii});
      freqhiall_audMSpN_fftaddbn{ii}=ft_freqbaseline(cfg,freqhiall_audMSpN_fftadd{ii});      
      freqloall_audPtac_comb1bn{ii}=ft_freqbaseline(cfg,freqloall_audPtac_comb1{ii});
      freqloall_audMSpN_comb1bn{ii}=ft_freqbaseline(cfg,freqloall_audMSpN_comb1{ii});
      freqhiall_audPtac_comb1bn{ii}=ft_freqbaseline(cfg,freqhiall_audPtac_comb1{ii});
      freqhiall_audMSpN_comb1bn{ii}=ft_freqbaseline(cfg,freqhiall_audMSpN_comb1{ii});

      %       baselinetime=[dsearchn(freqloall_tacPaud_fftadd{ii}.time',-1.2) dsearchn(freqloall_tacPaud_fftadd{ii}.time',-0.7)];
%       freqloall_tacPaud_fabn{ii}=freqloall_tacPaud_fftadd{ii};
%       freqloall_tacPaud_fabn{ii}.powspctrm=freqloall_tacPaud_fabn{ii}.powspctrm./repmat(nanmean(freqloall_tacPaud_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_fabn{ii}.time)]);
%       freqloall_tacMSpN_fabn{ii}=freqloall_tacMSpN_fftadd{ii};
%       freqloall_tacMSpN_fabn{ii}.powspctrm=freqloall_tacMSpN_fabn{ii}.powspctrm./repmat(nanmean(freqloall_tacMSpN_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_fabn{ii}.time)]);
%       freqloall_audPtac_fabn{ii}=freqloall_audPtac_fftadd{ii};
%       freqloall_audPtac_fabn{ii}.powspctrm=freqloall_audPtac_fabn{ii}.powspctrm./repmat(nanmean(freqloall_audPtac_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_fabn{ii}.time)]);
%       freqloall_audMSpN_fabn{ii}=freqloall_audMSpN_fftadd{ii};
%       freqloall_audMSpN_fabn{ii}.powspctrm=freqloall_audMSpN_fabn{ii}.powspctrm./repmat(nanmean(freqloall_audMSpN_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_fabn{ii}.time)]);
%       freqhiall_tacPaud_fabn{ii}=freqhiall_tacPaud_fftadd{ii};
%       freqhiall_tacPaud_fabn{ii}.powspctrm=freqhiall_tacPaud_fabn{ii}.powspctrm./repmat(nanmean(freqhiall_tacPaud_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_fabn{ii}.time)]);
%       freqhiall_tacMSpN_fabn{ii}=freqhiall_tacMSpN_fftadd{ii};
%       freqhiall_tacMSpN_fabn{ii}.powspctrm=freqhiall_tacMSpN_fabn{ii}.powspctrm./repmat(nanmean(freqhiall_tacMSpN_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_fabn{ii}.time)]);
%       freqhiall_audPtac_fabn{ii}=freqhiall_audPtac_fftadd{ii};
%       freqhiall_audPtac_fabn{ii}.powspctrm=freqhiall_audPtac_fabn{ii}.powspctrm./repmat(nanmean(freqhiall_audPtac_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_fabn{ii}.time)]);
%       freqhiall_audMSpN_fabn{ii}=freqhiall_audMSpN_fftadd{ii};
%       freqhiall_audMSpN_fabn{ii}.powspctrm=freqhiall_audMSpN_fabn{ii}.powspctrm./repmat(nanmean(freqhiall_audMSpN_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_fabn{ii}.time)]);
%       baselinetime=[dsearchn(freqloall_tacPaud_comb1{ii}.time',-1.2) dsearchn(freqloall_tacPaud_comb1{ii}.time',-0.7)];
%       freqloall_tacPaud_cbbn{ii}=freqloall_tacPaud_comb1{ii};
%       freqloall_tacPaud_cbbn{ii}.powspctrm=freqloall_tacPaud_cbbn{ii}.powspctrm./repmat(nanmean(freqloall_tacPaud_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_cbbn{ii}.time)]);
%       freqloall_tacMSpN_cbbn{ii}=freqloall_tacMSpN_comb1{ii};
%       freqloall_tacMSpN_cbbn{ii}.powspctrm=freqloall_tacMSpN_cbbn{ii}.powspctrm./repmat(nanmean(freqloall_tacMSpN_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_cbbn{ii}.time)]);
%       freqloall_audPtac_cbbn{ii}=freqloall_audPtac_comb1{ii};
%       freqloall_audPtac_cbbn{ii}.powspctrm=freqloall_audPtac_cbbn{ii}.powspctrm./repmat(nanmean(freqloall_audPtac_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_cbbn{ii}.time)]);
%       freqloall_audMSpN_cbbn{ii}=freqloall_audMSpN_comb1{ii};
%       freqloall_audMSpN_cbbn{ii}.powspctrm=freqloall_audMSpN_cbbn{ii}.powspctrm./repmat(nanmean(freqloall_audMSpN_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_cbbn{ii}.time)]);
%       freqhiall_tacPaud_cbbn{ii}=freqhiall_tacPaud_comb1{ii};
%       freqhiall_tacPaud_cbbn{ii}.powspctrm=freqhiall_tacPaud_cbbn{ii}.powspctrm./repmat(nanmean(freqhiall_tacPaud_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_cbbn{ii}.time)]);
%       freqhiall_tacMSpN_cbbn{ii}=freqhiall_tacMSpN_comb1{ii};
%       freqhiall_tacMSpN_cbbn{ii}.powspctrm=freqhiall_tacMSpN_cbbn{ii}.powspctrm./repmat(nanmean(freqhiall_tacMSpN_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_cbbn{ii}.time)]);
%       freqhiall_audPtac_cbbn{ii}=freqhiall_audPtac_comb1{ii};
%       freqhiall_audPtac_cbbn{ii}.powspctrm=freqhiall_audPtac_cbbn{ii}.powspctrm./repmat(nanmean(freqhiall_audPtac_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_cbbn{ii}.time)]);
%       freqhiall_audMSpN_cbbn{ii}=freqhiall_audMSpN_comb1{ii};
%       freqhiall_audMSpN_cbbn{ii}.powspctrm=freqhiall_audMSpN_cbbn{ii}.powspctrm./repmat(nanmean(freqhiall_audMSpN_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_cbbn{ii}.time)]);
    end
    
    cfg=[];
    cfg.keepindividual='yes';
    grindlo_tacPaud_fftadd=ft_freqgrandaverage(cfg,freqloall_tacPaud_fftadd{:});
    grindlo_tacMSpN_fftadd=ft_freqgrandaverage(cfg,freqloall_tacMSpN_fftadd{:});
    grindlo_audPtac_fftadd=ft_freqgrandaverage(cfg,freqloall_audPtac_fftadd{:});
    grindlo_audMSpN_fftadd=ft_freqgrandaverage(cfg,freqloall_audMSpN_fftadd{:});
    grindhi_tacPaud_fftadd=ft_freqgrandaverage(cfg,freqhiall_tacPaud_fftadd{:});
    grindhi_audPtac_fftadd=ft_freqgrandaverage(cfg,freqhiall_audPtac_fftadd{:});
    grindhi_tacMSpN_fftadd=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_fftadd{:});
    grindhi_audMSpN_fftadd=ft_freqgrandaverage(cfg,freqhiall_audMSpN_fftadd{:});
    
    grindlo_tacPaud_comb1=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb1{:});
    grindlo_tacMSpN_comb1=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb1{:});
    grindlo_audPtac_comb1=ft_freqgrandaverage(cfg,freqloall_audPtac_comb1{:});
    grindlo_audMSpN_comb1=ft_freqgrandaverage(cfg,freqloall_audMSpN_comb1{:});
    grindhi_tacPaud_comb1=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb1{:});
    grindhi_audPtac_comb1=ft_freqgrandaverage(cfg,freqhiall_audPtac_comb1{:});
    grindhi_tacMSpN_comb1=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb1{:});
    grindhi_audMSpN_comb1=ft_freqgrandaverage(cfg,freqhiall_audMSpN_comb1{:});
    
    grindlo_tacPaud_comb2=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb2{:});
    grindlo_tacMSpN_comb2=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb2{:});
    grindlo_audPtac_comb2=ft_freqgrandaverage(cfg,freqloall_audPtac_comb2{:});
    grindlo_audMSpN_comb2=ft_freqgrandaverage(cfg,freqloall_audMSpN_comb2{:});
    grindhi_tacPaud_comb2=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb2{:});
    grindhi_audPtac_comb2=ft_freqgrandaverage(cfg,freqhiall_audPtac_comb2{:});
    grindhi_tacMSpN_comb2=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb2{:});
    grindhi_audMSpN_comb2=ft_freqgrandaverage(cfg,freqhiall_audMSpN_comb2{:});

    grindlo_tacPaud_fftaddbn=ft_freqgrandaverage(cfg,freqloall_tacPaud_fftaddbn{:});
    grindlo_tacMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqloall_tacMSpN_fftaddbn{:});
    grindlo_audPtac_fftaddbn=ft_freqgrandaverage(cfg,freqloall_audPtac_fftaddbn{:});
    grindlo_audMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqloall_audMSpN_fftaddbn{:});
    grindhi_tacPaud_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_tacPaud_fftaddbn{:});
    grindhi_audPtac_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_audPtac_fftaddbn{:});
    grindhi_tacMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_fftaddbn{:});
    grindhi_audMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_audMSpN_fftaddbn{:});
    
    grindlo_tacPaud_comb1bn=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb1bn{:});
    grindlo_tacMSpN_comb1bn=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb1bn{:});
    grindlo_audPtac_comb1bn=ft_freqgrandaverage(cfg,freqloall_audPtac_comb1bn{:});
    grindlo_audMSpN_comb1bn=ft_freqgrandaverage(cfg,freqloall_audMSpN_comb1bn{:});
    grindhi_tacPaud_comb1bn=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb1bn{:});
    grindhi_audPtac_comb1bn=ft_freqgrandaverage(cfg,freqhiall_audPtac_comb1bn{:});
    grindhi_tacMSpN_comb1bn=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb1bn{:});
    grindhi_audMSpN_comb1bn=ft_freqgrandaverage(cfg,freqhiall_audMSpN_comb1bn{:});

    cfg=[];
    gravelo_tacPaud_fftadd=ft_freqgrandaverage(cfg,freqloall_tacPaud_fftadd{:});
    gravelo_tacMSpN_fftadd=ft_freqgrandaverage(cfg,freqloall_tacMSpN_fftadd{:});
    gravelo_audPtac_fftadd=ft_freqgrandaverage(cfg,freqloall_audPtac_fftadd{:});
    gravelo_audMSpN_fftadd=ft_freqgrandaverage(cfg,freqloall_audMSpN_fftadd{:});
    gravehi_tacPaud_fftadd=ft_freqgrandaverage(cfg,freqhiall_tacPaud_fftadd{:});
    gravehi_audPtac_fftadd=ft_freqgrandaverage(cfg,freqhiall_audPtac_fftadd{:});
    gravehi_tacMSpN_fftadd=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_fftadd{:});
    gravehi_audMSpN_fftadd=ft_freqgrandaverage(cfg,freqhiall_audMSpN_fftadd{:});
    
    gravelo_tacPaud_comb1=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb1{:});
    gravelo_tacMSpN_comb1=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb1{:});
    gravelo_audPtac_comb1=ft_freqgrandaverage(cfg,freqloall_audPtac_comb1{:});
    gravelo_audMSpN_comb1=ft_freqgrandaverage(cfg,freqloall_audMSpN_comb1{:});
    gravehi_tacPaud_comb1=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb1{:});
    gravehi_audPtac_comb1=ft_freqgrandaverage(cfg,freqhiall_audPtac_comb1{:});
    gravehi_tacMSpN_comb1=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb1{:});
    gravehi_audMSpN_comb1=ft_freqgrandaverage(cfg,freqhiall_audMSpN_comb1{:});
    
    gravelo_tacPaud_comb2=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb2{:});
    gravelo_tacMSpN_comb2=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb2{:});
    gravelo_audPtac_comb2=ft_freqgrandaverage(cfg,freqloall_audPtac_comb2{:});
    gravelo_audMSpN_comb2=ft_freqgrandaverage(cfg,freqloall_audMSpN_comb2{:});
    gravehi_tacPaud_comb2=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb2{:});
    gravehi_audPtac_comb2=ft_freqgrandaverage(cfg,freqhiall_audPtac_comb2{:});
    gravehi_tacMSpN_comb2=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb2{:});
    gravehi_audMSpN_comb2=ft_freqgrandaverage(cfg,freqhiall_audMSpN_comb2{:});

    gravelo_tacPaud_fftaddbn=ft_freqgrandaverage(cfg,freqloall_tacPaud_fftaddbn{:});
    gravelo_tacMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqloall_tacMSpN_fftaddbn{:});
    gravelo_audPtac_fftaddbn=ft_freqgrandaverage(cfg,freqloall_audPtac_fftaddbn{:});
    gravelo_audMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqloall_audMSpN_fftaddbn{:});
    gravehi_tacPaud_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_tacPaud_fftaddbn{:});
    gravehi_audPtac_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_audPtac_fftaddbn{:});
    gravehi_tacMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_fftaddbn{:});
    gravehi_audMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_audMSpN_fftaddbn{:});
    
    gravelo_tacPaud_comb1bn=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb1bn{:});
    gravelo_tacMSpN_comb1bn=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb1bn{:});
    gravelo_audPtac_comb1bn=ft_freqgrandaverage(cfg,freqloall_audPtac_comb1bn{:});
    gravelo_audMSpN_comb1bn=ft_freqgrandaverage(cfg,freqloall_audMSpN_comb1bn{:});
    gravehi_tacPaud_comb1bn=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb1bn{:});
    gravehi_audPtac_comb1bn=ft_freqgrandaverage(cfg,freqhiall_audPtac_comb1bn{:});
    gravehi_tacMSpN_comb1bn=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb1bn{:});
    gravehi_audMSpN_comb1bn=ft_freqgrandaverage(cfg,freqhiall_audMSpN_comb1bn{:});

    
    topoplotTFR_highlight(11,gravelo_tacPaud_fftadd,[],ylimlo(2,:),[-1.2 tacbasemax(ll)],[]);
    topoplotTFR_highlight(12,gravelo_tacMSpN_fftadd,[],ylimlo(2,:),[-1.2 tacbasemax(ll)],[]);
    topoplotTFR_highlight(13,gravehi_tacPaud_fftadd,[],ylimhi(2,:),[-1.2 tacbasemax(ll)],[]);
    topoplotTFR_highlight(14,gravehi_tacMSpN_fftadd,[],ylimhi(2,:),[-1.2 tacbasemax(ll)],[]);
    topoplotTFR_highlight(15,gravelo_audPtac_fftadd,[],ylimlo(2,:),[-1.2 audbasemax(ll)],[]);
    topoplotTFR_highlight(16,gravelo_audMSpN_fftadd,[],ylimlo(2,:),[-1.2 audbasemax(ll)],[]);
    topoplotTFR_highlight(17,gravehi_audPtac_fftadd,[],ylimhi(2,:),[-1.2 audbasemax(ll)],[]);
    topoplotTFR_highlight(18,gravehi_audMSpN_fftadd,[],ylimhi(2,:),[-1.2 audbasemax(ll)],[]);

    topoplotTFR_highlight(21,gravelo_tacPaud_comb1,[],ylimlo(2,:),[-1.2 tacbasemax(ll)],[]);
    topoplotTFR_highlight(22,gravelo_tacMSpN_comb1,[],ylimlo(2,:),[-1.2 tacbasemax(ll)],[]);
    topoplotTFR_highlight(23,gravehi_tacPaud_comb1,[],ylimhi(2,:),[-1.2 tacbasemax(ll)],[]);
    topoplotTFR_highlight(24,gravehi_tacMSpN_comb1,[],ylimhi(2,:),[-1.2 tacbasemax(ll)],[]);
    topoplotTFR_highlight(25,gravelo_audPtac_comb1,[],ylimlo(2,:),[-1.2 audbasemax(ll)],[]);
    topoplotTFR_highlight(26,gravelo_audMSpN_comb1,[],ylimlo(2,:),[-1.2 audbasemax(ll)],[]);
    topoplotTFR_highlight(27,gravehi_audPtac_comb1,[],ylimhi(2,:),[-1.2 audbasemax(ll)],[]);
    topoplotTFR_highlight(28,gravehi_audMSpN_comb1,[],ylimhi(2,:),[-1.2 audbasemax(ll)],[]);
    
    topoplotTFR_highlight(31,gravelo_tacPaud_comb2,[],ylimlo(2,:),[-1.2 tacbasemax(ll)],[]);
    topoplotTFR_highlight(32,gravelo_tacMSpN_comb2,[],ylimlo(2,:),[-1.2 tacbasemax(ll)],[]);
    topoplotTFR_highlight(33,gravehi_tacPaud_comb2,[],ylimhi(2,:),[-1.2 tacbasemax(ll)],[]);
    topoplotTFR_highlight(34,gravehi_tacMSpN_comb2,[],ylimhi(2,:),[-1.2 tacbasemax(ll)],[]);
    topoplotTFR_highlight(35,gravelo_audPtac_comb2,[],ylimlo(2,:),[-1.2 audbasemax(ll)],[]);
    topoplotTFR_highlight(36,gravelo_audMSpN_comb2,[],ylimlo(2,:),[-1.2 audbasemax(ll)],[]);
    topoplotTFR_highlight(37,gravehi_audPtac_comb2,[],ylimhi(2,:),[-1.2 audbasemax(ll)],[]);
    topoplotTFR_highlight(38,gravehi_audMSpN_comb2,[],ylimhi(2,:),[-1.2 audbasemax(ll)],[]);
    
    topoplotTFR_highlight(41,gravelo_tacPaud_fftaddbn,[],ylimlo(2,:),[-1.2 tacbasemax(ll)],[]);
    topoplotTFR_highlight(42,gravelo_tacMSpN_fftaddbn,[],ylimlo(2,:),[-1.2 tacbasemax(ll)],[]);
    topoplotTFR_highlight(43,gravehi_tacPaud_fftaddbn,[],ylimhi(2,:),[-1.2 tacbasemax(ll)],[]);
    topoplotTFR_highlight(44,gravehi_tacMSpN_fftaddbn,[],ylimhi(2,:),[-1.2 tacbasemax(ll)],[]);
    topoplotTFR_highlight(45,gravelo_audPtac_fftaddbn,[],ylimlo(2,:),[-1.2 audbasemax(ll)],[]);
    topoplotTFR_highlight(46,gravelo_audMSpN_fftaddbn,[],ylimlo(2,:),[-1.2 audbasemax(ll)],[]);
    topoplotTFR_highlight(47,gravehi_audPtac_fftaddbn,[],ylimhi(2,:),[-1.2 audbasemax(ll)],[]);
    topoplotTFR_highlight(48,gravehi_audMSpN_fftaddbn,[],ylimhi(2,:),[-1.2 audbasemax(ll)],[]);
    
    topoplotTFR_highlight(51,gravelo_tacPaud_comb1bn,[],ylimlo(2,:),[-1.2 tacbasemax(ll)],[]);
    topoplotTFR_highlight(52,gravelo_tacMSpN_comb1bn,[],ylimlo(2,:),[-1.2 tacbasemax(ll)],[]);
    topoplotTFR_highlight(53,gravehi_tacPaud_comb1bn,[],ylimhi(2,:),[-1.2 tacbasemax(ll)],[]);
    topoplotTFR_highlight(54,gravehi_tacMSpN_comb1bn,[],ylimhi(2,:),[-1.2 tacbasemax(ll)],[]);
    topoplotTFR_highlight(55,gravelo_audPtac_comb1bn,[],ylimlo(2,:),[-1.2 audbasemax(ll)],[]);
    topoplotTFR_highlight(56,gravelo_audMSpN_comb1bn,[],ylimlo(2,:),[-1.2 audbasemax(ll)],[]);
    topoplotTFR_highlight(57,gravehi_audPtac_comb1bn,[],ylimhi(2,:),[-1.2 audbasemax(ll)],[]);
    topoplotTFR_highlight(58,gravehi_audMSpN_comb1bn,[],ylimhi(2,:),[-1.2 audbasemax(ll)],[]);
    
    
    for ii=1:length(subuse)
      cfg=[];
      cfg.operation='subtract';
      cfg.parameter='powspctrm';
      freqloall_TPA_MSPN_fftadd{ii}=ft_math(cfg,freqloall_tacPaud_fftadd{ii},freqloall_tacMSpN_fftadd{ii});
      freqloall_APT_MSPN_fftadd{ii}=ft_math(cfg,freqloall_audPtac_fftadd{ii},freqloall_audMSpN_fftadd{ii});
      freqhiall_TPA_MSPN_fftadd{ii}=ft_math(cfg,freqhiall_tacPaud_fftadd{ii},freqhiall_tacMSpN_fftadd{ii});
      freqhiall_APT_MSPN_fftadd{ii}=ft_math(cfg,freqhiall_audPtac_fftadd{ii},freqhiall_audMSpN_fftadd{ii});
      freqloall_TPA_MSPN_comb1{ii}=ft_math(cfg,freqloall_tacPaud_comb1{ii},freqloall_tacMSpN_comb1{ii});
      freqloall_APT_MSPN_comb1{ii}=ft_math(cfg,freqloall_audPtac_comb1{ii},freqloall_audMSpN_comb1{ii});
      freqhiall_TPA_MSPN_comb1{ii}=ft_math(cfg,freqhiall_tacPaud_comb1{ii},freqhiall_tacMSpN_comb1{ii});
      freqhiall_APT_MSPN_comb1{ii}=ft_math(cfg,freqhiall_audPtac_comb1{ii},freqhiall_audMSpN_comb1{ii});
      freqloall_TPA_MSPN_comb2{ii}=ft_math(cfg,freqloall_tacPaud_comb2{ii},freqloall_tacMSpN_comb2{ii});
      freqloall_APT_MSPN_comb2{ii}=ft_math(cfg,freqloall_audPtac_comb2{ii},freqloall_audMSpN_comb2{ii});
      freqhiall_TPA_MSPN_comb2{ii}=ft_math(cfg,freqhiall_tacPaud_comb2{ii},freqhiall_tacMSpN_comb2{ii});
      freqhiall_APT_MSPN_comb2{ii}=ft_math(cfg,freqhiall_audPtac_comb2{ii},freqhiall_audMSpN_comb2{ii});
      freqloall_TPA_MSPN_fftaddbn{ii}=ft_math(cfg,freqloall_tacPaud_fftaddbn{ii},freqloall_tacMSpN_fftaddbn{ii});
      freqloall_APT_MSPN_fftaddbn{ii}=ft_math(cfg,freqloall_audPtac_fftaddbn{ii},freqloall_audMSpN_fftaddbn{ii});
      freqhiall_TPA_MSPN_fftaddbn{ii}=ft_math(cfg,freqhiall_tacPaud_fftaddbn{ii},freqhiall_tacMSpN_fftaddbn{ii});
      freqhiall_APT_MSPN_fftaddbn{ii}=ft_math(cfg,freqhiall_audPtac_fftaddbn{ii},freqhiall_audMSpN_fftaddbn{ii});
      freqloall_TPA_MSPN_comb1bn{ii}=ft_math(cfg,freqloall_tacPaud_comb1bn{ii},freqloall_tacMSpN_comb1bn{ii});
      freqloall_APT_MSPN_comb1bn{ii}=ft_math(cfg,freqloall_audPtac_comb1bn{ii},freqloall_audMSpN_comb1bn{ii});
      freqhiall_TPA_MSPN_comb1bn{ii}=ft_math(cfg,freqhiall_tacPaud_comb1bn{ii},freqhiall_tacMSpN_comb1bn{ii});
      freqhiall_APT_MSPN_comb1bn{ii}=ft_math(cfg,freqhiall_audPtac_comb1bn{ii},freqhiall_audMSpN_comb1bn{ii});
    end
    cfg=[];
    gravelo_TPA_MSPN_fftadd=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_fftadd{:});
    gravelo_APT_MSPN_fftadd=ft_freqgrandaverage(cfg,freqloall_APT_MSPN_fftadd{:});
    gravehi_TPA_MSPN_fftadd=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_fftadd{:});
    gravehi_APT_MSPN_fftadd=ft_freqgrandaverage(cfg,freqhiall_APT_MSPN_fftadd{:});
    gravelo_TPA_MSPN_comb1=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb1{:});
    gravelo_APT_MSPN_comb1=ft_freqgrandaverage(cfg,freqloall_APT_MSPN_comb1{:});
    gravehi_TPA_MSPN_comb1=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb1{:});
    gravehi_APT_MSPN_comb1=ft_freqgrandaverage(cfg,freqhiall_APT_MSPN_comb1{:});
    gravelo_TPA_MSPN_comb2=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb2{:});
    gravelo_APT_MSPN_comb2=ft_freqgrandaverage(cfg,freqloall_APT_MSPN_comb2{:});
    gravehi_TPA_MSPN_comb2=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb2{:});
    gravehi_APT_MSPN_comb2=ft_freqgrandaverage(cfg,freqhiall_APT_MSPN_comb2{:});
    gravelo_TPA_MSPN_fftaddbn=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_fftaddbn{:});
    gravelo_APT_MSPN_fftaddbn=ft_freqgrandaverage(cfg,freqloall_APT_MSPN_fftaddbn{:});
    gravehi_TPA_MSPN_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_fftaddbn{:});
    gravehi_APT_MSPN_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_APT_MSPN_fftaddbn{:});
    gravelo_TPA_MSPN_comb1bn=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb1bn{:});
    gravelo_APT_MSPN_comb1bn=ft_freqgrandaverage(cfg,freqloall_APT_MSPN_comb1bn{:});
    gravehi_TPA_MSPN_comb1bn=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb1bn{:});
    gravehi_APT_MSPN_comb1bn=ft_freqgrandaverage(cfg,freqhiall_APT_MSPN_comb1bn{:});
    
    
    nsub=length(freqloall_audMSpN_fftadd);
    
    cfg=[];
    cfg.latency=[-.8 0.8];
    cfg.neighbours=neighbours;
    cfg.parameter='powspctrm';
    cfg.method='montecarlo';
    cfg.numrandomization=2000;
    % cfg.correctm='holm';
    cfg.correctm='cluster';
    cfg.clusteralpha = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan = 2;
    cfg.statistic='depsamplesT';
    cfg.design=zeros(2,2*nsub);
    cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
    cfg.design(2,:)=[1:nsub 1:nsub];
    cfg.ivar=1;
    cfg.uvar=2;
    stattl_mc_fftadd=ft_freqstatistics(cfg, grindlo_tacPaud_fftadd, grindlo_tacMSpN_fftadd);
    statal_mc_fftadd=ft_freqstatistics(cfg, grindlo_audPtac_fftadd, grindlo_audMSpN_fftadd);
    statth_mc_fftadd=ft_freqstatistics(cfg, grindhi_tacPaud_fftadd, grindhi_tacMSpN_fftadd);
    statah_mc_fftadd=ft_freqstatistics(cfg, grindhi_audPtac_fftadd, grindhi_audMSpN_fftadd);
    stattl_mc_comb1=ft_freqstatistics(cfg, grindlo_tacPaud_comb1, grindlo_tacMSpN_comb1);
    statal_mc_comb1=ft_freqstatistics(cfg, grindlo_audPtac_comb1, grindlo_audMSpN_comb1);
    statth_mc_comb1=ft_freqstatistics(cfg, grindhi_tacPaud_comb1, grindhi_tacMSpN_comb1);
    statah_mc_comb1=ft_freqstatistics(cfg, grindhi_audPtac_comb1, grindhi_audMSpN_comb1);
    stattl_mc_comb2=ft_freqstatistics(cfg, grindlo_tacPaud_comb2, grindlo_tacMSpN_comb2);
    statal_mc_comb2=ft_freqstatistics(cfg, grindlo_audPtac_comb2, grindlo_audMSpN_comb2);
    statth_mc_comb2=ft_freqstatistics(cfg, grindhi_tacPaud_comb2, grindhi_tacMSpN_comb2);
    statah_mc_comb2=ft_freqstatistics(cfg, grindhi_audPtac_comb2, grindhi_audMSpN_comb2);
    stattl_mc_fftaddbn=ft_freqstatistics(cfg, grindlo_tacPaud_fftaddbn, grindlo_tacMSpN_fftaddbn);
    statal_mc_fftaddbn=ft_freqstatistics(cfg, grindlo_audPtac_fftaddbn, grindlo_audMSpN_fftaddbn);
    statth_mc_fftaddbn=ft_freqstatistics(cfg, grindhi_tacPaud_fftaddbn, grindhi_tacMSpN_fftaddbn);
    statah_mc_fftaddbn=ft_freqstatistics(cfg, grindhi_audPtac_fftaddbn, grindhi_audMSpN_fftaddbn);
    stattl_mc_comb1bn=ft_freqstatistics(cfg, grindlo_tacPaud_comb1bn, grindlo_tacMSpN_comb1bn);
    statal_mc_comb1bn=ft_freqstatistics(cfg, grindlo_audPtac_comb1bn, grindlo_audMSpN_comb1bn);
    statth_mc_comb1bn=ft_freqstatistics(cfg, grindhi_tacPaud_comb1bn, grindhi_tacMSpN_comb1bn);
    statah_mc_comb1bn=ft_freqstatistics(cfg, grindhi_audPtac_comb1bn, grindhi_audMSpN_comb1bn);
    
    
    timwin=[];
    base=[];
    for ss=1:size(ylimlo,1)
      ylim=ylimlo(ss,:);
      
      try
        fig=140+ss;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_fftadd,timwin,ylim,base,stattl_mc_fftadd);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_ta_fftadd_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_ta_fftadd_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a.png'],'-dpng')
        end
      catch
      end
      try
        fig=240+ss;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb1,timwin,ylim,base,stattl_mc_comb1);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_ta_comb1_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_ta_comb1_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a.png'],'-dpng')
        end
      catch
      end
      try
        fig=340+ss;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb2,timwin,ylim,base,stattl_mc_comb2);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_ta_comb2_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_ta_comb2_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a.png'],'-dpng')
        end
      catch
      end
      try
        fig=440+ss;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_fftaddbn,timwin,ylim,base,stattl_mc_fftaddbn);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_ta_fftaddbn_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25_bn.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_ta_fftaddbn_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a_bn.png'],'-dpng')
        end
      catch
      end
      try
        fig=540+ss;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb1bn,timwin,ylim,base,stattl_mc_comb1bn);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_ta_comb1bn_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25_bn.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_ta_comb1bn_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a_bn.png'],'-dpng')
        end
      catch
      end
      
      try
        fig=160+ss;
        topoplotTFR_highlight(fig,gravelo_APT_MSPN_fftadd,timwin,ylim,base,statal_mc_fftadd);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_at_fftadd_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_at_fftadd_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a.png'],'-dpng')
        end
      catch
      end
      try
        fig=260+ss;
        topoplotTFR_highlight(fig,gravelo_APT_MSPN_comb1,timwin,ylim,base,statal_mc_comb1);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_at_comb1_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_at_comb1_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a.png'],'-dpng')
        end
      catch
      end
      try
        fig=360+ss;
        topoplotTFR_highlight(fig,gravelo_APT_MSPN_comb2,timwin,ylim,base,statal_mc_comb2);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_at_comb2_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_at_comb2_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a.png'],'-dpng')
        end
      catch
      end
      try
        fig=460+ss;
        topoplotTFR_highlight(fig,gravelo_APT_MSPN_fftaddbn,timwin,ylim,base,statal_mc_fftaddbn);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_at_fftaddbn_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25_bn.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_at_fftaddbn_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a_bn.png'],'-dpng')
        end
      catch
      end
      try
        fig=560+ss;
        topoplotTFR_highlight(fig,gravelo_APT_MSPN_comb1bn,timwin,ylim,base,statal_mc_comb1bn);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_at_comb1bn_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25_bn.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdifflo_topoOverTime_mc_at_comb1bn_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a_bn.png'],'-dpng')
        end
      catch
      end
      
    end
    
    for ss=1:size(ylimhi,1)
      ylim=ylimhi(ss,:);
      try
        fig=150+ss;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_fftadd,timwin,ylim,base,statth_mc_fftadd);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_ta_fftadd_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_ta_fftadd_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a.png'],'-dpng')
        end
      catch
      end      
      try
        fig=250+ss;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb1,timwin,ylim,base,statth_mc_comb1);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_ta_comb1_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_ta_comb1_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a.png'],'-dpng')
        end
      catch
      end      
      try
        fig=350+ss;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb2,timwin,ylim,base,statth_mc_comb2);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_ta_comb2_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_ta_comb2_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a.png'],'-dpng')
        end
      catch
      end      
      try
        fig=450+ss;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_fftaddbn,timwin,ylim,base,statth_mc_fftaddbn);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_ta_fftaddbn_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25_bn.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_ta_fftaddbn_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a_bn.png'],'-dpng')
        end
      catch
      end
      try
        fig=550+ss;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb1bn,timwin,ylim,base,statth_mc_comb1bn);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_ta_comb1bn_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25_bn.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_ta_comb1bn_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a_bn.png'],'-dpng')
        end
      catch
      end
      
      try
        fig=170+ss;
        topoplotTFR_highlight(fig,gravehi_APT_MSPN_fftadd,timwin,ylim,base,statah_mc_fftadd);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_at_fftadd_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_at_fftadd_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a.png'],'-dpng')
        end
      catch
      end
      try
        fig=270+ss;
        topoplotTFR_highlight(fig,gravehi_APT_MSPN_comb1,timwin,ylim,base,statah_mc_comb1);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_at_comb1_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_at_comb1_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a.png'],'-dpng')
        end
      catch
      end
      try
        fig=370+ss;
        topoplotTFR_highlight(fig,gravehi_APT_MSPN_comb2,timwin,ylim,base,statah_mc_comb2);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_at_comb2_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_at_comb2_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a.png'],'-dpng')
        end
      catch
      end
      try
        fig=470+ss;
        topoplotTFR_highlight(fig,gravehi_APT_MSPN_fftaddbn,timwin,ylim,base,statah_mc_fftaddbn);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_at_fftaddbn_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25_bn.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_at_fftaddbn_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a_bn.png'],'-dpng')
        end
      catch
      end
      try
        fig=570+ss;
        topoplotTFR_highlight(fig,gravehi_APT_MSPN_comb1bn,timwin,ylim,base,statah_mc_comb1bn);
        if allcond_sameN
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_at_comb1bn_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a25_bn.png'],'-dpng')
        else
          print(fig,['D:\audtac\figs\grdiffhi_topoOverTime_mc_at_comb1bn_cond' num2str(ll) num2str(tt) num2str(ylim(1)) num2str(ylim(2)) 'a_bn.png'],'-dpng')
        end
      catch
      end
      
    end
    
    
    
    
    
    
    
    
  end
end






