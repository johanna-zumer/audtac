function [tlock_tac_s0,tlock_aud_s0,tlock_nul_s0,tr]=eeg_legomagic_trialSelection_freqwide(ii)
% function [tlock_tac_s0,tlock_aud_s0,tlock_nul_s0]=eeg_legomagic_trialSelection(ii)
% Trial selection: equal trial numbers per condition and removing NaN at ends

if ispc
  edir='D:\audtac\eeg_data\';
  ddir='D:\audtac\legomagic\diaries\';
else
  edir='/mnt/hgfs/D/audtac/eeg_data/';
  ddir='/mnt/hgfs/D/audtac/legomagic/diaries/';
end
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

cd([edir sub{ii} ])
%%

% load(['raw_each_rej_' sub{ii}],'raw_*'); % see eeg_legomagic_epoching.m for difference of these files
load(['raw1_each_rej_' sub{ii}],'raw_*'); % data has 0.2Hz highpass filter applied before epoching

cfg=[];
cfg.channel= {'all', '-ECG'};
cfg.reref='yes';
cfg.refchannel={'all', '-ECG'};
data_tac_ref=ft_preprocessing(cfg,raw_tac_rej);
data_aud_ref=ft_preprocessing(cfg,raw_aud_rej);
data_nul_ref=ft_preprocessing(cfg,raw_nul_rej);
clear raw_*
if ii<8
  soalist=[3 4 5 6 7];
else
  soalist=[1 3 4 5 6 7 9];
end

% MUst do shifting first before freqanalysis.
% include only trials during awake segments (in case participant went into N1 during sitting up portion)
for ll=1:length(soalist)
  for tt=3:6 % different thresholds of inclusion of SOA variance around desired SOA (see concat_stim_files.m)
    cfg=[];
    cfg.vartrllength=2;
    cfg.keeptrials='yes';
    cfg.trials=[data_tac_ref.trialinfo(:,tt)==soalist(ll) & data_tac_ref.trialinfo(:,2)==0 & ~isnan(data_tac_ref.trialinfo(:,10))];
    tlock_tac_s0{soalist(ll),tt-2}=ft_timelockanalysis(cfg,data_tac_ref); % only 'awake' state
    %   cfg.trials=data_tac_ref.trialinfo(:,3)==soalist(ll);
    %   tlock_tac_sall{soalist(ll)}=ft_timelockanalysis(cfg,data_tac_ref); % all trials, irrepective if fell asleep

    cfg=[];
    cfg.vartrllength=2;
    cfg.keeptrials='yes';
    cfg.trials=[data_aud_ref.trialinfo(:,tt)==soalist(ll) & data_aud_ref.trialinfo(:,2)==0 & ~isnan(data_aud_ref.trialinfo(:,10))];
    tlock_aud_s0{soalist(ll),tt-2}=ft_timelockanalysis(cfg,data_aud_ref);
    %   cfg.trials=data_aud_ref.trialinfo(:,3)==soalist(ll);
    %   tlock_aud_sall{soalist(ll)}=ft_timelockanalysis(cfg,data_aud_ref);
  end
end
cfg=[]; % tac alone 38 (-2)
cfg.vartrllength=2;
cfg.keeptrials='yes';
cfg.trials=[data_tac_ref.trialinfo(:,3)==-2 & data_tac_ref.trialinfo(:,2)==0];
tlock_tac_s0{10,1}=ft_timelockanalysis(cfg,data_tac_ref);
% cfg.trials=data_tac_ref.trialinfo(:,3)==-2;
% tlock_tac_sall{10}=ft_timelockanalysis(cfg,data_tac_ref);
cfg=[]; % aud alone 39 (-1)
cfg.vartrllength=2;
cfg.keeptrials='yes';
cfg.trials=[data_aud_ref.trialinfo(:,3)==-1 & data_aud_ref.trialinfo(:,2)==0];
tlock_aud_s0{10,1}=ft_timelockanalysis(cfg,data_aud_ref);
% cfg.trials=data_aud_ref.trialinfo(:,3)==-1;
% tlock_aud_sall{10}=ft_timelockanalysis(cfg,data_aud_ref);
cfg=[];
cfg.vartrllength=2;
cfg.keeptrials='yes';
cfg.trials=data_nul_ref.trialinfo(:,2)==0;
tlock_nul_s0{10,1}=ft_timelockanalysis(cfg,data_nul_ref);
% cfg.trials='all';
% tlock_nul_sall{10}=ft_timelockanalysis(cfg,data_nul_ref_ref);

clear data*

soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];
ttimwin=[min(-.5,soades-.5); max(0.8,soades+.8)]';
ttimwin(10,:)=[min(ttimwin(1:9,1)) max(ttimwin(1:9,2))];
atimwin=[min(-.5,fliplr(soades)-.5); max(0.8,fliplr(soades)+.8)]';
atimwin(10,:)=ttimwin(10,:);
ntimwin(10,:)=ttimwin(10,:);

chanuse=match_str(tlock_tac_s0{3}.label,'Fz');


for ll=[soalist 10] % this loop is to remove trials with NaN in (e.g. if artifact rejection cut through middle of a trial)
  for tt=1:4
    if ll==10 && tt>1
      continue
    else
      cfg=[];
%       cfg.keeptrials='yes';
      cfg.trials=find(~any(isnan(squeeze(tlock_tac_s0{ll,tt}.trial(:,chanuse,dsearchn(tlock_tac_s0{ll,tt}.time',ttimwin(ll,1)):dsearchn(tlock_tac_s0{ll,tt}.time',ttimwin(ll,2)) ))),2));
%       tlock_tac_s0{ll,tt}=ft_timelockanalysis(cfg,tlock_tac_s0{ll,tt});
      tlock_tac_s0{ll,tt}=ft_selectdata(cfg,tlock_tac_s0{ll,tt});
      
      tlock_tac_s0{ll,tt}=trim_nans(tlock_tac_s0{ll,tt});
            
      
      cfg=[];
%       cfg.keeptrials='yes';
%       cfg.trials=find(sum(squeeze(isnan(tlock_aud_s0{ll,tt}.trial(:,10,:)))')<.05*size(tlock_aud_s0{ll,tt}.trial,3));
      cfg.trials=find(~any(isnan(squeeze(tlock_aud_s0{ll,tt}.trial(:,chanuse,dsearchn(tlock_aud_s0{ll,tt}.time',atimwin(ll,1)):dsearchn(tlock_aud_s0{ll,tt}.time',atimwin(ll,2)) ))),2));
%       tlock_aud_s0{ll,tt}=ft_timelockanalysis(cfg,tlock_aud_s0{ll,tt});
      tlock_aud_s0{ll,tt}=ft_selectdata(cfg,tlock_aud_s0{ll,tt});
%       tmp=squeeze(mean(tlock_aud_s0{ll,tt}.trial,1));
%       if sum(isnan(tmp(1,:)))>.07*length(tmp(1,:))
%         error('this averaging didnt work');
%       end
      tlock_aud_s0{ll,tt}=trim_nans(tlock_aud_s0{ll,tt});
      
    end
    if ll==10 && tt==1
      cfg=[];
%       cfg.keeptrials='yes';
%       cfg.trials=find(sum(squeeze(isnan(tlock_nul_s0{ll}.trial(:,10,:)))')<.05*size(tlock_nul_s0{ll}.trial,3));
      cfg.trials=find(~any(isnan(squeeze(tlock_nul_s0{ll,tt}.trial(:,chanuse,dsearchn(tlock_nul_s0{ll,tt}.time',ntimwin(ll,1)):dsearchn(tlock_nul_s0{ll,tt}.time',ntimwin(ll,2)) ))),2));
%       tlock_nul_s0{ll}=ft_timelockanalysis(cfg,tlock_nul_s0{ll});
      tlock_nul_s0{ll}=ft_selectdata(cfg,tlock_nul_s0{ll});
%       tmp=squeeze(mean(tlock_nul_s0{ll,tt}.trial,1));
%       if sum(isnan(tmp(1,:)))>.07*length(tmp(1,:))
%         error('this averaging didnt work');
%       end
      tlock_nul_s0{ll,tt}=trim_nans(tlock_nul_s0{ll,tt});
      
    end
    if ll<10
      % make sure the multisensory auditory-locked and ms tactile-locked have same trials; sometimes not exactly!
      trcom=intersect(tlock_tac_s0{ll,tt}.trialinfo(:,11),tlock_aud_s0{ll,tt}.trialinfo(:,11));
      [~,tbb,~]=intersect(tlock_tac_s0{ll,tt}.trialinfo(:,11),trcom);
      [~,abb,~]=intersect(tlock_aud_s0{ll,tt}.trialinfo(:,11),trcom);      
      cfg=[];
      cfg.trials=tbb;
      tlock_tac_s0{ll,tt}=ft_selectdata(cfg,tlock_tac_s0{ll,tt});
      cfg=[];
      cfg.trials=abb;
      tlock_aud_s0{ll,tt}=ft_selectdata(cfg,tlock_aud_s0{ll,tt});
    end
  end
end

if 0  %why had we been loading this?
load([ddir sub{ii} '_audtac.mat']);
end

fsample=1/diff(tlock_tac_s0{10,1}.time(1:2)); % sampling rate


for ll=soalist
  for tt=1:4 % limits of plus/minus:  [no-limit .01 .005 .002]
    tnumtr(ll,tt)=min([size(tlock_tac_s0{ll,tt}.trial,1) size(tlock_tac_s0{10,1}.trial,1)]);
    anumtr(ll,tt)=min([size(tlock_aud_s0{ll,tt}.trial,1) size(tlock_aud_s0{10,1}.trial,1)]);
  end
end


% shifting unisensory to match for later sum of unisensory to match a multisensory condition
for ll=soalist
  for tt=1:4
    
    tmrandtr=Shuffle(1:size(tlock_tac_s0{ll,tt}.trial,1)); 
    turandtr=Shuffle(1:size(tlock_tac_s0{10}.trial,1)); 
    
    tmtrialkept{ll,tt}=sort(tmrandtr(1:tnumtr(ll,tt)));
    tutrialkept{ll,tt}=sort(turandtr(1:tnumtr(ll,tt)));

    cfg=[];
    % I have thought it through, and no minus sign needed here.
    cfg.offset=round(fsample*(tlock_tac_s0{ll,tt}.trialinfo(tmtrialkept{ll,tt},10)-soades(ll) -soades(ll))); % offset is in samples not seconds; -2*soades required for tactile not auditory
    cfg.trials=tutrialkept{ll,tt};
    tlock_tac_s0{ll+40,tt}=ft_redefinetrial(cfg,tlock_tac_s0{10}); % Talone with shift in time matching SOA jitter

    cfg=[];
    cfg.trials=tmtrialkept{ll,tt}; % only keep ms trials that have a paired shifted corresponding unisensory trial
    tlock_tac_s0{ll,tt}=ft_redefinetrial(cfg,tlock_tac_s0{ll,tt});
    
    amrandtr=Shuffle(1:size(tlock_aud_s0{ll,tt}.trial,1)); 
    aurandtr=Shuffle(1:size(tlock_aud_s0{10}.trial,1)); 
    
    amtrialkept{ll,tt}=sort(amrandtr(1:anumtr(ll,tt)));
    autrialkept{ll,tt}=sort(aurandtr(1:anumtr(ll,tt)));

    cfg=[];
    % I have thought it through, and no minus sign needed here.
    cfg.offset=round(fsample*(tlock_aud_s0{ll,tt}.trialinfo(amtrialkept{ll,tt},10)-soades(ll) +soades(ll))); % offset is in samples not seconds; -2*soades required for tactile not auditory
    cfg.trials=autrialkept{ll,tt};
    tlock_aud_s0{ll+40,tt}=ft_redefinetrial(cfg,tlock_aud_s0{10}); % Talone with shift in time matching SOA jitter

         
    cfg=[];
    cfg.trials=amtrialkept{ll,tt};
    tlock_aud_s0{ll,tt}=ft_redefinetrial(cfg,tlock_aud_s0{ll,tt});
    
  end
end

% timwin40=[timwin(1)+0.5 timwin(2)-0.5];

for ll=[soalist+40]
  for tt=1:4
% need to use atimwin to get right range
    cfg=[];
    cfg.trials=find(~any(isnan(squeeze(tlock_tac_s0{ll,tt}.trial(:,chanuse,dsearchn(tlock_tac_s0{ll,tt}.time',atimwin(ll-40,1)):dsearchn(tlock_tac_s0{ll,tt}.time',atimwin(ll-40,2)) ))),2));
    tlock_tac_s0{ll,tt}=ft_selectdata(cfg,tlock_tac_s0{ll,tt});
%     tlock_tac_s0{ll,tt}=trim_nans(tlock_tac_s0{ll,tt});

% need to use ttimwin to get right range
    cfg=[];
    cfg.trials=find(~any(isnan(squeeze(tlock_aud_s0{ll,tt}.trial(:,chanuse,dsearchn(tlock_aud_s0{ll,tt}.time',ttimwin(ll-40,1)):dsearchn(tlock_aud_s0{ll,tt}.time',ttimwin(ll-40,2)) ))),2));
    tlock_aud_s0{ll,tt}=ft_selectdata(cfg,tlock_aud_s0{ll,tt});
%     tlock_aud_s0{ll,tt}=trim_nans(tlock_aud_s0{ll,tt});
    
  end
end

for ll=soalist
  for tt=1:4
    numcondtfinal(ll,tt)=min([size(tlock_tac_s0{10}.trial,1) size(tlock_aud_s0{ll+40,tt}.trial,1) size(tlock_tac_s0{ll,tt}.trial,1) size(tlock_nul_s0{10}.trial,1)]);
    numcondafinal(ll,tt)=min([size(tlock_aud_s0{10}.trial,1) size(tlock_tac_s0{ll+40,tt}.trial,1) size(tlock_aud_s0{ll,tt}.trial,1) size(tlock_nul_s0{10}.trial,1)]);

    t10randtr=Shuffle(1:size(tlock_tac_s0{10}.trial,1));
    all40randtr=Shuffle(1:size(tlock_aud_s0{ll+40,tt}.trial,1));
    tllrandtr=Shuffle(1:size(tlock_tac_s0{ll,tt}.trial,1));
    nllrandtr=Shuffle(1:size(tlock_nul_s0{10}.trial,1));
    
    tr.t10trialkept{ll,tt}=sort(t10randtr(1:numcondtfinal(ll,tt)));
    tr.all40trialkept{ll,tt}=sort(all40randtr(1:numcondtfinal(ll,tt)));
    tr.tlltrialkept{ll,tt}=sort(tllrandtr(1:numcondtfinal(ll,tt)));
    tr.nllttrialkept{ll,tt}=sort(nllrandtr(1:numcondtfinal(ll,tt)));
    
    a10randtr=Shuffle(1:size(tlock_aud_s0{10}.trial,1));
    tll40randtr=Shuffle(1:size(tlock_tac_s0{ll+40,tt}.trial,1));
    allrandtr=Shuffle(1:size(tlock_aud_s0{ll,tt}.trial,1));
    nllrandtr=Shuffle(1:size(tlock_nul_s0{10}.trial,1));
    
    tr.a10trialkept{ll,tt}=sort(a10randtr(1:numcondafinal(ll,tt)));
    tr.tll40trialkept{ll,tt}=sort(tll40randtr(1:numcondafinal(ll,tt)));
    tr.alltrialkept{ll,tt}=sort(allrandtr(1:numcondafinal(ll,tt)));
    tr.nllatrialkept{ll,tt}=sort(nllrandtr(1:numcondafinal(ll,tt)));
    
    cfg=[];
    cfg.trials=tr.all40trialkept{ll,tt};
    tlock_aud_s0{ll+40,tt}=ft_selectdata(cfg,tlock_aud_s0{ll+40,tt}); 
    cfg=[];
    cfg.trials=tr.tlltrialkept{ll,tt};
    tlock_tac_s0{ll,tt}=ft_selectdata(cfg,tlock_tac_s0{ll,tt}); 
    if 0 % too memory intensive to do here
      cfg=[];
      cfg.trials=tr.t10trialkept{ll,tt};
      tlock_tac_s0{ll+20,tt}=ft_selectdata(cfg,tlock_tac_s0{10}); % create +20 here as 10 will be different for every ll,tt
      cfg=[];
      cfg.trials=tr.nllttrialkept{ll,tt};
      tlock_nul_s0{ll+50,tt}=ft_selectdata(cfg,tlock_nul_s0{10});
    end
    
    cfg=[];
    cfg.trials=tr.tll40trialkept{ll,tt};
    tlock_tac_s0{ll+40,tt}=ft_selectdata(cfg,tlock_tac_s0{ll+40,tt}); 
    cfg=[];
    cfg.trials=tr.alltrialkept{ll,tt};
    tlock_aud_s0{ll,tt}=ft_selectdata(cfg,tlock_aud_s0{ll,tt}); 
    if 0 % too memory intensive to do here
      cfg=[];
      cfg.trials=tr.a10trialkept{ll,tt};
      tlock_aud_s0{ll+20,tt}=ft_selectdata(cfg,tlock_aud_s0{10}); % create +20 here as 10 will be different for every ll,tt
      cfg=[];
      cfg.trials=tr.nllatrialkept{ll,tt};
      tlock_nul_s0{ll+60,tt}=ft_selectdata(cfg,tlock_nul_s0{10});
    end
  end
end

save('trialkept.mat','tr','numcond*')





