function [tlock_tac_s0,tlock_aud_s0,tlock_nul_s0]=eeg_legomagic_trialSelection(ii)
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
  end
end
for ll=1:length(soalist)
  for tt=3:6 % different thresholds of inclusion of SOA variance around desired SOA (see concat_stim_files.m)
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

timwin=[-0.7-0.5 0.8+0.5]; % ideal window is -0.7 to 0.5 but this will be shifted by up to 500ms either way 

for ll=[soalist 10] % this loop is to remove trials with NaN in (e.g. if artifact rejection cut through middle of a trial)
  for tt=1:4
    if ll==10 && tt>1
      continue
    else
      cfg=[];
      cfg.keeptrials='yes';
%       cfg.trials=find(sum(squeeze(isnan(tlock_tac_s0{ll,tt}.trial(:,10,:)))')<.05*size(tlock_tac_s0{ll,tt}.trial,3));
      cfg.trials=find(~any(isnan(squeeze(tlock_tac_s0{ll,tt}.trial(:,17,dsearchn(tlock_tac_s0{ll,tt}.time',timwin(1)):dsearchn(tlock_tac_s0{ll,tt}.time',timwin(2)) ))),2));
      tlock_tac_s0{ll,tt}=ft_timelockanalysis(cfg,tlock_tac_s0{ll,tt});
%       tmp=squeeze(mean(tlock_tac_s0{ll,tt}.trial,1));
%       if sum(isnan(tmp(1,:)))>.07*length(tmp(1,:))
%         error('this averaging didnt work');
%       end
      begsample=1;
      endsample=size(tlock_tac_s0{ll,tt}.trial,3);
      zerosample=dsearchn(tlock_tac_s0{ll,tt}.time',0);
      for rr=1:size(tlock_tac_s0{ll,tt}.trial,1)
        tmp=find(isnan(squeeze(tlock_tac_s0{ll,tt}.trial(rr,1,:))));
        if ~isempty(tmp)
          tmpr=max(tmp(tmp<zerosample));
          tmpo=min(tmp(tmp>zerosample));
          if begsample<tmpr
            begsample=tmpr;
          end
          if endsample>tmpo
            endsample=tmpo;
          end
        end
      end
      cfg=[];
      cfg.latency=[tlock_tac_s0{ll,tt}.time(begsample+1) tlock_tac_s0{ll,tt}.time(endsample-1)];
      tlock_tac_s0{ll,tt}=ft_selectdata(cfg,tlock_tac_s0{ll,tt});
      
      
      cfg=[];
      cfg.keeptrials='yes';
%       cfg.trials=find(sum(squeeze(isnan(tlock_aud_s0{ll,tt}.trial(:,10,:)))')<.05*size(tlock_aud_s0{ll,tt}.trial,3));
      cfg.trials=find(~any(isnan(squeeze(tlock_aud_s0{ll,tt}.trial(:,17,dsearchn(tlock_aud_s0{ll,tt}.time',timwin(1)):dsearchn(tlock_aud_s0{ll,tt}.time',timwin(2)) ))),2));
      tlock_aud_s0{ll,tt}=ft_timelockanalysis(cfg,tlock_aud_s0{ll,tt});
%       tmp=squeeze(mean(tlock_aud_s0{ll,tt}.trial,1));
%       if sum(isnan(tmp(1,:)))>.07*length(tmp(1,:))
%         error('this averaging didnt work');
%       end
      begsample=1;
      endsample=size(tlock_aud_s0{ll,tt}.trial,3);
      zerosample=dsearchn(tlock_aud_s0{ll,tt}.time',0);
      for rr=1:size(tlock_aud_s0{ll,tt}.trial,1)
        tmp=find(isnan(squeeze(tlock_aud_s0{ll,tt}.trial(rr,1,:))));
        if ~isempty(tmp)
          tmpr=max(tmp(tmp<zerosample));
          tmpo=min(tmp(tmp>zerosample));
          if begsample<tmpr
            begsample=tmpr;
          end
          if endsample>tmpo
            endsample=tmpo;
          end
        end
      end
      cfg=[];
      cfg.latency=[tlock_aud_s0{ll,tt}.time(begsample+1) tlock_aud_s0{ll,tt}.time(endsample-1)];
      tlock_aud_s0{ll,tt}=ft_selectdata(cfg,tlock_aud_s0{ll,tt});
      
    end
    if ll==10 && tt==1
      cfg=[];
      cfg.keeptrials='yes';
%       cfg.trials=find(sum(squeeze(isnan(tlock_nul_s0{ll}.trial(:,10,:)))')<.05*size(tlock_nul_s0{ll}.trial,3));
      cfg.trials=find(~any(isnan(squeeze(tlock_nul_s0{ll,tt}.trial(:,17,dsearchn(tlock_nul_s0{ll,tt}.time',timwin(1)):dsearchn(tlock_nul_s0{ll,tt}.time',timwin(2)) ))),2));
      tlock_nul_s0{ll}=ft_timelockanalysis(cfg,tlock_nul_s0{ll});
%       tmp=squeeze(mean(tlock_nul_s0{ll,tt}.trial,1));
%       if sum(isnan(tmp(1,:)))>.07*length(tmp(1,:))
%         error('this averaging didnt work');
%       end
      begsample=1;
      endsample=size(tlock_nul_s0{ll,tt}.trial,3);
      zerosample=dsearchn(tlock_nul_s0{ll,tt}.time',0);
      for rr=1:size(tlock_nul_s0{ll,tt}.trial,1)
        tmp=find(isnan(squeeze(tlock_nul_s0{ll,tt}.trial(rr,1,:))));
        if ~isempty(tmp)
          tmpr=max(tmp(tmp<zerosample));
          tmpo=min(tmp(tmp>zerosample));
          if begsample<tmpr
            begsample=tmpr;
          end
          if endsample>tmpo
            endsample=tmpo;
          end
        end
      end
      cfg=[];
      cfg.latency=[tlock_nul_s0{ll,tt}.time(begsample) tlock_nul_s0{ll,tt}.time(endsample)];
      tlock_nul_s0{ll,tt}=ft_selectdata(cfg,tlock_nul_s0{ll,tt});
      
    end
  end
end

plotflag=1;
load([ddir sub{ii} '_audtac.mat']);
soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];
chanuse=match_str(tlock_tac_s0{3}.label,'Fz');
fsample=1/diff(tlock_tac_s0{10,1}.time(1:2)); % sampling rate

% determine trial numbers per condition
numuni=min([size(tlock_tac_s0{10,1}.trialinfo,1),size(tlock_aud_s0{10,1}.trialinfo,1),size(tlock_nul_s0{10,1}.trialinfo,1)]);
for ll=soalist
  for tt=1:4 % limits of plus/minus:  [no-limit .01 .005 .002]
    trcom=intersect(tlock_tac_s0{ll,tt}.trialinfo(:,11),tlock_aud_s0{ll,tt}.trialinfo(:,11));
    nummul(ll,tt)=length(trcom);
    numcond(ll,tt)=min(nummul(ll,tt),numuni);
  end
end
% numcondall=min(numcond(soalist,:)); % min over all conditions for a given SOA variance threshold

% load(['tlock_diffs_' sub{ii} '.mat'],'*trialkept')

if 0 % with redo of raw1_* we want new selection of trials
try
  load(['trialkept.mat'])
  tk=1;
catch
  tk=0;
end
else
  tk=0;
end

% shifting unisensory to match for later sum of unisensory to match a multisensory condition
for ll=soalist
  for tt=1:4
    % make sure the multisensory auditory-locked and ms tactile-locked have same trials; sometimes not exactly!
    trcom=intersect(tlock_tac_s0{ll,tt}.trialinfo(:,11),tlock_aud_s0{ll,tt}.trialinfo(:,11));
    [~,tbb,~]=intersect(tlock_tac_s0{ll,tt}.trialinfo(:,11),trcom);
    [~,abb,~]=intersect(tlock_aud_s0{ll,tt}.trialinfo(:,11),trcom);
    
    if tk==0 % if not exist already from previous call to this function
      mrandtr=Shuffle(1:length(trcom)); % random trials to keep
      %     mrandtr=Shuffle(1:size(tlock_tac_s0{ll,tt}.trialinfo,1)); % random trials to keep
      trandtr=Shuffle(1:size(tlock_tac_s0{10}.trialinfo,1)); % random trials to keep
      arandtr=Shuffle(1:size(tlock_aud_s0{10}.trialinfo,1)); % random trials to keep
      nrandtr=Shuffle(1:size(tlock_nul_s0{10}.trialinfo,1)); % random trials to keep
      
      mtrialkept{ll,tt}=sort(mrandtr(1:numcond(ll,tt))); % only keep the appropriate number
      ttrialkept{ll,tt}=sort(trandtr(1:numcond(ll,tt))); % only keep the appropriate number
      atrialkept{ll,tt}=sort(arandtr(1:numcond(ll,tt))); % only keep the appropriate number
      ntrialkept{ll,tt}=sort(nrandtr(1:numcond(ll,tt))); % only keep the appropriate number
      
      % must also reduce trials for multisensory condition if they have more
      tarandtr=Shuffle(1:size(tlock_tac_s0{ll,tt}.trialinfo,1)); % random trials to keep
      atrandtr=Shuffle(1:size(tlock_aud_s0{ll,tt}.trialinfo,1)); % random trials to keep
      tatrialkept{ll,tt}=sort(tarandtr(1:numcond(ll,tt))); % only keep the appropriate number
      attrialkept{ll,tt}=sort(atrandtr(1:numcond(ll,tt))); % only keep the appropriate number
    end
    
    cfg=[];
    % I have thought it through, and no minus sign needed here.
    cfg.offset=round(fsample*(tlock_tac_s0{ll,tt}.trialinfo(tbb(mtrialkept{ll,tt}),10)-soades(ll) -soades(ll))); % offset is in samples not seconds; -2*soades required for tactile not auditory
    cfg.trials=ttrialkept{ll,tt};
    tlock_tac_s0{ll+40,tt}=ft_redefinetrial(cfg,tlock_tac_s0{10}); % Talone with shift in time matching SOA jitter
    
    cfg=[];
    % I have thought it through, and no minus sign needed here.
    cfg.offset=round(fsample*(tlock_aud_s0{ll,tt}.trialinfo(abb(mtrialkept{ll,tt}),10)-soades(ll) +soades(ll))); % tlock_aud_s0{ll,tt}.trialinfo(trialkept,10) is identical to tlock_tac_s0{ll,tt}.trialinfo(trialkept,10)
    cfg.trials=atrialkept{ll,tt};
    tlock_aud_s0{ll+40,tt}=ft_redefinetrial(cfg,tlock_aud_s0{10}); % Aalone with shift in time matching SOA jitter
    
    cfg=[];
    cfg.trials=ntrialkept{ll,tt};
    tlock_nul_s0{ll+40,tt}=ft_redefinetrial(cfg,tlock_nul_s0{10});
    
    cfg=[];
    cfg.trials=tatrialkept{ll,tt};
    tlock_tac_s0{ll,tt}=ft_redefinetrial(cfg,tlock_tac_s0{ll,tt});
    cfg=[];
    cfg.trials=attrialkept{ll,tt};
    tlock_aud_s0{ll,tt}=ft_redefinetrial(cfg,tlock_aud_s0{ll,tt});
    
    
  end
end

timwin40=[timwin(1)+0.5 timwin(2)-0.5];

for ll=[soalist+40]
  for tt=1:4
    %       cfg=[];
    %       cfg.keeptrials='yes';
    %       cfg.trials=find(sum(squeeze(isnan(tlock_tac_s0{ll,tt}.trial(:,10,:)))')<.05*size(tlock_tac_s0{ll,tt}.trial,3));
    %       tlock_tac_s0{ll,tt}=ft_timelockanalysis(cfg,tlock_tac_s0{ll,tt});
    %       tmp=squeeze(mean(tlock_tac_s0{ll,tt}.trial,1));
    %       if sum(isnan(tmp(1,:)))>.05*length(tmp(1,:))
    %         error('this averaging didnt work');
    %       end
%     begsample=1;
%     endsample=size(tlock_tac_s0{ll,tt}.trial,3);
%     zerosample=dsearchn(tlock_tac_s0{ll,tt}.time',0);
%     for rr=1:size(tlock_tac_s0{ll,tt}.trial,1)
%       tmp=find(isnan(squeeze(tlock_tac_s0{ll,tt}.trial(rr,1,:))));
%       if ~isempty(tmp)
%         tmpr=max(tmp(tmp<zerosample));
%         tmpo=min(tmp(tmp>zerosample));
%         if begsample<tmpr
%           begsample=tmpr;
%         end
%         if endsample>tmpo
%           endsample=tmpo;
%         end
%       end
%     end
%     cfg=[];
%     cfg.latency=[tlock_tac_s0{ll,tt}.time(begsample+1) tlock_tac_s0{ll,tt}.time(endsample-1)];
%     if diff(cfg.latency)<1.2
%       error('not long enough trials')
%     end
%     tlock_tac_s0{ll,tt}=ft_selectdata(cfg,tlock_tac_s0{ll,tt});
    
    cfg=[];
      cfg.trials=find(~any(isnan(squeeze(tlock_tac_s0{ll,tt}.trial(:,17,dsearchn(tlock_tac_s0{ll,tt}.time',timwin40(1)):dsearchn(tlock_tac_s0{ll,tt}.time',timwin40(2)) ))),2));
    tlock_tac_s0{ll,tt}=ft_selectdata(cfg,tlock_tac_s0{ll,tt});
    
    %       cfg=[];
    %       cfg.keeptrials='yes';
    %       cfg.trials=find(sum(squeeze(isnan(tlock_aud_s0{ll,tt}.trial(:,10,:)))')<.05*size(tlock_aud_s0{ll,tt}.trial,3));
    %       tlock_aud_s0{ll,tt}=ft_timelockanalysis(cfg,tlock_aud_s0{ll,tt});
    %       tmp=squeeze(mean(tlock_aud_s0{ll,tt}.trial,1));
    %       if sum(isnan(tmp(1,:)))>.05*length(tmp(1,:))
    %         error('this averaging didnt work');
    %       end
%     begsample=1;
%     endsample=size(tlock_aud_s0{ll,tt}.trial,3);
%     zerosample=dsearchn(tlock_aud_s0{ll,tt}.time',0);
%     for rr=1:size(tlock_aud_s0{ll,tt}.trial,1)
%       tmp=find(isnan(squeeze(tlock_aud_s0{ll,tt}.trial(rr,1,:))));
%       if ~isempty(tmp)
%         tmpr=max(tmp(tmp<zerosample));
%         tmpo=min(tmp(tmp>zerosample));
%         if begsample<tmpr
%           begsample=tmpr;
%         end
%         if endsample>tmpo
%           endsample=tmpo;
%         end
%       end
%     end
%     cfg=[];
%     cfg.latency=[tlock_aud_s0{ll,tt}.time(begsample+1) tlock_aud_s0{ll,tt}.time(endsample-1)];
%     if diff(cfg.latency)<1.2
%       error('not long enough trials')
%     end
    cfg=[];
      cfg.trials=find(~any(isnan(squeeze(tlock_aud_s0{ll,tt}.trial(:,17,dsearchn(tlock_aud_s0{ll,tt}.time',timwin40(1)):dsearchn(tlock_aud_s0{ll,tt}.time',timwin40(2)) ))),2));
    tlock_aud_s0{ll,tt}=ft_selectdata(cfg,tlock_aud_s0{ll,tt});
    

      cfg=[];
      cfg.trials=find(~any(isnan(squeeze(tlock_nul_s0{ll,tt}.trial(:,17,dsearchn(tlock_nul_s0{ll,tt}.time',timwin40(1)):dsearchn(tlock_nul_s0{ll,tt}.time',timwin40(2)) ))),2));
    tlock_nul_s0{ll,tt}=ft_selectdata(cfg,tlock_nul_s0{ll,tt});

  end
end

for ll=soalist
  for tt=1:4
    numcondtfinal(ll,tt)=min([size(tlock_tac_s0{10}.trial,1) size(tlock_aud_s0{ll+40,tt}.trial,1) size(tlock_tac_s0{ll,tt}.trial,1) size(tlock_nul_s0{ll+40,tt}.trial,1)]);
    numcondafinal(ll,tt)=min([size(tlock_aud_s0{10}.trial,1) size(tlock_tac_s0{ll+40,tt}.trial,1) size(tlock_aud_s0{ll,tt}.trial,1) size(tlock_nul_s0{ll+40,tt}.trial,1)]);
    keyboard

    t10randtr=Shuffle(1:size(tlock_tac_s0{10}.trial,1));
    all40randtr=Shuffle(1:size(tlock_aud_s0{ll+40,tt}.trial,1));
    tllrandtr=Shuffle(1:size(tlock_tac_s0{ll,tt}.trial,1));
    nll40randtr=Shuffle(1:size(tlock_nul_s0{ll+40,tt}.trial,1));
    
    t10trialkept{ll,tt}=sort(t10randtr(1:numcondtfinal(ll,tt)));
    all40trialkept{ll,tt}=sort(all40randtr(1:numcondtfinal(ll,tt)));
    tlltrialkept{ll,tt}=sort(tllrandtr(1:numcondtfinal(ll,tt)));
    nll40ttrialkept{ll,tt}=sort(nll40randtr(1:numcondtfinal(ll,tt)));
    
    a10randtr=Shuffle(1:size(tlock_aud_s0{10}.trial,1));
    tll40randtr=Shuffle(1:size(tlock_tac_s0{ll+40,tt}.trial,1));
    allrandtr=Shuffle(1:size(tlock_aud_s0{ll,tt}.trial,1));
    nll40randtr=Shuffle(1:size(tlock_nul_s0{ll+40,tt}.trial,1));
    
    a10trialkept{ll,tt}=sort(a10randtr(1:numcondafinal(ll,tt)));
    tll40trialkept{ll,tt}=sort(tll40randtr(1:numcondafinal(ll,tt)));
    alltrialkept{ll,tt}=sort(allrandtr(1:numcondafinal(ll,tt)));
    nll40atrialkept{ll,tt}=sort(nll40randtr(1:numcondafinal(ll,tt)));
    
    cfg=[];
    cfg.trials=t10trialkept{ll,tt};
    tlock_tac_s0{ll+20,tt}=ft_selectdata(cfg,tlock_tac_s0{10}); % create +20 here as it will be different for every ll,tt
    cfg=[];
    cfg.trials=all40trialkept{ll,tt};
    tlock_aud_s0{ll+40,tt}=ft_selectdata(cfg,tlock_aud_s0{ll+40,tt}); 
    cfg=[];
    cfg.trials=tlltrialkept{ll,tt};
    tlock_tac_s0{ll,tt}=ft_selectdata(cfg,tlock_tac_s0{ll,tt}); 
    cfg=[];
    cfg.trials=nll40ttrialkept{ll,tt};
    tlock_nul_s0{ll+50,tt}=ft_selectdata(cfg,tlock_nul_s0{ll+40,tt}); 

    cfg=[];
    cfg.trials=a10trialkept{ll,tt};
    tlock_aud_s0{ll+20,tt}=ft_selectdata(cfg,tlock_aud_s0{10}); % create +20 here as it will be different for every ll,tt
    cfg=[];
    cfg.trials=tll40trialkept{ll,tt};
    tlock_tac_s0{ll+40,tt}=ft_selectdata(cfg,tlock_tac_s0{ll+40,tt}); 
    cfg=[];
    cfg.trials=alltrialkept{ll,tt};
    tlock_aud_s0{ll,tt}=ft_selectdata(cfg,tlock_aud_s0{ll,tt}); 
    cfg=[];
    cfg.trials=nll40atrialkept{ll,tt};
    tlock_nul_s0{ll+60,tt}=ft_selectdata(cfg,tlock_nul_s0{ll+40,tt}); 
  end
end

save('trialkept.mat','*trialkept')




