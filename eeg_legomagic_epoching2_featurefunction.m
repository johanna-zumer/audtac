function [raw_tac, raw_aud, raw_nul, artflag]=eeg_legomagic_epoching2_featurefunction(ii,sleep)

warning('this function is now deprecated; see eeg_legomagic_epoching2');
% return

% preprocessing of EEG data from Hills, 64ch MR-cap
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

cd([edir sub{ii} ])

%%


% see eeg_legomagic_viewFeatures.m (these only exist now for iiuse)
if sleep
  try
    load('artfctdef_bed_automan.mat')
    artfctdef=artfctdef_new;clear artfctdef_new
    artflag=1;
  catch
    load('artfctdef_bed_auto.mat')
    artflag=0;
  end
  files=dir('cc*b.mat');
else
  try
    load('artfctdef_sit_automan.mat')
    artfctdef=artfctdef_new;clear artfctdef_new
    artflag=1;
  catch
    load('artfctdef_sit_auto.mat')
    artflag=0;
  end
  files=dir('cc*s.mat');
end


%   if 0
%     load([cfg.dataset(1:end-4) '_event.mat']); % loads event saved out during ft_trialfun_general_sleep.m
%     eventS10=ft_filter_event(event,'value','S 10');
%     eventS2=ft_filter_event(event,'value','S  2');
%     eventS1=ft_filter_event(event,'value','S  1');
%     [eventS10.type]=deal('A');
%     [eventS2.type]=deal('T');
%     [eventS1.type]=deal('N');
%     %   [eventS2.offset]=deal(0);
%     for ll=1:length(eventS2),eventS2(ll).sample=eventS2(ll).sample+200;end
%     eventStim=[eventS2 eventS10 eventS1];
%
%     cfg=[];
%     cfg.layout='EEG1010.lay';
%     cfg.viewmode='vertical';
%     cfg.blocksize=30;
%     cfg.plotlabels='yes';
%     cfg.channel={'Fz' 'Cz' 'Pz' 'Oz'};
%     cfg.event=eventStim;
%     cfg.ploteventlabels='colorvalue';
%     dcfg=ft_databrowser(cfg,raw_all_ica_rej);
%   end



load([ddir sub{ii} '_audtac.mat']);
saudtrials=find(saudseq);
baudtrials=find(baudseq);
% stactrials=find(stime_touch>0);
% btactrials=find(btime_touch>0);
if sleep
  numtac=bnumtac;
  numaud=bnumaud;
  numnul=bnumnul;
  audtrials=baudtrials;
  nulseq=bnulseq;
  time_touch=btime_touch;
  try
    load time_sleep.mat
  catch
    cfg=[];
    cfg.dataset=files.name;
    raw_all=ft_preprocessing(cfg);
    time=raw_all.time;
    save time_sleep.mat time
    clear raw_all
  end
else
  numtac=snumtac;
  numaud=snumaud;
  numnul=snumnul;
  audtrials=saudtrials;
  nulseq=snulseq;
  time_touch=stime_touch;
  load time.mat
%   switch ii
%     case 8
%       time_touch=time_touch(1:end-1);
%     otherwise
%   end
end

% only care to have at least one channel in data to carry through the .trialinfo
cfg=[];
cfg.dataset=files.name;
cfg.demean='yes';
%   cfg.bsfilter='yes';
%   cfg.bsfreq=[49 51; 99 101; 149 151];
%   cfg.hpfilter='yes';
%   cfg.hpfiltord=3;
%   cfg.hpfreq=0.2;  % Could do higher for awake data but not want higher for sleep.
%   cfg.channel={'all' '-ECG' '-VEOG' '-HEOG' '-EMG'};
cfg.channel={'Fz' 'Cz' 'Oz'};
%   cfg.reref='yes';
%   cfg.refchannel='all'; % used linked-mastoids for sleep staging but use average reference for ERPs and TFRs and source localisation
raw_hpf=ft_preprocessing(cfg);



prestim=5;  % seconds to evaluate overall
poststim=5;
prenar=1; % more narrow focus
postnar=1;



% All trials
cfg=[];
%     if sleep
%       cfg.dataset=fileb.name;
%     else
cfg.dataset=files.name;
%     end
cfg.trialfun='ft_trialfun_general_sleep';
cfg.trialdef.eventtype  = 'Stimulus';
cfg.trialdef.eventvalue = 'S 20'; % corresponds to lj.prepareStrobe(20)
cfg.trialdef.prestim = prenar; % these shouldn't be same prestim/poststim as for stimulus-locked definetrial
cfg.trialdef.poststim = postnar; % these shouldn't be same prestim/poststim as for stimulus-locked definetrial
cfg.time=time;
cfgtr=ft_definetrial(cfg);


if sleep
  bsoarealextra=[bsoareal(1,[1 end])];
  soause=setdiff(1:size(bsoareal,2),cfgtr.trlselexcl); % sometimes first or last trial may be excluded if too close to dataset boundary as function of prestim/poststim length

  %   if bsoarealextra(1) % this is weird. if it's 0 (Nul), don't get rid of it.
%   bsoareal=bsoareal(:,2:end);
%   bsoaeff=bsoaeff(2:end);
%   %   end
%   %   if bsoarealextra(end)
%   bsoareal=bsoareal(:,1:end-1);
%   bsoaeff=bsoaeff(1:end-1);
  %   end
  cfgtr.trl(:,7)=bsoareal(1,soause); % <---- which index 1,2,3,4 to use????
  cfgtr.trl(:,8)=bsoareal(2,soause); % <---- which index 1,2,3,4 to use????
  cfgtr.trl(:,9)=bsoareal(3,soause); % <---- which index 1,2,3,4 to use????
  cfgtr.trl(:,10)=bsoareal(4,soause); % <---- which index 1,2,3,4 to use????
  cfgtr.trl(:,13)=bsoaeff(soause);
else
  ssoarealextra=[ssoareal(1,[1 end])];
  soause=setdiff(1:size(ssoareal,2),cfgtr.trlselexcl);
  %   if ssoarealextra(1)
%   switch ii
%     case [8 ]
%     otherwise
%       ssoareal=ssoareal(:,2:end);
%       ssoaeff=ssoaeff(2:end);
%   end
%   %   end
%   %   if ssoarealextra(end)
%   switch ii
%     case [5 8]
%     otherwise
%       ssoareal=ssoareal(:,1:end-1);
%       ssoaeff=ssoaeff(1:end-1);
%   end
%   switch ii
%     case 8
%       cfgtr.trl=cfgtr.trl(1:end-1,:);
%     otherwise
%   end
  %   end
  cfgtr.trl(:,7)=ssoareal(1,soause); % <---- which index 1,2,3,4 to use????
  cfgtr.trl(:,8)=ssoareal(2,soause); % <---- which index 1,2,3,4 to use????
  cfgtr.trl(:,9)=ssoareal(3,soause); % <---- which index 1,2,3,4 to use????
  cfgtr.trl(:,10)=ssoareal(4,soause); % <---- which index 1,2,3,4 to use????
  cfgtr.trl(:,13)=ssoaeff(soause);
end
cfgtr.trl(:,11)=cfgtr.trl(:,6)>0 | cfgtr.trl(:,6)==-2; % tactile trials
cfgtr.trl(:,12)=cfgtr.trl(:,6)>0 | cfgtr.trl(:,6)==-1; % auditory trials
cfgtr.trl(:,14)=1:size(cfgtr.trl,1);

% discovered a code '214' for ii=[8b 10b 12s]
if find(cfgtr.trl(:,6)<-2 | cfgtr.trl(:,6)>9)
  disp('this is weird')
  cfgtr.trl(find(cfgtr.trl(:,6)<-2 | cfgtr.trl(:,6)>9),1:6)
  keyboard
  cfgtr.trl(find(cfgtr.trl(:,6)<-2 | cfgtr.trl(:,6)>9),:)=[];
end

ns=nan(1,length(cfgtr.event));
for ll=1:length(cfgtr.event),ns(ll)=strcmp(cfgtr.event(ll).type,'New Segment');end
nsind=find(ns);

if sleep==0 && ii==20 % Must have started/stopped within same acquistion
  nsind=nsind([2 3 5 6]);
elseif sleep==1 && ii==8
  nsind=nsind([1 3:6]);
elseif sleep==0 && ii==21
  nsind=nsind([2 3 4 5]);
elseif sleep==0 && ii==32
  nsind=nsind([1:6 8]);
end


% Tactile
% need to get timing for each trial centered on stimulus correctly
cfg=[];
cfg.dataset=files.name;
cfg.trialfun='ft_trialfun_general_sleep';
cfg.trialdef.eventtype  = 'Stimulus';
cfg.trialdef.eventvalue = 'S  2'; % corresponds to lj.prepareStrobe(2)
cfg.trialdef.prestim = prestim;
cfg.trialdef.poststim = poststim;
cfg.time=time;  %% <- why is this here?
cfgtactr=ft_definetrial(cfg);

for ll=1:length(cfgtr.event),tac(ll)=strcmp(cfgtr.event(ll).value,'S  2');end
tacind=find(tac);
for ll=1:length(nsind)-1
  cfgtacdiff(ll)=length(find(tac(nsind(ll):nsind(ll+1))))-numtac(ll);
end
cfgtacdiff(length(nsind))=length(find(tac(nsind(ll+1):end)))-numtac(ll+1);
ss=0;
trluse=[];
for ll=1:length(nsind)
  %     trluse=[trluse [ss+1]:[ss+numtac(ll)-cfgtacdiff(ll)]];
  trluse=[trluse [ss+1]:[ss+numtac(ll)]];
  ss=ss+numtac(ll)+cfgtacdiff(ll);
end
if length(trluse)~=length(time_touch) % should be the same
  error('num tac wrong');  %% FIX systematicaly not case by case!!!
end

tactr=find(cfgtr.trl(:,11));
if length(trluse)~=size(cfgtactr.trl,1) % due to prenar prestim there might be differences
  if length(trluse)-size(cfgtactr.trl,1) ==2 % safe enought to assume chop off a trial at either end
    trluse=trluse(2:end-1);
  elseif length(trluse)-size(cfgtactr.trl,1) ==1 % figure out if its beginning or end to chop off
    if any(cfgtactr.trl(1:5,6)-cfgtr.trl(tactr(trluse(1:5)),6)) && ~any(cfgtactr.trl(end-4:end,6)-cfgtr.trl(tactr(trluse(end-4:end)),6))
      trluse=trluse(2:end); % remove first trial
    elseif any(cfgtactr.trl(end-4:end,6)-cfgtr.trl(tactr(trluse(end-4:end)),6)) && ~any(cfgtactr.trl(1:5,6)-cfgtr.trl(tactr(trluse(1:5)),6))
      trluse=trluse(1:end-1);
    else
      disp('not sure what to throw out')
      keyboard
    end
  else
    disp('why number of trial mismatch?')
    keyboard
  end
end
% if sleep
%   if bsoarealextra(1)==-2 || bsoarealextra(1)>0
%     trluse=trluse(2:end)-1;
%   end
%   if bsoarealextra(end)==-2 || bsoarealextra(end)>0
%     trluse=trluse(1:end-1);
%   end
% else
%   switch ii
%     case [8 ]
%     otherwise
%       if ssoarealextra(1)==-2 || ssoarealextra(1)>0
%         trluse=trluse(2:end)-1;
%       end
%   end
%   
%   switch ii
%     case []
%     otherwise
%       if ssoarealextra(end)==-2 || ssoarealextra(end)>0
%         trluse=trluse(1:end-1);
%       end
%   end
%   switch ii
%     case 8
%       cfgtactr.trl=cfgtactr.trl(1:end-1,:);
%     otherwise
%   end
% end
tactrl=[cfgtactr.trl cfgtr.trl(tactr(trluse),7:end)];
% tactrl=cfgtr.trl(tactr(trluse),[1:3 5:end]); % 4th is irrelavent
% numtac=min(size(cfgtr.trl,1),length(time_touch));
% cfgtr.trl=cfgtr.trl(1:numtac,:);
% cfgtr.trl(:,3)=cfgtr.trl(:,3)-round(time_touch(1:numtac)*raw_all_ica_rej{ff}.fsample)';
% cfgtr.trl(:,3)=cfgtr.trl(:,3)-round(time_touch/diff(time{1}(1:2)))'; %using the difference between time as the sampling rate
tactrl(:,3)=tactrl(:,3)-round(time_touch/diff(time{1}(1:2)))'; %using the difference between time as the sampling rate

% Some time_touch are NaN which means we don't know trial onset.  Throw all those out.
tactrl(find(isnan(tactrl(:,3))),:)=[];

find(abs(diff(time_touch))>.01)

%  what's going on here?
%   load(['raw_all_ica_rej_' sub{ii}])
%   try
%     raw_all_ica_rej=raw_all_rej;clear raw_all_rej
%   end

cfg=[];
cfg.trl=tactrl;
raw_tac=ft_redefinetrial(cfg,raw_hpf);
clear cfgtactr

% Auditory
% need to get timing for each trial centered on stimulus correctly
cfg=[];
% cfg.dataset=files(ff).name;
% cfg.dataset=['cc_cc_spm8_' sub{ii} '_01s.mat'];
cfg.dataset=files.name;
cfg.trialfun='ft_trialfun_general_sleep';
cfg.trialdef.eventtype  = 'Stimulus';
cfg.trialdef.eventvalue = 'S 10'; % corresponds to lj.prepareStrobe(10)
cfg.trialdef.prestim = prestim;
cfg.trialdef.poststim = poststim;
cfg.time=time;
cfgaudtr=ft_definetrial(cfg);

for ll=1:length(cfgtr.event),aud(ll)=strcmp(cfgtr.event(ll).value,'S 10');end
audind=find(aud);
for ll=1:length(nsind)-1
  cfgauddiff(ll)=length(find(aud(nsind(ll):nsind(ll+1))))-numaud(ll);
end
cfgauddiff(length(nsind))=length(find(aud(nsind(ll+1):end)))-numaud(ll+1);
ss=0;
trluse=[];
for ll=1:length(nsind)
  %     trluse=[trluse [ss+1]:[ss+numtac(ll)-cfgtacdiff(ll)]];
  trluse=[trluse [ss+1]:[ss+numaud(ll)]];
  ss=ss+numaud(ll)+cfgauddiff(ll);
end
if length(trluse)~=length(audtrials) % should be the same
  error('num aud wrong');
end

audtr1=find(cfgtr.trl(:,12));
if length(trluse)~=size(cfgaudtr.trl,1) % due to prenar prestim there might be differences
  if length(trluse)-size(cfgaudtr.trl,1) ==2 % safe enought to assume chop off a trial at either end
    trluse=trluse(2:end-1);
  elseif length(trluse)-size(cfgaudtr.trl,1) ==1 % figure out if its beginning or end to chop off
    if any(cfgaudtr.trl(1:5,6)-cfgtr.trl(audtr1(trluse(1:5)),6)) && ~any(cfgaudtr.trl(end-4:end,6)-cfgtr.trl(audtr1(trluse(end-4:end)),6))
          trluse=trluse(2:end); % remove first trial
    elseif any(cfgaudtr.trl(end-4:end,6)-cfgtr.trl(audtr1(trluse(end-4:end)),6)) && ~any(cfgaudtr.trl(1:5,6)-cfgtr.trl(audtr1(trluse(1:5)),6))
      trluse=trluse(1:end-1);
    else
      disp('not sure what to throw out')
      keyboard
    end      
  else
    disp('why number of trial mismatch?')
    keyboard
  end  
end
% if sleep
%   if bsoarealextra(1)==-1 || bsoarealextra(1)>0
%     trluse=trluse(2:end)-1;
%   end
%   if bsoarealextra(end)==-1 || bsoarealextra(end)>0
%     trluse=trluse(1:end-1);
%   end
% else
%   switch ii
%     case [8 ]
%     otherwise
%       if ssoarealextra(1)==-1 || ssoarealextra(1)>0
%         trluse=trluse(2:end)-1;
%       end
%   end
%   switch ii
%     case [5 8]
%     otherwise
%       if ssoarealextra(end)==-1 || ssoarealextra(end)>0
%         trluse=trluse(1:end-1);
%       end
%   end
% end
audtrl=[cfgaudtr.trl cfgtr.trl(audtr1(trluse),7:end)];
% audtrl=cfgtr.trl(audtr1(trluse),[1:3 5:end]); % 4th is irrelavent
% cfgtr.trl(:,3)=cfgtr.trl(:,3)-round(fin.audDelay/1000*raw_all_ica_rej{ff}.fsample);
% cfgtr.trl(:,3)=cfgtr.trl(:,3)-round(.001*audDelay/diff(time{1}(1:2)))'; %using the difference between time as the sampling rate
audtrl(:,3)=audtrl(:,3)-round(.001*audDelay/diff(time{1}(1:2)))'; %using the difference between time as the sampling rate
cfg=[];
cfg.trl=audtrl;
raw_aud=ft_redefinetrial(cfg,raw_hpf);
clear cfgaudtr

% Null events
cfg=[];
% cfg.dataset=files(ff).name;
% cfg.dataset=['cc_cc_spm8_' sub{ii} '_01s.mat'];
cfg.dataset=files.name;
cfg.trialfun='ft_trialfun_general_sleep';
cfg.trialdef.eventtype  = 'Stimulus';
cfg.trialdef.eventvalue = 'S  1'; % corresponds to lj.prepareStrobe(10)
cfg.trialdef.prestim = prestim;
cfg.trialdef.poststim = poststim;
cfg.time=time;
cfgnultr=ft_definetrial(cfg);

for ll=1:length(cfgtr.event),nul(ll)=strcmp(cfgtr.event(ll).value,'S  1');end
nulind=find(nul);
for ll=1:length(nsind)-1
  cfgnuldiff(ll)=length(find(nul(nsind(ll):nsind(ll+1))))-numnul(ll);
end
cfgnuldiff(length(nsind))=length(find(nul(nsind(ll+1):end)))-numnul(ll+1);
ss=0;
trluse=[];
for ll=1:length(nsind)
  %     trluse=[trluse [ss+1]:[ss+numtac(ll)-cfgtacdiff(ll)]];
  trluse=[trluse [ss+1]:[ss+numnul(ll)]];
  ss=ss+numnul(ll)+cfgnuldiff(ll);
end
if length(trluse)~=length(find(nulseq)) % should be the same
  error('num nul wrong');
end
nultr=find(~cfgtr.trl(:,6));
if length(trluse)~=size(cfgnultr.trl,1) % due to prenar prestim there might be differences
  if length(trluse)-size(cfgnultr.trl,1) ==2 % safe enought to assume chop off a trial at either end
    trluse=trluse(2:end-1);
  elseif length(trluse)-size(cfgnultr.trl,1) ==1 % figure out if its beginning or end to chop off
    if any(cfgnultr.trl(1:5,6)-cfgtr.trl(nultr(trluse(1:5)),6)) && ~any(cfgnultr.trl(end-4:end,6)-cfgtr.trl(nultr(trluse(end-4:end)),6))
      trluse=trluse(2:end); % remove first trial
    elseif any(cfgnultr.trl(end-4:end,6)-cfgtr.trl(nultr(trluse(end-4:end)),6)) && ~any(cfgnultr.trl(1:5,6)-cfgtr.trl(nultr(trluse(1:5)),6))
      trluse=trluse(1:end-1);
    else
      disp('not sure what to throw out')
      keyboard
    end
  else
    disp('why number of trial mismatch?')
    keyboard
  end
end

% if sleep
%   if bsoarealextra(1)==0
%     trluse=trluse(2:end)-1;
%     cfgnultr.trl=cfgnultr.trl(2:end,:);
%   end
%   if bsoarealextra(end)==0
%     trluse=trluse(1:end-1);
%     cfgnultr.trl=cfgnultr.trl(1:end-1,:);
%   end
% else
%   switch ii
%     case [8 ]
%     otherwise
%       if ssoarealextra(1)==0
%         trluse=trluse(2:end)-1;
%         cfgnultr.trl=cfgnultr.trl(2:end,:);
%       end
%   end
%   switch ii
%     case [8 ]
%     otherwise
%       if ssoarealextra(end)==0
%         trluse=trluse(1:end-1);
%         cfgnultr.trl=cfgnultr.trl(1:end-1,:);
%       end
%   end
% end
nultrl=[cfgnultr.trl cfgtr.trl(nultr(trluse),7:end)];
% nultrl=cfgtr.trl(nultr(trluse),[1:3 5:end]); % 4th is irrelavent
% cfgtr.trl=cfgtr.trl(trluse,:);
% [size(cfgtr.trl,1) length(find(snulseq))] % should be the same

% use the auditory delay; it doesn't really matter, as it's random anyway
cfgtr.trl(:,3)=cfgtr.trl(:,3)-round(.001*audDelay/diff(time{1}(1:2)))'; %using the difference between time as the sampling rate
cfg=[];
cfg.trl=nultrl;
raw_nul=ft_redefinetrial(cfg,raw_hpf);
clear cfgnultr

N=raw_hpf.sampleinfo(end);
clear raw_hpf


if length(intersect(find(diff(raw_tac.trialinfo(:,3:4),[],2)) , find(raw_tac.trialinfo(:,4)~=-3)))>10
  disp('trl misaligned?'); keyboard
end
if length(intersect(find(diff(raw_aud.trialinfo(:,3:4),[],2)) , find(raw_aud.trialinfo(:,4)~=-3)))>10
  disp('trl misaligned?'); keyboard
end
if length(intersect(find(diff(raw_nul.trialinfo(:,3:4),[],2)) , find(raw_nul.trialinfo(:,4)~=-3)))>10
  disp('trl misaligned?'); keyboard
end
  
tmp=[corr(raw_tac.trialinfo(raw_tac.trialinfo(:,4)~=-3,3:4)) corr(raw_aud.trialinfo(raw_aud.trialinfo(:,4)~=-3,3:4)) corr(raw_nul.trialinfo(raw_nul.trialinfo(:,4)~=-3,3:4))];
if any(tmp(1,[2 4 6])<.99)
  disp('trl misaligned?'); keyboard
end

%%
% Add sleep, artifact, Kc, spindle tags

artnames=fieldnames(artfctdef);
samp.tmp=[];
raw_tac.trialinfo(:,12:169)=0; % put prestim in 12:19 and poststim in 22:29
raw_aud.trialinfo(:,12:169)=0;
raw_nul.trialinfo(:,12:169)=0;
for aa=1:length(artnames)
  switch artnames{aa}
    case 'muscle'
      trialcol=12;
    case 'eog'
      trialcol=13;
    case 'kc'
      trialcol=14;
    case 'delta'
      trialcol=15;
    case 'sw'
      trialcol=16;
    case 'sp_fast'
      trialcol=17;
    case 'sp_slow'
      trialcol=18;
    case 'sleepblink'
      trialcol=19;
    case {'n2' 'n1' 'wake' 'n3' 'r' 'visual'}
      trialcol=nan;
    otherwise
      error('what did I miss?')
  end
  
  % 12-19: percentage of samples(time) during pre-stim (5s) involved in this feature
  % 22-29: percentage of samples(time) during post-stim (5s) involved in this feature
  % 32-39: percentage of samples(time) during pre-stim (1s) involved in this feature
  % 42-49: percentage of samples(time) during post-stim (1s) involved in this feature
  
  % 50-59:   [amplitude negmax(s)] pairs in prestim of Kc (up to 5 pairs alllowed)
  % 60-69:   [amplitude negmax(s)] pairs in poststim of Kc (up to 5 pairs alllowed)
  % 70-79:   [amplitude negmax(s)] pairs in prestim of delta (up to 5 pairs alllowed)
  % 80-89:   [amplitude negmax(s)] pairs in poststim of delta (up to 5 pairs alllowed)
  % 90-99:   [amplitude negmax(s)] pairs in prestim of sw (up to 5 pairs alllowed)
  % 100-109: [amplitude negmax(s)] pairs in poststim of sw (up to 5 pairs alllowed)
  
  if isfield(artfctdef.(artnames{aa}),'artifact') && ~isnan(trialcol)
    samp.(artnames{aa})=binarise_artifact_begendpoints(artfctdef.(artnames{aa}).artifact,N);
    for ll=1:size(raw_tac.sampleinfo,1)
      %       raw_tac.trialinfo(ll,trialcol)   =any(samp.(artnames{aa})(raw_tac.sampleinfo(ll,1):raw_tac.sampleinfo(ll,1)+dsearchn(raw_tac.time{ll}',0)));
      %       raw_tac.trialinfo(ll,trialcol+10)=any(samp.(artnames{aa})(raw_tac.sampleinfo(ll,1)+dsearchn(raw_tac.time{ll}',0):raw_tac.sampleinfo(ll,2)));
      raw_tac.trialinfo(ll,trialcol)   =mean(samp.(artnames{aa})(raw_tac.sampleinfo(ll,1):raw_tac.sampleinfo(ll,1)+dsearchn(raw_tac.time{ll}',0)));
      raw_tac.trialinfo(ll,trialcol+10)=mean(samp.(artnames{aa})(raw_tac.sampleinfo(ll,1)+dsearchn(raw_tac.time{ll}',0):raw_tac.sampleinfo(ll,2)));
      raw_tac.trialinfo(ll,trialcol+20)=mean(samp.(artnames{aa})(dsearchn(raw_tac.time{ll}',0)-prenar*1000: dsearchn(raw_tac.time{ll}',0)));
      raw_tac.trialinfo(ll,trialcol+30)=mean(samp.(artnames{aa})(dsearchn(raw_tac.time{ll}',0):dsearchn(raw_tac.time{ll}',0)+prenar*1000));
    end
    for ll=1:size(raw_aud.sampleinfo,1)
      raw_aud.trialinfo(ll,trialcol)   =mean(samp.(artnames{aa})(raw_aud.sampleinfo(ll,1):raw_aud.sampleinfo(ll,1)+dsearchn(raw_aud.time{ll}',0)));
      raw_aud.trialinfo(ll,trialcol+10)=mean(samp.(artnames{aa})(raw_aud.sampleinfo(ll,1)+dsearchn(raw_aud.time{ll}',0):raw_aud.sampleinfo(ll,2)));
      %       raw_aud.trialinfo(ll,trialcol)   =any(samp.(artnames{aa})(raw_aud.sampleinfo(ll,1):raw_aud.sampleinfo(ll,1)+dsearchn(raw_aud.time{ll}',0)));
      %       raw_aud.trialinfo(ll,trialcol+10)=any(samp.(artnames{aa})(raw_aud.sampleinfo(ll,1)+dsearchn(raw_aud.time{ll}',0):raw_aud.sampleinfo(ll,2)));
      %       raw_aud.trialinfo(ll,trialcol)=any(samp.(artnames{aa})(raw_aud.sampleinfo(ll,1):raw_aud.sampleinfo(ll,2)));
      raw_aud.trialinfo(ll,trialcol+20)=mean(samp.(artnames{aa})(dsearchn(raw_aud.time{ll}',0)-prenar*1000: dsearchn(raw_aud.time{ll}',0)));
      raw_aud.trialinfo(ll,trialcol+30)=mean(samp.(artnames{aa})(dsearchn(raw_aud.time{ll}',0):dsearchn(raw_aud.time{ll}',0)+prenar*1000));
    end
    for ll=1:size(raw_nul.sampleinfo,1)
      raw_nul.trialinfo(ll,trialcol)   =mean(samp.(artnames{aa})(raw_nul.sampleinfo(ll,1):raw_nul.sampleinfo(ll,1)+dsearchn(raw_nul.time{ll}',0)));
      raw_nul.trialinfo(ll,trialcol+10)=mean(samp.(artnames{aa})(raw_nul.sampleinfo(ll,1)+dsearchn(raw_nul.time{ll}',0):raw_nul.sampleinfo(ll,2)));
      %       raw_nul.trialinfo(ll,trialcol)   =any(samp.(artnames{aa})(raw_nul.sampleinfo(ll,1):raw_nul.sampleinfo(ll,1)+dsearchn(raw_nul.time{ll}',0)));
      %       raw_nul.trialinfo(ll,trialcol+10)=any(samp.(artnames{aa})(raw_nul.sampleinfo(ll,1)+dsearchn(raw_nul.time{ll}',0):raw_nul.sampleinfo(ll,2)));
      %       raw_nul.trialinfo(ll,trialcol)=any(samp.(artnames{aa})(raw_nul.sampleinfo(ll,1):raw_nul.sampleinfo(ll,2)));
      raw_nul.trialinfo(ll,trialcol+20)=mean(samp.(artnames{aa})(dsearchn(raw_nul.time{ll}',0)-prenar*1000: dsearchn(raw_nul.time{ll}',0)));
      raw_nul.trialinfo(ll,trialcol+30)=mean(samp.(artnames{aa})(dsearchn(raw_nul.time{ll}',0):dsearchn(raw_nul.time{ll}',0)+prenar*1000));
    end
    samp=rmfield(samp,artnames{aa}); % clearing as this can get large
    
    % save out BigWave negpeak time and max amplitude
    tacusepre=zeros(size(raw_tac.trialinfo,1),1);
    audusepre=zeros(size(raw_aud.trialinfo,1),1);
    nulusepre=zeros(size(raw_nul.trialinfo,1),1);
    tacusepost=zeros(size(raw_tac.trialinfo,1),1);
    audusepost=zeros(size(raw_aud.trialinfo,1),1);
    nulusepost=zeros(size(raw_nul.trialinfo,1),1);
    if strcmp(artnames{aa},'kc') || strcmp(artnames{aa},'delta') || strcmp(artnames{aa},'sw')
      if strcmp(artnames{aa},'kc')
        colampstartpre=50;
        colnegstartpre=60;
        colampstartpost=70;        
        colnegstartpost=80;        
      elseif strcmp(artnames{aa},'delta')
        colampstartpre=90;
        colnegstartpre=100;
        colampstartpost=110;        
        colnegstartpost=120;        
      elseif strcmp(artnames{aa},'sw')
        colampstartpre=130;
        colnegstartpre=140;
        colampstartpost=150;        
        colnegstartpost=160;        
      end
      for kk=1:length(artfctdef.(artnames{aa}).negmax)
        if ~isempty( find(raw_tac.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_tac.sampleinfo(:,2)))
          ind=find(raw_tac.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_tac.sampleinfo(:,2));
          for nn=1:length(ind)
            if raw_tac.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_tac.sampleinfo(ind(nn),1)) < 0  % prestim
              if tacusepre(ind(nn))>9
                error('too many BigWaves on this trial')
              end
              raw_tac.trialinfo(ind(nn),colampstartpre+tacusepre(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
              raw_tac.trialinfo(ind(nn),colnegstartpre+tacusepre(ind(nn)))=raw_tac.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_tac.sampleinfo(ind(nn),1));
              tacusepre(ind(nn))=tacusepre(ind(nn))+1;
            else % poststim
              if tacusepost(ind(nn))>9
                error('too many BigWaves on this trial')
              end
              raw_tac.trialinfo(ind(nn),colampstartpost+tacusepost(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
              raw_tac.trialinfo(ind(nn),colnegstartpost+tacusepost(ind(nn)))=raw_tac.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_tac.sampleinfo(ind(nn),1));
%               raw_tac.trialinfo(ind(nn),colstartpost+tacusepost(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
%               raw_tac.trialinfo(ind(nn),colstartpost+1+tacusepost(ind(nn)))=raw_tac.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_tac.sampleinfo(ind(nn),1));
              tacusepost(ind(nn))=tacusepost(ind(nn))+1;
            end
          end
        end
        if ~isempty( find(raw_aud.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_aud.sampleinfo(:,2)))
          ind=find(raw_aud.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_aud.sampleinfo(:,2));
          for nn=1:length(ind)
            if raw_aud.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_aud.sampleinfo(ind(nn),1)) < 0  % prestim
              if audusepre(ind(nn))>9
                error('too many BigWaves on this trial')
              end
              raw_aud.trialinfo(ind(nn),colampstartpre+audusepre(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
              raw_aud.trialinfo(ind(nn),colnegstartpre+audusepre(ind(nn)))=raw_aud.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_aud.sampleinfo(ind(nn),1));
%               raw_aud.trialinfo(ind(nn),colstartpre+audusepre(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
%               raw_aud.trialinfo(ind(nn),colstartpre+1+audusepre(ind(nn)))=raw_aud.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_aud.sampleinfo(ind(nn),1));
              audusepre(ind(nn))=audusepre(ind(nn))+1;
            else % poststim
              if audusepost(ind(nn))>9
                error('too many BigWaves on this trial')
              end
              raw_aud.trialinfo(ind(nn),colampstartpost+audusepost(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
              raw_aud.trialinfo(ind(nn),colnegstartpost+audusepost(ind(nn)))=raw_aud.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_aud.sampleinfo(ind(nn),1));
%               raw_aud.trialinfo(ind(nn),colstartpost+audusepost(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
%               raw_aud.trialinfo(ind(nn),colstartpost+1+audusepost(ind(nn)))=raw_aud.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_aud.sampleinfo(ind(nn),1));
              audusepost(ind(nn))=audusepost(ind(nn))+1;
            end
          end
        end
        if ~isempty( find(raw_nul.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_nul.sampleinfo(:,2)))
          ind=find(raw_nul.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_nul.sampleinfo(:,2));
          for nn=1:length(ind)
            if raw_nul.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_nul.sampleinfo(ind(nn),1)) < 0  % prestim
              if nulusepre(ind(nn))>9
                error('too many BigWaves on this trial')
              end
              raw_nul.trialinfo(ind(nn),colampstartpre+nulusepre(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
              raw_nul.trialinfo(ind(nn),colnegstartpre+nulusepre(ind(nn)))=raw_nul.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_nul.sampleinfo(ind(nn),1));
%               raw_nul.trialinfo(ind(nn),colstartpre+nulusepre(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
%               raw_nul.trialinfo(ind(nn),colstartpre+1+nulusepre(ind(nn)))=raw_nul.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_nul.sampleinfo(ind(nn),1));
              nulusepre(ind(nn))=nulusepre(ind(nn))+1;
            else % poststim
              if nulusepost(ind(nn))>9
                error('too many BigWaves on this trial')
              end
              raw_nul.trialinfo(ind(nn),colampstartpost+nulusepost(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
              raw_nul.trialinfo(ind(nn),colnegstartpost+nulusepost(ind(nn)))=raw_nul.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_nul.sampleinfo(ind(nn),1));
%               raw_nul.trialinfo(ind(nn),colstartpost+nulusepost(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
%               raw_nul.trialinfo(ind(nn),colstartpost+1+nulusepost(ind(nn)))=raw_nul.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_nul.sampleinfo(ind(nn),1));
              nulusepost(ind(nn))=nulusepost(ind(nn))+1;
            end
          end
        end
        
%         if ~isempty( find(raw_tac.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)< raw_tac.sampleinfo(:,1)+prestim*1000))
%           ind=find(raw_tac.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)< raw_tac.sampleinfo(:,1)+prestim*1000);
% %           ind=find(raw_tac.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_tac.sampleinfo(:,2));
%           for nn=1:length(ind)
%             if tacusepre(ind(nn))>9
%               error('too many BigWaves on this trial')
%             end
%             raw_tac.trialinfo(ind(nn),colstartpre+tacusepre(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
%             raw_tac.trialinfo(ind(nn),colstartpre+1+tacusepre(ind(nn)))=raw_tac.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_tac.sampleinfo(ind(nn),1));
%             tacusepre(ind(nn))=tacusepre(ind(nn))+2;
% %             raw_tac.trialinfo(ind(nn),trialcol+20)=raw_tac.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_tac.sampleinfo(ind(nn),1));
% %             raw_tac.trialinfo(ind(nn),trialcol+30)=artfctdef.(artnames{aa}).amplitude(kk);
%           end
%         end
%         if ~isempty( find( raw_tac.sampleinfo(:,1)+prestim*1000 <artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_tac.sampleinfo(:,2)  ))
%           ind=find( raw_tac.sampleinfo(:,1)+prestim*1000 <artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_tac.sampleinfo(:,2) );
%           for nn=1:length(ind)
%             if ind(nn)==440
%               keyboard
%             end
%             if tacusepost(ind(nn))>9
%               error('too many BigWaves on this trial')
%             end
%             raw_tac.trialinfo(ind(nn),colstartpost+tacusepost(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
%             raw_tac.trialinfo(ind(nn),colstartpost+1+tacusepost(ind(nn)))=raw_tac.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_tac.sampleinfo(ind(nn),1));
%             tacusepost(ind(nn))=tacusepost(ind(nn))+2;
%           end
%         end
%         if ~isempty( find(raw_aud.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)< raw_aud.sampleinfo(:,1)+prestim*1000))
%           ind=find(raw_aud.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)< raw_aud.sampleinfo(:,1)+prestim*1000);
%           for nn=1:length(ind)
%             if audusepre(ind(nn))>9
%               error('too many BigWaves on this trial')
%             end
%             raw_aud.trialinfo(ind(nn),colstartpre+audusepre(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
%             raw_aud.trialinfo(ind(nn),colstartpre+1+audusepre(ind(nn)))=raw_aud.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_aud.sampleinfo(ind(nn),1));
%             audusepre(ind(nn))=audusepre(ind(nn))+2;
%           end
%         end
%         if ~isempty( find( raw_aud.sampleinfo(:,1)+prestim*1000 <artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_aud.sampleinfo(:,2)  ))
%           ind=find( raw_aud.sampleinfo(:,1)+prestim*1000 <artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_aud.sampleinfo(:,2) );
%           for nn=1:length(ind)
%             if audusepost(ind(nn))>9
%               error('too many BigWaves on this trial')
%             end
%             raw_aud.trialinfo(ind(nn),colstartpost+audusepost(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
%             raw_aud.trialinfo(ind(nn),colstartpost+1+audusepost(ind(nn)))=raw_aud.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_aud.sampleinfo(ind(nn),1));
%             audusepost(ind(nn))=audusepost(ind(nn))+2;
%           end
%         end
%         if ~isempty( find(raw_nul.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)< raw_nul.sampleinfo(:,1)+prestim*1000))
%           ind=find(raw_nul.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)< raw_nul.sampleinfo(:,1)+prestim*1000);
%           for nn=1:length(ind)
%             if nulusepre(ind(nn))>9
%               error('too many BigWaves on this trial')
%             end
%             raw_nul.trialinfo(ind(nn),colstartpre+nulusepre(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
%             raw_nul.trialinfo(ind(nn),colstartpre+1+nulusepre(ind(nn)))=raw_nul.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_nul.sampleinfo(ind(nn),1));
%             nulusepre(ind(nn))=nulusepre(ind(nn))+2;
%           end
%         end
%         if ~isempty( find( raw_nul.sampleinfo(:,1)+prestim*1000 <artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_nul.sampleinfo(:,2)  ))
%           ind=find( raw_nul.sampleinfo(:,1)+prestim*1000 <artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_nul.sampleinfo(:,2) );
%           for nn=1:length(ind)
%             if nulusepost(ind(nn))>9
%               error('too many BigWaves on this trial')
%             end
%             raw_nul.trialinfo(ind(nn),colstartpost+nulusepost(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
%             raw_nul.trialinfo(ind(nn),colstartpost+1+nulusepost(ind(nn)))=raw_nul.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_nul.sampleinfo(ind(nn),1));
%             nulusepost(ind(nn))=nulusepost(ind(nn))+2;
%           end
%         end

%         if ~isempty( find(raw_tac.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_tac.sampleinfo(:,2)))
%           ind=find(raw_tac.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_tac.sampleinfo(:,2));
%           for nn=1:length(ind)
%             raw_tac.trialinfo(ind(nn),trialcol+20)=raw_tac.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_tac.sampleinfo(ind(nn),1));
%             raw_tac.trialinfo(ind(nn),trialcol+30)=artfctdef.(artnames{aa}).amplitude(kk);
%           end
%         end        
%         if ~isempty( find(raw_aud.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_aud.sampleinfo(:,2)))
%           ind=find(raw_aud.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_aud.sampleinfo(:,2));
%           for nn=1:length(ind)
%             raw_aud.trialinfo(ind(nn),trialcol+20)=raw_aud.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_aud.sampleinfo(ind(nn),1));
%             raw_aud.trialinfo(ind(nn),trialcol+30)=artfctdef.(artnames{aa}).amplitude(kk);
%           end
%         end
%         if ~isempty( find(raw_nul.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_nul.sampleinfo(:,2)))
%           ind=find(raw_nul.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_nul.sampleinfo(:,2));
%           for nn=1:length(ind)
%             raw_nul.trialinfo(ind(nn),trialcol+20)=raw_nul.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_nul.sampleinfo(ind(nn),1));
%             raw_nul.trialinfo(ind(nn),trialcol+30)=artfctdef.(artnames{aa}).amplitude(kk);
%           end
%         end  
      end
    end
    
  elseif ~isnan(trialcol) % no artifact of this type; fill in column of zeros
    raw_tac.trialinfo(:,trialcol)=0;
    raw_aud.trialinfo(:,trialcol)=0;
    raw_nul.trialinfo(:,trialcol)=0;
  end
end % aa

clear cfgtr


% some trials have been discarded by only setting to NaN.  Here we
% explicitly reject them from the data structure.
keeptrials=ones(1,length(raw_tac.time));
for ll=1:length(raw_tac.time)
  if isnan(raw_tac.time{ll}(1))
    keeptrials(ll)=0;
  end
end
cfg=[];
cfg.trials=find(keeptrials);
raw_tac=ft_redefinetrial(cfg,raw_tac);

keeptrials=ones(1,length(raw_aud.time));
for ll=1:length(raw_aud.time)
  if isnan(raw_aud.time{ll}(1))
    keeptrials(ll)=0;
  end
end
cfg=[];
cfg.trials=find(keeptrials);
raw_aud=ft_redefinetrial(cfg,raw_aud);

keeptrials=ones(1,length(raw_nul.time));
for ll=1:length(raw_nul.time)
  if isnan(raw_nul.time{ll}(1))
    keeptrials(ll)=0;
  end
end
cfg=[];
cfg.trials=find(keeptrials);
raw_nul=ft_redefinetrial(cfg,raw_nul);


%   if 0 % not needed now with muscle and EOG identification ?
%     cfg=[];
%     cfg.method='channel';
%     cfg.alim=100;
%     raw_tac_rej=ft_rejectvisual(cfg,raw_tac);
%
%     tacart=raw_tac_rej.cfg.artfctdef;
%
%     cfg=[];
%     cfg.artfctdef=tacart;
%     raw_aud_rej=ft_rejectartifact(cfg,raw_aud);
%
%     cfg=[];
%     cfg.method='channel';
%     cfg.alim=100;
%     raw_aud_rej=ft_rejectvisual(cfg,raw_aud_rej);
%
%     audart=raw_aud_rej.cfg.artfctdef;
%
%     % if auditory was thrown out here that tactile still in, then it will later
%     % be thrown out in eeg_legomagic_trialSelection1
%
%     raw_nul_rej=ft_rejectvisual(cfg,raw_nul);
%     nulart=raw_nul_rej.cfg.artfctdef;
%   end

%   save(['raw2_each_rej_' sub{ii} '_sleep' num2str(sleep)],'raw_*','artflag','-v7.3'); % with 0.2Hz highpass before epoching
%   clear raw_aud raw_nul raw_tac


return
