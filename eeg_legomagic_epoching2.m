function [raw_tac, raw_aud, raw_nul, artflag]=eeg_legomagic_epoching2(ii,sleep,featfull,saveflag)
% function [raw_tac, raw_aud, raw_nul, artflag]=eeg_legomagic_epoching2(ii,sleep,featfull,saveflag)
%
% ii:       subject index
% sleep:    0 sitting, 1 bed
% featfull: 0 features with minimal channels, 1 for channels but minimal features
% saveflag: 0 not save to disk raw* results, 1 do save to disk raw* results

% preprocessing of EEG data from Hills, 64ch MR-cap
tomflag=0; % set to 1 if you are Tom (0 for Johanna)
if ispc
  if tomflag
    edir='G:\jz_sleep_data\data\';
    ddir='E:\JohannaProject\diaries\';
  else
    edir='D:\audtac\eeg_data\';
    ddir='D:\audtac\legomagic\diaries\';
    bdir='D:\audtac\behav_data\';
    sdir='D:\audtac\spss_stuff\';
  end
else
  [~,hostname]=system('hostname');
  if ~isempty(strfind(hostname,'les')) | ~isempty(strfind(hostname,'LES')) % either COLLES-151401 or LES-LINUX_FS3
    edir='~/audtac/eeg_data/';
    ddir='~/audtac/legomagic/diaries/';
    bdir='~/audtac/behav_data/';
    sdir='~/audtac/spss_stuff/';
  else % assume on VM linux of psychl-132432
    edir='/mnt/hgfs/D/audtac/eeg_data/';
    ddir='/mnt/hgfs/D/audtac/legomagic/diaries/';
    bdir='/mnt/hgfs/D/audtac/behav_data/';
    sdir='/mnt/hgfs/D/audtac/spss_stuff/';
  end
end
cd(edir)

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

% see eeg_legomagic_viewFeatures.m
if sleep
  try
    load('artfctdef_bed_automan.mat')
    artfctdef=artfctdef_new;clear artfctdef_new
    artflag=1;
  catch
    disp('shame, need to do manual sorting of EOG here')
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
    disp('shame, need to do manual sorting of EOG here')
    load('artfctdef_sit_auto.mat')
    artflag=0;
  end
  files=dir('cc*s.mat');
end


%     if 0
%       load([cfg.dataset(1:end-4) '_event.mat']); % loads event saved out during ft_trialfun_general_sleep.m
%       eventS10=ft_filter_event(event,'value','S 10');
%       eventS2=ft_filter_event(event,'value','S  2');
%       eventS1=ft_filter_event(event,'value','S  1');
%       [eventS10.type]=deal('A');
%       [eventS2.type]=deal('T');
%       [eventS1.type]=deal('N');
%       %   [eventS2.offset]=deal(0);
%       for ll=1:length(eventS2),eventS2(ll).sample=eventS2(ll).sample+200;end
%       eventStim=[eventS2 eventS10 eventS1];
%
%       cfg=[];
%       cfg.layout='EEG1010.lay';
%       cfg.viewmode='vertical';
%       cfg.blocksize=30;
%       cfg.plotlabels='yes';
%       cfg.channel={'Fz' 'Cz' 'Pz' 'Oz'};
%       cfg.event=eventStim;
%       cfg.ploteventlabels='colorvalue';
%       dcfg=ft_databrowser(cfg,raw_all_ica_rej);
%     end

%     % create trl with 30s segments
%     for ll=1:ceil( raw_hpf.sampleinfo(2)/30000)
%       if ll==ceil( raw_hpf.sampleinfo(2)/30000)
%         trl(ll,:)=[(30000*(ll-1))+1 raw_hpf.sampleinfo(2) 0];
%       elseif ll==1
%         trl(ll,:)=[200 (30000*(ll-1))+30000 0];
%       else
%         trl(ll,:)=[(30000*(ll-1))+1 (30000*(ll-1))+30000 0];
%       end
%     end

%     if 0
%       % Kcomplex
%       cfg=[];
%       cfg.trl=trl; % 30s segments
%       cfg.continuous='yes';
%       cfg.artfctdef.zvalue.channel={'Pz'};
%       cfg.artfctdef.zvalue.cutoff=4;
%       cfg.artfctdef.zvalue.interactive = 'yes';
%       cfg.artfctdef.zvalue.lpfilter='yes';
%       cfg.artfctdef.zvalue.lpfreq=8;
%       cfg.artfctdef.zvalue.lpfiltord=3;
%       [cfg, artifact] = ft_artifact_zvalue(cfg, raw_all_ica_rej)
%
%       cfg=[];
%       cfg.layout='EEG1010.lay';
%       cfg.viewmode='vertical';
%       cfg.blocksize=30;
%       cfg.plotlabels='yes';
%       cfg.channel={'Fz' 'Cz' 'Pz' 'Oz'};
%       cfg.artfctdef.zvalue.artifact=artifact;
%       cfg.preproc.lpfilter='yes';
%       cfg.preproc.lpfreq=30;
%       dcfg=ft_databrowser(cfg,raw_all_ica_rej);
%
%       % muscle
%       cfg=[];
%       cfg.trl=trl(2:end-1,:); % 30s segments
%       cfg.continuous='yes';
%       cfg.artfctdef.muscle.channel={'FT7' 'FT8' 'TP7' 'TP8'};
%       cfg.artfctdef.muscle.bpfilter    = 'yes';
%       cfg.artfctdef.muscle.bpfreq      = [80 140];
%       cfg.artfctdef.muscle.bpfiltord   = 8;
%       cfg.artfctdef.muscle.bpfilttype  = 'but';
%       cfg.artfctdef.muscle.hilbert     = 'yes';
%       cfg.artfctdef.muscle.cutoff = 2;
%       cfg.artfctdef.muscle.interactive = 'yes';
%       %   cfg.artfctdef.muscle.trlpadding  = 1;
%       %   cfg.artfctdef.muscle.fltpadding  = 1;
%       cfg.artfctdef.muscle.artpadding  = 3;
%       [cfg,martifact]=ft_artifact_muscle(cfg,raw_all_ica_rej);
%
%       cfg=[];
%       cfg.layout='EEG1010.lay';
%       cfg.viewmode='vertical';
%       cfg.blocksize=30;
%       cfg.plotlabels='yes';
%       cfg.channel={'Fz' 'Cz' 'Pz' 'Oz'};
%       cfg.artfctdef.zvalue.artifact=martifact;
%       %   cfg.preproc.lpfilter='yes';
%       %   cfg.preproc.lpfreq=30;
%       dcfg=ft_databrowser(cfg,raw_all_ica_rej);
%
%     end
%


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
end

% this is to ensure enough memory to call ft_preprocessing next
if ispc
  S=memory; % memory doesn't exist on linux
  gbavail=S.MemAvailableAllArrays/(1024^3);
else
  [aa,bb]=system('head /proc/meminfo');
  if ~isempty(strfind(hostname,'les')) % les-linux-fs3
    gbavail=str2num(bb(strfind(bb,'MemFree')+8:strfind(bb,'MemFree')+23))/(1024^2);
  elseif ~isempty(strfind(hostname,'LES')) || ~isempty(strfind(hostname,'psychl-132432')) % Linux VM of psychl-132432-1  % COLLES-151401
    gb1=str2num(bb(strfind(bb,'Inactive:')+9:strfind(bb,'Inactive:')+23))/(1024^2);
    gb2=str2num(bb(strfind(bb,'MemFree')+8:strfind(bb,'MemFree')+23))/(1024^2);
    gbavail=max(gb1,gb2);
  end
end
tic
while gbavail < 1 % don't push things if waiting a bit could help
  if toc<(10*60) % 10 minutes
    pause(10)
  else
    pause(120)
    disp('waiting here for memory')
  end
end

if featfull==1
  % file=dir(['cc_cc_spm8_' sub{ii} '*.mat']);
  cfg=[];
  cfg.dataset=files.name;
  cfg.demean='yes';
  cfg.bsfilter='yes';
  cfg.bsfreq=[49 51; 99 101; 149 151];
  cfg.hpfilter='yes';
  cfg.hpfiltord=3;
  cfg.hpfreq=0.2;  % Could do higher for awake data but not want higher for sleep.
  cfg.channel={'all' '-ECG' '-VEOG' '-HEOG' '-EMG'};
  cfg.reref='yes';
  cfg.refchannel='all'; % used linked-mastoids for sleep staging but use average reference for ERPs and TFRs and source localisation
  raw_hpf=ft_preprocessing(cfg);
  
  prestim=1.5;  % seconds to evaluate overall
  poststim=2.3;
  prenar=1; % more narrow focus
  postnar=1;
elseif featfull==0
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
  
  prestim=4;  % seconds to evaluate overall;  % changed from 5 to 4 Oct 2018
  poststim=4;
  prenar=1.2; % more narrow focus
  postnar=1.2; % changed from 1.0 to 1.2 on 22 Oct, 2018
end


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
cfg.trialdef.prestim = prenar; % kept this long from ii==12 and afterwards
cfg.trialdef.poststim = postnar;
cfg.time=time;
cfgtr=ft_definetrial(cfg);

if sleep
%   bsoarealextra=[bsoareal(1,[1 end])];
  soause=setdiff(1:size(bsoareal,2),cfgtr.trlselexcl); % sometimes first or last trial may be excluded if too close to dataset boundary as function of prestim/poststim length
  cfgtr.trl(:,7)=bsoareal(1,:); % 
  cfgtr.trl(:,8)=bsoareal(2,:); % <---- Final decision to use this (awake paper)
  cfgtr.trl(:,9)=bsoareal(3,:); % 
  cfgtr.trl(:,10)=bsoareal(4,:); % 
  cfgtr.trl(:,13)=bsoaeff;
else
%   ssoarealextra=[ssoareal(1,[1 end])];
  soause=setdiff(1:size(ssoareal,2),cfgtr.trlselexcl);
  cfgtr.trl(:,7)=ssoareal(1,:); % 
  cfgtr.trl(:,8)=ssoareal(2,:); % <---- Final decision to use this (awake paper)
  cfgtr.trl(:,9)=ssoareal(3,:); % 
  cfgtr.trl(:,10)=ssoareal(4,:); %
  cfgtr.trl(:,13)=ssoaeff;
end
cfgtr.trl(:,11)=cfgtr.trl(:,6)>0 | cfgtr.trl(:,6)==-2; % tactile trials
cfgtr.trl(:,12)=cfgtr.trl(:,6)>0 | cfgtr.trl(:,6)==-1; % auditory trials
cfgtr.trl(:,14)=1:size(cfgtr.trl,1);

% discovered a code '214' for some ii
if find(cfgtr.trl(:,6)<-2 | cfgtr.trl(:,6)>9)
  if [sleep==1 && any(find(ii==[8 9 10 11 15 16 17 18 20 21 22 23 24 25 26 27 28 29 30 31 32]))] || [sleep==0 && any(find(ii==[12 17 24 31 32]))]
    % known cases where this happens
  else
    disp('this is weird')
    cfgtr.trl(find(cfgtr.trl(:,6)<-2 | cfgtr.trl(:,6)>9),1:6)
    keyboard
  end
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
  error('num tac wrong');
end

tactr=find(cfgtr.trl(:,11));
if length(trluse)~=size(cfgtactr.trl,1) % due to prenar prestim there might be differences
  if length(trluse)-size(cfgtactr.trl,1) ==2 % safe enought to assume chop off a trial at either end
    trluse=trluse(2:end-1);
    time_touch=time_touch(2:end-1);
  elseif length(trluse)-size(cfgtactr.trl,1) ==1 % figure out if its beginning or end to chop off
    if any(cfgtactr.trl(1:5,6)-cfgtr.trl(tactr(trluse(1:5)),6)) && ~any(cfgtactr.trl(end-4:end,6)-cfgtr.trl(tactr(trluse(end-4:end)),6))
      trluse=trluse(2:end); % remove first trial
      time_touch=time_touch(2:end); % remove first trial
    elseif any(cfgtactr.trl(end-4:end,6)-cfgtr.trl(tactr(trluse(end-4:end)),6)) && ~any(cfgtactr.trl(1:5,6)-cfgtr.trl(tactr(trluse(1:5)),6))
      trluse=trluse(1:end-1);
      time_touch=time_touch(1:end-1);
    else
      disp('not sure what to throw out')
      keyboard
    end
  else
    disp('why number of trial mismatch?')
    keyboard
  end
end
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
nultr=find(cfgtr.trl(:,6)==0);

if length(trluse)~=size(cfgnultr.trl,1) % due to prenar prestim there might be differences
  if length(trluse)-size(cfgnultr.trl,1) ==2 % safe enought to assume chop off a trial at either end
    trluse=trluse(2:end-1);
  elseif length(trluse)-size(cfgnultr.trl,1) ==1 % figure out if its beginning or end to chop off
    if any(cfgnultr.trl(1:5,6)-cfgtr.trl(nultr(trluse(1:5)),6)) && ~any(cfgnultr.trl(end-4:end,6)-cfgtr.trl(nultr(trluse(end-4:end)),6))
      trluse=trluse(2:end); % remove first trial
    elseif any(cfgnultr.trl(end-4:end,6)-cfgtr.trl(nultr(trluse(end-4:end)),6)) && ~any(cfgnultr.trl(1:5,6)-cfgtr.trl(nultr(trluse(1:5)),6))
      trluse=trluse(1:end-1);
    elseif ii==17 && sleep==0 % determined manually for this case
      trluse=trluse(2:end); % remove first trial
    else
      disp('not sure what to throw out')
      keyboard
    end
  else
    disp('why number of trial mismatch?')
    keyboard
  end
end
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
clear cfgtr

N=raw_hpf.sampleinfo(end);
clear raw_hpf

if length(intersect(find(diff(raw_tac.trialinfo(:,3:4),[],2)) , find(raw_tac.trialinfo(:,4)~=-3)))>10
  if [sleep==1 && any(find(ii==[11 25 27 28]))] || [sleep==0 && any(find(ii==[10 17 25 30]))]
    % known cases where this happens
  else
    disp('trl misaligned?'); keyboard
  end
end
if length(intersect(find(diff(raw_aud.trialinfo(:,3:4),[],2)) , find(raw_aud.trialinfo(:,4)~=-3)))>10
  if [sleep==1 && any(find(ii==[11 25 27 28]))] || [sleep==0 && any(find(ii==[10 17 25 30]))]
    % known cases where this happens
  else
    disp('trl misaligned?'); keyboard
  end
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
bigwavefeats={'down' 'up' 'negmax' 'negmax_tp' 'posmax' 'posmax_tp' 'upstart' 'upend' 'uponset'};
% spindlefeats={'amplitude' 'duration' 'maximum'};
spindlefeats={'artifact' 'amplitude' 'duration'};

% 12-19: percentage of samples(time) during pre-stim (prestim) involved in this feature
% 22-29: percentage of samples(time) during post-stim (poststim) involved in this feature
% 32-39: percentage of samples(time) during pre-stim (prenar) involved in this feature; only relevant for featfull=0
% 42-49: percentage of samples(time) during post-stim (postnar) involved in this feature; only relevant for featfull=0
if featfull==1
  raw_tac.trialinfo(:,12:29)=0;
  raw_aud.trialinfo(:,12:29)=0;
  raw_nul.trialinfo(:,12:29)=0;
elseif featfull==0
  raw_tac.trialinfo(:,12:249)=0;
  raw_aud.trialinfo(:,12:249)=0;
  raw_nul.trialinfo(:,12:249)=0;
end
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
  if isfield(artfctdef.(artnames{aa}),'artifact') && ~isnan(trialcol)
    samp.(artnames{aa})=binarise_artifact_begendpoints(artfctdef.(artnames{aa}).artifact,N);
    for ll=1:size(raw_tac.sampleinfo,1)
      %           raw_tac.trialinfo(ll,trialcol)=any(samp.(artnames{aa})(raw_tac.sampleinfo(ll,1):raw_tac.sampleinfo(ll,2)));
      raw_tac.trialinfo(ll,trialcol)   =mean(samp.(artnames{aa})(raw_tac.sampleinfo(ll,1):raw_tac.sampleinfo(ll,1)+dsearchn(raw_tac.time{ll}',0)));
      raw_tac.trialinfo(ll,trialcol+10)=mean(samp.(artnames{aa})(raw_tac.sampleinfo(ll,1)+dsearchn(raw_tac.time{ll}',0):raw_tac.sampleinfo(ll,2)));
      if featfull==0
        %               keyboard % should '0' as defintion of prestim/poststim
        %               boundary be changed/dependent on ll??  No: we don't have
        %               trial types 'll' yet.
        raw_tac.trialinfo(ll,trialcol+20)=mean(samp.(artnames{aa})(dsearchn(raw_tac.time{ll}',0)-prenar*1000: dsearchn(raw_tac.time{ll}',0)));
        raw_tac.trialinfo(ll,trialcol+30)=mean(samp.(artnames{aa})(dsearchn(raw_tac.time{ll}',0):dsearchn(raw_tac.time{ll}',0)+prenar*1000));
      end
    end
    for ll=1:size(raw_aud.sampleinfo,1)
      %           raw_aud.trialinfo(ll,trialcol)=any(samp.(artnames{aa})(raw_aud.sampleinfo(ll,1):raw_aud.sampleinfo(ll,2)));
      raw_aud.trialinfo(ll,trialcol)   =mean(samp.(artnames{aa})(raw_aud.sampleinfo(ll,1):raw_aud.sampleinfo(ll,1)+dsearchn(raw_aud.time{ll}',0)));
      raw_aud.trialinfo(ll,trialcol+10)=mean(samp.(artnames{aa})(raw_aud.sampleinfo(ll,1)+dsearchn(raw_aud.time{ll}',0):raw_aud.sampleinfo(ll,2)));
      if featfull==0
        raw_aud.trialinfo(ll,trialcol+20)=mean(samp.(artnames{aa})(dsearchn(raw_aud.time{ll}',0)-prenar*1000: dsearchn(raw_aud.time{ll}',0)));
        raw_aud.trialinfo(ll,trialcol+30)=mean(samp.(artnames{aa})(dsearchn(raw_aud.time{ll}',0):dsearchn(raw_aud.time{ll}',0)+prenar*1000));
      end
    end
    for ll=1:size(raw_nul.sampleinfo,1)
      %           raw_nul.trialinfo(ll,trialcol)   =mean(samp.(artnames{aa})(raw_nul.sampleinfo(ll,1):raw_nul.sampleinfo(ll,1)+dsearchn(raw_nul.time{ll}',0)));
      raw_nul.trialinfo(ll,trialcol+10)=mean(samp.(artnames{aa})(raw_nul.sampleinfo(ll,1)+dsearchn(raw_nul.time{ll}',0):raw_nul.sampleinfo(ll,2)));
      raw_nul.trialinfo(ll,trialcol)=any(samp.(artnames{aa})(raw_nul.sampleinfo(ll,1):raw_nul.sampleinfo(ll,2)));
      if featfull==0
        raw_nul.trialinfo(ll,trialcol+20)=mean(samp.(artnames{aa})(dsearchn(raw_nul.time{ll}',0)-prenar*1000: dsearchn(raw_nul.time{ll}',0)));
        raw_nul.trialinfo(ll,trialcol+30)=mean(samp.(artnames{aa})(dsearchn(raw_nul.time{ll}',0):dsearchn(raw_nul.time{ll}',0)+prenar*1000));
      end
    end
    samp=rmfield(samp,artnames{aa}); % clearing as this can get large
    
    if featfull==0
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
        raw_tac.(artnames{aa})=repmat({[]},1,size(raw_tac.trialinfo,1));
        raw_aud.(artnames{aa})=repmat({[]},1,size(raw_aud.trialinfo,1));
        raw_nul.(artnames{aa})=repmat({[]},1,size(raw_nul.trialinfo,1));
        for kk=1:length(artfctdef.(artnames{aa}).negmax)
%           if ~isempty( find(raw_tac.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_tac.sampleinfo(:,2)))
          if ~isempty( find(raw_tac.sampleinfo(:,1)<artfctdef.(artnames{aa}).struct(kk).down & artfctdef.(artnames{aa}).struct(kk).upend<raw_tac.sampleinfo(:,2) & artfctdef.(artnames{aa}).struct(kk).posmax<raw_tac.sampleinfo(:,2)))
%             ind=find(raw_tac.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_tac.sampleinfo(:,2));
            ind=find(raw_tac.sampleinfo(:,1)<artfctdef.(artnames{aa}).struct(kk).down & artfctdef.(artnames{aa}).struct(kk).upend<raw_tac.sampleinfo(:,2) & artfctdef.(artnames{aa}).struct(kk).posmax<raw_tac.sampleinfo(:,2));

            % assign Kc to only one trial: added Oct 2018
            if length(ind)>1
              Kcdown=[artfctdef.(artnames{aa}).struct(kk).down-mean(raw_tac.sampleinfo(ind,1:2),2)]*diff(time{1}(1:2));
              Kcupend=[artfctdef.(artnames{aa}).struct(kk).upend-mean(raw_tac.sampleinfo(ind,1:2),2)]*diff(time{1}(1:2));              
              if length(ind)>1 && Kcdown(1)>1.2
                ind=ind(2:end);
                Kcdown=Kcdown(2:end);
                Kcupend=Kcupend(2:end);
              end
              if length(ind)>1 && Kcupend(end)<-0.7
                ind=ind(1:end-1);
                Kcdown=Kcdown(1:end-1);
                Kcupend=Kcupend(1:end-1);
              end
              if length(ind)==2 && Kcdown(1)>-0.4 && Kcupend(2)<0
                ind=ind(1);
                Kcdown=Kcdown(1);
                Kcupend=Kcupend(1);
              elseif length(ind)==2 && Kcdown(1)>.7 && Kcdown(2)<.2 % ignore first if it's far in post-stim and second is *during* stim onset
                ind=ind(2);
                Kcdown=Kcdown(2);
                Kcupend=Kcupend(2);
              elseif length(ind)==1
                % carry on
              else % carry on
%                 disp('help!')
%                 keyboard
              end
            end
            for nn=1:length(ind)
              %               keyboard % set boundary as 250ms after onset of first stim,
              %               which will be definded using stimcode -2:9
              if raw_tac.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_tac.sampleinfo(ind(nn),1)) < 0  % prestim
                if tacusepre(ind(nn))>9
                  error('too many BigWaves on this trial')
                end
                raw_tac.trialinfo(ind(nn),colampstartpre+tacusepre(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
                raw_tac.trialinfo(ind(nn),colnegstartpre+tacusepre(ind(nn)))=raw_tac.time{ind(nn)}(artfctdef.(artnames{aa}).negmax(kk)-raw_tac.sampleinfo(ind(nn),1));
                %                 keyboard % save out negslope as well
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
              % and save out whole structure/information about BigWave; might be more than one BigWave per trial
              if ~isempty(raw_tac.(artnames{aa}){ind(nn)})
                old_numKc=length(raw_tac.(artnames{aa}){ind(nn)});
                if isfield(raw_tac.(artnames{aa}){ind(nn)}(old_numKc),'upend') && raw_tac.(artnames{aa}){ind(nn)}(old_numKc).upend<0
                  raw_tac.(artnames{aa}){ind(nn)}(old_numKc)=[];
                  old_numKc=old_numKc-1;
                  new_numKc=1;
                elseif isfield(raw_tac.(artnames{aa}){ind(nn)}(old_numKc),'upend') && Kcdown(nn)<1
                  new_numKc=1;                  
                elseif isfield(raw_tac.(artnames{aa}){ind(nn)}(old_numKc),'upend')
                  new_numKc=0;
                  disp('not including later post-stim Kc')
                  continue
                end
              else
                old_numKc=0;
                new_numKc=1;
              end
              numKc=old_numKc+new_numKc;
              raw_tac.(artnames{aa}){ind(nn)}=[raw_tac.(artnames{aa}){ind(nn)}  artfctdef.(artnames{aa}).struct(kk)];
              for bg=1:length(bigwavefeats)
                raw_tac.(artnames{aa}){ind(nn)}(numKc).(bigwavefeats{bg})=raw_tac.time{ind(nn)}(raw_tac.(artnames{aa}){ind(nn)}(numKc).(bigwavefeats{bg})-raw_tac.sampleinfo(ind(nn),1));
              end
            end
          end
%           if ~isempty( find(raw_aud.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_aud.sampleinfo(:,2)))
%             ind=find(raw_aud.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_aud.sampleinfo(:,2));
          if ~isempty( find(raw_aud.sampleinfo(:,1)<artfctdef.(artnames{aa}).struct(kk).down & artfctdef.(artnames{aa}).struct(kk).upend<raw_aud.sampleinfo(:,2) & artfctdef.(artnames{aa}).struct(kk).posmax<raw_aud.sampleinfo(:,2)))
            ind=find(raw_aud.sampleinfo(:,1)<artfctdef.(artnames{aa}).struct(kk).down & artfctdef.(artnames{aa}).struct(kk).upend<raw_aud.sampleinfo(:,2) & artfctdef.(artnames{aa}).struct(kk).posmax<raw_aud.sampleinfo(:,2));
            % assign Kc to only one trial: added Oct 2018
            if length(ind)>1
              Kcdown=[artfctdef.(artnames{aa}).struct(kk).down-mean(raw_aud.sampleinfo(ind,1:2),2)]*diff(time{1}(1:2));
              Kcupend=[artfctdef.(artnames{aa}).struct(kk).upend-mean(raw_aud.sampleinfo(ind,1:2),2)]*diff(time{1}(1:2));
              if length(ind)>1 && Kcdown(1)>1.2
                ind=ind(2:end);
                Kcdown=Kcdown(2:end);
                Kcupend=Kcupend(2:end);
              end
              if length(ind)>1 && Kcupend(end)<-0.7
                ind=ind(1:end-1);
                Kcdown=Kcdown(1:end-1);
                Kcupend=Kcupend(1:end-1);
              end
              if length(ind)==2 && Kcdown(1)>-0.4 && Kcupend(2)<0  % keep post-stim (evoked) and ignore prestim if it ends before time=0
                ind=ind(1);
                Kcdown=Kcdown(1);
                Kcupend=Kcupend(1);
              elseif length(ind)==2 && Kcdown(1)>.7 && Kcdown(2)<.2 % ignore first if it's far in post-stim and second is *during* stim onset
                ind=ind(2);
                Kcdown=Kcdown(2);
                Kcupend=Kcupend(2);
              elseif length(ind)==1  % carry on
              else % carry on
%                 disp('help!')
%                 keyboard
              end
            end
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
              % and save out whole structure/information about BigWave; might be more than one BigWave per trial
              if ~isempty(raw_aud.(artnames{aa}){ind(nn)})
                old_numKc=length(raw_aud.(artnames{aa}){ind(nn)});
                if isfield(raw_aud.(artnames{aa}){ind(nn)}(old_numKc),'upend') && raw_aud.(artnames{aa}){ind(nn)}(old_numKc).upend<0
                  raw_aud.(artnames{aa}){ind(nn)}(old_numKc)=[];
                  old_numKc=old_numKc-1;
                  new_numKc=1;
                elseif isfield(raw_aud.(artnames{aa}){ind(nn)}(old_numKc),'upend') && Kcdown(nn)<1
                  new_numKc=1;                  
                elseif isfield(raw_aud.(artnames{aa}){ind(nn)}(old_numKc),'upend')
                  new_numKc=0;
                  disp('not including later post-stim Kc')
                  continue
                end
              else
                old_numKc=0;
                new_numKc=1;
              end
              numKc=old_numKc+new_numKc;
              raw_aud.(artnames{aa}){ind(nn)}=[raw_aud.(artnames{aa}){ind(nn)}  artfctdef.(artnames{aa}).struct(kk)];
              for bg=1:length(bigwavefeats)
                raw_aud.(artnames{aa}){ind(nn)}(numKc).(bigwavefeats{bg})=raw_aud.time{ind(nn)}(raw_aud.(artnames{aa}){ind(nn)}(numKc).(bigwavefeats{bg})-raw_aud.sampleinfo(ind(nn),1));
              end
            end % nn
          end
%           if ~isempty( find(raw_nul.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_nul.sampleinfo(:,2)))
%             ind=find(raw_nul.sampleinfo(:,1)<artfctdef.(artnames{aa}).negmax(kk) & artfctdef.(artnames{aa}).negmax(kk)<raw_nul.sampleinfo(:,2));
          if ~isempty( find(raw_nul.sampleinfo(:,1)<artfctdef.(artnames{aa}).struct(kk).down & artfctdef.(artnames{aa}).struct(kk).upend<raw_nul.sampleinfo(:,2) & artfctdef.(artnames{aa}).struct(kk).posmax<raw_nul.sampleinfo(:,2)))
            ind=find(raw_nul.sampleinfo(:,1)<artfctdef.(artnames{aa}).struct(kk).down & artfctdef.(artnames{aa}).struct(kk).upend<raw_nul.sampleinfo(:,2) & artfctdef.(artnames{aa}).struct(kk).posmax<raw_nul.sampleinfo(:,2));
            % try to assign Kc to only one trial: added Oct 2018
            if length(ind)>1
              Kcdown=[artfctdef.(artnames{aa}).struct(kk).down-mean(raw_nul.sampleinfo(ind,1:2),2)]*diff(time{1}(1:2));
              Kcupend=[artfctdef.(artnames{aa}).struct(kk).upend-mean(raw_nul.sampleinfo(ind,1:2),2)]*diff(time{1}(1:2));
              if length(ind)>1 && Kcdown(1)>1.2
                ind=ind(2:end);
                Kcdown=Kcdown(2:end);
                Kcupend=Kcupend(2:end);
              end
              if length(ind)>1 && Kcupend(end)<-0.7
                ind=ind(1:end-1);
                Kcdown=Kcdown(1:end-1);
                Kcupend=Kcupend(1:end-1);
              end
              if length(ind)==2 && Kcdown(1)>-0.4 && Kcupend(2)<0
                ind=ind(1);
                Kcdown=Kcdown(1);
                Kcupend=Kcupend(1);
              elseif length(ind)==2 && Kcdown(1)>.7 && Kcdown(2)<.2 % ignore first if it's far in post-stim and second is *during* stim onset
                ind=ind(2);
                Kcdown=Kcdown(2);
                Kcupend=Kcupend(2);
              elseif length(ind)==1
                % carry on
              else % carry on
%                 disp('help!')
%                 keyboard
              end
            end
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
              % and save out whole structure/information about BigWave; might be more than one BigWave per trial
              if ~isempty(raw_nul.(artnames{aa}){ind(nn)})
                old_numKc=length(raw_nul.(artnames{aa}){ind(nn)});
                if isfield(raw_nul.(artnames{aa}){ind(nn)}(old_numKc),'upend') && raw_nul.(artnames{aa}){ind(nn)}(old_numKc).upend<0
                  raw_nul.(artnames{aa}){ind(nn)}(old_numKc)=[];
                  old_numKc=old_numKc-1;
                  new_numKc=1;
                elseif isfield(raw_nul.(artnames{aa}){ind(nn)}(old_numKc),'upend') && Kcdown(nn)<1
                  new_numKc=1;                  
                elseif isfield(raw_nul.(artnames{aa}){ind(nn)}(old_numKc),'upend')
                  new_numKc=0;
                  disp('not including later post-stim Kc')
                  continue
                end
              else
                old_numKc=0;
                new_numKc=1;
              end
              numKc=old_numKc+new_numKc;
              raw_nul.(artnames{aa}){ind(nn)}=[raw_nul.(artnames{aa}){ind(nn)}  artfctdef.(artnames{aa}).struct(kk)];
              for bg=1:length(bigwavefeats)  %changing sample index to peri-stimulus time
                raw_nul.(artnames{aa}){ind(nn)}(numKc).(bigwavefeats{bg})=raw_nul.time{ind(nn)}(raw_nul.(artnames{aa}){ind(nn)}(numKc).(bigwavefeats{bg})-raw_nul.sampleinfo(ind(nn),1));
              end
            end
          end
        end % kk
      end % if artnames BigWaves
      tacusepre=zeros(size(raw_tac.trialinfo,1),1);
      audusepre=zeros(size(raw_aud.trialinfo,1),1);
      nulusepre=zeros(size(raw_nul.trialinfo,1),1);
      tacusepost=zeros(size(raw_tac.trialinfo,1),1);
      audusepost=zeros(size(raw_aud.trialinfo,1),1);
      nulusepost=zeros(size(raw_nul.trialinfo,1),1);
      if strcmp(artnames{aa},'sp_fast') || strcmp(artnames{aa},'sp_slow')
        % add spindle features here!
        if strcmp(artnames{aa},'sp_fast')
          colampstartpre=170;
          coldurstartpre=180;
          colampstartpost=190;
          coldurstartpost=200;
        elseif strcmp(artnames{aa},'sp_slow')
          colampstartpre=210;
          coldurstartpre=220;
          colampstartpost=230;
          coldurstartpost=240;
        end
        raw_tac.(artnames{aa})=repmat({[]},1,size(raw_tac.trialinfo,1));
        raw_aud.(artnames{aa})=repmat({[]},1,size(raw_aud.trialinfo,1));
        raw_nul.(artnames{aa})=repmat({[]},1,size(raw_nul.trialinfo,1));
        if isfield(artfctdef.(artnames{aa}),'artifact') && ~isfield(artfctdef.(artnames{aa}),'amplitude')
          error('rerun *spindles.m to get amplitude and duration')
        end
        for kk=1:length(artfctdef.(artnames{aa}).amplitude)
          artmid=round(mean(artfctdef.(artnames{aa}).artifact(kk,:))); % must round as must be integer sample location
          %           if ~isempty( find(raw_tac.sampleinfo(:,1)<artmid & artmid<raw_tac.sampleinfo(:,2)))
          %             ind=find(raw_tac.sampleinfo(:,1)<artmid & artmid<raw_tac.sampleinfo(:,2));
          if ~isempty(raw_tac.sampleinfo(:,1)<artfctdef.(artnames{aa}).artifact(kk,1) & artfctdef.(artnames{aa}).artifact(kk,2)<raw_tac.sampleinfo(:,2));
            ind=find(raw_tac.sampleinfo(:,1)<artfctdef.(artnames{aa}).artifact(kk,1) & artfctdef.(artnames{aa}).artifact(kk,2)<raw_tac.sampleinfo(:,2));
            for nn=1:length(ind)
              if raw_tac.time{ind(nn)}(artmid-raw_tac.sampleinfo(ind(nn),1)) < 0  % prestim
                if tacusepre(ind(nn))>9
                  error('too many spindles on this trial')
                end
                raw_tac.trialinfo(ind(nn),colampstartpre+tacusepre(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
                raw_tac.trialinfo(ind(nn),coldurstartpre+tacusepre(ind(nn)))=artfctdef.(artnames{aa}).duration(kk)/1000;
                tacusepre(ind(nn))=tacusepre(ind(nn))+1;
              else % poststim
                if tacusepost(ind(nn))>9
                  error('too many spindles on this trial')
                end
                raw_tac.trialinfo(ind(nn),colampstartpost+tacusepost(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
                raw_tac.trialinfo(ind(nn),coldurstartpost+tacusepost(ind(nn)))=artfctdef.(artnames{aa}).duration(kk)/1000;
                tacusepost(ind(nn))=tacusepost(ind(nn))+1;
              end
              % and save out whole structure/information about Spindle; might be more than one Spindle per trial
              if isfield(raw_tac.(artnames{aa}){ind(nn)},'endtime')
                numSp=length(raw_tac.(artnames{aa}){ind(nn)});
              else
                numSp=0;
              end
              if numSp && raw_tac.(artnames{aa}){ind(nn)}(numSp).endtime<-0.3
                raw_tac.(artnames{aa}){ind(nn)}(numSp)=[];
                numSp=numSp-1;
              else
                % nothing
              end
              numSp=numSp+1;
              for bg=1:length(spindlefeats)
                if bg==1
                  raw_tac.(artnames{aa}){ind(nn)}(numSp).starttime=raw_tac.time{ind(nn)}(artfctdef.(artnames{aa}).(spindlefeats{bg})(kk,1)-raw_tac.sampleinfo(ind(nn),1));
                  raw_tac.(artnames{aa}){ind(nn)}(numSp).endtime  =raw_tac.time{ind(nn)}(artfctdef.(artnames{aa}).(spindlefeats{bg})(kk,2)-raw_tac.sampleinfo(ind(nn),1));
                elseif bg==2
                  raw_tac.(artnames{aa}){ind(nn)}(numSp).(spindlefeats{bg})=artfctdef.(artnames{aa}).(spindlefeats{bg})(kk);
                elseif bg==3
                  raw_tac.(artnames{aa}){ind(nn)}(numSp).(spindlefeats{bg})=artfctdef.(artnames{aa}).(spindlefeats{bg})(kk)/1000;
                end
              end
            end %nn
          end % if
%           if ~isempty( find(raw_aud.sampleinfo(:,1)<artmid & artmid<raw_aud.sampleinfo(:,2)))
%             ind=find(raw_aud.sampleinfo(:,1)<artmid & artmid<raw_aud.sampleinfo(:,2));
          if ~isempty(raw_aud.sampleinfo(:,1)<artfctdef.(artnames{aa}).artifact(kk,1) & artfctdef.(artnames{aa}).artifact(kk,2)<raw_aud.sampleinfo(:,2));
            ind=find(raw_aud.sampleinfo(:,1)<artfctdef.(artnames{aa}).artifact(kk,1) & artfctdef.(artnames{aa}).artifact(kk,2)<raw_aud.sampleinfo(:,2));
            for nn=1:length(ind)
              if raw_aud.time{ind(nn)}(artmid-raw_aud.sampleinfo(ind(nn),1)) < 0  % prestim
                if audusepre(ind(nn))>9
                  error('too many spindles on this trial')
                end
                raw_aud.trialinfo(ind(nn),colampstartpre+audusepre(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
                raw_aud.trialinfo(ind(nn),coldurstartpre+audusepre(ind(nn)))=artfctdef.(artnames{aa}).duration(kk)/1000;
                audusepre(ind(nn))=audusepre(ind(nn))+1;
              else % poststim
                if audusepost(ind(nn))>9
                  error('too many spindles on this trial')
                end
                raw_aud.trialinfo(ind(nn),colampstartpost+audusepost(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
                raw_aud.trialinfo(ind(nn),coldurstartpost+audusepost(ind(nn)))=artfctdef.(artnames{aa}).duration(kk)/1000;
                audusepost(ind(nn))=audusepost(ind(nn))+1;
              end
              if isfield(raw_aud.(artnames{aa}){ind(nn)},'endtime')
                numSp=length(raw_aud.(artnames{aa}){ind(nn)});
              else
                numSp=0;
              end
              if numSp && raw_aud.(artnames{aa}){ind(nn)}(numSp).endtime<-0.3
                raw_aud.(artnames{aa}){ind(nn)}(numSp)=[];
                numSp=numSp-1;
              else
                % nothing
              end
              numSp=numSp+1;
              for bg=1:length(spindlefeats)
                if bg==1
                  raw_aud.(artnames{aa}){ind(nn)}(numSp).starttime=raw_aud.time{ind(nn)}(artfctdef.(artnames{aa}).(spindlefeats{bg})(kk,1)-raw_aud.sampleinfo(ind(nn),1));
                  raw_aud.(artnames{aa}){ind(nn)}(numSp).endtime  =raw_aud.time{ind(nn)}(artfctdef.(artnames{aa}).(spindlefeats{bg})(kk,2)-raw_aud.sampleinfo(ind(nn),1));
                elseif bg==2
                  raw_aud.(artnames{aa}){ind(nn)}(numSp).(spindlefeats{bg})=artfctdef.(artnames{aa}).(spindlefeats{bg})(kk);
                elseif bg==3
                  raw_aud.(artnames{aa}){ind(nn)}(numSp).(spindlefeats{bg})=artfctdef.(artnames{aa}).(spindlefeats{bg})(kk)/1000;
                end
              end
            end %nn
          end % if
%           if ~isempty( find(raw_nul.sampleinfo(:,1)<artmid & artmid<raw_nul.sampleinfo(:,2)))
%             ind=find(raw_nul.sampleinfo(:,1)<artmid & artmid<raw_nul.sampleinfo(:,2));
          if ~isempty(raw_nul.sampleinfo(:,1)<artfctdef.(artnames{aa}).artifact(kk,1) & artfctdef.(artnames{aa}).artifact(kk,2)<raw_nul.sampleinfo(:,2));
            ind=find(raw_nul.sampleinfo(:,1)<artfctdef.(artnames{aa}).artifact(kk,1) & artfctdef.(artnames{aa}).artifact(kk,2)<raw_nul.sampleinfo(:,2));
            for nn=1:length(ind)
              if raw_nul.time{ind(nn)}(artmid-raw_nul.sampleinfo(ind(nn),1)) < 0  % prestim
                if nulusepre(ind(nn))>9
                  error('too many spindles on this trial')
                end
                raw_nul.trialinfo(ind(nn),colampstartpre+nulusepre(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
                raw_nul.trialinfo(ind(nn),coldurstartpre+nulusepre(ind(nn)))=artfctdef.(artnames{aa}).duration(kk)/1000;
                nulusepre(ind(nn))=nulusepre(ind(nn))+1;
              else % poststim
                if nulusepost(ind(nn))>9
                  error('too many spindles on this trial')
                end
                raw_nul.trialinfo(ind(nn),colampstartpost+nulusepost(ind(nn)))=artfctdef.(artnames{aa}).amplitude(kk);
                raw_nul.trialinfo(ind(nn),coldurstartpost+nulusepost(ind(nn)))=artfctdef.(artnames{aa}).duration(kk)/1000;
                nulusepost(ind(nn))=nulusepost(ind(nn))+1;
              end
              if isfield(raw_nul.(artnames{aa}){ind(nn)},'endtime')
                numSp=length(raw_nul.(artnames{aa}){ind(nn)});
              else
                numSp=0;
              end
              if numSp && raw_nul.(artnames{aa}){ind(nn)}(numSp).endtime<-0.3
                raw_nul.(artnames{aa}){ind(nn)}(numSp)=[];
                numSp=numSp-1;
              else
                % nothing
              end
              numSp=numSp+1;
              for bg=1:length(spindlefeats)
                if bg==1
                  raw_nul.(artnames{aa}){ind(nn)}(numSp).starttime=raw_nul.time{ind(nn)}(artfctdef.(artnames{aa}).(spindlefeats{bg})(kk,1)-raw_nul.sampleinfo(ind(nn),1));
                  raw_nul.(artnames{aa}){ind(nn)}(numSp).endtime  =raw_nul.time{ind(nn)}(artfctdef.(artnames{aa}).(spindlefeats{bg})(kk,2)-raw_nul.sampleinfo(ind(nn),1));
                elseif bg==2
                  raw_nul.(artnames{aa}){ind(nn)}(numSp).(spindlefeats{bg})=artfctdef.(artnames{aa}).(spindlefeats{bg})(kk);
                elseif bg==3
                  raw_nul.(artnames{aa}){ind(nn)}(numSp).(spindlefeats{bg})=artfctdef.(artnames{aa}).(spindlefeats{bg})(kk)/1000;
                end
              end
            end % nn
          end % if
        end % kk
        
      end % if artnames Spindles
    end % featfull
    
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
if featfull==0
  raw_tac.delta={raw_tac.delta{find(keeptrials)}};
  raw_tac.kc={raw_tac.kc{find(keeptrials)}};
  raw_tac.sw={raw_tac.sw{find(keeptrials)}};
end


keeptrials=ones(1,length(raw_aud.time));
for ll=1:length(raw_aud.time)
  if isnan(raw_aud.time{ll}(1))
    keeptrials(ll)=0;
  end
end
cfg=[];
cfg.trials=find(keeptrials);
raw_aud=ft_redefinetrial(cfg,raw_aud);
if featfull==0
  raw_aud.delta={raw_aud.delta{find(keeptrials)}};
  raw_aud.kc={raw_aud.kc{find(keeptrials)}};
  raw_aud.sw={raw_aud.sw{find(keeptrials)}};
end

keeptrials=ones(1,length(raw_nul.time));
for ll=1:length(raw_nul.time)
  if isnan(raw_nul.time{ll}(1))
    keeptrials(ll)=0;
  end
end
cfg=[];
cfg.trials=find(keeptrials);
raw_nul=ft_redefinetrial(cfg,raw_nul);
if featfull==0
  raw_nul.delta={raw_nul.delta{find(keeptrials)}};
  raw_nul.kc={raw_nul.kc{find(keeptrials)}};
  raw_nul.sw={raw_nul.sw{find(keeptrials)}};
end

%   cfg=[];
%   cfg.hpfilter='yes';
%   cfg.hpfreq=1;
%   raw_tac_hp=ft_preprocessing(cfg,raw_tac);
%   raw_aud_hp=ft_preprocessing(cfg,raw_aud);
%   raw_nul_hp=ft_preprocessing(cfg,raw_nul);

% final trial/channel rejection
%   cfg=[];
%   cfg.method='summary';
%   % raw_all_rej=ft_rejectvisual(cfg,raw_all);
% %   raw_tac_rej=ft_rejectvisual(cfg,raw_tac_hp);
% %   raw_aud_rej=ft_rejectvisual(cfg,raw_aud_hp);
% %   raw_nul_rej=ft_rejectvisual(cfg,raw_nul_hp);
%   raw_tac_rej=ft_rejectvisual(cfg,raw_tac);
%   raw_aud_rej=ft_rejectvisual(cfg,raw_aud);
%   raw_nul_rej=ft_rejectvisual(cfg,raw_nul);


if 0 % not needed now with muscle and EOG identification ?
  cfg=[];
  cfg.method='channel';
  cfg.alim=100;
  raw_tac_rej=ft_rejectvisual(cfg,raw_tac);
  
  tacart=raw_tac_rej.cfg.artfctdef;
  
  cfg=[];
  cfg.artfctdef=tacart;
  raw_aud_rej=ft_rejectartifact(cfg,raw_aud);
  
  cfg=[];
  cfg.method='channel';
  cfg.alim=100;
  raw_aud_rej=ft_rejectvisual(cfg,raw_aud_rej);
  
  audart=raw_aud_rej.cfg.artfctdef;
  
  % if auditory was thrown out here that tactile still in, then it will later
  % be thrown out in eeg_legomagic_trialSelection1
  
  raw_nul_rej=ft_rejectvisual(cfg,raw_nul);
  nulart=raw_nul_rej.cfg.artfctdef;
end

if saveflag
  % save(['raw_each_rej_' sub{ii}],'raw_*rej'); % with 1Hz highpass after epoching
  %     save(['raw1_each_rej_' sub{ii}],'raw_*rej','*art'); % with 0.2Hz highpass before epoching
  %     save(['raw2_each_rej_' sub{ii} '_sleep' num2str(sleep)],'raw_*','artfctdef','-v7.3'); % with 0.2Hz highpass before epoching
  save(['raw2_each_rej_' sub{ii} '_sleep' num2str(sleep)],'raw_*','artflag','-v7.3'); % with 0.2Hz highpass before epoching
end

% clear raw_aud raw_nul raw_tac

return

%%
%     musclesample=zeros(1,raw_hpf.sampleinfo(end));
%     if isfield(artfctdef,'muscle')
%     for ll=1:size(artfctdef.muscle.artifact,1)
%       musclesample(artfctdef.muscle.artifact(ll,1):artfctdef.muscle.artifact(ll,2))=1;
%     end
%     end
%     for ll=1:size(raw_tac.sampleinfo,1)
%       raw_tac.trialinfo(ll,12)=any(musclesample(raw_tac.sampleinfo(ll,1):raw_tac.sampleinfo(ll,2)));
%     end
%     for ll=1:size(raw_aud.sampleinfo,1)
%       raw_aud.trialinfo(ll,12)=any(musclesample(raw_aud.sampleinfo(ll,1):raw_aud.sampleinfo(ll,2)));
%     end
%     for ll=1:size(raw_nul.sampleinfo,1)
%       raw_nul.trialinfo(ll,12)=any(musclesample(raw_nul.sampleinfo(ll,1):raw_nul.sampleinfo(ll,2)));
%     end
%     clear musclesample
%
%     eogsample=zeros(1,raw_hpf.sampleinfo(end));
%     if isfield(artfctdef,'eog')
%     for ll=1:size(artfctdef.eog.artifact,1)
%       eogsample(artfctdef.eog.artifact(ll,1):artfctdef.eog.artifact(ll,2))=1;
%     end
%     end
%     for ll=1:size(raw_tac.sampleinfo,1)
%       raw_tac.trialinfo(ll,13)=any(eogsample(raw_tac.sampleinfo(ll,1):raw_tac.sampleinfo(ll,2)));
%     end
%     for ll=1:size(raw_aud.sampleinfo,1)
%       raw_aud.trialinfo(ll,13)=any(eogsample(raw_aud.sampleinfo(ll,1):raw_aud.sampleinfo(ll,2)));
%     end
%     for ll=1:size(raw_nul.sampleinfo,1)
%       raw_nul.trialinfo(ll,13)=any(eogsample(raw_nul.sampleinfo(ll,1):raw_nul.sampleinfo(ll,2)));
%     end
%     clear eogsample
%
%     kcsample=zeros(1,raw_hpf.sampleinfo(end));
%     if isfield(artfctdef,'kc')
%       for ll=1:size(artfctdef.kc.artifact,1)
%         kcsample(artfctdef.kc.artifact(ll,1):artfctdef.kc.artifact(ll,2))=1;
%       end
%     end
%     for ll=1:size(raw_tac.sampleinfo,1)
%       raw_tac.trialinfo(ll,14)=any(kcsample(raw_tac.sampleinfo(ll,1):raw_tac.sampleinfo(ll,2)));
%     end
%     for ll=1:size(raw_aud.sampleinfo,1)
%       raw_aud.trialinfo(ll,14)=any(kcsample(raw_aud.sampleinfo(ll,1):raw_aud.sampleinfo(ll,2)));
%     end
%     for ll=1:size(raw_nul.sampleinfo,1)
%       raw_nul.trialinfo(ll,14)=any(kcsample(raw_nul.sampleinfo(ll,1):raw_nul.sampleinfo(ll,2)));
%     end
%     clear kcsample
%
%     deltasample=zeros(1,raw_hpf.sampleinfo(end));
%     if isfield(artfctdef,'delta')
%     for ll=1:size(artfctdef.delta.artifact,1)
%       deltasample(artfctdef.delta.artifact(ll,1):artfctdef.delta.artifact(ll,2))=1;
%     end
%     end
%     for ll=1:size(raw_tac.sampleinfo,1)
%       raw_tac.trialinfo(ll,15)=any(deltasample(raw_tac.sampleinfo(ll,1):raw_tac.sampleinfo(ll,2)));
%     end
%     for ll=1:size(raw_aud.sampleinfo,1)
%       raw_aud.trialinfo(ll,15)=any(deltasample(raw_aud.sampleinfo(ll,1):raw_aud.sampleinfo(ll,2)));
%     end
%     for ll=1:size(raw_nul.sampleinfo,1)
%       raw_nul.trialinfo(ll,15)=any(deltasample(raw_nul.sampleinfo(ll,1):raw_nul.sampleinfo(ll,2)));
%     end
%     clear deltasample
%
%     swsample=zeros(1,raw_hpf.sampleinfo(end));
%     if isfield(artfctdef,'sw')
%     for ll=1:size(artfctdef.sw.artifact,1)
%       swsample(artfctdef.sw.artifact(ll,1):artfctdef.sw.artifact(ll,2))=1;
%     end
%     end
%     for ll=1:size(raw_tac.sampleinfo,1)
%       raw_tac.trialinfo(ll,16)=any(swsample(raw_tac.sampleinfo(ll,1):raw_tac.sampleinfo(ll,2)));
%     end
%     for ll=1:size(raw_aud.sampleinfo,1)
%       raw_aud.trialinfo(ll,16)=any(swsample(raw_aud.sampleinfo(ll,1):raw_aud.sampleinfo(ll,2)));
%     end
%     for ll=1:size(raw_nul.sampleinfo,1)
%       raw_nul.trialinfo(ll,16)=any(swsample(raw_nul.sampleinfo(ll,1):raw_nul.sampleinfo(ll,2)));
%     end
%     clear swsample
%
%     sp_fastsample=zeros(1,raw_hpf.sampleinfo(end));
%     if isfield(artfctdef,'sp_fast')
%     for ll=1:size(artfctdef.sp_fast.artifact,1)
%       sp_fastsample(artfctdef.sp_fast.artifact(ll,1):artfctdef.sp_fast.artifact(ll,2))=1;
%     end
%     end
%     for ll=1:size(raw_tac.sampleinfo,1)
%       raw_tac.trialinfo(ll,17)=any(sp_fastsample(raw_tac.sampleinfo(ll,1):raw_tac.sampleinfo(ll,2)));
%     end
%     for ll=1:size(raw_aud.sampleinfo,1)
%       raw_aud.trialinfo(ll,17)=any(sp_fastsample(raw_aud.sampleinfo(ll,1):raw_aud.sampleinfo(ll,2)));
%     end
%     for ll=1:size(raw_nul.sampleinfo,1)
%       raw_nul.trialinfo(ll,17)=any(sp_fastsample(raw_nul.sampleinfo(ll,1):raw_nul.sampleinfo(ll,2)));
%     end
%     clear sp_fastsample
%
%     sp_slowsample=zeros(1,raw_hpf.sampleinfo(end));
%     if isfield(artfctdef,'sp_slow')
%     for ll=1:size(artfctdef.sp_slow.artifact,1)
%       sp_slowsample(artfctdef.sp_slow.artifact(ll,1):artfctdef.sp_slow.artifact(ll,2))=1;
%     end
%     end
%     for ll=1:size(raw_tac.sampleinfo,1)
%       raw_tac.trialinfo(ll,18)=any(sp_slowsample(raw_tac.sampleinfo(ll,1):raw_tac.sampleinfo(ll,2)));
%     end
%     for ll=1:size(raw_aud.sampleinfo,1)
%       raw_aud.trialinfo(ll,18)=any(sp_slowsample(raw_aud.sampleinfo(ll,1):raw_aud.sampleinfo(ll,2)));
%     end
%     for ll=1:size(raw_nul.sampleinfo,1)
%       raw_nul.trialinfo(ll,18)=any(sp_slowsample(raw_nul.sampleinfo(ll,1):raw_nul.sampleinfo(ll,2)));
%     end
%     clear sp_slowsample

%%

% cfg=[];
% cfg.reref='yes';
% cfg.refchannel='all';
% raw_all_ref=ft_preprocessing(cfg,raw_all);
%
% cfg=[];
% cfg.reref='yes';
% cfg.refchannel='all';
% raw_all_deref=ft_preprocessing(cfg,raw_all_demean);
%
% figure;cfg=[];cfg.xlim=[30 34];cfg.channel='Fz';cfg.ylim=[-10 40];ft_singleplotER(cfg,raw_all);
% figure;cfg=[];cfg.xlim=[30 34];cfg.channel='Fz';cfg.ylim=[-10 40];ft_singleplotER(cfg,raw_all_demean);
%
% figure;cfg=[];cfg.xlim=[30 34];cfg.channel='Fz';cfg.ylim=[-10 40];ft_singleplotER(cfg,raw_all_ref);
% figure;cfg=[];cfg.xlim=[30 34];cfg.channel='Fz';cfg.ylim=[-10 40];ft_singleplotER(cfg,raw_all_deref);

%%

% this is incorporated into eeg_legomagic_trialSelection. Use that instead
%
% % ERP filtering
%
% load(['raw_each_rej_' sub{ii}],'raw_*');
%
% if 0
%   cfg=[];
%   cfg.bpfilter='yes';
%   cfg.bpfreq=[1 40];
%   cfg.demean='yes';
%   cfg.baselinewindow=[-.6 -.1];
%   cfg.channel= {'all', '-ECG'}
%   data_tac_filt=ft_preprocessing(cfg,raw_tac_rej);
%   data_aud_filt=ft_preprocessing(cfg,raw_aud_rej);
%   data_nul_filt=ft_preprocessing(cfg,raw_nul_rej);
%
%   % % Solve this: what should be reference channel for scalp ERP? fixme
%   cfg=[];
%   cfg.reref='yes';
%   cfg.refchannel='all'; % necessary if doing source localisation
%   cfg.refchannel={'FT9', 'FT10'}; % sort of like linked mastoids?
%   data_tac_filt_ref=ft_preprocessing(cfg,data_tac_filt);
%   data_aud_filt_ref=ft_preprocessing(cfg,data_aud_filt);
%   data_nul_filt_ref=ft_preprocessing(cfg,data_nul_filt);
% else
%   % % Solve this: what should be reference channel for scalp ERP? fixme
%   cfg=[];
%   cfg.reref='yes';
%   cfg.refchannel='all'; % necessary if doing source localisation
% %   cfg.refchannel={'FT9', 'FT10'}; % sort of like linked mastoids?
%   data_tac_ref=ft_preprocessing(cfg,raw_tac_rej);
%   data_aud_ref=ft_preprocessing(cfg,raw_aud_rej);
%   data_nul_ref=ft_preprocessing(cfg,raw_nul_rej);
%
%   cfg=[];
%   cfg.bpfilter='yes';
%   cfg.bpfreq=[1 40];
%   cfg.demean='yes';
%   cfg.baselinewindow=[-.6 -.1];
%   cfg.channel= {'all', '-ECG'}
%   data_tac_ref_filt=ft_preprocessing(cfg,data_tac_ref);
%   data_aud_ref_filt=ft_preprocessing(cfg,data_aud_ref);
%   data_nul_ref_filt=ft_preprocessing(cfg,data_nul_ref);
% end
%
%
% clear raw_*
%
% if ii>7
%   soalist=[1 3 4 5 6 7 9];
% else
%   soalist=[3 4 5 6 7];
% end
% % include only trials during awake segments (in case participant went into N1 during sitting up portion)
% for ll=1:length(soalist)
%   for tt=3:6 % different thresholds of inclusion of SOA variance around desired SOA (see concat_stim_files.m)
%     cfg=[];
%     cfg.vartrllength=2;
%     cfg.keeptrials='yes';
%     cfg.trials=[data_tac_filt_ref.trialinfo(:,tt)==soalist(ll) & data_tac_filt_ref.trialinfo(:,2)==0 & ~isnan(data_tac_filt_ref.trialinfo(:,10))];
%     tlock_tac_s0{soalist(ll),tt-2}=ft_timelockanalysis(cfg,data_tac_filt_ref); % only 'awake' state
%     %   cfg.trials=data_tac_filt_ref.trialinfo(:,3)==soalist(ll);
%     %   tlock_tac_sall{soalist(ll)}=ft_timelockanalysis(cfg,data_tac_filt_ref); % all trials, irrepective if fell asleep
%   end
% end
% for ll=1:length(soalist)
%   for tt=3:6 % different thresholds of inclusion of SOA variance around desired SOA (see concat_stim_files.m)
%     cfg=[];
%     cfg.vartrllength=2;
%     cfg.keeptrials='yes';
%     cfg.trials=[data_aud_filt_ref.trialinfo(:,tt)==soalist(ll) & data_aud_filt.trialinfo(:,2)==0 & ~isnan(data_aud_filt_ref.trialinfo(:,10))];
%     tlock_aud_s0{soalist(ll),tt-2}=ft_timelockanalysis(cfg,data_aud_filt_ref);
%     %   cfg.trials=data_aud_filt_ref.trialinfo(:,3)==soalist(ll);
%     %   tlock_aud_sall{soalist(ll)}=ft_timelockanalysis(cfg,data_aud_filt_ref);
%   end
% end
% cfg=[]; % tac alone 38 (-2)
% cfg.vartrllength=2;
% cfg.keeptrials='yes';
% cfg.trials=[data_tac_filt_ref.trialinfo(:,3)==-2 & data_tac_filt_ref.trialinfo(:,2)==0];
% tlock_tac_s0{10,1}=ft_timelockanalysis(cfg,data_tac_filt_ref);
% % cfg.trials=data_tac_filt_ref.trialinfo(:,3)==-2;
% % tlock_tac_sall{10}=ft_timelockanalysis(cfg,data_tac_filt_ref);
% cfg=[]; % aud alone 39 (-1)
% cfg.vartrllength=2;
% cfg.keeptrials='yes';
% cfg.trials=[data_aud_filt_ref.trialinfo(:,3)==-1 & data_aud_filt_ref.trialinfo(:,2)==0];
% tlock_aud_s0{10,1}=ft_timelockanalysis(cfg,data_aud_filt_ref);
% % cfg.trials=data_aud_filt_ref.trialinfo(:,3)==-1;
% % tlock_aud_sall{10}=ft_timelockanalysis(cfg,data_aud_filt_ref);
% cfg=[];
% cfg.vartrllength=2;
% cfg.keeptrials='yes';
% cfg.trials=data_nul_filt_ref.trialinfo(:,2)==0;
% tlock_nul_s0{10,1}=ft_timelockanalysis(cfg,data_nul_filt_ref);
% % cfg.trials='all';
% % tlock_nul_sall{10}=ft_timelockanalysis(cfg,data_nul_filt_ref);
%
% clear data*
%
% for ll=[soalist 10] % this loop is to remove trials with NaN in (e.g. if artifact rejection cut through middle of a trial)
%   for tt=1:4
%     if ll==10 && tt>1
%       continue
%     else
%     cfg=[];
%     cfg.keeptrials='yes';
%     cfg.trials=find(sum(squeeze(isnan(tlock_tac_s0{ll,tt}.trial(:,10,:)))')<.05*size(tlock_tac_s0{ll,tt}.trial,3));
%     tlock_tac_s0{ll,tt}=ft_timelockanalysis(cfg,tlock_tac_s0{ll,tt});
%     %   cfg.trials=find(sum(squeeze(isnan(tlock_tac_sall{ll}.trial(:,10,:)))')<.05*size(tlock_tac_sall{ll}.trial,3));
%     %   tlock_tac_sall{ll}=ft_timelockanalysis(cfg,tlock_tac_sall{ll});
%     cfg=[];
%     cfg.keeptrials='yes';
%     cfg.trials=find(sum(squeeze(isnan(tlock_aud_s0{ll,tt}.trial(:,10,:)))')<.05*size(tlock_aud_s0{ll,tt}.trial,3));
%     tlock_aud_s0{ll,tt}=ft_timelockanalysis(cfg,tlock_aud_s0{ll,tt});
%     %   cfg.trials=find(sum(squeeze(isnan(tlock_aud_sall{ll}.trial(:,10,:)))')<.05*size(tlock_aud_sall{ll}.trial,3));
%     %   tlock_aud_sall{ll}=ft_timelockanalysis(cfg,tlock_aud_sall{ll});
%     end
%     if ll==10 && tt==1
%       cfg=[];
%       cfg.keeptrials='yes';
%       cfg.trials=find(sum(squeeze(isnan(tlock_nul_s0{ll}.trial(:,10,:)))')<.05*size(tlock_nul_s0{ll}.trial,3));
%       tlock_nul_s0{ll}=ft_timelockanalysis(cfg,tlock_nul_s0{ll});
%       %     cfg.trials=find(sum(squeeze(isnan(tlock_nul_sall{ll}.trial(:,10,:)))')<.05*size(tlock_nul_sall{ll}.trial,3));
%       %     tlock_nul_sall{ll}=ft_timelockanalysis(cfg,tlock_nul_sall{ll});
%     end
%   end
% end
%
% % % collapsing over all SOA conditions
% % clear keeptrial
% % for ll=1:size(data_tac_filt_ref.trial,2),
% %   keeptrial(ll)=~any(any(isnan(data_tac_filt_ref.trial{ll})));
% % end
% % cfg=[];
% % cfg.vartrllength=2;
% % cfg.trials=find(keeptrial);
% % cfg.keeptrials='yes';
% % tlock_tac{11}=ft_timelockanalysis(cfg,data_tac_filt_ref);
% %
% % % collapsing over all SOA conditions
% % clear keeptrial
% % for ll=1:size(data_aud_filt_ref.trial,2),
% %   keeptrial(ll)=~any(any(isnan(data_aud_filt_ref.trial{ll})));
% % end
% % cfg=[];
% % cfg.vartrllength=2;
% % cfg.trials=find(keeptrial);
% % cfg.keeptrials='yes';
% % tlock_aud{11}=ft_timelockanalysis(cfg,data_aud_filt_ref);
%
%
% if 0
%   cfg=[];
%   cfg.interactive='yes';
%   cfg.layout='EEG1010.lay';
%   ft_multiplotER(cfg,tlock_tac_s0{10,1});
%   %   ft_multiplotER(cfg,tlock_tac_sall{10});
%   ft_multiplotER(cfg,tlock_aud_s0{10,1});
%   %   ft_multiplotER(cfg,tlock_aud_sall{10});
%   ft_multiplotER(cfg,tlock_nul_s0{10,1});
%   %   ft_multiplotER(cfg,tlock_nul_sall{10});
%
%   % ft_multiplotER(cfg,tlock_tac{11});
%   % ft_multiplotER(cfg,tlock_aud{11});
%
%   ft_multiplotER(cfg,tlock_tac_s0{5,1});
%   ft_multiplotER(cfg,tlock_tac_s0{5,4});
%   ft_multiplotER(cfg,tlock_aud_s0{5,1});
%   ft_multiplotER(cfg,tlock_aud_s0{5,4});
% end
%
% %  contrasting conditions
%
% plotflag=1;
% load([ddir sub{ii} '_audtac.mat']);
% soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];
% chanuse=match_str(tlock_tac_s0{3}.label,'Fz');
%
% % determine trial numbers per condition
% numuni=min([size(tlock_tac_s0{10,1}.trialinfo,1),size(tlock_aud_s0{10,1}.trialinfo,1),size(tlock_nul_s0{10,1}.trialinfo,1)]);
% for ll=soalist
%   for tt=1:4 % limits of plus/minus:  [no-limit .01 .005 .002]
%     trcom=intersect(tlock_tac_s0{ll,tt}.trialinfo(:,11),tlock_aud_s0{ll,tt}.trialinfo(:,11));
%     nummul(ll,tt)=length(trcom);
%     numcond(ll,tt)=min(nummul(ll,tt),numuni);
%   end
% end
% numcondall=min(numcond(soalist,:)); % min over all conditions for a given SOA variance threshold
%
% fsample=1/diff(tlock_tac_s0{ll,tt}.time(1:2)); % sampling rate
%
% try close(10);end
% try close(11);end
% try close(19);end
% try close(20);end
% % create sum of unisensory, specific to the exact SOAs and actual jitters, for each multisensory condition
% for ll=soalist
%   for tt=1:4
%
%     % make sure the multisensory auditory-locked and ms tactile-locked have same trials; sometimes not exactly!
%     trcom=intersect(tlock_tac_s0{ll,tt}.trialinfo(:,11),tlock_aud_s0{ll,tt}.trialinfo(:,11));
%     [~,tbb,~]=intersect(tlock_tac_s0{ll,tt}.trialinfo(:,11),trcom);
%     [~,abb,~]=intersect(tlock_aud_s0{ll,tt}.trialinfo(:,11),trcom);
%
%     mrandtr=Shuffle(1:length(trcom)); % random trials to keep
% %     mrandtr=Shuffle(1:size(tlock_tac_s0{ll,tt}.trialinfo,1)); % random trials to keep
%     trandtr=Shuffle(1:size(tlock_tac_s0{10}.trialinfo,1)); % random trials to keep
%     arandtr=Shuffle(1:size(tlock_aud_s0{10}.trialinfo,1)); % random trials to keep
%     nrandtr=Shuffle(1:size(tlock_nul_s0{10}.trialinfo,1)); % random trials to keep
%     % Decide later: ???
%     if 1
%       mtrialkept{ll,tt}=sort(mrandtr(1:numcond(ll,tt))); % only keep the appropriate number
%       ttrialkept{ll,tt}=sort(trandtr(1:numcond(ll,tt))); % only keep the appropriate number
%       atrialkept{ll,tt}=sort(arandtr(1:numcond(ll,tt))); % only keep the appropriate number
%       ntrialkept{ll,tt}=sort(nrandtr(1:numcond(ll,tt))); % only keep the appropriate number
%     else
%       mtrialkept{ll,tt}=sort(mrandtr(1:numcondall(tt))); % only keep the appropriate number
%       ttrialkept{ll,tt}=sort(trandtr(1:numcondall(tt))); % only keep the appropriate number
%       atrialkept{ll,tt}=sort(arandtr(1:numcondall(tt))); % only keep the appropriate number
%       ntrialkept{ll,tt}=sort(nrandtr(1:numcondall(tt))); % only keep the appropriate number
%     end
%
%     cfg=[];
%     % I have thought it through, and no minus sign needed here.
%     %     cfg.offset=round(fsample*(tlock_tac_s0{ll,tt}.trialinfo(trialkept,10)-soades(ll))); % offset is in samples not seconds
%     cfg.offset=round(fsample*(tlock_tac_s0{ll,tt}.trialinfo(tbb(mtrialkept{ll,tt}),10)-soades(ll) -soades(ll))); % offset is in samples not seconds; -2*soades required for tactile not auditory
%     cfg.trials=ttrialkept{ll,tt};
%     tlock_tac_s0{ll+40,tt}=ft_redefinetrial(cfg,tlock_tac_s0{10}); % Talone with shift in time matching SOA jitter
%     cfg=[];
%     cfg.vartrllength=2;
%     tlock_tac_s0tlock{ll+40,tt}=ft_timelockanalysis(cfg,tlock_tac_s0{ll+40,tt}); % Talone with shift in time matching SOA jitter
%
%     cfg=[];
%     % I have thought it through, and no minus sign needed here.
%     %     cfg.offset=round(fsample*(tlock_aud_s0{ll,tt}.trialinfo(trialkept,10)-soades(ll))); % tlock_aud_s0{ll,tt}.trialinfo(trialkept,10) is identical to tlock_tac_s0{ll,tt}.trialinfo(trialkept,10)
%     cfg.offset=round(fsample*(tlock_aud_s0{ll,tt}.trialinfo(abb(mtrialkept{ll,tt}),10)-soades(ll) +soades(ll))); % tlock_aud_s0{ll,tt}.trialinfo(trialkept,10) is identical to tlock_tac_s0{ll,tt}.trialinfo(trialkept,10)
%     cfg.trials=atrialkept{ll,tt};
%     tlock_aud_s0{ll+40,tt}=ft_redefinetrial(cfg,tlock_aud_s0{10}); % Aalone with shift in time matching SOA jitter
%     cfg=[];
%     cfg.vartrllength=2;
%     tlock_aud_s0tlock{ll+40,tt}=ft_timelockanalysis(cfg,tlock_aud_s0{ll+40,tt}); % Aalone with shift in time matching SOA jitter
%
%     cfg=[];
%     cfg.trials=ntrialkept{ll,tt};
%     tlock_nul_s0{ll+40,tt}=ft_redefinetrial(cfg,tlock_nul_s0{10});
%
%     % create sum of unisensory conditions
%     cfg=[];
%     cfg.operation='add';
%     cfg.parameter='avg';
%     % the call to ft_timelockanalysis is b/c .avg field not present after ft_redefinetrial
%     tlock_tacPaud_s0{ll,tt}=ft_math(cfg,tlock_tac_s0{10},tlock_aud_s0tlock{ll+40,tt}); % keeping tac centred but jitter aud (to match tlock_tac_s0{X,tt})
%     tlock_audPtac_s0{ll,tt}=ft_math(cfg,tlock_aud_s0{10},tlock_tac_s0tlock{ll+40,tt}); % keeping aud centred but jitter tac (to match tlock_aud_s0{X,tt})
%
%     %     timets0=dsearchn(tlock_tac_s0{10}.time',tlock_tacPaud_s0{ll,tt}.time(1)):dsearchn(tlock_tac_s0{10}.time',tlock_tacPaud_s0{ll,tt}.time(end));
%     %     timeas0=dsearchn(tlock_aud_s0{10}.time',tlock_audPtac_s0{ll,tt}.time(1)):dsearchn(tlock_aud_s0{10}.time',tlock_audPtac_s0{ll,tt}.time(end));
%
%     if plotflag
%       figure(19); % Tac-alone (black) and Aud-alone (red) with their subset-jittered + shifted(blue and magenta)
%       subplot(9,4,(ll-1)*4+tt);
%       plot(tlock_tac_s0{10}.time,tlock_tac_s0{10}.avg(chanuse,:),'k');hold on;
%       plot(tlock_tac_s0{ll+40,tt}.time,tlock_tac_s0tlock{ll+40,tt}.avg(chanuse,:),'b');hold on;
%       plot(tlock_aud_s0{10}.time,tlock_aud_s0{10}.avg(chanuse,:),'r');hold on;
%       plot(tlock_aud_s0{ll+40,tt}.time,tlock_aud_s0tlock{ll+40,tt}.avg(chanuse,:),'m');hold on;
%       axis([-1 2.5 -10 10])
%
%       figure(20); % subset-jitted t-alone and a-alone (blue and magenta), as well as tacPaud and audPtac (green and red)
%       subplot(9,4,(ll-1)*4+tt);
%       plot(tlock_tac_s0{ll+40,tt}.time,tlock_tac_s0tlock{ll+40,tt}.avg(chanuse,:),'b');hold on;
%       plot(tlock_aud_s0{ll+40,tt}.time,tlock_aud_s0tlock{ll+40,tt}.avg(chanuse,:),'m');hold on;
%       plot(tlock_tacPaud_s0{ll,tt}.time,tlock_tacPaud_s0{ll,tt}.avg(chanuse,:),'g');hold on
%       plot(tlock_audPtac_s0{ll,tt}.time,tlock_audPtac_s0{ll,tt}.avg(chanuse,:),'r');hold on
%       axis([-1 2.5 -10 10])
%     end
%
% %     % create sum of multisensory plus nul conditions
% %     cfg=[];
% %     cfg.operation='add';
% %     cfg.parameter='avg';
% %     tlock_tac_s0{ll+20,tt}=ft_math(cfg,tlock_tac_s0{ll,tt},ft_timelockanalysis([],tlock_nul_s0{ll+40,tt})); % TA+N
% %     %   tlock_tac_sall{ll+20,tt}=ft_math(cfg,tlock_tac_sall{ll,tt},tlock_nul_sall{10}); % TA+N
% %     tlock_aud_s0{ll+20,tt}=ft_math(cfg,tlock_aud_s0{ll,tt},ft_timelockanalysis([],tlock_nul_s0{ll+40,tt})); % AT+N
% %     %   tlock_aud_sall{ll+20,tt}=ft_math(cfg,tlock_aud_sall{ll,tt},tlock_nul_sall{10}); % AT+N
%
%     % create sum of multisensory plus nul conditions
%     cfg=[];
%     cfg.operation='add';
%     cfg.parameter='avg';
%     tlock_tacMSpN_s0{ll,tt}=ft_math(cfg,tlock_tac_s0{ll,tt},ft_timelockanalysis([],tlock_nul_s0{ll+40,tt})); % TA+N
%     tlock_audMSpN_s0{ll,tt}=ft_math(cfg,tlock_aud_s0{ll,tt},ft_timelockanalysis([],tlock_nul_s0{ll+40,tt})); % AT+N
%
%     %   cfg=[];
%     %   cfg.operation='subtract';
%     %   cfg.parameter='avg';
%     % %   tlock_tac_s0{ll+10}=ft_math(cfg,ft_timelockanalysis([],tlock_tac_s0{ll}),ft_timelockanalysis([],tlock_tac_s0{10})); % TA-T = 'aud' isolated from multisensory
%     % %   tlock_tac_sall{ll+10}=ft_math(cfg,ft_timelockanalysis([],tlock_tac_sall{ll}),ft_timelockanalysis([],tlock_tac_sall{10})); % TA-T = 'aud' isolated from multisensory
%     % %   tlock_aud_s0{ll+10}=ft_math(cfg,ft_timelockanalysis([],tlock_aud_s0{ll}),ft_timelockanalysis([],tlock_aud_s0{10})); % AT-A = 'tac' isolated from mulitsensory
%     % %   tlock_aud_sall{ll+10}=ft_math(cfg,ft_timelockanalysis([],tlock_aud_sall{ll}),ft_timelockanalysis([],tlock_aud_sall{10})); % AT-A = 'tac' isolated from mulitsensory
%     %   tlock_tac_s0{ll+10}=ft_math(cfg,tlock_tac_s0{ll},tlock_tac_s0{10}); % TA-T = 'aud' isolated from multisensory
%     %   tlock_tac_sall{ll+10}=ft_math(cfg,tlock_tac_sall{ll},tlock_tac_sall{10}); % TA-T = 'aud' isolated from multisensory
%     %   tlock_aud_s0{ll+10}=ft_math(cfg,tlock_aud_s0{ll},tlock_aud_s0{10}); % AT-A = 'tac' isolated from mulitsensory
%     %   tlock_aud_sall{ll+10}=ft_math(cfg,tlock_aud_sall{ll},tlock_aud_sall{10}); % AT-A = 'tac' isolated from mulitsensory
%     %
%     %   if plotflag
%     %     figure(10);
%     %     subplot(2,5,ll-2);plot(tlock_tac_sallavg{30+ll}.time,tlock_tac_sallavg{30+ll}.avg(chanuse,:),'k');
%     %     hold on;plot(tlock_aud_sall{10+ll}.time,tlock_aud_sall{10+ll}.avg(chanuse,:),'b');axis([-.2 .5 -inf inf])
%     %     legend('tactile-alone shifted','AT-A: residual tactile')
%     %     subplot(2,5,ll-2+5);plot(tlock_aud_sallavg{30+ll}.time,tlock_aud_sallavg{30+ll}.avg(chanuse,:),'k');
%     %     hold on;plot(tlock_tac_sall{10+ll}.time,tlock_tac_sall{10+ll}.avg(chanuse,:),'b');axis([-.2 .5 -inf inf])
%     %     legend('auditory-alone shifted','TA-T: residual auditory')
%     %   end
%
%
%     cfg=[];
%     cfg.operation='subtract';
%     cfg.parameter='avg';
%     %   tlock_tpa_mtamn_s0{ll}=ft_math(cfg,ft_timelockanalysis([],tlock_tacPaud_s0{ll}),ft_timelockanalysis([],tlock_tac_s0{ll+20})); % (T + As) - (TA + N)
%     %   tlock_tpa_mtamn_sall{ll}=ft_math(cfg,ft_timelockanalysis([],tlock_tacPaud_sall{ll}),ft_timelockanalysis([],tlock_tac_sall{ll+20})); % (T + As) - (TA + N)
%     %   tlock_apt_matmn_s0{ll}=ft_math(cfg,ft_timelockanalysis([],tlock_audPtac_s0{ll}),ft_timelockanalysis([],tlock_aud_s0{ll+20})); % (A + Ts) - (AT + N)
%     %   tlock_apt_matmn_sall{ll}=ft_math(cfg,ft_timelockanalysis([],tlock_audPtac_sall{ll}),ft_timelockanalysis([],tlock_aud_sall{ll+20})); % (A + Ts) - (AT + N)
%     tlock_tpa_mtamn_s0{ll,tt}=ft_math(cfg,tlock_tacPaud_s0{ll,tt},tlock_tacMSpN_s0{ll,tt}); % (T + As) - (TA + N)
% %     tlock_tpa_mtamn_sall{ll}=ft_math(cfg,tlock_tacPaud_sall{ll},tlock_tac_sall{ll+20}); % (T + As) - (TA + N)
%     tlock_apt_matmn_s0{ll,tt}=ft_math(cfg,tlock_audPtac_s0{ll,tt},tlock_audMSpN_s0{ll,tt}); % (A + Ts) - (AT + N)
% %     tlock_apt_matmn_sall{ll}=ft_math(cfg,tlock_audPtac_sall{ll},tlock_aud_sall{ll+20}); % (A + Ts) - (AT + N)
%
%     if plotflag && tt==2
%       figure(11); % (Sum of Unisensory) minus (Multisensory plus Nul)
%       if ii>7
%         if ll==9
%           lll=7;
%         elseif ll>1
%           lll=ll-1;
%         else
%           lll=ll;
%         end
%         subplot(2,7,lll);plot(tlock_tpa_mtamn_s0{ll,tt}.time,tlock_tpa_mtamn_s0{ll,tt}.avg(chanuse,:),'k');axis([-.2 .5 -7 7]);
%         legend('(T+A)-(TA-N), time0 is tactile')
%         subplot(2,7,lll+7);plot(tlock_apt_matmn_s0{ll,tt}.time,tlock_apt_matmn_s0{ll,tt}.avg(chanuse,:),'k');axis([-.2 .5 -7 7])
%         legend('(A+T)-(AT-N), time0 is auditory')
%       else
%         subplot(2,5,ll-2);plot(tlock_tpa_mtamn_s0{ll,tt}.time,tlock_tpa_mtamn_s0{ll,tt}.avg(chanuse,:),'k');axis([-.2 .5 -7 7]);
%         legend('(T+A)-(TA-N), time0 is tactile')
%         subplot(2,5,ll-2+5);plot(tlock_apt_matmn_s0{ll,tt}.time,tlock_apt_matmn_s0{ll,tt}.avg(chanuse,:),'k');axis([-.2 .5 -7 7])
%         legend('(A+T)-(AT-N), time0 is auditory')
%       end
%     end
%
%   end
% end
%
%
% % cfg=[];
% % cfg.operation='add';
% % cfg.parameter='avg';
% % tlock_tacPaud_s0{10}=ft_math(cfg,tlock_aud_s0{10},tlock_tac_s0{10}); % Talone + Aalone
% % tlock_tacPaud_sall{10}=ft_math(cfg,tlock_aud_sall{10},tlock_tac_sall{10}); % Talone + Aalone
% % timets0=dsearchn(tlock_tac_s0{10}.time',tlock_tacPaud_s0{10}.time(1)):dsearchn(tlock_tac_s0{10}.time',tlock_tacPaud_s0{10}.time(end));
% % timetsall=dsearchn(tlock_tac_sall{10}.time',tlock_tacPaud_sall{10}.time(1)):dsearchn(tlock_tac_sall{10}.time',tlock_tacPaud_sall{10}.time(end));
% % timeas0=dsearchn(tlock_aud_s0{10}.time',tlock_tacPaud_s0{10}.time(1)):dsearchn(tlock_aud_s0{10}.time',tlock_tacPaud_s0{10}.time(end));
% % timeasall=dsearchn(tlock_aud_sall{10}.time',tlock_tacPaud_sall{10}.time(1)):dsearchn(tlock_aud_sall{10}.time',tlock_tacPaud_sall{10}.time(end));
% % % red and cyan line should be identical.
% % figure;plot(tlock_tac_s0{10}.time(timets0),[tlock_tac_s0{10}.avg(17,timets0); tlock_aud_s0{10}.avg(17,timeas0); tlock_tacPaud_s0{10}.avg(17,:); tlock_tac_s0{10}.avg(17,timets0)+tlock_aud_s0{10}.avg(17,timeas0)]')
% % figure;plot(tlock_tac_sall{10}.time(timetsall),[tlock_tac_sall{10}.avg(17,timetsall); tlock_aud_sall{10}.avg(17,timeasall); tlock_tacPaud_sall{10}.avg(17,:); tlock_tac_sall{10}.avg(17,timetsall)+tlock_aud_sall{10}.avg(17,timeasall)]')
%
%
%
% % for ll=1:length(soalist)
%   %   % this shifts tactile-alone to match where it would be for SOA conditions of soalist
%   %   cfg=[];
%   %   cfg.offset=-round(soatimes(soalist(ll))/1000/diff(tlock_aud_sall{3}.time(1:2))); % Note, the negative sign is on purpose here only, not for auditory
%   %   tlock_tac_s0{soalist(ll)+30}=ft_redefinetrial(cfg,tlock_tac_s0{10}); % Talone plus shift in time
%   %   tlock_tac_sall{soalist(ll)+30}=ft_redefinetrial(cfg,tlock_tac_sall{10}); % Talone plus shift in time
%   %   tlock_tac_s0avg{soalist(ll)+30}=ft_timelockanalysis([],tlock_tac_s0{soalist(ll)+30})
%   %   tlock_tac_sallavg{soalist(ll)+30}=ft_timelockanalysis([],tlock_tac_sall{soalist(ll)+30})
%   %   % tlock_tac_s0{33} is with tac starting 70ms later.
%   %
%   %   % this shifts auditory-alone to match where it would be for SOA conditions of soalist
%   %   cfg=[];
%   %   cfg.offset=round(soatimes(soalist(ll))/1000/diff(tlock_aud_sall{3}.time(1:2)));
%   %   tlock_aud_s0{soalist(ll)+30}=ft_redefinetrial(cfg,tlock_aud_s0{10}); % Aalone plus shift in time
%   %   tlock_aud_sall{soalist(ll)+30}=ft_redefinetrial(cfg,tlock_aud_sall{10}); % Aalone plus shift in time
%   %   tlock_aud_s0avg{soalist(ll)+30}=ft_timelockanalysis([],tlock_aud_s0{soalist(ll)+30})
%   %   tlock_aud_sallavg{soalist(ll)+30}=ft_timelockanalysis([],tlock_aud_sall{soalist(ll)+30})
%   %   % tlock_aud_s0{33} is with aud starting 70ms sooner.
%
%   %   cfg=[];
%   %   cfg.operation='add';
%   %   cfg.parameter='avg';
%   % %   tlock_tacPaud_s0{soalist(ll)}=ft_math(cfg,ft_timelockanalysis([],tlock_tac_s0{10}),ft_timelockanalysis([],tlock_aud_s0{soalist(ll)+30})); % Talone + Aalone_shifted
%   % %   tlock_tacPaud_sall{soalist(ll)}=ft_math(cfg,ft_timelockanalysis([],tlock_tac_sall{10}),ft_timelockanalysis([],tlock_aud_sall{soalist(ll)+30})); % Talone + Aalone_shifted
%   %   tlock_tacPaud_s0{soalist(ll)}=ft_math(cfg,tlock_tac_s0{10},tlock_aud_s0avg{soalist(ll)+30}); % Talone + Aalone_shifted
%   %   tlock_tacPaud_sall{soalist(ll)}=ft_math(cfg,tlock_tac_sall{10},tlock_aud_sallavg{soalist(ll)+30}); % Talone + Aalone_shifted
%   %
%   %   cfg=[];
%   %   cfg.operation='add';
%   %   cfg.parameter='avg';
%   % %   tlock_audPtac_s0{soalist(ll)}=ft_math(cfg,ft_timelockanalysis([],tlock_aud_s0{10}),ft_timelockanalysis([],tlock_tac_s0{soalist(ll)+30})); % Aalone + Talone_shifted
%   % %   tlock_audPtac_sall{soalist(ll)}=ft_math(cfg,ft_timelockanalysis([],tlock_aud_sall{10}),ft_timelockanalysis([],tlock_tac_sall{soalist(ll)+30})); % Aalone + Talone_shifted
%   %   tlock_audPtac_s0{soalist(ll)}=ft_math(cfg,tlock_aud_s0{10},tlock_tac_s0avg{soalist(ll)+30}); % Aalone + Talone_shifted
%   %   tlock_audPtac_sall{soalist(ll)}=ft_math(cfg,tlock_aud_sall{10},tlock_tac_sallavg{soalist(ll)+30}); % Aalone + Talone_shifted
% % end
%
% if plotflag
%   try close(9);end
%   figure(9);
%   subplot(2,3,1);
%   plot(tlock_tac_s0{10}.time,tlock_tac_s0{10}.avg(chanuse,:),'k');
%   hold on;plot(tlock_aud_s0{10}.time,tlock_aud_s0{10}.avg(chanuse,:),'b');
%   hold on;plot(tlock_tacPaud_s0{5}.time,tlock_tacPaud_s0{5,1}.avg(chanuse,:),'g');
%   hold on;plot(tlock_audPtac_s0{5}.time,tlock_audPtac_s0{5,1}.avg(chanuse,:),'r');axis([-.2 .5 -inf inf])
%   legend('tactile alone','auditory alone','sum unisensory, tactile at 0','sum unisensory, auditory at 0')
%
%
%   for ll=3:7
%     subplot(2,3,ll-1);
%     plot(tlock_tac_s0{ll}.time,tlock_tac_s0{ll}.avg(chanuse,:),'k');
%     hold on;plot(tlock_aud_s0{ll}.time,tlock_aud_s0{ll}.avg(chanuse,:),'b');
%     hold on;plot(tlock_tacPaud_s0{ll,1}.time,tlock_tacPaud_s0{ll,1}.avg(chanuse,:),'g');
%     hold on;plot(tlock_audPtac_s0{ll,1}.time,tlock_audPtac_s0{ll,1}.avg(chanuse,:),'r');axis([-.2 .5 -inf inf])
%     legend('multisensory, tactile at 0','multisensory, auditory at 0','sum unisensory, tactile at 0','sum unisensory, auditory at 0')
%   end
%
% end
%
%
%
% % save(['tlock_diffs_' sub{ii} '.mat'],'tlock_*_m*','*trialkept')
% save(['tlock_diffs_' sub{ii} '.mat'],'*trialkept','tlock_*P*','tlock_*N*')
%
% % figure;
% % channel={'Fz','Cz'};
% % for ll=1:length(soalist)
% %   cfg=[];
% %   cfg.channel=channel;
% %   cfg.ylim=[-6 6];
% %   cfg.xlim=[-.1 .4];
% %   cfg.layout='EEG1010.lay';
% %
% %   subplot(2,5,ll)
% %   ft_singleplotER(cfg,tlock_tpa_mtamn_s0{soalist(ll)})
% %   ft_singleplotER(cfg,tlock_tpa_mtamn_sall{soalist(ll)})
% %   subplot(2,5,ll+5)
% %   ft_singleplotER(cfg,tlock_apt_matmn_s0{soalist(ll)})
% %   ft_singleplotER(cfg,tlock_apt_matmn_sall{soalist(ll)})
% % end
% % figure;
% % channel={'Fz'};
% % for ll=1:length(soalist)
% %   cfg=[];
% %   cfg.channel=channel;
% %   cfg.ylim=[-4 4];
% %   cfg.xlim=[-.1 .4];
% %   cfg.layout='EEG1010.lay';
% %
% %   subplot(2,5,ll)
% %   ft_singleplotER(cfg,tlock_tpa_mtamn_s0{soalist(ll)})
% %   ft_singleplotER(cfg,tlock_tpa_mtamn_sall{soalist(ll)})
% %   subplot(2,5,ll+5)
% %   ft_singleplotER(cfg,tlock_apt_matmn_s0{soalist(ll)})
% %   ft_singleplotER(cfg,tlock_apt_matmn_sall{soalist(ll)})
% % end
%
% fzpeaks0=nan(2,soalist(end));
% fzpeaksall=nan(2,soalist(end));
% for ll=soalist
%   fzpeaks0(1,ll)=tlock_tpa_mtamn_s0{ll}.avg(match_str(tlock_tac_s0{ll}.label,'Fz'),dsearchn(tlock_tpa_mtamn_s0{ll}.time',.18));
%   fzpeaks0(2,ll)=tlock_apt_matmn_s0{ll}.avg(match_str(tlock_tac_s0{ll}.label,'Fz'),dsearchn(tlock_apt_matmn_s0{ll}.time',.18));
%   fzpeaks0(1,ll)=tlock_tpa_mtamn_s0{ll}.avg(match_str(tlock_tac_s0{ll}.label,'Fz'),dsearchn(tlock_tpa_mtamn_s0{ll}.time',.1));
%   fzpeaks0(2,ll)=tlock_apt_matmn_s0{ll}.avg(match_str(tlock_tac_s0{ll}.label,'Fz'),dsearchn(tlock_apt_matmn_s0{ll}.time',.1));
% %   fzpeaksall(1,ll)=tlock_tpa_mtamn_sall{ll}.avg(match_str(tlock_tac_sall{ll}.label,'Fz'),dsearchn(tlock_tpa_mtamn_sall{ll}.time',.18));
% %   fzpeaksall(2,ll)=tlock_apt_matmn_sall{ll}.avg(match_str(tlock_tac_sall{ll}.label,'Fz'),dsearchn(tlock_apt_matmn_sall{ll}.time',.18));
% end
%
% if 0
%   cfg=[];
%   cfg.interactive='yes';
%   cfg.layout='EEG1010.lay';
%   ft_multiplotER(cfg,tlock_tac_sall{15});
%   ft_multiplotER(cfg,tlock_tac{chanuse});
%   ft_multiplotER(cfg,tlock_tpa_mtamn_sall{5});
%   ft_multiplotER(cfg,tlock_tpa_matmn{5});
% end
%
