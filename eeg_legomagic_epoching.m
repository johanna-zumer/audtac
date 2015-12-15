% preprocessing of EEG data from Hills, 64ch MR-cap
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


component{2}=[9 3 25]; % ii=2, ff=3; 9=heart, 3=eyeblink, 25 artifact
component{3}=[2 3 6 8 10 17 20 28]; % ii=3, ff=7;  =heart, 2,3,6,8,10,17,20,28=eyeblink,  artifact
component{4}=[1 2 20 4 5 6 7 8 10 ]; % ii=4, sitting_concat;  =heart, 1,2,20=eyeblink,  4,5,6,7,8,10=artifact
component{5}=[1 2]; % ii=05, 's' 1,2 eyeblink
component{6}=[1 2 5 6 29]; % ii=06, 's' 1,2 eyeblink, 5,6 heartbeat, 29 artifact
component{7}=[4 5]; %ii=07, 's' 4,5 eyemovement
component{8}=[1 5 23]; %ii=08, 's' 1,5 eye;  23 artifact
component{9}=[1]; %iii=09;  's' 1 eye
component{10}=[]; %ii=10; none fit obvious eye or artifact (all blend)
component{11}=[2 3]; %ii=11; 's' 2,3 eye
component{12}=[]; %ii=12; none fit obvious eye or artifact (all blend)
component{13}=[3 4]; %ii=13; 's' 3 4 eye
component{14}=[7 8]; %ii=14; 's' 7,8 eye
component{15}=[1 2]; %ii=15; 's' 1,2 eye
component{16}=[4 6 23 27]; %ii=16; 's'  4 6eye   23 27 artifact
component{17}=[1 2 16 17 27]; %ii=17; 's' 1 2  eye  16 17 27 artifact
component{18}=[1 19]; %ii=18; 's' 1  eye  19 artifact
component{19}=[1 ]; %ii=19; 's' 1  eye  
component{20}=[1 2 3 30]; %ii=20; 's' 1,2,3  eye  30 artifact
component{21}=[1 3 23 29]; %ii=21; 's' 1,,3  eye  23, 29 artifact
component{22}=[1 6]; %ii=22; 's' 1,6  eye   artifact
component{23}=[2  16]; %ii=23; 's' 2 eye  16 artifact
component{24}=[2 3]; %ii=24; 's' 2,3 eye   artifact
component{25}=[1 4 9 20]; %ii=25; 's' 1,4,9 eye  20 artifact
component{26}=[1 4 30]; %ii=26; 's' 1,4 eye 30  artifact
component{27}=[1 2 3 5 8 9 10 17 18 22 28]; %ii=27; 's' 1,3 eye 2,5, 8, 9,10,17 ,18,22,28 artifact
component{28}=[2 3 5 15 26]; %ii=28; 's' 2,3,5 eye  15, 26 artifact
component{29}=[1 2 5 6 10 17 24 28 30]; %ii=29; 's' 1 2 5 6 10 eye 17 24 28 30  artifact
component{30}=[1 2 6 21 28 29 30]; %ii=30; 's' 1 2 6 eye 21 28 29 30  artifact
component{31}=[1 2 12 18 20 24 27 28]; %ii=31; 's' 1 2 eye 12 18 20 artifact
component{32}=[1 16 18 19 20 22 29]; %ii=32; 's' 1 eye  16 18 19 20 22 29 artifact

%%

sleep=1;
rej_comp=0;
rej_artfct=0; % change this to 1

%%

for ii=5:32
  
  clearvars -except ii sub component edir ddir sleep
  
  cd([edir sub{ii} ])
  
  if sleep
    files=dir('cc*b.mat');
  else
    files=dir('cc*s.mat');
  end
  cfg=[];
  cfg.dataset=files.name;
  raw_all=ft_preprocessing(cfg);
  cfg=[];
  cfg.demean='yes';
  cfg.channel={'all' '-ECG' '-EMG' '-EOG'};
  cfg.bsfilter='yes';
  cfg.bsfreq=[49 51; 99 101; 149 151];
  raw_all_demean=ft_preprocessing(cfg,raw_all);
  clear raw_all
  
  % apply filter to unepoched data
  cfg=[];
  cfg.hpfilter='yes';
  cfg.hpfiltord=3;
  cfg.hpfreq=0.2;  % Could do higher for awake data but not want higher for sleep.
  %   raw_all_ica_rej=ft_preprocessing(cfg,raw_all_ica_rej);
  tmp2=ft_preprocessing(cfg,raw_all_demean);
  clear raw_all_demean
    
  if rej_comp    
    load comp30.mat    
    cfg=[];
    cfg.component=component{ii};
    raw_all_ica=ft_rejectcomponent(cfg,comp30,tmp2);
    clear comp30
  else
    raw_all_ica=tmp2;
  end
  clear tmp2
  
  if rej_artfct && ~sleep
    load tmp.mat
    cfg=[];
    cfg.artfctdef=artfctdef;
    cfg.artfctdef.reject='partial';
    raw_all_ica_rej=ft_rejectartifact(cfg,raw_all_ica);
  elseif rej_artfct && sleep
    cfg=[];
    cfg.layout='EEG1010.lay';
    cfg=ft_databrowser(cfg,raw_all_ica);     
    cfg.artfctdef.reject='partial';
    raw_all_ica_rej=ft_rejectartifact(cfg,raw_all_ica);
  else
    raw_all_ica_rej=raw_all_ica;
  end
  clear raw_all_ica
  
  % create trl with 30s segments
  for ll=1:ceil( raw_all_ica_rej.sampleinfo(2)/30000)
    if ll==ceil( raw_all_ica_rej.sampleinfo(2)/30000)
      trl(ll,:)=[(30000*(ll-1))+1 raw_all_ica_rej.sampleinfo(2) 0];
    elseif ll==1
      trl(ll,:)=[200 (30000*(ll-1))+30000 0];
    else
      trl(ll,:)=[(30000*(ll-1))+1 (30000*(ll-1))+30000 0];
    end
  end
  
  if 0
    load([cfg.dataset(1:end-4) '_event.mat']); % loads event saved out during ft_trialfun_general_sleep.m
    eventS10=ft_filter_event(event,'value','S 10');
    eventS2=ft_filter_event(event,'value','S  2');
    eventS1=ft_filter_event(event,'value','S  1');
    [eventS10.type]=deal('A');
    [eventS2.type]=deal('T');
    [eventS1.type]=deal('N');
    %   [eventS2.offset]=deal(0);
    for ll=1:length(eventS2),eventS2(ll).sample=eventS2(ll).sample+200;end
    eventStim=[eventS2 eventS10 eventS1];
    
    cfg=[];
    cfg.layout='EEG1010.lay';
    cfg.viewmode='vertical';
    cfg.blocksize=30;
    cfg.plotlabels='yes';
    cfg.channel={'Fz' 'Cz' 'Pz' 'Oz'};
    cfg.event=eventStim;
    cfg.ploteventlabels='colorvalue';
    dcfg=ft_databrowser(cfg,raw_all_ica_rej);
  end
  
  if 0
    % Kcomplex
    cfg=[];
    cfg.trl=trl; % 30s segments
    cfg.continuous='yes';
    cfg.artfctdef.zvalue.channel={'Pz'};
    cfg.artfctdef.zvalue.cutoff=4;
    cfg.artfctdef.zvalue.interactive = 'yes';
    cfg.artfctdef.zvalue.lpfilter='yes';
    cfg.artfctdef.zvalue.lpfreq=8;
    cfg.artfctdef.zvalue.lpfiltord=3;
    [cfg, artifact] = ft_artifact_zvalue(cfg, raw_all_ica_rej)
  
  cfg=[];
  cfg.layout='EEG1010.lay';
  cfg.viewmode='vertical';
  cfg.blocksize=30;
  cfg.plotlabels='yes';
  cfg.channel={'Fz' 'Cz' 'Pz' 'Oz'};
  cfg.artfctdef.zvalue.artifact=artifact;
  cfg.preproc.lpfilter='yes';
  cfg.preproc.lpfreq=30;
  dcfg=ft_databrowser(cfg,raw_all_ica_rej);

  % muscle
  cfg=[];
  cfg.trl=trl(2:end-1,:); % 30s segments
  cfg.continuous='yes';
  cfg.artfctdef.muscle.channel={'FT7' 'FT8' 'TP7' 'TP8'};
  cfg.artfctdef.muscle.bpfilter    = 'yes';
  cfg.artfctdef.muscle.bpfreq      = [80 140];
  cfg.artfctdef.muscle.bpfiltord   = 8;
  cfg.artfctdef.muscle.bpfilttype  = 'but';
  cfg.artfctdef.muscle.hilbert     = 'yes';
  cfg.artfctdef.muscle.cutoff = 2;
  cfg.artfctdef.muscle.interactive = 'yes';
  %   cfg.artfctdef.muscle.trlpadding  = 1;
  %   cfg.artfctdef.muscle.fltpadding  = 1;
  cfg.artfctdef.muscle.artpadding  = 3;
  [cfg,martifact]=ft_artifact_muscle(cfg,raw_all_ica_rej);
  
  cfg=[];
  cfg.layout='EEG1010.lay';
  cfg.viewmode='vertical';
  cfg.blocksize=30;
  cfg.plotlabels='yes';
  cfg.channel={'Fz' 'Cz' 'Pz' 'Oz'};
  cfg.artfctdef.zvalue.artifact=martifact;
%   cfg.preproc.lpfilter='yes';
%   cfg.preproc.lpfreq=30;
  dcfg=ft_databrowser(cfg,raw_all_ica_rej);
  
  end
  
  
  load([ddir sub{ii} '_audtac.mat']);
  saudtrials=find(saudseq);
  baudtrials=find(baudseq);
  % stactrials=find(stime_touch>0);
  % btactrials=find(btime_touch>0);
  if sleep
    load time_sleep.mat
  else
    load time.mat
  end
  
  % file=dir(['cc_cc_spm8_' sub{ii} '*.mat']);
  % cfg=[];
  % cfg.dataset=file.name;
  % raw_all=ft_preprocessing(cfg);
  
  % All trials
  cfg=[];
  if sleep
    cfg.dataset=fileb.name;
  else
    cfg.dataset=files.name;
  end
  cfg.trialfun='ft_trialfun_general_sleep';
  cfg.trialdef.eventtype  = 'Stimulus';
  cfg.trialdef.eventvalue = 'S 20'; % corresponds to lj.prepareStrobe(20)
  cfg.trialdef.prestim = 1.5; % kept this long from ii==12 and afterwards
  cfg.trialdef.poststim = 2.3;
  cfg.time=time;
  cfgtr=ft_definetrial(cfg);
  
  cfgtr.trl(:,7)=ssoareal(1,:); % <---- which index 1,2,3,4 to use????
  cfgtr.trl(:,8)=ssoareal(2,:); % <---- which index 1,2,3,4 to use????
  cfgtr.trl(:,9)=ssoareal(3,:); % <---- which index 1,2,3,4 to use????
  cfgtr.trl(:,10)=ssoareal(4,:); % <---- which index 1,2,3,4 to use????
  cfgtr.trl(:,11)=cfgtr.trl(:,6)>0 | cfgtr.trl(:,6)==-2; % tactile trials
  cfgtr.trl(:,12)=cfgtr.trl(:,6)>0 | cfgtr.trl(:,6)==-1; % auditory trials
  cfgtr.trl(:,13)=ssoaeff;
  cfgtr.trl(:,14)=1:size(cfgtr.trl,1);
  
  ns=nan(1,length(cfgtr.event));
  for ll=1:length(cfgtr.event),ns(ll)=strcmp(cfgtr.event(ll).type,'New Segment');end
  nsind=find(ns);
  
  if ii==20 % Must have started/stopped within same acquistion
    nsind=nsind([2 3 5 6]);
  elseif ii==21
    nsind=nsind([2 3 4 5]);
  elseif ii==32
    nsind=nsind([1:6 8]);
  end
  
  % Tactile
  % need to get timing for each trial centered on stimulus correctly
  cfg=[];
  cfg.dataset=files.name;
  cfg.trialfun='ft_trialfun_general_sleep';
  cfg.trialdef.eventtype  = 'Stimulus';
  cfg.trialdef.eventvalue = 'S  2'; % corresponds to lj.prepareStrobe(2)
  cfg.trialdef.prestim = 1.5;
  cfg.trialdef.poststim = 2.3;
  cfg.time=time;  %% <- why is this here?
  cfgtactr=ft_definetrial(cfg);
  
  for ll=1:length(cfgtr.event),tac(ll)=strcmp(cfgtr.event(ll).value,'S  2');end
  tacind=find(tac);
  for ll=1:length(nsind)-1
    cfgtacdiff(ll)=length(find(tac(nsind(ll):nsind(ll+1))))-snumtac(ll);
  end
  cfgtacdiff(length(nsind))=length(find(tac(nsind(ll+1):end)))-snumtac(ll+1);
  ss=0;
  trluse=[];
  for ll=1:length(nsind)
    %     trluse=[trluse [ss+1]:[ss+snumtac(ll)-cfgtacdiff(ll)]];
    trluse=[trluse [ss+1]:[ss+snumtac(ll)]];
    ss=ss+snumtac(ll)+cfgtacdiff(ll);
  end
  if length(trluse)~=length(stime_touch) % should be the same
    error('num tac wrong');
  end
  
  tactr=find(cfgtr.trl(:,11));
  tactrl=[cfgtactr.trl cfgtr.trl(tactr(trluse),7:end)];
  % tactrl=cfgtr.trl(tactr(trluse),[1:3 5:end]); % 4th is irrelavent
  % numtac=min(size(cfgtr.trl,1),length(time_touch));
  % cfgtr.trl=cfgtr.trl(1:numtac,:);
  % cfgtr.trl(:,3)=cfgtr.trl(:,3)-round(time_touch(1:numtac)*raw_all_ica_rej{ff}.fsample)';
  % cfgtr.trl(:,3)=cfgtr.trl(:,3)-round(stime_touch/diff(time{1}(1:2)))'; %using the difference between time as the sampling rate
  tactrl(:,3)=tactrl(:,3)-round(stime_touch/diff(time{1}(1:2)))'; %using the difference between time as the sampling rate
  
  find(abs(diff(stime_touch))>.01)
  
%  what's going on here?
%   load(['raw_all_ica_rej_' sub{ii}])
%   try
%     raw_all_ica_rej=raw_all_rej;clear raw_all_rej
%   end
  
  cfg=[];
  cfg.trl=tactrl;
  raw_tac=ft_redefinetrial(cfg,raw_all_ica_rej);
  
  % Auditory
  % need to get timing for each trial centered on stimulus correctly
  cfg=[];
  % cfg.dataset=files(ff).name;
  % cfg.dataset=['cc_cc_spm8_' sub{ii} '_01s.mat'];
  cfg.dataset=files.name;
  cfg.trialfun='ft_trialfun_general_sleep';
  cfg.trialdef.eventtype  = 'Stimulus';
  cfg.trialdef.eventvalue = 'S 10'; % corresponds to lj.prepareStrobe(10)
  cfg.trialdef.prestim = 1.5;
  cfg.trialdef.poststim = 2.3;
  cfg.time=time;
  cfgaudtr=ft_definetrial(cfg);
  
  for ll=1:length(cfgtr.event),aud(ll)=strcmp(cfgtr.event(ll).value,'S 10');end
  audind=find(aud);
  for ll=1:length(nsind)-1
    cfgauddiff(ll)=length(find(aud(nsind(ll):nsind(ll+1))))-snumaud(ll);
  end
  cfgauddiff(length(nsind))=length(find(aud(nsind(ll+1):end)))-snumaud(ll+1);
  ss=0;
  trluse=[];
  for ll=1:length(nsind)
    %     trluse=[trluse [ss+1]:[ss+snumtac(ll)-cfgtacdiff(ll)]];
    trluse=[trluse [ss+1]:[ss+snumaud(ll)]];
    ss=ss+snumaud(ll)+cfgauddiff(ll);
  end
  if length(trluse)~=length(saudtrials) % should be the same
    error('num aud wrong');
  end
  audtr1=find(cfgtr.trl(:,12));
  audtrl=[cfgaudtr.trl cfgtr.trl(audtr1(trluse),7:end)];
  % audtrl=cfgtr.trl(audtr1(trluse),[1:3 5:end]); % 4th is irrelavent
  % cfgtr.trl(:,3)=cfgtr.trl(:,3)-round(fin.audDelay/1000*raw_all_ica_rej{ff}.fsample);
  % cfgtr.trl(:,3)=cfgtr.trl(:,3)-round(.001*audDelay/diff(time{1}(1:2)))'; %using the difference between time as the sampling rate
  audtrl(:,3)=audtrl(:,3)-round(.001*audDelay/diff(time{1}(1:2)))'; %using the difference between time as the sampling rate
  cfg=[];
  cfg.trl=audtrl;
  raw_aud=ft_redefinetrial(cfg,raw_all_ica_rej);
  
  % Null events
  cfg=[];
  % cfg.dataset=files(ff).name;
  % cfg.dataset=['cc_cc_spm8_' sub{ii} '_01s.mat'];
  cfg.dataset=files.name;
  cfg.trialfun='ft_trialfun_general_sleep';
  cfg.trialdef.eventtype  = 'Stimulus';
  cfg.trialdef.eventvalue = 'S  1'; % corresponds to lj.prepareStrobe(10)
  cfg.trialdef.prestim = 1.5;
  cfg.trialdef.poststim = 2.3;
  cfg.time=time;
  cfgnultr=ft_definetrial(cfg);
  
  for ll=1:length(cfgtr.event),nul(ll)=strcmp(cfgtr.event(ll).value,'S  1');end
  nulind=find(nul);
  for ll=1:length(nsind)-1
    cfgnuldiff(ll)=length(find(nul(nsind(ll):nsind(ll+1))))-snumnul(ll);
  end
  cfgnuldiff(length(nsind))=length(find(nul(nsind(ll+1):end)))-snumnul(ll+1);
  ss=0;
  trluse=[];
  for ll=1:length(nsind)
    %     trluse=[trluse [ss+1]:[ss+snumtac(ll)-cfgtacdiff(ll)]];
    trluse=[trluse [ss+1]:[ss+snumnul(ll)]];
    ss=ss+snumnul(ll)+cfgnuldiff(ll);
  end
  if length(trluse)~=length(find(snulseq)) % should be the same
    error('num nul wrong');
  end
  nultr=find(~cfgtr.trl(:,6));
  nultrl=[cfgnultr.trl cfgtr.trl(nultr(trluse),7:end)];
  % nultrl=cfgtr.trl(nultr(trluse),[1:3 5:end]); % 4th is irrelavent
  % cfgtr.trl=cfgtr.trl(trluse,:);
  % [size(cfgtr.trl,1) length(find(snulseq))] % should be the same
  
  % use the auditory delay; it doesn't really matter, as it's random anyway
  cfgtr.trl(:,3)=cfgtr.trl(:,3)-round(.001*audDelay/diff(time{1}(1:2)))'; %using the difference between time as the sampling rate
  cfg=[];
  cfg.trl=nultrl;
  raw_nul=ft_redefinetrial(cfg,raw_all_ica_rej);
  
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

  clear raw_all_ica_rej
  % save(['raw_each_rej_' sub{ii}],'raw_*rej'); % with 1Hz highpass after epoching
  save(['raw1_each_rej_' sub{ii}],'raw_*rej','*art'); % with 0.2Hz highpass before epoching
  clear raw_aud raw_nul raw_tac
  
  
end

return



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
