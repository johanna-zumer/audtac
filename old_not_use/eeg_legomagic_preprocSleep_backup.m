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
sub{16}='e16';% t.p.  24/07/14

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

%% Load data

% data preprocessed by spm_preproc.m

% there should only be 1 file for each if successfully removed intermediates above
% files=dir('cc*s.mat');
fileb=dir('cc*b.mat');

% if any(ii==[3 ])
%   file=dir(['cc_cc_spm8_' sub{ii} '*.mat']);
% elseif any(ii==[4 ])
%   file=dir(['cc_spm8_' sub{ii} '*.mat']);  
% end

cfg=[];
cfg.dataset=fileb.name;
raw_all=ft_preprocessing(cfg);
time=raw_all.time;
save time.mat time
% % size of sleep-score chunk in samples
% chunksamp=raw_all.hdr.orig.other.CRC.score{3}*raw_all.fsample;
% % create sample-by-sample vector of sleep score
% for ll=1:length(raw_all.hdr.orig.other.CRC.score{1})
%     samplesleep((ll-1)*chunksamp+1:(ll*chunksamp),1)=raw_all.hdr.orig.other.CRC.score{1}(ll);
% end
% % add as new field to events structure the sleep score
% for ll=1:length(raw_all.hdr.orig.trials.events)
%     raw_all.hdr.orig.trials.events(ll).sleep= samplesleep(dsearchn(raw_all.time{1}',raw_all.hdr.orig.trials.events(ll).time));
% end
%
% trials.events=raw_all.hdr.orig.trials.events;
% save trialevents.mat trials
% % make sure this 'events' is getting carried through, or save as separate
% % struct to load in later?

cfg=[];
cfg.demean='yes';
cfg.channel={'all' '-ECG'};
cfg.bsfilter='yes';
cfg.bsfreq=[49 51; 99 101; 149 151];
raw_all_demean=ft_preprocessing(cfg,raw_all);
cfg=[];cfg.layout='EEG1010.lay';
cfg=ft_databrowser(cfg,raw_all_demean);  %% make sure to get all segments marked here!!

cfg.artfctdef.reject='partial';
tmp=ft_rejectartifact(cfg,raw_all_demean);

cfg=[];cfg.method='trial';
raw_all_rej=ft_rejectvisual(cfg,tmp);

save raw_all_rej.mat raw_all_rej

if 0
  cfg=[];
  cfg.numcomponent=30;
  cfg.method='fastica';
  cfg.randomseed=17;
else
  cfg=[];
  cfg.numcomponent=30;
  cfg.method='runica';
end
comp30=ft_componentanalysis(cfg,raw_all_rej);
% comp20=ft_componentanalysis(cfg,raw_all_rej);
save comp30 comp30
cfg=[];cfg.layout='EEG1010.lay';
cfg=ft_databrowser(cfg,comp30);
artfctdef=cfg.artfctdef;
save tmp.mat artfctdef

cfg=[];
cfg.component=[9 3 25]; % ii=2, ff=3; 9=heart, 3=eyeblink, 25 artifact
cfg.component=[2 3 6 8 10 17 20 28]; % ii=3, ff=7;  =heart, 2,3,6,8,10,17,20,28=eyeblink,  artifact
cfg.component=[1 2 20 4 5 6 7 8 10 ]; % ii=3, sitting_concat;  =heart, 1,2,20=eyeblink,  4,5,6,7,8,10=artifact
cfg.component=[1 2]; % ii=5, 's' 1,2 eyeblink
raw_all_ica=ft_rejectcomponent(cfg,comp30,raw_all_rej);
cfg=[];
cfg.artfctdef=artfctdef;
cfg.artfctdef.reject='partial';
raw_all_ica_rej=ft_rejectartifact(cfg,raw_all_ica);

cfg=[];cfg.layout='EEG1010.lay';cfg=ft_databrowser(cfg,raw_all_ica_rej);

cfg.artfctdef=artfctdef;
cfg.artfctdef.reject='partial';
raw_all_ica_rej=ft_rejectartifact(cfg,raw_all_ica_rej);


save(['raw_all_ica_rej_' sub{ii}],'raw_all_ica_rej','-v7.3');

%%  Trial sorting

concat_stim_files(ii);
cd([edir sub{ii} ])


% there should only be 1 file for each if successfully removed intermediates above
% files=dir('cc*s.mat');
fileb=dir('cc*b.mat');


load([ddir sub{ii} '_audtac.mat']);
saudtrials=find(saudseq);
baudtrials=find(baudseq);
% stactrials=find(stime_touch>0);
% btactrials=find(btime_touch>0);
load time.mat

if ii<3
  for tt=1:length(tactrials)
    indx5=find(info.lightsensor{tactrials(tt)}>.5*[max(info.lightsensor{tactrials(tt)})-min(info.lightsensor{tactrials(tt)})]+min(info.lightsensor{tactrials(tt)}));
    indx8=find(info.lightsensor{tactrials(tt)}>.8*[max(info.lightsensor{tactrials(tt)})-min(info.lightsensor{tactrials(tt)})]+min(info.lightsensor{tactrials(tt)}));
    try
      mint5(tt)=indx5(1);
      mint8(tt)=indx8(1);
    end
  end
  % turns out that 50% is better criterion (change it in main code) but then
  % for timing to be consistent, add back 4ms.
  time_touch=.001*(mint5+4-1); % add 4ms for diff of 50% vs 80%, and subtract 1 for the time relative to info.lighttime{tt}(1)
  if median(time_touch-info.time_touch(tactrials))>.001
    error('time_touch not right?')
  end
else
  %     stime_touch=stime_touch(stactrials);
  %     btime_touch=btime_touch(btactrials);
end

% file=dir(['cc_cc_spm8_' sub{ii} '*.mat']);
% cfg=[];
% cfg.dataset=file.name;
% raw_all=ft_preprocessing(cfg);


% Do tactile
cfg=[];
cfg.dataset=fileb.name;
cfg.trialfun='ft_trialfun_general_sleep';
cfg.trialdef.eventtype  = 'Stimulus';
cfg.trialdef.eventvalue = 'S  2'; % corresponds to lj.prepareStrobe(2)
cfg.trialdef.prestim = 0.5;
cfg.trialdef.poststim = 2.3;
cfg.time=time;
cfgtr=ft_definetrial(cfg);

for ll=1:length(cfgtr.event),ns(ll)=strcmp(cfgtr.event(ll).type,'New Segment');end
nsind=find(ns);
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
cfgtr.trl=cfgtr.trl(trluse,:);
[size(cfgtr.trl,1) length(stime_touch)] % should be the same

% numtac=min(size(cfgtr.trl,1),length(time_touch));
% cfgtr.trl=cfgtr.trl(1:numtac,:);
% cfgtr.trl(:,3)=cfgtr.trl(:,3)-round(time_touch(1:numtac)*raw_all_ica_rej{ff}.fsample)';
cfgtr.trl(:,3)=cfgtr.trl(:,3)-round(stime_touch/diff(time{1}(1:2)))'; %using the difference between time as the sampling rate

load(['raw_all_ica_rej_' sub{ii}])
try
  raw_all_ica_rej=raw_all_rej;clear raw_all_rej
end

cfg=[];
cfg.trl=cfgtr.trl;
raw_tac=ft_redefinetrial(cfg,raw_all_ica_rej);

% Do auditory
cfg=[];
% cfg.dataset=files(ff).name;
% cfg.dataset=['cc_cc_spm8_' sub{ii} '_01s.mat'];
cfg.dataset=files.name;
cfg.trialfun='ft_trialfun_general_sleep';
cfg.trialdef.eventtype  = 'Stimulus';
cfg.trialdef.eventvalue = 'S 10'; % corresponds to lj.prepareStrobe(10)
cfg.trialdef.prestim = 0.7;
cfg.trialdef.poststim = 2.1;
cfg.time=time;
cfgtr=ft_definetrial(cfg);

for ll=1:length(cfgtr.event),ns(ll)=strcmp(cfgtr.event(ll).type,'New Segment');end
nsind=find(ns);
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
cfgtr.trl=cfgtr.trl(trluse,:);
[size(cfgtr.trl,1) length(saudtrials)] % should be the same

% cfgtr.trl(:,3)=cfgtr.trl(:,3)-round(fin.audDelay/1000*raw_all_ica_rej{ff}.fsample);
cfgtr.trl(:,3)=cfgtr.trl(:,3)-round(.001*audDelay/diff(time{1}(1:2)))'; %using the difference between time as the sampling rate
cfg=[];
cfg.trl=cfgtr.trl;
raw_aud=ft_redefinetrial(cfg,raw_all_ica_rej);


% Do null events
cfg=[];
% cfg.dataset=files(ff).name;
% cfg.dataset=['cc_cc_spm8_' sub{ii} '_01s.mat'];
cfg.dataset=files.name;
cfg.trialfun='ft_trialfun_general_sleep';
cfg.trialdef.eventtype  = 'Stimulus';
cfg.trialdef.eventvalue = 'S  1'; % corresponds to lj.prepareStrobe(10)
cfg.trialdef.prestim = 0.7;
cfg.trialdef.poststim = 2.1;
cfg.time=time;
cfgtr=ft_definetrial(cfg);

for ll=1:length(cfgtr.event),ns(ll)=strcmp(cfgtr.event(ll).type,'New Segment');end
nsind=find(ns);
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
cfgtr.trl=cfgtr.trl(trluse,:);
[size(cfgtr.trl,1) length(find(snulseq))] % should be the same

% use the auditory delay; it doesn't really matter, as it's random anyway
cfgtr.trl(:,3)=cfgtr.trl(:,3)-round(.001*audDelay/diff(time{1}(1:2)))'; %using the difference between time as the sampling rate
cfg=[];
cfg.trl=cfgtr.trl;
raw_nul=ft_redefinetrial(cfg,raw_all_ica_rej);


% final trial/channel rejection
cfg=[];
cfg.method='summary';
% raw_all_rej=ft_rejectvisual(cfg,raw_all);
raw_tac_rej=ft_rejectvisual(cfg,raw_tac);
raw_aud_rej=ft_rejectvisual(cfg,raw_aud);
raw_nul_rej=ft_rejectvisual(cfg,raw_nul);

clear raw_all_ica_rej
save(['raw_each_rej_' sub{ii}],'raw_*rej');
clear raw_aud raw_nul raw_tac

%% temp test 50hz reref e04



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



%% ERP filtering

load(['raw_each_rej_' sub{ii}],'raw_*');

cfg=[];
cfg.bpfilter='yes';
cfg.bpfreq=[1 40];
cfg.demean='yes';
cfg.baselinewindow=[-.6 -.1];
cfg.channel= {'all', '-ECG'}
data_tac_filt=ft_preprocessing(cfg,raw_tac_rej);
data_aud_filt=ft_preprocessing(cfg,raw_aud_rej);
data_nul_filt=ft_preprocessing(cfg,raw_nul_rej);

% % Solve this: what should be reference channel for scalp ERP?
cfg=[];
cfg.reref='yes';
cfg.refchannel='all'; % necessary if doing source localisation
cfg.refchannel={'FT9', 'FT10'}; % sort of like linked mastoids?
data_tac_filt_ref=ft_preprocessing(cfg,data_tac_filt);
data_aud_filt_ref=ft_preprocessing(cfg,data_aud_filt);
data_nul_filt_ref=ft_preprocessing(cfg,data_nul_filt);

clear raw_*


% soalist=[1 3 4 5 6 7 9];
soalist=[3 4 5 6 7];
for ll=1:length(soalist)
  cfg=[];
  cfg.vartrllength=2;
  cfg.keeptrials='yes';
  cfg.trials=[data_tac_filt_ref.trialinfo(:,3)==soalist(ll) & data_tac_filt_ref.trialinfo(:,2)==0];
  tlock_tac_s0{soalist(ll)}=ft_timelockanalysis(cfg,data_tac_filt_ref); % only 'awake' state
  cfg.trials=data_tac_filt_ref.trialinfo(:,3)==soalist(ll);
  tlock_tac_sall{soalist(ll)}=ft_timelockanalysis(cfg,data_tac_filt_ref); % all trials, irrepective if fell asleep
end
for ll=1:length(soalist)
  cfg=[];
  cfg.vartrllength=2;
  cfg.keeptrials='yes';
  cfg.trials=[data_aud_filt_ref.trialinfo(:,3)==soalist(ll) & data_aud_filt.trialinfo(:,2)==0];
  tlock_aud_s0{soalist(ll)}=ft_timelockanalysis(cfg,data_aud_filt_ref);
  cfg.trials=data_aud_filt_ref.trialinfo(:,3)==soalist(ll);
  tlock_aud_sall{soalist(ll)}=ft_timelockanalysis(cfg,data_aud_filt_ref);
end
cfg=[]; % tac alone 38 (-2)
cfg.vartrllength=2;
cfg.keeptrials='yes';
cfg.trials=[data_tac_filt_ref.trialinfo(:,3)==-2 & data_tac_filt_ref.trialinfo(:,2)==0];
tlock_tac_s0{10}=ft_timelockanalysis(cfg,data_tac_filt_ref);
cfg.trials=data_tac_filt_ref.trialinfo(:,3)==-2;
tlock_tac_sall{10}=ft_timelockanalysis(cfg,data_tac_filt_ref);
cfg=[]; % aud alone 39 (-1)
cfg.vartrllength=2;
cfg.keeptrials='yes';
cfg.trials=[data_aud_filt_ref.trialinfo(:,3)==-1 & data_aud_filt_ref.trialinfo(:,2)==0];
tlock_aud_s0{10}=ft_timelockanalysis(cfg,data_aud_filt_ref);
cfg.trials=data_aud_filt_ref.trialinfo(:,3)==-1;
tlock_aud_sall{10}=ft_timelockanalysis(cfg,data_aud_filt_ref);
cfg=[];
cfg.vartrllength=2;
cfg.keeptrials='yes';
cfg.trials=data_nul_filt_ref.trialinfo(:,2)==0;
tlock_nul_s0{10}=ft_timelockanalysis(cfg,data_nul_filt_ref);
cfg.trials='all';
tlock_nul_sall{10}=ft_timelockanalysis(cfg,data_nul_filt_ref);

for ll=[soalist 10] % this loop is to remove trials with NaN in (e.g. if artifact rejection cut through middle of a trial)
  cfg=[];
  cfg.keeptrials='yes';
  cfg.trials=find(sum(squeeze(isnan(tlock_tac_s0{ll}.trial(:,10,:)))')<.05*size(tlock_tac_s0{ll}.trial,3));
  tlock_tac_s0{ll}=ft_timelockanalysis(cfg,tlock_tac_s0{ll});
  cfg.trials=find(sum(squeeze(isnan(tlock_tac_sall{ll}.trial(:,10,:)))')<.05*size(tlock_tac_sall{ll}.trial,3));
  tlock_tac_sall{ll}=ft_timelockanalysis(cfg,tlock_tac_sall{ll});
  cfg=[];
  cfg.keeptrials='yes';
  cfg.trials=find(sum(squeeze(isnan(tlock_aud_s0{ll}.trial(:,10,:)))')<.05*size(tlock_aud_s0{ll}.trial,3));
  tlock_aud_s0{ll}=ft_timelockanalysis(cfg,tlock_aud_s0{ll});
  cfg.trials=find(sum(squeeze(isnan(tlock_aud_sall{ll}.trial(:,10,:)))')<.05*size(tlock_aud_sall{ll}.trial,3));
  tlock_aud_sall{ll}=ft_timelockanalysis(cfg,tlock_aud_sall{ll});
  if ll==10
    cfg=[];
    cfg.keeptrials='yes';
    cfg.trials=find(sum(squeeze(isnan(tlock_nul_s0{ll}.trial(:,10,:)))')<.05*size(tlock_nul_s0{ll}.trial,3));
    tlock_nul_s0{ll}=ft_timelockanalysis(cfg,tlock_nul_s0{ll});
    cfg.trials=find(sum(squeeze(isnan(tlock_nul_sall{ll}.trial(:,10,:)))')<.05*size(tlock_nul_sall{ll}.trial,3));
    tlock_nul_sall{ll}=ft_timelockanalysis(cfg,tlock_nul_sall{ll});
  end
end


% % collapsing over all SOA conditions
% clear keeptrial
% for ll=1:size(data_tac_filt_ref.trial,2),
%   keeptrial(ll)=~any(any(isnan(data_tac_filt_ref.trial{ll})));
% end
% cfg=[];
% cfg.vartrllength=2;
% cfg.trials=find(keeptrial);
% cfg.keeptrials='yes';
% tlock_tac{11}=ft_timelockanalysis(cfg,data_tac_filt_ref);
% 
% % collapsing over all SOA conditions
% clear keeptrial
% for ll=1:size(data_aud_filt_ref.trial,2),
%   keeptrial(ll)=~any(any(isnan(data_aud_filt_ref.trial{ll})));
% end
% cfg=[];
% cfg.vartrllength=2;
% cfg.trials=find(keeptrial);
% cfg.keeptrials='yes';
% tlock_aud{11}=ft_timelockanalysis(cfg,data_aud_filt_ref);


if 0
  cfg=[];
  cfg.interactive='yes';
  cfg.layout='EEG1010.lay';
  ft_multiplotER(cfg,tlock_tac_s0{10});
  ft_multiplotER(cfg,tlock_tac_sall{10});
  ft_multiplotER(cfg,tlock_aud_s0{10});
  ft_multiplotER(cfg,tlock_aud_sall{10});
  ft_multiplotER(cfg,tlock_nul_s0{10});
  ft_multiplotER(cfg,tlock_nul_sall{10});
  
  % ft_multiplotER(cfg,tlock_tac{11});
  % ft_multiplotER(cfg,tlock_aud{11});
  
  ft_multiplotER(cfg,tlock_tac_sall{5});
  ft_multiplotER(cfg,tlock_aud_sall{5});
end

%%  contrasting conditions

plotflag=1;
load([ddir sub{ii} '_audtac.mat']);
cfg=[];
cfg.operation='add';
cfg.parameter='avg';
tlock_tacPaud_s0{10}=ft_math(cfg,tlock_aud_s0{10},tlock_tac_s0{10}); % Talone + Aalone
tlock_tacPaud_sall{10}=ft_math(cfg,tlock_aud_sall{10},tlock_tac_sall{10}); % Talone + Aalone


timets0=dsearchn(tlock_tac_s0{10}.time',tlock_tacPaud_s0{10}.time(1)):dsearchn(tlock_tac_s0{10}.time',tlock_tacPaud_s0{10}.time(end));
timetsall=dsearchn(tlock_tac_sall{10}.time',tlock_tacPaud_sall{10}.time(1)):dsearchn(tlock_tac_sall{10}.time',tlock_tacPaud_sall{10}.time(end));
timeas0=dsearchn(tlock_aud_s0{10}.time',tlock_tacPaud_s0{10}.time(1)):dsearchn(tlock_aud_s0{10}.time',tlock_tacPaud_s0{10}.time(end));
timeasall=dsearchn(tlock_aud_sall{10}.time',tlock_tacPaud_sall{10}.time(1)):dsearchn(tlock_aud_sall{10}.time',tlock_tacPaud_sall{10}.time(end));

% red and cyan line should be identical.
figure;plot(tlock_tac_s0{10}.time(timets0),[tlock_tac_s0{10}.avg(17,timets0); tlock_aud_s0{10}.avg(17,timeas0); tlock_tacPaud_s0{10}.avg(17,:); tlock_tac_s0{10}.avg(17,timets0)+tlock_aud_s0{10}.avg(17,timeas0)]')
figure;plot(tlock_tac_sall{10}.time(timetsall),[tlock_tac_sall{10}.avg(17,timetsall); tlock_aud_sall{10}.avg(17,timeasall); tlock_tacPaud_sall{10}.avg(17,:); tlock_tac_sall{10}.avg(17,timetsall)+tlock_aud_sall{10}.avg(17,timeasall)]')


chanuse=match_str(tlock_tac_sall{3}.label,'Fz');

for ll=1:length(soalist)
  % this shifts tactile-alone to match where it would be for SOA conditions of soalist
  cfg=[];
  cfg.offset=-round(soatimes(soalist(ll))/1000/diff(tlock_aud_sall{3}.time(1:2))); % Note, the negative sign is on purpose here only, not for auditory
  tlock_tac_s0{soalist(ll)+30}=ft_redefinetrial(cfg,tlock_tac_s0{10}); % Talone plus shift in time
  tlock_tac_sall{soalist(ll)+30}=ft_redefinetrial(cfg,tlock_tac_sall{10}); % Talone plus shift in time
  tlock_tac_s0avg{soalist(ll)+30}=ft_timelockanalysis([],tlock_tac_s0{soalist(ll)+30})
  tlock_tac_sallavg{soalist(ll)+30}=ft_timelockanalysis([],tlock_tac_sall{soalist(ll)+30})
  % tlock_tac_s0{33} is with tac starting 70ms later.
  
  % this shifts auditory-alone to match where it would be for SOA conditions of soalist
  cfg=[];
  cfg.offset=round(soatimes(soalist(ll))/1000/diff(tlock_aud_sall{3}.time(1:2)));
  tlock_aud_s0{soalist(ll)+30}=ft_redefinetrial(cfg,tlock_aud_s0{10}); % Aalone plus shift in time
  tlock_aud_sall{soalist(ll)+30}=ft_redefinetrial(cfg,tlock_aud_sall{10}); % Aalone plus shift in time
  tlock_aud_s0avg{soalist(ll)+30}=ft_timelockanalysis([],tlock_aud_s0{soalist(ll)+30})
  tlock_aud_sallavg{soalist(ll)+30}=ft_timelockanalysis([],tlock_aud_sall{soalist(ll)+30})
  % tlock_aud_s0{33} is with aud starting 70ms sooner.
  
  cfg=[];
  cfg.operation='add';
  cfg.parameter='avg';
%   tlock_tacPaud_s0{soalist(ll)}=ft_math(cfg,ft_timelockanalysis([],tlock_tac_s0{10}),ft_timelockanalysis([],tlock_aud_s0{soalist(ll)+30})); % Talone + Aalone_shifted
%   tlock_tacPaud_sall{soalist(ll)}=ft_math(cfg,ft_timelockanalysis([],tlock_tac_sall{10}),ft_timelockanalysis([],tlock_aud_sall{soalist(ll)+30})); % Talone + Aalone_shifted
  tlock_tacPaud_s0{soalist(ll)}=ft_math(cfg,tlock_tac_s0{10},tlock_aud_s0avg{soalist(ll)+30}); % Talone + Aalone_shifted
  tlock_tacPaud_sall{soalist(ll)}=ft_math(cfg,tlock_tac_sall{10},tlock_aud_sallavg{soalist(ll)+30}); % Talone + Aalone_shifted
  
  cfg=[];
  cfg.operation='add';
  cfg.parameter='avg';
%   tlock_audPtac_s0{soalist(ll)}=ft_math(cfg,ft_timelockanalysis([],tlock_aud_s0{10}),ft_timelockanalysis([],tlock_tac_s0{soalist(ll)+30})); % Aalone + Talone_shifted
%   tlock_audPtac_sall{soalist(ll)}=ft_math(cfg,ft_timelockanalysis([],tlock_aud_sall{10}),ft_timelockanalysis([],tlock_tac_sall{soalist(ll)+30})); % Aalone + Talone_shifted
  tlock_audPtac_s0{soalist(ll)}=ft_math(cfg,tlock_aud_s0{10},tlock_tac_s0avg{soalist(ll)+30}); % Aalone + Talone_shifted
  tlock_audPtac_sall{soalist(ll)}=ft_math(cfg,tlock_aud_sall{10},tlock_tac_sallavg{soalist(ll)+30}); % Aalone + Talone_shifted
end

if plotflag
  try close(9);end
  figure(9);
  subplot(2,3,1);plot(tlock_tac_sall{10}.time,tlock_tac_sall{10}.avg(chanuse,:),'k');
  hold on;plot(tlock_aud_sall{10}.time,tlock_aud_sall{10}.avg(chanuse,:),'b');
  hold on;plot(tlock_tacPaud_sall{5}.time,tlock_tacPaud_sall{5}.avg(chanuse,:),'g');
  hold on;plot(tlock_audPtac_sall{5}.time,tlock_audPtac_sall{5}.avg(chanuse,:),'r');axis([-.2 .5 -inf inf])
  legend('tactile alone','auditory alone','sum unisensory, tactile at 0','sum unisensory, auditory at 0')
  
  
  for ll=3:7
    subplot(2,3,ll-1);plot(tlock_tac_sall{ll}.time,tlock_tac_sall{ll}.avg(chanuse,:),'k');
    hold on;plot(tlock_aud_sall{ll}.time,tlock_aud_sall{ll}.avg(chanuse,:),'b');
    hold on;plot(tlock_tacPaud_sall{ll}.time,tlock_tacPaud_sall{ll}.avg(chanuse,:),'g');
    hold on;plot(tlock_audPtac_sall{ll}.time,tlock_audPtac_sall{ll}.avg(chanuse,:),'r');axis([-.2 .5 -inf inf])
  legend('multisensory, tactile at 0','multisensory, auditory at 0','sum unisensory, tactile at 0','sum unisensory, auditory at 0')
  end
  
end


try close(10);end
try close(11);end
for ll=soalist
  cfg=[];
  cfg.operation='subtract';
  cfg.parameter='avg';
%   tlock_tac_s0{ll+10}=ft_math(cfg,ft_timelockanalysis([],tlock_tac_s0{ll}),ft_timelockanalysis([],tlock_tac_s0{10})); % TA-T = 'aud' isolated from multisensory
%   tlock_tac_sall{ll+10}=ft_math(cfg,ft_timelockanalysis([],tlock_tac_sall{ll}),ft_timelockanalysis([],tlock_tac_sall{10})); % TA-T = 'aud' isolated from multisensory
%   tlock_aud_s0{ll+10}=ft_math(cfg,ft_timelockanalysis([],tlock_aud_s0{ll}),ft_timelockanalysis([],tlock_aud_s0{10})); % AT-A = 'tac' isolated from mulitsensory
%   tlock_aud_sall{ll+10}=ft_math(cfg,ft_timelockanalysis([],tlock_aud_sall{ll}),ft_timelockanalysis([],tlock_aud_sall{10})); % AT-A = 'tac' isolated from mulitsensory
  tlock_tac_s0{ll+10}=ft_math(cfg,tlock_tac_s0{ll},tlock_tac_s0{10}); % TA-T = 'aud' isolated from multisensory
  tlock_tac_sall{ll+10}=ft_math(cfg,tlock_tac_sall{ll},tlock_tac_sall{10}); % TA-T = 'aud' isolated from multisensory
  tlock_aud_s0{ll+10}=ft_math(cfg,tlock_aud_s0{ll},tlock_aud_s0{10}); % AT-A = 'tac' isolated from mulitsensory
  tlock_aud_sall{ll+10}=ft_math(cfg,tlock_aud_sall{ll},tlock_aud_sall{10}); % AT-A = 'tac' isolated from mulitsensory
  
  if plotflag
    figure(10); 
    subplot(2,5,ll-2);plot(tlock_tac_sallavg{30+ll}.time,tlock_tac_sallavg{30+ll}.avg(chanuse,:),'k');
    hold on;plot(tlock_aud_sall{10+ll}.time,tlock_aud_sall{10+ll}.avg(chanuse,:),'b');axis([-.2 .5 -inf inf])
    legend('tactile-alone shifted','AT-A: residual tactile')
    subplot(2,5,ll-2+5);plot(tlock_aud_sallavg{30+ll}.time,tlock_aud_sallavg{30+ll}.avg(chanuse,:),'k');
    hold on;plot(tlock_tac_sall{10+ll}.time,tlock_tac_sall{10+ll}.avg(chanuse,:),'b');axis([-.2 .5 -inf inf])
    legend('auditory-alone shifted','TA-T: residual auditory')
  end
  
  cfg=[];
  cfg.operation='add';
  cfg.parameter='avg';
  tlock_tac_s0{ll+20}=ft_math(cfg,tlock_tac_s0{ll},tlock_nul_s0{10}); % TA+N
  tlock_tac_sall{ll+20}=ft_math(cfg,tlock_tac_sall{ll},tlock_nul_sall{10}); % TA+N
  tlock_aud_s0{ll+20}=ft_math(cfg,tlock_aud_s0{ll},tlock_nul_s0{10}); % AT+N
  tlock_aud_sall{ll+20}=ft_math(cfg,tlock_aud_sall{ll},tlock_nul_sall{10}); % AT+N
  
  cfg=[];
  cfg.operation='subtract';
  cfg.parameter='avg';
%   tlock_tpa_mtamn_s0{ll}=ft_math(cfg,ft_timelockanalysis([],tlock_tacPaud_s0{ll}),ft_timelockanalysis([],tlock_tac_s0{ll+20})); % (T + As) - (TA + N)
%   tlock_tpa_mtamn_sall{ll}=ft_math(cfg,ft_timelockanalysis([],tlock_tacPaud_sall{ll}),ft_timelockanalysis([],tlock_tac_sall{ll+20})); % (T + As) - (TA + N)
%   tlock_apt_matmn_s0{ll}=ft_math(cfg,ft_timelockanalysis([],tlock_audPtac_s0{ll}),ft_timelockanalysis([],tlock_aud_s0{ll+20})); % (A + Ts) - (AT + N)
%   tlock_apt_matmn_sall{ll}=ft_math(cfg,ft_timelockanalysis([],tlock_audPtac_sall{ll}),ft_timelockanalysis([],tlock_aud_sall{ll+20})); % (A + Ts) - (AT + N)
  tlock_tpa_mtamn_s0{ll}=ft_math(cfg,tlock_tacPaud_s0{ll},tlock_tac_s0{ll+20}); % (T + As) - (TA + N)
  tlock_tpa_mtamn_sall{ll}=ft_math(cfg,tlock_tacPaud_sall{ll},tlock_tac_sall{ll+20}); % (T + As) - (TA + N)
  tlock_apt_matmn_s0{ll}=ft_math(cfg,tlock_audPtac_s0{ll},tlock_aud_s0{ll+20}); % (A + Ts) - (AT + N)
  tlock_apt_matmn_sall{ll}=ft_math(cfg,tlock_audPtac_sall{ll},tlock_aud_sall{ll+20}); % (A + Ts) - (AT + N)
    
  if plotflag
    figure(11); % (Sum of Unisensory) minus (Multisensory plus Nul)
    subplot(2,5,ll-2);plot(tlock_tpa_mtamn_sall{ll}.time,tlock_tpa_mtamn_sall{ll}.avg(chanuse,:),'k');axis([-.2 .5 -7 7]);
    legend('(T+A)-(TA-N), time0 is tactile')
    subplot(2,5,ll-2+5);plot(tlock_apt_matmn_sall{ll}.time,tlock_apt_matmn_sall{ll}.avg(chanuse,:),'k');axis([-.2 .5 -7 7])
    legend('(A+T)-(AT-N), time0 is auditory')
  end
  
end

save(['tlock_diffs_' sub{ii} '.mat'],'tlock_*_m*')

% figure;
% channel={'Fz','Cz'};
% for ll=1:length(soalist)
%   cfg=[];
%   cfg.channel=channel;
%   cfg.ylim=[-6 6];
%   cfg.xlim=[-.1 .4];
%   cfg.layout='EEG1010.lay';
%   
%   subplot(2,5,ll)
%   ft_singleplotER(cfg,tlock_tpa_mtamn_s0{soalist(ll)})
%   ft_singleplotER(cfg,tlock_tpa_mtamn_sall{soalist(ll)})
%   subplot(2,5,ll+5)
%   ft_singleplotER(cfg,tlock_apt_matmn_s0{soalist(ll)})
%   ft_singleplotER(cfg,tlock_apt_matmn_sall{soalist(ll)})
% end
% figure;
% channel={'Fz'};
% for ll=1:length(soalist)
%   cfg=[];
%   cfg.channel=channel;
%   cfg.ylim=[-4 4];
%   cfg.xlim=[-.1 .4];
%   cfg.layout='EEG1010.lay';
%   
%   subplot(2,5,ll)
%   ft_singleplotER(cfg,tlock_tpa_mtamn_s0{soalist(ll)})
%   ft_singleplotER(cfg,tlock_tpa_mtamn_sall{soalist(ll)})
%   subplot(2,5,ll+5)
%   ft_singleplotER(cfg,tlock_apt_matmn_s0{soalist(ll)})
%   ft_singleplotER(cfg,tlock_apt_matmn_sall{soalist(ll)})
% end

fzpeaks0=nan(2,soalist(end));
fzpeaksall=nan(2,soalist(end));
for ll=soalist
  fzpeaks0(1,ll)=tlock_tpa_mtamn_s0{ll}.avg(match_str(tlock_tac_s0{ll}.label,'Fz'),dsearchn(tlock_tpa_mtamn_s0{ll}.time',.18));
  fzpeaks0(2,ll)=tlock_apt_matmn_s0{ll}.avg(match_str(tlock_tac_s0{ll}.label,'Fz'),dsearchn(tlock_apt_matmn_s0{ll}.time',.18));
  fzpeaksall(1,ll)=tlock_tpa_mtamn_sall{ll}.avg(match_str(tlock_tac_sall{ll}.label,'Fz'),dsearchn(tlock_tpa_mtamn_sall{ll}.time',.18));
  fzpeaksall(2,ll)=tlock_apt_matmn_sall{ll}.avg(match_str(tlock_tac_sall{ll}.label,'Fz'),dsearchn(tlock_apt_matmn_sall{ll}.time',.18));
end

if 0
  cfg=[];
  cfg.interactive='yes';
  cfg.layout='EEG1010.lay';
  ft_multiplotER(cfg,tlock_tac_sall{15});
  ft_multiplotER(cfg,tlock_tac{chanuse});
  ft_multiplotER(cfg,tlock_tpa_mtamn_sall{5});
  ft_multiplotER(cfg,tlock_tpa_matmn{5});
end


%% Statistics

% contrast of conditions in GLM: Aalone + Talone_shifted - ATmult - N

load eeg1010_neighb
clear stat*
cfg=[];
cfg.latency=[0 .3];
% cfg.parameter='avg';
% cfg.method='montecarlo';
cfg.method='analytic';
% cfg.statistic='indepsamplesregrT';
cfg.statistic='indepsamplesT';
cfg.ivar=1;
% cfg.wvar=2;
cfg.cvar=2;
cfg.correctm='holm';
% cfg.correctm='cluster';
cfg.numrandomization=100;
cfg.neighbours=neighbours;


for ll=soalist
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
    
