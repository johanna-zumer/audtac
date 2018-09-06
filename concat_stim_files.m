function concat_stim_files(ii,plotflag)
% concatenate stim files

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
sub{101}='e01'; 
sub{102}='e02'; 
sub{103}='e03'; 
sub{104}='e04'; 
%%%  above is pilot, below is real
sub{105}='e05'; 
sub{106}='e06'; 
sub{107}='e07'; 
sub{108}='e08'; 
sub{109}='e09';
sub{110}='e10';
sub{111}='e11';
sub{112}='e12';
sub{113}='e13';
sub{114}='e14';
sub{115}='e15';
sub{116}='e16';
sub{117}='e17';
sub{118}='e18';
sub{119}='e19';
sub{120}='e20';
sub{121}='e21';
sub{122}='e22';
sub{123}='e23';
sub{124}='e24';
sub{125}='e25';
sub{126}='e26';
sub{127}='e27';
sub{128}='e28';
sub{129}='e29';
sub{130}='e30';
sub{131}='e31';
sub{132}='e32';


% ii=103;
if ~exist('plotflag')
  plotflag=0;
end

% thresh for max(:)-min(:) within one trial (if too small, then probably not touch)
maxminthresh=140;
medthresh=.03; % +/- 30ms from median of data.

cd([ddir ])
%%

% note, this works for ii>101

sfiles=dir([ddir sub{ii} '_*s.mat']);
bfiles=dir([ddir sub{ii} '_*b.mat']);
rfiles=dir([ddir sub{ii} '_*r.mat']);
dosfiles=0;
dobfiles=0;

stime_touch=[];
btime_touch=[];
rtime_touch=[];
saudseq=[];
baudseq=[];
raudseq=[];
snulseq=[];
bnulseq=[];
rnulseq=[];

snumtac=[];
snumaud=[];
snumnul=[];
snumtr=[];
bnumtac=[];
bnumaud=[];
bnumnul=[];
bnumtr=[];
rnumtac=[];
rnumaud=[];
rnumnul=[];
rnumtr=[];

ssoaeff=[];
bsoaeff=[];
rsoaeff=[];
ssoaseq=[];
bsoaseq=[];
rsoaseq=[];

rrts=[]; % RTs for response-condition only

figind=1;

if dosfiles
for ll=1:length(sfiles)
  sstim{ll}=load(sfiles(ll).name);
  
  % get time_touch0, the time of touch for all trials (=0 for non-tactile trials)
  if ii<103
    tactrials=find(sstim{ll}.info.time_touch>0);
    mint5=[];mint8=[];
    for tt=1:length(tactrials)
      indx5=find(sstim{ll}.info.lightsensor{tactrials(tt)}>.5*[max(sstim{ll}.info.lightsensor{tactrials(tt)})-min(sstim{ll}.info.lightsensor{tactrials(tt)})]+min(sstim{ll}.info.lightsensor{tactrials(tt)}));
      indx8=find(sstim{ll}.info.lightsensor{tactrials(tt)}>.8*[max(sstim{ll}.info.lightsensor{tactrials(tt)})-min(sstim{ll}.info.lightsensor{tactrials(tt)})]+min(sstim{ll}.info.lightsensor{tactrials(tt)}));
      try
        mint5(tt)=indx5(1);
        mint8(tt)=indx8(1);
      end
    end
    % turns out that 50% is better criterion (change it in main code) but then
    % for timing to be consistent, add back 4ms.
    time_touch=.001*(mint5+4-1); % add 4ms for diff of 50% vs 80%, and subtract 1 for the time relative to info.lighttime{tt}(1)
    if median(time_touch-sstim{ll}.info.time_touch(tactrials))>.001
      error('time_touch not right?')
    end
    origtac=find(sstim{ll}.info.valves_seq);
    time_touch0(origtac(1:length(time_touch)))=time_touch;
    if length(time_touch0)>length(sstim{ll}.info.lightsensor)
      time_touch0=time_touch0(1:length(sstim{ll}.info.lightsensor));
    end
    clear tactrials
  else
    time_touch0=sstim{ll}.info.time_touch(1:length(sstim{ll}.info.lightsensor));
  end
  
  tactr=find(time_touch0>0); % tactile trials
  time_touch1=time_touch0(tactr);
  
  if plotflag
    figure(figind);subplot(2,1,1);plot(time_touch1);axis([-inf inf .17 .25])
    figind=figind+1;
  end
  
  tacrej=find(abs(time_touch1-median(time_touch1))>medthresh); % reject if deviates from median by more than 20ms (good for catching sudden movements but bad for slow drifts)
  time_touch1(tacrej)=nan; % reject those
  tacrej=[tacrej find(abs(diff(time_touch1))>.01)+1]; % additionally reject those with change of 10ms from one trial to next
  
  % if the deviation of the light plot from its baseline is not enough, probably means it wasn't touching the skin
  for nn=1:length(tactr)
    % Note, that the threshold of 120 is arbitrary but works for subject 06
    if max(sstim{ll}.info.lightsensor{tactr(nn)})-min(sstim{ll}.info.lightsensor{tactr(nn)})<maxminthresh
      tacrej=[tacrej nn];
    end
  end
  tacrej=unique(tacrej);
  
  time_touch1(tacrej)=nan;
  time_touch0(tactr(tacrej))=nan;
  
  cfg=[];
  cfg.dataset=[edir sub{ii} '/' sfiles(ll).name(1:end-3) 'eeg'];
  if ii==108 & ll==2
    cfg.dataset=[edir sub{ii} '/' sfiles(ll).name(1:4) '2s.eeg'];
  end
  cfg.trialfun='ft_trialfun_general';
  cfg.trialdef.eventtype  = 'Stimulus';
  cfg.trialdef.eventvalue = 'S  2'; % corresponds to lj.prepareStrobe(2)
  cfg.trialdef.prestim = 0.5;
  cfg.trialdef.poststim = 2.3;
  cfgtr=ft_definetrial(cfg);
  
  % If number of trials in the data is greater than in the sstim{ll}.info, then have
  % remainder trials of info just be NaN and they will be excluded later.
  % If number of trials in the info is greater than in the data, then cut
  % it out now.
  
  numtacdat=size(cfgtr.trl,1);
  numtactap=length(time_touch1);
  snumtac=[snumtac numtacdat];
  time_touch2=nan(1,numtacdat); % ensure that length is at least numtacdat.
  time_touch2(1:numtactap)=time_touch1; %if this is shorter than numtacdat, than remainder is Nan; else it's too long
  time_touch2=time_touch2(1:numtacdat);  % so then cut it short just in case
  stime_touch=[stime_touch time_touch2]; % concatenate over blocks
  
  audseq=sstim{ll}.info.auditory_seq(1,1:length(sstim{ll}.info.lightsensor));
  audtr=find(audseq);
  cfg=[];
  cfg.dataset=[edir sub{ii} '/' sfiles(ll).name(1:end-3) 'eeg'];
  if ii==108 & ll==2
    cfg.dataset=[edir sub{ii} '/' sfiles(ll).name(1:4) '2s.eeg'];
  end
  cfg.trialfun='ft_trialfun_general';
  cfg.trialdef.eventtype  = 'Stimulus';
  cfg.trialdef.eventvalue = 'S 10'; % corresponds to lj.prepareStrobe(10)
  cfg.trialdef.prestim = 0.7;
  cfg.trialdef.poststim = 2.1;
  cfgtr=ft_definetrial(cfg);
  
  numauddat=size(cfgtr.trl,1);
  numaudtap=length(audtr);
  snumaud=[snumaud numauddat];
  audsequse=nan(1,numauddat);
  audsequse(audtr(1:numaudtap))=1;
  audsequse=audsequse(1:numauddat);
  saudseq=[saudseq audsequse];
  
  nultr=find(~isnan(sstim{ll}.info.null_time));
  cfg=[];
  cfg.dataset=[edir sub{ii} '/' sfiles(ll).name(1:end-3) 'eeg'];
  if ii==108 & ll==2
    cfg.dataset=[edir sub{ii} '/' sfiles(ll).name(1:4) '2s.eeg'];
  end
  cfg.trialfun='ft_trialfun_general';
  cfg.trialdef.eventtype  = 'Stimulus';
  cfg.trialdef.eventvalue = 'S  1'; % corresponds to lj.prepareStrobe(1)
  cfg.trialdef.prestim = 0.7;
  cfg.trialdef.poststim = 2.1;
  cfgtr=ft_definetrial(cfg);
  
  numnuldat=size(cfgtr.trl,1);
  numnultap=length(nultr);
  snumnul=[snumnul numnuldat];
  nulsequse=nan(1,numnuldat);
  nulsequse(nultr(1:numnultap))=1;
  nulsequse=nulsequse(1:numnuldat);
  snulseq=[snulseq nulsequse];
  
  % total trials
  cfg=[];
  cfg.dataset=[edir sub{ii} '/' sfiles(ll).name(1:end-3) 'eeg'];
  if ii==108 & ll==2
    cfg.dataset=[edir sub{ii} '/' sfiles(ll).name(1:4) '2s.eeg'];
  end
  cfg.trialfun='ft_trialfun_general';
  cfg.trialdef.eventtype  = 'Stimulus';
  cfg.trialdef.eventvalue = 'S 20';
  cfg.trialdef.prestim = 0.7;
  cfg.trialdef.poststim = 2.1;
  cfgtr=ft_definetrial(cfg);
  
  
  cfg.trialdef.prestim = 1;
  cfg.trialdef.poststim = 1;
  cfgtr1=ft_definetrial(cfg);
  
%   keyboard % this commented out 4 lines used to work...why not now?
%   cfg.trialfun='ft_trialfun_general_sleep';
%   cfg.trialdef.prestim = 1;
%   cfg.trialdef.poststim = 1;
%   cfgtr2=ft_definetrial(cfg);
  
  
  numtrtap=length(time_touch0);
  numtrdat=size(cfgtr.trl,1);
  snumtr=[snumtr numtrdat];
  soaeff=nan(1,numtrdat);
  %   soaeff(1:numtrtap)=sstim{ll}.info.soa_effective2(1:numtrtap);
  soaeff(1:numtrtap)=-(sstim{ll}.info.valve_time - sstim{ll}.info.audio_time + time_touch0 - sstim{ll}.fin.audDelay/1000); %using time_touch0, we have eliminated bad tactile trials if photodiode messed up.
  soaeff=soaeff(1:numtrdat);
  ssoaeff=[ssoaeff soaeff];
  soaseq=nan(1,numtrdat);
  soaseq(1:numtrtap)=sstim{ll}.info.soa_seq(1:numtrtap);
  soaseq=soaseq(1:numtrdat);
  ssoaseq=[ssoaseq soaseq];
  %   soaeff=sstim{ll}.info.soa_effective2;
  %   ssoaeff=[ssoaeff soaeff];
  %   soaseq=sstim{ll}.info.soa_seq;
  %   ssoaseq=[ssoaseq soaseq];
  
  numplots=16; % should be a square number
  if plotflag
    figure(figind);
    figind=figind+1;
    for nn=1:numplots
      ind=round((numtrdat/numplots)*nn);
      %     ind=min(100,numtrdat); % no point in starting at 100 if number of trials is much fewer!
      while ind>length(sstim{ll}.info.lightsensor) | sstim{ll}.info.valves_seq(ind)~=1
        ind=ind-1;
      end
      subplot(sqrt(numplots),sqrt(numplots),nn);    plot(sstim{ll}.info.lighttime{ind},sstim{ll}.info.lightsensor{ind})
      axis([-inf inf 150 800])
    end
    %     figure;plot(sstim{ll}.info.time_touch(sstim{ll}.info.time_touch>0));axis([-inf inf .17 .25])
    figure(figind-2);subplot(2,1,2);plot(time_touch1);axis([-inf inf .17 .25])
  end
end
end

%% Sleep data
if dobfiles
for ll=1:length(bfiles)
  bstim{ll}=load(bfiles(ll).name);
  
  if ii<103
    tactrials=find(bstim{ll}.info.time_touch>0);
    mint5=[];mint8=[];
    for tt=1:length(tactrials)
      indx5=find(bstim{ll}.info.lightsensor{tactrials(tt)}>.5*[max(bstim{ll}.info.lightsensor{tactrials(tt)})-min(bstim{ll}.info.lightsensor{tactrials(tt)})]+min(bstim{ll}.info.lightsensor{tactrials(tt)}));
      indx8=find(bstim{ll}.info.lightsensor{tactrials(tt)}>.8*[max(bstim{ll}.info.lightsensor{tactrials(tt)})-min(bstim{ll}.info.lightsensor{tactrials(tt)})]+min(bstim{ll}.info.lightsensor{tactrials(tt)}));
      try
        mint5(tt)=indx5(1);
        mint8(tt)=indx8(1);
      end
    end
    % turns out that 50% is better criterion (change it in main code) but then
    % for timing to be consistent, add back 4ms.
    time_touch=.001*(mint5+4-1); % add 4ms for diff of 50% vs 80%, and subtract 1 for the time relative to info.lighttime{tt}(1)
    if median(time_touch-bstim{ll}.info.time_touch(tactrials))>.001
      error('time_touch not right?')
    end
    origtac=find(bstim{ll}.info.valves_seq);
    time_touch0(origtac(1:length(time_touch)))=time_touch;
    if length(time_touch0)>length(bstim{ll}.info.lightsensor)
      time_touch0=time_touch0(1:length(bstim{ll}.info.lightsensor));
    end
    clear tactrials
  else
    time_touch0=bstim{ll}.info.time_touch(1:length(bstim{ll}.info.lightsensor));
  end
  
  if ii==107 && ll==3
    keyboard % this doesn't seem to match with later... figure out why if we care about this participant 7
    time_touch0=time_touch0(1:79);
  end
  if ii==107 && ll==4
    time_touch0=time_touch0(1:25);
  end
  
  %   time_touch1=time_touch0(time_touch0>0);
  
  tactr=find(time_touch0>0);
  time_touch1=time_touch0(tactr);
  
  if plotflag
    figure(figind);subplot(2,1,1);plot(time_touch1);axis([-inf inf .17 .25])
    figind=figind+1;
  end
  
  tacrej=find(abs(time_touch1-median(time_touch1))>medthresh);
  time_touch1(tacrej)=nan;
  tacrej=[tacrej find(abs(diff(time_touch1))>.01)+1];
  
  % if the deviation of the light plot from its baseline is not enough, probably means it wasn't touching the skin
  for nn=1:length(tactr)
    if max(bstim{ll}.info.lightsensor{tactr(nn)})-min(bstim{ll}.info.lightsensor{tactr(nn)})<maxminthresh % the number 110 was chosen after seeing several subjects
      tacrej=[tacrej nn];
    end
  end
  tacrej=unique(tacrej);
  
  time_touch1(tacrej)=nan;
  time_touch0(tactr(tacrej))=nan;
  
  cfg=[];
  cfg.dataset=[edir sub{ii} '/' bfiles(ll).name(1:end-3) 'eeg'];
  cfg.trialfun='ft_trialfun_general';
  cfg.trialdef.eventtype  = 'Stimulus';
  cfg.trialdef.eventvalue = 'S  2'; % corresponds to lj.prepareStrobe(2)
  cfg.trialdef.prestim = 0.5;
  cfg.trialdef.poststim = 2.3;
  try
    cfgtr=ft_definetrial(cfg);
    
    numtacdat=size(cfgtr.trl,1);
    numtactap=length(time_touch1);
    bnumtac=[bnumtac numtacdat];
    time_touch2=nan(1,numtacdat); % ensure that length is at least numtacdat.
    time_touch2(1:numtactap)=time_touch1; %if this is shorter than numtacdat, than remainder is Nan; else it's too long
    time_touch2=time_touch2(1:numtacdat);  % so then cut it short just in case
    btime_touch=[btime_touch time_touch2]; % concatenate over blocks
    %     numtac=min(size(cfgtr.trl,1),length(time_touch1));
    %     btime_touch=[btime_touch time_touch1(1:numtac)];
    %     bnumtac=[bnumtac length(time_touch1)];
    
    audseq=bstim{ll}.info.auditory_seq(1,1:length(bstim{ll}.info.lightsensor));
    audtr=find(audseq);
    cfg=[];
    cfg.dataset=[edir sub{ii} '/' bfiles(ll).name(1:end-3) 'eeg'];
    cfg.trialfun='ft_trialfun_general';
    cfg.trialdef.eventtype  = 'Stimulus';
    cfg.trialdef.eventvalue = 'S 10'; % corresponds to lj.prepareStrobe(10)
    cfg.trialdef.prestim = 0.7;
    cfg.trialdef.poststim = 2.1;
    cfgtr=ft_definetrial(cfg);
    
    numauddat=size(cfgtr.trl,1);
    numaudtap=length(audtr);
    bnumaud=[bnumaud numauddat];
    audsequse=nan(1,numauddat);
    audsequse(audtr(1:numaudtap))=1;
    audsequse=audsequse(1:numauddat);
    baudseq=[baudseq audsequse];
    %     numaud=min(size(cfgtr.trl,1),length(audtr));
    %     audsequse=[];
    %     audsequse(audtr(1:numaud))=1;
    %     baudseq=[baudseq audsequse];
    %     bnumaud=[bnumaud length(audtr)];
    
    nultr=find(~isnan(bstim{ll}.info.null_time));
    cfg=[];
    cfg.dataset=[edir sub{ii} '/' bfiles(ll).name(1:end-3) 'eeg'];
    cfg.trialfun='ft_trialfun_general';
    cfg.trialdef.eventtype  = 'Stimulus';
    cfg.trialdef.eventvalue = 'S  1'; % corresponds to lj.prepareStrobe(1)
    cfg.trialdef.prestim = 0.7;
    cfg.trialdef.poststim = 2.1;
    cfgtr=ft_definetrial(cfg);
    
    numnuldat=size(cfgtr.trl,1);
    numnultap=length(nultr);
    bnumnul=[bnumnul numnuldat];
    nulsequse=nan(1,numnuldat);
    nulsequse(nultr(1:numnultap))=1;
    nulsequse=nulsequse(1:numnuldat);
    bnulseq=[bnulseq nulsequse];
    %     numnul=min(size(cfgtr.trl,1),length(nultr));
    %     nulsequse=[];
    %     nulsequse(nultr(1:numnul))=1;
    %     bnulseq=[bnulseq nulsequse];
    %     bnumnul=[bnumnul length(nultr)];
    
    % total number of trials is:
    cfg=[];
    cfg.dataset=[edir sub{ii} '/' bfiles(ll).name(1:end-3) 'eeg'];
    cfg.trialfun='ft_trialfun_general';
    cfg.trialdef.eventtype  = 'Stimulus';
    cfg.trialdef.eventvalue = 'S 20'; % corresponds to lj.prepareStrobe(1)
    cfg.trialdef.prestim = 0.7;
    cfg.trialdef.poststim = 2.1;
    cfgtr=ft_definetrial(cfg);
    
    numtrtap=length(time_touch0);
    numtrdat=size(cfgtr.trl,1);
    bnumtr=[bnumtr numtrdat];
    soaeff=nan(1,numtrdat);
    %   soaeff(1:numtrtap)=bstim{ll}.info.soa_effective2(1:numtrtap);
    soaeff(1:numtrtap)=-(bstim{ll}.info.valve_time - bstim{ll}.info.audio_time + time_touch0 - bstim{ll}.fin.audDelay/1000); %using time_touch0, we have eliminated bad tactile trials if photodiode messed up.
    soaeff=soaeff(1:numtrdat);
    bsoaeff=[bsoaeff soaeff];
    soaseq=nan(1,numtrdat);
    soaseq(1:numtrtap)=bstim{ll}.info.soa_seq(1:numtrtap);
    soaseq=soaseq(1:numtrdat);
    bsoaseq=[bsoaseq soaseq];
    %     numtr=length(time_touch0);
    %     soaeff=bstim{ll}.info.soa_effective2(1:numtr);
    %     bsoaeff=[bsoaeff soaeff];
    %     soaseq=bstim{ll}.info.soa_seq(1:numtr);
    %     bsoaseq=[bsoaseq soaseq];
    
    %   soaeff=bstim{ll}.info.soa_effective2;
    % %     soaeff(isnan(time_touch0))=nan;
    %     bsoaeff=[bsoaeff soaeff];
    %     soaseq=bstim{ll}.info.soa_seq;
    %     bsoaseq=[bsoaseq soaseq];
    
    numplots=16; % should be a square number
    if plotflag
      figure(figind);
      figind=figind+1;
      for nn=1:numplots
        ind=round((numtrdat/numplots)*nn);
        %     ind=min(100,numtrdat); % no point in starting at 100 if number of trials is much fewer!
        while ind>length(bstim{ll}.info.lightsensor) | bstim{ll}.info.valves_seq(ind)~=1
          ind=ind-1;
        end
        subplot(sqrt(numplots),sqrt(numplots),nn);    plot(bstim{ll}.info.lighttime{ind},bstim{ll}.info.lightsensor{ind})
        axis([-inf inf 150 800])
      end
      %     figure;plot(sstim{ll}.info.time_touch(sstim{ll}.info.time_touch>0));axis([-inf inf .17 .25])
      figure(figind-2);subplot(2,1,2);plot(time_touch1);axis([-inf inf .17 .25])
    end
    
  catch
    disp(cfg.dataset) % this sometimes happened when I restarted the stim but left the EEG acqu running
  end
end
end

%% Sitting responding no-EEG data
for ll=1:length(rfiles)
  rstim{ll}=load(rfiles(ll).name);
  
  nt=length(rstim{ll}.info.lightsensor);
  if ii==117 && ll==2
    nt=148; % she stopped to ask a question
  end
  
  if nt<5
    continue
  end
  
  % get time_touch0, the time of touch for all trials (=0 for non-tactile trials)
%   if ii<103
%     tactrials=find(rstim{ll}.info.time_touch>0);
%     mint5=[];mint8=[];
%     for tt=1:length(tactrials)
%       indx5=find(sstim{ll}.info.lightsensor{tactrials(tt)}>.5*[max(sstim{ll}.info.lightsensor{tactrials(tt)})-min(sstim{ll}.info.lightsensor{tactrials(tt)})]+min(sstim{ll}.info.lightsensor{tactrials(tt)}));
%       indx8=find(sstim{ll}.info.lightsensor{tactrials(tt)}>.8*[max(sstim{ll}.info.lightsensor{tactrials(tt)})-min(sstim{ll}.info.lightsensor{tactrials(tt)})]+min(sstim{ll}.info.lightsensor{tactrials(tt)}));
%       try
%         mint5(tt)=indx5(1);
%         mint8(tt)=indx8(1);
%       end
%     end
%     % turns out that 50% is better criterion (change it in main code) but then
%     % for timing to be consistent, add back 4ms.
%     time_touch=.001*(mint5+4-1); % add 4ms for diff of 50% vs 80%, and subtract 1 for the time relative to info.lighttime{tt}(1)
%     if median(time_touch-sstim{ll}.info.time_touch(tactrials))>.001
%       error('time_touch not right?')
%     end
%     origtac=find(sstim{ll}.info.valves_seq);
%     time_touch0(origtac(1:length(time_touch)))=time_touch;
%     if length(time_touch0)>length(sstim{ll}.info.lightsensor)
%       time_touch0=time_touch0(1:length(sstim{ll}.info.lightsensor));
%     end
%     clear tactrials
%   else
    time_touch0=rstim{ll}.info.time_touch(1:nt);
%   end
  
  tactr=find(time_touch0>0); % tactile trials
  time_touch1=time_touch0(tactr);
  
  if plotflag
    figure(figind);subplot(2,1,1);plot(time_touch1);axis([-inf inf .17 .25])
    figind=figind+1;
  end
  
  tacrej=find(abs(time_touch1-median(time_touch1))>medthresh); % reject if deviates from median by more than 20ms (good for catching sudden movements but bad for slow drifts)
  time_touch1(tacrej)=nan; % reject those
  tacrej=[tacrej find(abs(diff(time_touch1))>.01)+1]; % additionally reject those with change of 10ms from one trial to next
  
  % if the deviation of the light plot from its baseline is not enough, probably means it wasn't touching the skin
  for nn=1:length(tactr)
    % Note, that the threshold of 120 is arbitrary but works for subject 06
    if max(rstim{ll}.info.lightsensor{tactr(nn)})-min(rstim{ll}.info.lightsensor{tactr(nn)})<maxminthresh
      tacrej=[tacrej nn];
    end
  end
  tacrej=unique(tacrej);
  
  time_touch1(tacrej)=nan;
  time_touch0(tactr(tacrej))=nan;
  
%   cfg=[];
%   cfg.dataset=[edir sub{ii} '/' sfiles(ll).name(1:end-3) 'eeg'];
%   if ii==108 & ll==2
%     cfg.dataset=[edir sub{ii} '/' sfiles(ll).name(1:4) '2s.eeg'];
%   end
%   cfg.trialfun='ft_trialfun_general';
%   cfg.trialdef.eventtype  = 'Stimulus';
%   cfg.trialdef.eventvalue = 'S  2'; % corresponds to lj.prepareStrobe(2)
%   cfg.trialdef.prestim = 0.5;
%   cfg.trialdef.poststim = 2.3;
%   cfgtr=ft_definetrial(cfg);
%   
%   % If number of trials in the data is greater than in the sstim{ll}.info, then have
%   % remainder trials of info just be NaN and they will be excluded later.
%   % If number of trials in the info is greater than in the data, then cut
%   % it out now.
%   
%   numtacdat=size(cfgtr.trl,1);
  numtactap=length(time_touch1);
%   snumtac=[snumtac numtacdat];
%   time_touch2=nan(1,numtacdat); % ensure that length is at least numtacdat.
  time_touch2(1:numtactap)=time_touch1; %if this is shorter than numtacdat, than remainder is Nan; else it's too long
%   time_touch2=time_touch2(1:numtacdat);  % so then cut it short just in case
  rtime_touch=[rtime_touch time_touch2]; % concatenate over blocks
  
  audseq=rstim{ll}.info.auditory_seq(1,1:nt);
  audtr=find(audseq);
%   cfg=[];
%   cfg.dataset=[edir sub{ii} '/' sfiles(ll).name(1:end-3) 'eeg'];
%   if ii==108 & ll==2
%     cfg.dataset=[edir sub{ii} '/' sfiles(ll).name(1:4) '2s.eeg'];
%   end
%   cfg.trialfun='ft_trialfun_general';
%   cfg.trialdef.eventtype  = 'Stimulus';
%   cfg.trialdef.eventvalue = 'S 10'; % corresponds to lj.prepareStrobe(10)
%   cfg.trialdef.prestim = 0.7;
%   cfg.trialdef.poststim = 2.1;
%   cfgtr=ft_definetrial(cfg);
%   
%   numauddat=size(cfgtr.trl,1);
  numaudtap=length(audtr);
%   snumaud=[snumaud numauddat];
%   audsequse=nan(1,numauddat);
  audsequse(audtr(1:numaudtap))=1;
%   audsequse=audsequse(1:numauddat);
  raudseq=[raudseq audsequse];
  
  nultr=find(~isnan(rstim{ll}.info.null_time(1:nt)));
%   cfg=[];
%   cfg.dataset=[edir sub{ii} '/' sfiles(ll).name(1:end-3) 'eeg'];
%   if ii==108 & ll==2
%     cfg.dataset=[edir sub{ii} '/' sfiles(ll).name(1:4) '2s.eeg'];
%   end
%   cfg.trialfun='ft_trialfun_general';
%   cfg.trialdef.eventtype  = 'Stimulus';
%   cfg.trialdef.eventvalue = 'S  1'; % corresponds to lj.prepareStrobe(1)
%   cfg.trialdef.prestim = 0.7;
%   cfg.trialdef.poststim = 2.1;
%   cfgtr=ft_definetrial(cfg);
%   
%   numnuldat=size(cfgtr.trl,1);
  numnultap=length(nultr);
%   snumnul=[snumnul numnuldat];
%   nulsequse=nan(1,numnuldat);
  nulsequse(nultr(1:numnultap))=1;
%   nulsequse=nulsequse(1:numnuldat);
  rnulseq=[rnulseq nulsequse];
  
%   % total trials
%   cfg=[];
%   cfg.dataset=[edir sub{ii} '/' sfiles(ll).name(1:end-3) 'eeg'];
%   if ii==108 & ll==2
%     cfg.dataset=[edir sub{ii} '/' sfiles(ll).name(1:4) '2s.eeg'];
%   end
%   cfg.trialfun='ft_trialfun_general';
%   cfg.trialdef.eventtype  = 'Stimulus';
%   cfg.trialdef.eventvalue = 'S 20';
%   cfg.trialdef.prestim = 0.7;
%   cfg.trialdef.poststim = 2.1;
%   cfgtr=ft_definetrial(cfg);
%   
%   
%   cfg.trialdef.prestim = 1;
%   cfg.trialdef.poststim = 1;
%   cfgtr1=ft_definetrial(cfg);
%   
%   keyboard % this commented out 4 lines used to work...why not now?
% %   cfg.trialfun='ft_trialfun_general_sleep';
% %   cfg.trialdef.prestim = 1;
% %   cfg.trialdef.poststim = 1;
% %   cfgtr2=ft_definetrial(cfg);
  
  
  numtrtap=length(time_touch0);
%   numtrdat=size(cfgtr.trl,1);
%   snumtr=[snumtr numtrdat];
%   soaeff=nan(1,numtrdat);
  soaeff=nan(1,numtrtap);
  soaeff(1:numtrtap)=-(rstim{ll}.info.valve_time(1:nt) - rstim{ll}.info.audio_time(1:nt) + time_touch0 - rstim{ll}.fin.audDelay/1000); %using time_touch0, we have eliminated bad tactile trials if photodiode messed up.
%   soaeff=soaeff(1:numtrdat);
  rsoaeff=[rsoaeff soaeff];
%   soaseq=nan(1,numtrdat);
  soaseq=nan(1,numtrtap);
  soaseq(1:numtrtap)=rstim{ll}.info.soa_seq(1:numtrtap);
%   soaseq=soaseq(1:numtrdat);
  rsoaseq=[rsoaseq soaseq];
  %   soaeff=sstim{ll}.info.soa_effective2;
  %   ssoaeff=[ssoaeff soaeff];
  %   soaseq=sstim{ll}.info.soa_seq;
  %   ssoaseq=[ssoaseq soaseq];
  
  numplots=16; % should be a square number
  if plotflag
    figure(figind);
    figind=figind+1;
    for nn=1:numplots
      if 0
      ind=round((numtrdat/numplots)*nn);
      else
      ind=round((numtrtap/numplots)*nn);
      end
      %     ind=min(100,numtrdat); % no point in starting at 100 if number of trials is much fewer!
      while ind>nt | rstim{ll}.info.valves_seq(ind)~=1
        ind=ind-1;
      end
      subplot(sqrt(numplots),sqrt(numplots),nn);    plot(rstim{ll}.info.lighttime{ind},rstim{ll}.info.lightsensor{ind})
      axis([-inf inf 150 800])
    end
    %     figure;plot(sstim{ll}.info.time_touch(sstim{ll}.info.time_touch>0));axis([-inf inf .17 .25])
    figure(figind-2);subplot(2,1,2);plot(time_touch1);axis([-inf inf .17 .25])
  end
  
  auddelay=rstim{ll}.fin.audDelay/1000;
  resptime=rstim{ll}.info.prsout(1:nt)-rstim{ll}.info.start_time;
  firststimtime=min(rstim{ll}.info.audio_time(1:nt)+auddelay-rstim{ll}.info.start_time, rstim{ll}.info.valve_time(1:nt)+time_touch0-rstim{ll}.info.start_time);
  rts=resptime-firststimtime;  
  keyboard
  rrts = [rrts rts];
  
end
keyboard

%%
audDelay=sstim{1}.fin.audDelay;
soatimes=sstim{1}.fin.soa_desired;

% exclude trials for which SOA was incorrect due to time_touch errors

lim=[.01 .005 .002 .0005];
for ll=1:length(lim)
  ssoareal(ll,:)=ssoaseq;
  for nn=1:length(ssoaseq)
    if ssoaeff(nn)>-.5-lim(ll) && ssoaeff(nn)<-.5+lim(ll) && ssoaseq(nn)>0 % -500
      ssoareal(ll,nn)=1;
    elseif ssoaeff(nn)>-.07-lim(ll) && ssoaeff(nn)<-.07+lim(ll) && ssoaseq(nn)>0 % -70
      ssoareal(ll,nn)=3;
    elseif ssoaeff(nn)>-.02-lim(ll) && ssoaeff(nn)<-.02+lim(ll) && ssoaseq(nn)>0 % -20
      ssoareal(ll,nn)=4;
    elseif ssoaeff(nn)>0-lim(ll) && ssoaeff(nn)<0+lim(ll) && ssoaseq(nn)>0 % 0
      ssoareal(ll,nn)=5;
    elseif ssoaeff(nn)>.02-lim(ll) && ssoaeff(nn)<.02+lim(ll) && ssoaseq(nn)>0 % 20
      ssoareal(ll,nn)=6;
    elseif ssoaeff(nn)>.07-lim(ll) && ssoaeff(nn)<.07+lim(ll) && ssoaseq(nn)>0 % 70
      ssoareal(ll,nn)=7;
    elseif ssoaeff(nn)>.5-lim(ll) && ssoaeff(nn)<.5+lim(ll) && ssoaseq(nn)>0 % 500
      ssoareal(ll,nn)=9;
    elseif any(ssoaseq(nn)==[-2 -1 0])
    else
      ssoareal(ll,nn)=-3; % not use
    end
  end
  %   disp(['This number was total available in awake: ' num2str(length(find(~isnan(ssoaeff))))])
  %   disp(['This number was changed based on actual time delay in awake: ' num2str(length(find(ssoareal(ll,:)-ssoaseq)))])
  %   disp(['Number not used and not a NaN in awake for limit ' num2str(lim(ll)) ': ' num2str(length(find(ssoareal(ll,:)-ssoaseq~=0 & ~isnan(ssoaeff))) )])
  spctexcl(ll)=100*length(find(ssoareal(ll,:)-ssoaseq~=0 & ~isnan(ssoaeff)))/length(find(~isnan(ssoaeff)));
end
spctexcl % sitting percent excluded

lim=[.01 .005 .002 .0005];
for ll=1:length(lim)
  bsoareal(ll,:)=bsoaseq;
  for nn=1:length(bsoaseq)
    if bsoaeff(nn)>-.5-lim(ll) && bsoaeff(nn)<-.5+lim(ll) && bsoaseq(nn)>0 % -500
      bsoareal(ll,nn)=1;
    elseif bsoaeff(nn)>-.07-lim(ll) && bsoaeff(nn)<-.07+lim(ll) && bsoaseq(nn)>0 % -70
      bsoareal(ll,nn)=3;
    elseif bsoaeff(nn)>-.02-lim(ll) && bsoaeff(nn)<-.02+lim(ll) && bsoaseq(nn)>0 % -20
      bsoareal(ll,nn)=4;
    elseif bsoaeff(nn)>0-lim(ll) && bsoaeff(nn)<0+lim(ll) && bsoaseq(nn)>0 % 0
      bsoareal(ll,nn)=5;
    elseif bsoaeff(nn)>.02-lim(ll) && bsoaeff(nn)<.02+lim(ll) && bsoaseq(nn)>0 % 20
      bsoareal(ll,nn)=6;
    elseif bsoaeff(nn)>.07-lim(ll) && bsoaeff(nn)<.07+lim(ll) && bsoaseq(nn)>0 % 70
      bsoareal(ll,nn)=7;
    elseif bsoaeff(nn)>.5-lim(ll) && bsoaeff(nn)<.5+lim(ll) && bsoaseq(nn)>0 % 500
      bsoareal(ll,nn)=9;
    elseif any(bsoaseq(nn)==[-2 -1 0])
    else
      bsoareal(ll,nn)=-3; % not use
    end
  end
  %   disp(['This number was total available in sleep: ' num2str(length(find(~isnan(bsoaeff))))])
  %   disp(['This number was changed based on actual time delay in awake: ' num2str(length(find(bsoareal(ll,:)-bsoaseq)))])
  %   disp(['Number not used and not a NaN in awake for limit ' num2str(lim(ll)) ': ' num2str(length(find(bsoareal(ll,:)-bsoaseq~=0 & ~isnan(bsoaeff))) )])
  bpctexcl(ll)=100*length(find(bsoareal(ll,:)-bsoaseq~=0 & ~isnan(bsoaeff)))/length(find(~isnan(bsoaeff)));
end
bpctexcl % bed percent excluded

lim=[.01 .005 .002 .0005];
for ll=1:length(lim)
  rsoareal(ll,:)=rsoaseq;
  for nn=1:length(rsoaseq)
    if rsoaeff(nn)>-.5-lim(ll) && rsoaeff(nn)<-.5+lim(ll) && rsoaseq(nn)>0 % -500
      rsoareal(ll,nn)=1;
    elseif rsoaeff(nn)>-.07-lim(ll) && rsoaeff(nn)<-.07+lim(ll) && rsoaseq(nn)>0 % -70
      rsoareal(ll,nn)=3;
    elseif rsoaeff(nn)>-.02-lim(ll) && rsoaeff(nn)<-.02+lim(ll) && rsoaseq(nn)>0 % -20
      rsoareal(ll,nn)=4;
    elseif rsoaeff(nn)>0-lim(ll) && rsoaeff(nn)<0+lim(ll) && rsoaseq(nn)>0 % 0
      rsoareal(ll,nn)=5;
    elseif rsoaeff(nn)>.02-lim(ll) && rsoaeff(nn)<.02+lim(ll) && rsoaseq(nn)>0 % 20
      rsoareal(ll,nn)=6;
    elseif rsoaeff(nn)>.07-lim(ll) && rsoaeff(nn)<.07+lim(ll) && rsoaseq(nn)>0 % 70
      rsoareal(ll,nn)=7;
    elseif rsoaeff(nn)>.5-lim(ll) && rsoaeff(nn)<.5+lim(ll) && rsoaseq(nn)>0 % 500
      rsoareal(ll,nn)=9;
    elseif any(rsoaseq(nn)==[-2 -1 0])
    else
      rsoareal(ll,nn)=-3; % not use
    end
  end
  %   disp(['This number was total available in awake: ' num2str(length(find(~isnan(ssoaeff))))])
  %   disp(['This number was changed based on actual time delay in awake: ' num2str(length(find(ssoareal(ll,:)-ssoaseq)))])
  %   disp(['Number not used and not a NaN in awake for limit ' num2str(lim(ll)) ': ' num2str(length(find(ssoareal(ll,:)-ssoaseq~=0 & ~isnan(ssoaeff))) )])
  rpctexcl(ll)=100*length(find(rsoareal(ll,:)-rsoaseq~=0 & ~isnan(rsoaeff)))/length(find(~isnan(rsoaeff)));
end
rpctexcl % sitting percent excluded


save([sub{ii} '_audtac.mat'],'*time*','*aud*','*numtac','soatimes','s*nul*','b*nul*','r*nul*','ssoa*','bsoa*','rsoa*','*excl','rrts')

return

%%
for ii=[105:128 130:132]
  load([sub{ii} '_audtac.mat']);
  spctexcli(ii,:)=spctexcl;
  bpctexcli(ii,:)=bpctexcl;
  
  for ll=1:4
    [saa,sbb]=sort(ssoareal(ll,:));
    [baa,bbb]=sort(bsoareal(ll,:));
    %     tmp=[{ssoaeff(bb(find(aa==1)))'} {ssoaeff(bb(find(aa==3)))'} {ssoaeff(bb(find(aa==4)))'} {ssoaeff(bb(find(aa==5)))'} {ssoaeff(bb(find(aa==6)))'} {ssoaeff(bb(find(aa==7)))'} {ssoaeff(bb(find(aa==9)))'}]
    %     figure;distributionPlot(tmp);
    for ss=[1 3 4 5 6 7 9]
      ssoavals{ll}{ii-100,ss}=ssoaeff(sbb(find(saa==ss)))';
      bsoavals{ll}{ii-100,ss}=bsoaeff(bbb(find(baa==ss)))';
    end
  end
  
end
spctexcli([105:128 130:132],:)
bpctexcli([105:128 130:132],:)

% yes, this distributionPlot extrapolates and estimates a histogram function, so it will be smooth at end cutoff points.
figure;distributionPlot(ssoavals{1})
figure;distributionPlot(ssoavals{2})
figure;distributionPlot(ssoavals{3})
figure;distributionPlot(ssoavals{4})
figure;distributionPlot(bsoavals{1})
figure;distributionPlot(bsoavals{2})
figure;distributionPlot(bsoavals{3})
figure;distributionPlot(bsoavals{4})

% figure;plotspread(ssoavals{2})


