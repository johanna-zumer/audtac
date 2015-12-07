% preprocessing of EEG data from Hills, 64ch MR-cap
clear all
if ispc
    edir='D:\audtac\eeg_data\';
else
    edir='/mnt/hgfs/D/audtac/eeg_data/';
    ddir='/mnt/hgfs/D/audtac/legomagic/diaries/';
end
cd(edir)
sub{1}='p01';
sub{2}='e01'; % a.m. 21/05/14
sub{3}='e02'; % m.a. 04/06/14



ii=2;
%% Load data

files=dir([edir sub{ii} '*.eeg']);
% read all data in first to run ICA over all conditions
for ff=1:length(files)
    cfg=[];
    cfg.dataset=files(ff).name;
    raw_all=ft_preprocessing(cfg);
    
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
    
    
    if 1
        cfg=[];
        cfg.numcomponent=30;
        cfg.method='runica';
        comp30=ft_componentanalysis(cfg,raw_all_rej);
        cfg=[];cfg.layout='EEG1010.lay';ft_databrowser(cfg,comp30)
        
        cfg=[];
        cfg.component=[9 2 12 14 25 ]; % ii=2, ff=3; 9=heart, 2=eyeblink, 12, 14 25 artifact
        cfg.component=[2 3 8 9 28]; % ii=2, ff=4; 28=heart, 2,3=eyeblink, 8,9 artifact
        raw_all_ica=ft_rejectcomponent(cfg,comp30,raw_all_rej);
        cfg=[];cfg.layout='EEG1010.lay';ft_databrowser(cfg,raw_all_ica)
        
    else
        cfg=[];
        cfg.numcomponent=30;
        cfg.method='fastica';
        cfg.randomseed=17;
        comp30=ft_componentanalysis(cfg,raw_all_demean);
        cfg=[];cfg.layout='EEG1010.lay';ft_databrowser(cfg,comp30)
        % fastica
        cfg=[];
        cfg.component=[3 4 5 13 14]; %
        raw_all_ica=ft_rejectcomponent(cfg,comp30,raw_all_demean);
        
    end
    
    cfg=[];
    cfg.demean='yes';
    cfg.reref='yes';
    cfg.channel={'all' '-ECG'};
    cfg.refchannel={'all' '-ECG'};
    %     raw_all_proc{ff}=ft_preprocessing(cfg,raw_all_ica);
    raw_all_ica_reref{ff}=ft_preprocessing(cfg,raw_all_ica);
    
end

save(['raw_all_ica_reref_' sub{ii}],'raw_all_ica_reref','-v7.3');

    clear raw_all
if 0 % they are often too big
    save(['raw_all_proc_' sub{ii}],'raw_all_proc','-v7.3');
    clear raw_all*
end

%% ICA for eyeblinks

load(['raw_all_proc_' sub{ii}],'raw_all_proc');

% cfg=[];
% cfg.randomseed=17;
% cfg.channel= {'all', '-ECG'};
% cfg.numcomponent=40;
% cfg.method='fastica';
% comp_all=ft_componentanalysis(cfg,raw_all_proc);
% save(['comp_all_' cfg.method '_' sub{ii}],'comp_all');

for ff=1:length(files)
    cfg=[];
    %     cfg.randomseed=17;
    cfg.channel= {'all', '-ECG'};
    cfg.runica.pca=30;
    cfg.numcomponent=30;
    cfg.method='runica';
    comp_all{ff}=ft_componentanalysis(cfg,raw_all_proc{ff});
    
    cfg=[];
    cfg.randomseed=17;
    cfg.channel= {'all', '-ECG'};
    cfg.numcomponent=30;
    cfg.method='fastica';
    comp_allf{ff}=ft_componentanalysis(cfg,raw_all_proc{ff});
    
end
save(['comp_all_' cfg.method '_' sub{ii}],'comp_all');

cfg=[];cfg.layout='EEG1010.lay';ft_databrowser(cfg,comp_all{ff})
cfg=[];cfg.layout='EEG1010.lay';ft_databrowser(cfg,comp_allf{ff})

% runica seed 17, 30 pca, all data
cfg=[];
cfg.component=[1 2 4 6 12]; % first 3 eyeblink; last two 50hz
raw_all_ica=ft_rejectcomponent(cfg,comp_all,raw_all_proc);
save(['raw_all_runica_' sub{ii}],'raw_all_ica');

% cfg=[];
% cfg.toilim=[0 500];
% raw_all_proc_1half=ft_redefinetrial(cfg,raw_all_proc);
%
% cfg=[];
% cfg.randomseed=17;
% cfg.channel= {'all', '-ECG'};
% cfg.numcomponent=20;
% cfg.runica.pca=20;
% cfg.method='runica';
% comp_1half=ft_componentanalysis(cfg,raw_all_proc_1half);
% save(['comp_1half_' cfg.method '_' sub{ii}],'comp_1half');
%
% cfg=[];cfg.layout='EEG1010.lay';ft_databrowser(cfg,comp_1half)


% % runica seed 17, 20 pca
% cfg=[];
% cfg.component=[1 2 3 6 9]; % first 3 eyeblink; last two 50hz
% raw_all_ica=ft_rejectcomponent(cfg,comp_1half,raw_all_proc);
% save(['raw_all_ica_' sub{ii}],'raw_all_ica');

% % fastica seed 17
% cfg=[];
% cfg.component=[6 7 20 35]; % first 3 eyeblink; 26 muscle
% raw_all_ica=ft_rejectcomponent(cfg,comp_all);

% cfg=[];
% cfg.component=[1 2 3 26]; % first 3 eyeblink; 26 muscle
% data_tac_ica=ft_rejectcomponent(cfg,comp_tac,data_tac);

%% Separate to trials
if 0
    load(['raw_all_runica_' sub{ii}],'raw_all_ica');
    
    % cfg=[];
    % % cfg.dataset=[edir sub{ii} '_passive_audtac_2.eeg'];
    % cfg.dataset=files(ff).name;
    % cfg.trialfun='ft_trialfun_general';
    % cfg.trialdef.eventtype  = 'Stimulus';
    % cfg.trialdef.eventvalue = 'S 20'; % corresponds to lj.prepareStrobe(20) (end of trial)
    % cfg.trialdef.prestim = 1;
    % cfg.trialdef.poststim = 4;
    % cfgtr=ft_definetrial(cfg);
    % cfg=[];
    % cfg.trl=cfgtr.trl;
    % % raw_tac=ft_redefinetrial(cfg,raw_all_ica);
    % raw_all=ft_redefinetrial(cfg,raw_all_ica_reref{ff});
    
    
    for trig=40+unique(info.soa_seq)
        cfg=[];
        % cfg.dataset=[edir sub{ii} '_passive_audtac_2.eeg'];
        cfg.dataset=files(ff).name;
        %     cfg.trialfun='ft_trialfun_general';
        cfg.trialfun='ft_trialfun_trialS20';
        cfg.trialdef.eventtype  = 'Stimulus';
        cfg.trialdef.eventvalue = ['S ' num2str(trig)]; % corresponds to lj.prepareStrobe
        cfg.trialdef.prestim = 'S 20';
        cfg.trialdef.poststim = {'S 20','S  8'};
        cfgtr=ft_definetrial(cfg);
        cfg=[];
        cfg.trl=cfgtr.trl;
        % raw_tac=ft_redefinetrial(cfg,raw_all_ica);
        raw_soa{trig}=ft_redefinetrial(cfg,raw_all_ica_reref{ff});
        for tt=1:length(raw_soa{trig}.trial)
            eventin=find([cfgtr.event.sample]>raw_soa{trig}.sampleinfo(tt,1) & [cfgtr.event.sample]<raw_soa{trig}.sampleinfo(tt,2));
            evind=1;
            while ~strcmp(cfgtr.event(eventin(evind)).value,'S 20')
                evind=evind+1;
            end
            evind=evind+1;
            while ~strcmp(cfgtr.event(eventin(evind)).value,'S 20')
                evind=evind+1;
            end
            
        end
        
    end
    
    cfg=[];
    % cfg.dataset=[edir sub{ii} '_passive_audtac_2.eeg'];
    cfg.dataset=files(ff).name;
    cfg.trialfun='ft_trialfun_general';
    cfg.trialdef.eventtype  = 'Stimulus';
    cfg.trialdef.eventvalue = 'S  2'; % corresponds to lj.prepareStrobe(2)
    cfg.trialdef.prestim = 0.5;
    cfg.trialdef.poststim = 2.7;
    cfgtr=ft_definetrial(cfg);
    cfgtr.trl(:,3)=cfgtr.trl(:,3)-round(time_touch*1000)';
    cfg=[];
    cfg.trl=cfgtr.trl;
    cfg.dataset=files(ff).name;
    raw_tac=ft_preprocessing(cfg);
%     raw_tac=ft_redefinetrial(cfg,raw_all_ica_reref{ff});
    
    cfg=[];
    % cfg.dataset=[edir sub{ii} '_passive_audtac_2.eeg'];
    cfg.dataset=files(ff).name;
    cfg.trialfun='ft_trialfun_general';
    cfg.trialdef.eventtype  = 'Stimulus';
    cfg.trialdef.eventvalue = 'S 10'; % corresponds to lj.prepareStrobe(10)
    cfg.trialdef.prestim = 0.7;
    cfg.trialdef.poststim = 1.3;
    cfgtr=ft_definetrial(cfg);
    cfg=[];
    cfg.trl=cfgtr.trl;
    raw_aud=ft_redefinetrial(cfg,raw_all_ica_reref{ff});
    
    cfg=[];
    % cfg.dataset=[edir sub{ii} '_passive_audtac_2.eeg'];
    cfg.dataset=files(ff).name;
    cfg.trialfun='ft_trialfun_general';
    cfg.trialdef.eventtype  = 'Stimulus';
    cfg.trialdef.eventvalue = 'S  1'; % corresponds to lj.prepareStrobe(1)
    cfg.trialdef.prestim = 0.7;
    cfg.trialdef.poststim = 1.3;
    cfgtr=ft_definetrial(cfg);
    cfg=[];
    cfg.trl=cfgtr.trl;
    raw_nul=ft_redefinetrial(cfg,raw_all_ica_reref{ff});
    
    % cfg=[];
    % cfg.dataset=[edir sub{ii} '_passive_audtac_2.eeg'];
    % % hdr=ft_read_header(cfg.dataset);
    % cfg.trialfun='ft_trialfun_general';
    % cfg.trialdef.eventtype  = 'Stimulus';
    % cfg.trialdef.eventvalue = 'S 16';
    % cfg.trialdef.prestim = 0.7;
    % cfg.trialdef.poststim = 1.3;
    % cfg=ft_definetrial(cfg);
    % data_tac=ft_preprocessing(cfg);
    
    % cfg=[];
    % cfg.dataset=[edir sub{ii} '_passive_audtac_2.eeg'];
    % cfg.trialfun='ft_trialfun_general';
    % cfg.trialdef.eventtype  = 'Stimulus';
    % cfg.trialdef.eventvalue = 'S 48';
    % cfg.trialdef.prestim = 0.7;
    % cfg.trialdef.poststim = 1.3;
    % cfg=ft_definetrial(cfg);
    % data_aud=ft_preprocessing(cfg);
    %
    % cfg=[];
    % cfg.dataset=[edir sub{ii} '_passive_audtac_2.eeg'];
    % % hdr=ft_read_header(cfg.dataset);
    % cfg.trialfun='ft_trialfun_general';
    % cfg.trialdef.eventtype  = 'Stimulus';
    % cfg.trialdef.eventvalue = 'S  1';
    % cfg.trialdef.prestim = 0.7;
    % cfg.trialdef.poststim = 1.3;
    % cfg=ft_definetrial(cfg);
    % data_nul=ft_preprocessing(cfg);
end

%% sorting by SOA

files=dir([edir sub{ii} '*.eeg']);
ff=3;
% soatimes=[-500 -100 -70 -20 0 20 70 100 500];
load([ddir files(ff).name(1:end-4) '.mat']);
soatimes=fin.soa_desired;

load raw_all_ica_reref_e01.mat

cfg=[];
% cfg.dataset=[edir sub{ii} '_passive_audtac_2.eeg'];
cfg.dataset=files(ff).name;
%     cfg.trialfun='ft_trialfun_general';
cfg.trialfun='ft_trialfun_trialS20';
cfg.trialdef.eventtype  = 'Stimulus';
cfg.trialdef.eventvalue = 'S 20'; % corresponds to lj.prepareStrobe
cfg.trialdef.prestim = 'S 20';
cfg.trialdef.poststim = {'S 20','S  8'};
cfgtr=ft_definetrial(cfg);
cfg=[];
cfg.trl=cfgtr.trl;
% raw_tac=ft_redefinetrial(cfg,raw_all_ica);
raw_all=ft_redefinetrial(cfg,raw_all_ica_reref{ff});

% tac_touch=info.time_touch(info.time_touch>0);
tactrials=find(info.time_touch>0);
for tt=1:length(tactrials)
    indx5=find(info.lightsensor{tactrials(tt)}>.5*[max(info.lightsensor{tactrials(tt)})-min(info.lightsensor{tactrials(tt)})]+min(info.lightsensor{tactrials(tt)}));
    indx8=find(info.lightsensor{tactrials(tt)}>.8*[max(info.lightsensor{tactrials(tt)})-min(info.lightsensor{tactrials(tt)})]+min(info.lightsensor{tactrials(tt)}));
    try
        mint5(tt)=indx5(1);
        mint8(tt)=indx8(1);
    end
end
% turns out that 50% is better criteria (change it in main code?) but then
% for timing to be consistent, add back 4ms.
time_touch=.001*(mint5+4-1); % add 4ms for diff of 50% vs 80%, and subtract 1 for the time relative to info.lighttime{tt}(1)
if median(time_touch-info.time_touch(tactrials))>.001
    error('time_touch not right?')
end

% if length(tac_touch)~=length(raw_tac.trial);
%     error('not same trials');
% end

% the actual sample when the tactile stim touched
raw_all.trialinfo(:,8)=nan;
raw_all.trialinfo(tactrials(raw_all.trialinfo(tactrials,4)==2),8) = raw_all.trialinfo(raw_all.trialinfo(:,4)==2,5)+round(time_touch(raw_all.trialinfo(tactrials,4)==2)*raw_all.fsample)';
raw_all.trialinfo(tactrials(raw_all.trialinfo(tactrials,6)==2),8) = raw_all.trialinfo(raw_all.trialinfo(:,6)==2,7)+round(time_touch(raw_all.trialinfo(tactrials,6)==2)*raw_all.fsample)';

% the actual sample when the auditory stim starts
raw_all.trialinfo(:,9)=nan;
raw_all.trialinfo(raw_all.trialinfo(:,4)==10,9) = raw_all.trialinfo(raw_all.trialinfo(:,4)==10,5) + round(.008*raw_all.fsample);
raw_all.trialinfo(raw_all.trialinfo(:,6)==10,9) = raw_all.trialinfo(raw_all.trialinfo(:,6)==10,7) + round(.008*raw_all.fsample);

% these should be similar.
[raw_all.trialinfo(:,9)-raw_all.trialinfo(:,8) round(1000*info.soa_effective2)']


% tactile trials
cfg=[];
cfg.trials=tactrials;
for tt=1:length(tactrials)
    cfg.begsample(tt)=max(1,raw_all.trialinfo(tactrials(tt),8)-1000);
    cfg.endsample(tt)=min(raw_all.trialinfo(tactrials(tt),8)+2000, size(raw_all.trial{tactrials(tt)},2));
end
tmp=ft_redefinetrial(cfg,raw_all);
cfg=[];
cfg.offset=-raw_all.trialinfo(tactrials,8);
raw_tac=ft_redefinetrial(cfg,tmp);
clear tmp;

% auditory trials
audtrials=find(info.auditory_seq(1,:));
cfg=[];
cfg.trials=audtrials;
for tt=1:length(audtrials)
    cfg.begsample(tt)=max(1,raw_all.trialinfo(audtrials(tt),9)-1000);
    cfg.endsample(tt)=min(raw_all.trialinfo(audtrials(tt),9)+2000, size(raw_all.trial{audtrials(tt)},2));
end
tmp=ft_redefinetrial(cfg,raw_all);
cfg=[];
cfg.offset=-raw_all.trialinfo(audtrials,9);
raw_aud=ft_redefinetrial(cfg,tmp);

% null trials
nultrials=find(raw_all.trialinfo(:,4)==1);
cfg=[];
cfg.trials=nultrials;
for tt=1:length(nultrials)
    cfg.begsample(tt)=max(1,raw_all.trialinfo(nultrials(tt),5)-1000);
    cfg.endsample(tt)=min(raw_all.trialinfo(nultrials(tt),5)+2000, size(raw_all.trial{nultrials(tt)},2));
end
tmp=ft_redefinetrial(cfg,raw_all);
cfg=[];
cfg.offset=-raw_all.trialinfo(nultrials,5);
raw_nul=ft_redefinetrial(cfg,tmp);


% % tactile trials
% raw_tac.trialinfo(:,2:5)=nan;
% for tt=1:length(raw_tac.trial)
%      eventin=find([cfgtr.event.sample]>raw_tac.sampleinfo(tt,1) & [cfgtr.event.sample]<raw_tac.sampleinfo(tt,2));
%      aud(1:2)=0;
%      tac(1:2)=0;
%      for ee=1:length(eventin)
%          switch cfgtr.event(eventin(ee)).value
%              case 'S  2'
%                  if ee==1
%                      tac(1)=1;
%                  elseif ee==2
%                      tac(2)=1;
%                  end
%              case 'S 10'
%                  if ee==1
%                      aud(1)=1;
%                  elseif ee==2
%                      aud(2)=1;
%                  end
%          end
%      end
% %      if length(eventin)==2 % old code
%      if length(eventin)==4 % two events for both aud and tac, plus two for trial type and end of trial
%          soa=cfgtr.event(eventin(2)).sample-cfgtr.event(eventin(1)).sample;
%          soams=soa/raw_tac.fsample*1000; % convert to ms
%          switch tac(1)
%              case 1 % if tactile is first
%                  raw_tac.trialinfo(tt,2)=soams-10; % -2 is code for tactile only.
%              case 0 % if tactile is second
%                  raw_tac.trialinfo(tt,2)=-soams-10; % -2 is code for tactile only.
%          end
%          switch aud(1)
%              case 1 % if aud is first
%                  raw_tac.trialinfo(tt,4)=soams+10; % -2 is code for tactile only.
%              case 0 % if aud is second
%                  raw_tac.trialinfo(tt,4)=-soams+10; % -2 is code for tactile only.
%          end
%      end
%      try
%          raw_tac.trialinfo(tt,3)=nearest(soatimes, raw_tac.trialinfo(tt,2));
%      catch
%          raw_tac.trialinfo(tt,3)=nan;
%      end
%      try
%          raw_tac.trialinfo(tt,5)=nearest(soatimes, raw_tac.trialinfo(tt,4));
%      catch
%          raw_tac.trialinfo(tt,5)=nan;
%      end
% %      raw_tac.soatype{tt}=num2str(raw_tac.trialinfo(tt,3));
% end
%
% % auditory trials
% raw_aud.trialinfo(:,2:5)=nan;
% for tt=1:length(raw_aud.trial)
%      eventin=find([cfgtr.event.sample]>raw_aud.sampleinfo(tt,1) & [cfgtr.event.sample]<raw_aud.sampleinfo(tt,2));
%      aud(1:2)=0;
%      tac(1:2)=0;
%      for ee=1:length(eventin)
%          switch cfgtr.event(eventin(ee)).value
%              case 'S 16'
%                  if ee==1
%                      tac(1)=1;
%                  elseif ee==2
%                      tac(2)=1;
%                  end
%              case 'S 48'
%                  if ee==1
%                      aud(1)=1;
%                  elseif ee==2
%                      aud(2)=1;
%                  end
%          end
%      end
%      if length(eventin)==2
%          soa=cfgtr.event(eventin(2)).sample-cfgtr.event(eventin(1)).sample;
%          soams=soa/raw_aud.fsample*1000; % convert to ms
%          switch tac(1)
%              case 1 % if tactile is first
% %                  soatype{tt}=num2str(soams-10); % and 10ms delay on trigger
%                  raw_aud.trialinfo(tt,2)=soams-10; % -2 is code for tactile only.
%              case 0 % if tactile is second
% %                  soatype{tt}=num2str(-soams-10); % convert to ms
%                  raw_aud.trialinfo(tt,2)=-soams-10; % -2 is code for tactile only.
%          end
%          switch aud(1)
%              case 1 % if aud is first
%                  raw_aud.trialinfo(tt,4)=soams+10; % -2 is code for tactile only.
%              case 0 % if aud is second
%                  raw_aud.trialinfo(tt,4)=-soams+10; % -2 is code for tactile only.
%          end
%      end
%      try
%          raw_aud.trialinfo(tt,3)=nearest(soatimes, raw_aud.trialinfo(tt,2));
%      catch
%          raw_aud.trialinfo(tt,3)=nan;
%      end
%      try
%          raw_aud.trialinfo(tt,5)=nearest(soatimes, raw_aud.trialinfo(tt,4));
%      catch
%          raw_aud.trialinfo(tt,5)=nan;
%      end
% %      raw_aud.soatype{tt}=num2str(raw_tac.trialinfo(tt,3));
% end

% thus, both raw_tac and raw_aud have .trialinfo(:,2) indicating the soa
% defined as the amount of time that tactile leads auditory (if positive
% number) and if negative means that auditory is first.
% the (:,3) column gives category as defined in call_audtac_SOA_variations

% the (:,4) column is reverse of that (i.e. the time that auditory leads
% tactile if positive).

%% ICA for eyeblinks

% cfg=[];
% cfg.randomseed=13;
% cfg.channel= {'all', '-ECG'};
% cfg.numcomponent=50;
% comp_tac=ft_componentanalysis(cfg,data_tac);
% cfg=[];cfg.layout='EEG1010.lay';ft_databrowser(cfg,comp_tac)
% cfg=[];
% cfg.component=[1 2 3 26]; % first 3 eyeblink; 26 muscle
% data_tac_ica=ft_rejectcomponent(cfg,comp_tac,data_tac);
%
% cfg=[];
% cfg.randomseed=13;
% cfg.channel= {'all', '-ECG'};
% cfg.numcomponent=50;
% comp_aud=ft_componentanalysis(cfg,data_aud);
% cfg=[];cfg.layout='EEG1010.lay';ft_databrowser(cfg,comp_tac)
% cfg=[];
% cfg.component=[1 2 3 26]; % first 3 eyeblink; 26 muscle
% data_tac_ica=ft_rejectcomponent(cfg,comp_tac,data_tac);

%% rejecting artifacts

cfg=[];
cfg.method='summary';
% raw_all_rej=ft_rejectvisual(cfg,raw_all);
raw_tac_rej=ft_rejectvisual(cfg,raw_tac);
raw_aud_rej=ft_rejectvisual(cfg,raw_aud);
raw_nul_rej=ft_rejectvisual(cfg,raw_nul);

clear raw_all_rej
save(['raw_each_rej_' sub{ii}],'raw_*rej');
clear raw_aud raw_nul raw_tac

