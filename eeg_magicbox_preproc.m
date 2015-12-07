% preprocessing of EEG data from Hills, 64ch MR-cap

if ispc
    edir='D:\audtac\eeg_data\';
else
    edir='/mnt/hgfs/D/audtac/eeg_data/';
end
cd(edir);

sub{1}='p01';
sub{2}='e01';



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
    raw_all_demean=ft_preprocessing(cfg,raw_all);

    cfg=[];
    cfg.numcomponent=30;
    cfg.method='runica';
    cfg.randomseed=17; % except that this doesn't work.
    comp30=ft_componentanalysis(cfg,raw_all_demean);
    cfg=[];cfg.layout='EEG1010.lay';ft_databrowser(cfg,comp30)

    % runica seed 17, 30 pca, all data
    cfg=[];
    cfg.component=[3 4 5 13 14 15 17 18]; % 
    raw_all_ica=ft_rejectcomponent(cfg,comp30,raw_all_demean);
    
    cfg=[];
    cfg.demean='yes';
    cfg.reref='yes';
    cfg.channel={'all' '-ECG'};
    cfg.refchannel={'all' '-ECG'};
%     raw_all_proc{ff}=ft_preprocessing(cfg,raw_all_ica);
    raw_all_ica_reref{ff}=ft_preprocessing(cfg,raw_all_ica);
    
end
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

load(['raw_all_runica_' sub{ii}],'raw_all_ica');

cfg=[];
% cfg.dataset=[edir sub{ii} '_passive_audtac_2.eeg'];
cfg.dataset=files(ff).name;
cfg.trialfun='ft_trialfun_general';
cfg.trialdef.eventtype  = 'Stimulus';
cfg.trialdef.eventvalue = 'S  2'; % corresponds to lj.prepareStrobe(2)
cfg.trialdef.prestim = 0.7;
cfg.trialdef.poststim = 1.3;
cfgtr=ft_definetrial(cfg);
cfg=[];
cfg.trl=cfgtr.trl;
% raw_tac=ft_redefinetrial(cfg,raw_all_ica);
raw_tac=ft_redefinetrial(cfg,raw_all_ica_reref{ff});

cfg=[];
cfg.dataset=[edir sub{ii} '_passive_audtac_2.eeg'];
cfg.trialfun='ft_trialfun_general';
cfg.trialdef.eventtype  = 'Stimulus';
cfg.trialdef.eventvalue = 'S 10'; % corresponds to lj.prepareStrobe(10)
cfg.trialdef.prestim = 0.7;
cfg.trialdef.poststim = 1.3;
cfgtr=ft_definetrial(cfg);
cfg=[];
cfg.trl=cfgtr.trl;
raw_aud=ft_redefinetrial(cfg,raw_all_ica);

cfg=[];
cfg.dataset=[edir sub{ii} '_passive_audtac_2.eeg'];
cfg.trialfun='ft_trialfun_general';
cfg.trialdef.eventtype  = 'Stimulus';
cfg.trialdef.eventvalue = 'S  1'; % corresponds to lj.prepareStrobe(1)
cfg.trialdef.prestim = 0.7;
cfg.trialdef.poststim = 1.3;
cfgtr=ft_definetrial(cfg);
cfg=[];
cfg.trl=cfgtr.trl;
raw_nul=ft_redefinetrial(cfg,raw_all_ica);

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

%% sorting by SOA

soatimes=[-500 -100 -70 -20 0 20 70 100 500];


% tactile trials
raw_tac.trialinfo(:,2:5)=nan;
for tt=1:length(raw_tac.trial)
     eventin=find([cfgtr.event.sample]>raw_tac.sampleinfo(tt,1) & [cfgtr.event.sample]<raw_tac.sampleinfo(tt,2));
     aud(1:2)=0;
     tac(1:2)=0;
     for ee=1:length(eventin)
         switch cfgtr.event(eventin(ee)).value
             case 'S 16' 
                 if ee==1
                     tac(1)=1;
                 elseif ee==2
                     tac(2)=1;
                 end
             case 'S 48' 
                 if ee==1
                     aud(1)=1;
                 elseif ee==2
                     aud(2)=1;
                 end
         end
     end
     if length(eventin)==2
         soa=cfgtr.event(eventin(2)).sample-cfgtr.event(eventin(1)).sample;
         soams=soa/raw_tac.fsample*1000; % convert to ms 
         switch tac(1)
             case 1 % if tactile is first
                 raw_tac.trialinfo(tt,2)=soams-10; % -2 is code for tactile only.
             case 0 % if tactile is second
                 raw_tac.trialinfo(tt,2)=-soams-10; % -2 is code for tactile only.
         end
         switch aud(1)
             case 1 % if aud is first
                 raw_tac.trialinfo(tt,4)=soams+10; % -2 is code for tactile only.
             case 0 % if aud is second
                 raw_tac.trialinfo(tt,4)=-soams+10; % -2 is code for tactile only.
         end
     end
     try
         raw_tac.trialinfo(tt,3)=nearest(soatimes, raw_tac.trialinfo(tt,2));
     catch
         raw_tac.trialinfo(tt,3)=nan;
     end
     try
         raw_tac.trialinfo(tt,5)=nearest(soatimes, raw_tac.trialinfo(tt,4));
     catch
         raw_tac.trialinfo(tt,5)=nan;
     end
%      raw_tac.soatype{tt}=num2str(raw_tac.trialinfo(tt,3));
end

% auditory trials
raw_aud.trialinfo(:,2:5)=nan;
for tt=1:length(raw_aud.trial)
     eventin=find([cfgtr.event.sample]>raw_aud.sampleinfo(tt,1) & [cfgtr.event.sample]<raw_aud.sampleinfo(tt,2));
     aud(1:2)=0;
     tac(1:2)=0;
     for ee=1:length(eventin)
         switch cfgtr.event(eventin(ee)).value
             case 'S 16'
                 if ee==1
                     tac(1)=1;
                 elseif ee==2
                     tac(2)=1;
                 end
             case 'S 48'
                 if ee==1
                     aud(1)=1;
                 elseif ee==2
                     aud(2)=1;
                 end
         end
     end
     if length(eventin)==2
         soa=cfgtr.event(eventin(2)).sample-cfgtr.event(eventin(1)).sample;
         soams=soa/raw_aud.fsample*1000; % convert to ms 
         switch tac(1)
             case 1 % if tactile is first
%                  soatype{tt}=num2str(soams-10); % and 10ms delay on trigger
                 raw_aud.trialinfo(tt,2)=soams-10; % -2 is code for tactile only.
             case 0 % if tactile is second
%                  soatype{tt}=num2str(-soams-10); % convert to ms
                 raw_aud.trialinfo(tt,2)=-soams-10; % -2 is code for tactile only.
         end
         switch aud(1)
             case 1 % if aud is first
                 raw_aud.trialinfo(tt,4)=soams+10; % -2 is code for tactile only.
             case 0 % if aud is second
                 raw_aud.trialinfo(tt,4)=-soams+10; % -2 is code for tactile only.
         end
     end
     try
         raw_aud.trialinfo(tt,3)=nearest(soatimes, raw_aud.trialinfo(tt,2));
     catch
         raw_aud.trialinfo(tt,3)=nan;
     end
     try
         raw_aud.trialinfo(tt,5)=nearest(soatimes, raw_aud.trialinfo(tt,4));
     catch
         raw_aud.trialinfo(tt,5)=nan;
     end
%      raw_aud.soatype{tt}=num2str(raw_tac.trialinfo(tt,3));
end

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
raw_tac_rej=ft_rejectvisual(cfg,raw_tac);
raw_aud_rej=ft_rejectvisual(cfg,raw_aud);
raw_nul_rej=ft_rejectvisual(cfg,raw_nul);

clear raw_aud raw_nul raw_tac

save(['raw_each_rej_' sub{ii}],'raw_*');
