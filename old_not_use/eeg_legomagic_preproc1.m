% preprocessing of EEG data from Hills, 64ch MR-cap
clear all
if ispc
    edir='D:\audtac\eeg_data\';
else
    edir='/mnt/hgfs/D/audtac/eeg_data/';
    ddir='/mnt/hgfs/D/audtac/legomagic/diaries/';
end
cd(edir)

sub{1}='p01'; % ma.a. 03/04/14
sub{2}='e01'; % ab.m. 21/05/14
sub{3}='e02'; % m.a. 04/06/14
sub{4}='e03'; % ag.m. 10/06/14


ii=3;

cd([edir sub{ii} ])

%% Load data

% files=dir([edir sub{ii} '*.eeg']);

file=dir(['cc_cc_spm8_' sub{ii} '*.mat']);
file=dir([sub{ii} '_12b.eeg']);

cfg=[];
cfg.dataset=file.name;
raw_all=ft_preprocessing(cfg);

% size of sleep-score chunk in samples
chunksamp=raw_all.hdr.orig.other.CRC.score{3}*raw_all.fsample;
% create sample-by-sample vector of sleep score
for ll=1:length(raw_all.hdr.orig.other.CRC.score{1})
    samplesleep((ll-1)*chunksamp+1:(ll*chunksamp),1)=raw_all.hdr.orig.other.CRC.score{1}(ll);
end
% add as new field to events structure the sleep score
for ll=1:length(raw_all.hdr.orig.trials.events)
    raw_all.hdr.orig.trials.events(ll).sleep= samplesleep(dsearchn(raw_all.time{1}',raw_all.hdr.orig.trials.events(ll).time));
end

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
clear tmp

%     cfg=[];
%     cfg.numcomponent=30;
%     cfg.method='fastica';
%     cfg.randomseed=17;
cfg=[];
cfg.numcomponent=30;
cfg.method='runica';
comp30=ft_componentanalysis(cfg,raw_all_rej);
save(['comp30' file.name(1:end-4)],'comp30')
cfg=[];cfg.layout='EEG1010.lay';ft_databrowser(cfg,comp30)

cfg=[];
cfg.component=[9 3 25]; % ii=2, ff=3; 9=heart, 3=eyeblink, 25 artifact
cfg.component=[2 3 6 8 10 17 20 28]; % ii=3, ff=7;  =heart, 2,3,6,8,10,17,20,28=eyeblink,  artifact
cfg.component=[1 2 20 4 5 6 7 8 10 ]; % ii=3, sitting_concat;  =heart, 1,2,20=eyeblink,  4,5,6,7,8,10=artifact
cfg.component=[2 6 9]; % ii=3, e02_12b;  =heart, =eyeblink,  2,6,9=artifact
raw_all_ica=ft_rejectcomponent(cfg,comp30,raw_all_rej);
cfg=[];cfg.layout='EEG1010.lay';cfg=ft_databrowser(cfg,raw_all_ica);

cfg.artfctdef.reject='partial';
raw_all_ica_rej=ft_rejectartifact(cfg,raw_all_ica);


save(['raw_all_ica_rej_' sub{ii} '_' file.name(1:end-4)],'raw_all_ica_rej','-v7.3');

%%  Trial sorting

for ff=1:length(files)
    load([ddir files(ff).name(1:end-4) '.mat']);
    load([ddir file.name(1:end-4) '.mat']);
    clear dia2
    try
        dia2=load([ddir files(ff).name(1:end-4) 'b.mat']);
        info.time_touch=[info.time_touch dia2.info.time_touch];
        info.auditory_seq=[info.auditory_seq dia2.info.auditory_seq];
    end
    try
        dia2=load([ddir file.name(1:end-4) 'b.mat']);
        info.time_touch=[info.time_touch dia2.info.time_touch];
        info.auditory_seq=[info.auditory_seq dia2.info.auditory_seq];
    end
    audtrials=find(info.auditory_seq(1,1:length(info.lightsensor)));
    tactrials=find(info.time_touch>0);
    soatimes=fin.soa_desired;
    
    
    if ii<3
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
    else
        time_touch=info.time_touch(tactrials);
    end
    
    
    
    cfg=[];
    try
        cfg.dataset=files(ff).name;
    catch
        cfg.dataset=file.name;
    end
    cfg.trialfun='ft_trialfun_general';
    cfg.trialdef.eventtype  = 'Stimulus';
    cfg.trialdef.eventvalue = 'S  2'; % corresponds to lj.prepareStrobe(2)
    cfg.trialdef.prestim = 0.5;
    cfg.trialdef.poststim = 2.3;
    cfgtr=ft_definetrial(cfg);
    numtac=min(size(cfgtr.trl,1),length(time_touch));
    cfgtr.trl=cfgtr.trl(1:numtac,:);
    cfgtr.trl(:,3)=cfgtr.trl(:,3)-round(time_touch(1:numtac)*raw_all_ica_rej.fsample)';
    cfg=[];
    cfg.trl=cfgtr.trl;
    raw_tac=ft_redefinetrial(cfg,raw_all_ica_rej);
    
    cfg=[];
    try
        cfg.dataset=files(ff).name;
    catch
        cfg.dataset=file.name;
    end
    cfg.trialfun='ft_trialfun_general';
    cfg.trialdef.eventtype  = 'Stimulus';
    cfg.trialdef.eventvalue = 'S 10'; % corresponds to lj.prepareStrobe(10)
    cfg.trialdef.prestim = 0.7;
    cfg.trialdef.poststim = 2.1;
    cfgtr=ft_definetrial(cfg);
    numaud=min(size(cfgtr.trl,1),length(audtrials));
    cfgtr.trl=cfgtr.trl(1:numaud,:);
    cfgtr.trl(:,3)=cfgtr.trl(:,3)-round(fin.audDelay/1000*raw_all_ica_rej.fsample);
    cfg=[];
    cfg.trl=cfgtr.trl;
    raw_aud=ft_redefinetrial(cfg,raw_all_ica_rej);
    
end


cfg=[];
cfg.method='summary';
% raw_all_rej=ft_rejectvisual(cfg,raw_all);
raw_tac_rej=ft_rejectvisual(cfg,raw_tac);
raw_aud_rej=ft_rejectvisual(cfg,raw_aud);
raw_nul_rej=ft_rejectvisual(cfg,raw_nul);

clear raw_all_rej
save(['raw_each_rej_' sub{ii} '_' file.name(1:end-4)],'raw_*rej');
clear raw_aud raw_nul raw_tac

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

cfg=[];
cfg.reref='yes';
cfg.refchannel='all';
data_tac_filt_ref=ft_preprocessing(cfg,data_tac_filt);
data_aud_filt_ref=ft_preprocessing(cfg,data_aud_filt);
data_nul_filt_ref=ft_preprocessing(cfg,data_nul_filt);

clear raw_*


% collapsing over all SOA conditions
for ll=1:size(data_tac_filt_ref.trial,2),
    keeptrial(ll)=~any(any(isnan(data_tac_filt_ref.trial{ll})));
end
cfg=[];
cfg.vartrllength=2;
cfg.trials=find(keeptrial);
% cfg.keeptrials='yes';
tlock_tac{11}=ft_timelockanalysis(cfg,data_tac_filt_ref);

% collapsing over all SOA conditions
for ll=1:size(data_aud_filt_ref.trial,2),
    keeptrial(ll)=~any(any(isnan(data_aud_filt_ref.trial{ll})));
end
cfg=[];
cfg.vartrllength=2;
cfg.trials=find(keeptrial);
% cfg.keeptrials='yes';
tlock_aud{11}=ft_timelockanalysis(cfg,data_aud_filt_ref);


cfg=[];
cfg.interactive='yes';
cfg.layout='EEG1010.lay';
ft_multiplotER(cfg,tlock_tac{11});
ft_multiplotER(cfg,tlock_aud{11});

