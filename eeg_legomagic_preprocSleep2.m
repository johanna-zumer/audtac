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

ii=3;

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
save time_sleep.mat time

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
cfg.channel={'all' '-ECG' '-EMG' '-EOG'};
cfg.bsfilter='yes';
cfg.bsfreq=[49 51; 99 101; 149 151];
raw_all_demean=ft_preprocessing(cfg,raw_all);

% No point artifact reject here yet, wait until after component rejection
% cfg=[];cfg.layout='EEG1010.lay';
% cfg=ft_databrowser(cfg,raw_all_demean);  %% make sure to get all segments marked here!!
% 
% cfg.artfctdef.reject='partial';
% tmp=ft_rejectartifact(cfg,raw_all_demean);
% 
% cfg=[];cfg.method='trial';
% raw_all_rej=ft_rejectvisual(cfg,tmp);

% save raw_all_rej_sleep.mat raw_all_rej

%% Skip this if using comp from 's'????
% cfg=[];
% cfg.numcomponent=30;
% cfg.method='runica';
% comp30=ft_componentanalysis(cfg,raw_all_rej);
% % comp20=ft_componentanalysis(cfg,raw_all_rej);
% save comp30_sleep comp30

%%

% % do this on Linux (crash on windows) see:
% % http://mailman.science.ru.nl/pipermail/fieldtrip/2014-May/007919.html 
% % set(gcf,'renderer','painters');
% load comp30_sleep
% cfg=[];cfg.layout='EEG1010.lay';
% cfg=ft_databrowser(cfg,comp30);
% artfctdef=cfg.artfctdef;
% save tmp_sleep.mat artfctdef
% 
% % Question:  Do I need to have 's' and 'b' concat for ICA so throwing out
% % same components?   I can always apply the rejection based on 's' on to
% % 'b', but it may not be optimal for 'b'.   hm....  But have to redo 's'
% % then?

%%
% back to windows
if 1 % load comp to throw out based on 's'.
  load comp30.mat
else
  load comp30_sleep.mat
  load raw_all_rej_sleep.mat
end

component1{2}=[9 3 25]; % ii=2, ff=3; 9=heart, 3=eyeblink, 25 artifact
component1{3}=[2 3 6 8 10 17 20 28]; % ii=3, ff=7;  =heart, 2,3,6,8,10,17,20,28=eyeblink,  artifact
component1{4}=[1 2 20 4 5 6 7 8 10 ]; % ii=3, sitting_concat;  =heart, 1,2,20=eyeblink,  4,5,6,7,8,10=artifact
component1{5}=[1 2]; % ii=05, 's' 1,2 eyeblink
component1{6}=[1 2 5 6 29]; % ii=06, 's' 1,2 eyeblink, 5,6 heartbeat, 29 artifact
component1{7}=[4 5]; %ii=07, 's' 4,5 eyemovement
component1{8}=[1 5 23]; %ii=08, 's' 1,5 eye;  23 artifact
component1{9}=[1]; %iii=09;  's' 1 eye
component1{10}=[]; %ii=10; none fit obvious eye or artifact (all blend)
component1{11}=[2 3]; %ii=11; 's' 2,3 eye
component1{12}=[]; %ii=12; none fit obvious eye or artifact (all blend)
component1{13}=[3 4]; %ii=13; 's' 3 4 eye
component1{14}=[7 8]; %ii=14; 's' 7,8 eye
component1{15}=[1 2]; %ii=15; 's' 1,2 eye
component1{16}=[4 6 23 27]; %ii=16; 's'  4 6eye   23 27 artifact
component1{17}=[1 2 16 17 27]; %ii=17; 's' 1 2  eye  16 17 27 artifact
component1{18}=[1 19]; %ii=18; 's' 1  eye  19 artifact
component1{19}=[1 ]; %ii=19; 's' 1  eye  
component1{20}=[1 2 3 30]; %ii=20; 's' 1,2,3  eye  30 artifact
component1{21}=[1 3 23 29]; %ii=21; 's' 1,,3  eye  23, 29 artifact
component1{22}=[1 6]; %ii=22; 's' 1,6  eye   artifact
component1{23}=[2  16]; %ii=23; 's' 2 eye  16 artifact
component1{24}=[2 3]; %ii=24; 's' 2,3 eye   artifact
component1{25}=[1 4 9 20]; %ii=25; 's' 1,4,9 eye  20 artifact
component1{26}=[1 4 30]; %ii=26; 's' 1,4 eye 30  artifact
component1{27}=[1 2 3 5 8 9 10 17 18 22 28]; %ii=27; 's' 1,3 eye 2,5, 8, 9,10,17 ,18,22,28 artifact
component1{28}=[2 3 5 15 26]; %ii=28; 's' 2,3,5 eye  15, 26 artifact
component1{29}=[1 2 5 6 10 17 24 28 30]; %ii=29; 's' 1 2 5 6 10 eye 17 24 28 30  artifact
component1{30}=[1 2 6 21 28 29 30]; %ii=30; 's' 1 2 6 eye 21 28 29 30  artifact
component1{31}=[1 2 12 18 20 24 27 28]; %ii=31; 's' 1 2 eye 12 18 20 artifact
component1{32}=[1 16 18 19 20 22 29]; %ii=32; 's' 1 eye  16 18 19 20 22 29 artifact

cfg=[];
cfg.component=component1{ii};
raw_all_ica=ft_rejectcomponent(cfg,comp30,raw_all_demean);
clear raw_all

cfg=[];cfg.layout='EEG1010.lay';
cfg.blocksize=30;
cfg=ft_databrowser(cfg,raw_all_ica);  %% make sure to get all segments marked here!!

cfg=[];cfg.method='trial';
tmp1=ft_rejectvisual(cfg,raw_all_ica);

cfg=[];cfg.method='channel';
tmp1=ft_rejectvisual(cfg,raw_all_ica);

cfg.artfctdef.reject='partial';
artfctdef=cfg.artfctdef;
tmp=ft_rejectartifact(cfg,raw_all_demean);

cfg=[];cfg.method='trial';
raw_all_rej=ft_rejectvisual(cfg,tmp);


% load tmp.mat 
cfg=[];
cfg.artfctdef=artfctdef;
cfg.artfctdef.reject='partial';
raw_all_ica_rej=ft_rejectartifact(cfg,raw_all_ica);

cfg=[];cfg.layout='EEG1010.lay';cfg=ft_databrowser(cfg,raw_all_ica_rej);

cfg.artfctdef.reject='partial';
raw_all_ica_rej=ft_rejectartifact(cfg,raw_all_ica_rej);

save(['raw_all_ica_rej_sleep_' sub{ii}],'raw_all_ica_rej','-v7.3');
delete raw_all_rej_sleep.mat 

%  Trial sorting

% % this runs for both awake and asleep, thus already run
% concat_stim_files(ii+100,1);
% keyboard
% close all

% See eeg_legomagic_epoching.m next

