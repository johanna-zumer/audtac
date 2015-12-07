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

ii=32;

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
files=dir('cc*s.mat');
fileb=dir('cc*b.mat');

% if any(ii==[3 ])
%   file=dir(['cc_cc_spm8_' sub{ii} '*.mat']);
% elseif any(ii==[4 ])
%   file=dir(['cc_spm8_' sub{ii} '*.mat']);
% end

cfg=[];
cfg.dataset=files.name;
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
cfg.channel={'all' '-ECG' '-EMG' '-EOG'};
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

%%
% do this on Linux (crash on windows) see:
% http://mailman.science.ru.nl/pipermail/fieldtrip/2014-May/007919.html 
% set(gcf,'renderer','painters');
load comp30
cfg=[];cfg.layout='EEG1010.lay';
cfg=ft_databrowser(cfg,comp30);
artfctdef=cfg.artfctdef;
save tmp.mat artfctdef

% Question:  Do I need to have 's' and 'b' concat for ICA so throwing out
% same components?   I can always apply the rejection based on 's' on to
% 'b', but it may not be optimal for 'b'.   hm....  But have to redo 's'
% then?

%%
% back to windows
load comp30.mat
load raw_all_rej.mat

cfg=[];
cfg.component=[9 3 25]; % ii=2, ff=3; 9=heart, 3=eyeblink, 25 artifact
cfg.component=[2 3 6 8 10 17 20 28]; % ii=3, ff=7;  =heart, 2,3,6,8,10,17,20,28=eyeblink,  artifact
cfg.component=[1 2 20 4 5 6 7 8 10 ]; % ii=3, sitting_concat;  =heart, 1,2,20=eyeblink,  4,5,6,7,8,10=artifact
cfg.component=[1 2]; % ii=05, 's' 1,2 eyeblink
cfg.component=[1 2 5 6 29]; % ii=06, 's' 1,2 eyeblink, 5,6 heartbeat, 29 artifact
cfg.component=[4 5]; %ii=07, 's' 4,5 eyemovement
cfg.component=[1 5 23]; %ii=08, 's' 1,5 eye;  23 artifact
cfg.component=[1]; %iii=09;  's' 1 eye
cfg.component=[]; %ii=10; none fit obvious eye or artifact (all blend)
cfg.component=[2 3]; %ii=11; 's' 2,3 eye
cfg.component=[]; %ii=12; none fit obvious eye or artifact (all blend)
cfg.component=[3 4]; %ii=13; 's' 3 4 eye
cfg.component=[7 8]; %ii=14; 's' 7,8 eye
cfg.component=[1 2]; %ii=15; 's' 1,2 eye
cfg.component=[4 6 23 27]; %ii=16; 's'  4 6eye   23 27 artifact
cfg.component=[1 2 16 17 27]; %ii=17; 's' 1 2  eye  16 17 27 artifact
cfg.component=[1 19]; %ii=18; 's' 1  eye  19 artifact
cfg.component=[1 ]; %ii=19; 's' 1  eye  
cfg.component=[1 2 3 30]; %ii=20; 's' 1,2,3  eye  30 artifact
cfg.component=[1 3 23 29]; %ii=21; 's' 1,,3  eye  23, 29 artifact
cfg.component=[1 6]; %ii=22; 's' 1,6  eye   artifact
cfg.component=[2  16]; %ii=23; 's' 2 eye  16 artifact
cfg.component=[2 3]; %ii=24; 's' 2,3 eye   artifact
cfg.component=[1 4 9 20]; %ii=25; 's' 1,4,9 eye  20 artifact
cfg.component=[1 4 30]; %ii=26; 's' 1,4 eye 30  artifact
cfg.component=[1 2 3 5 8 9 10 17 18 22 28]; %ii=27; 's' 1,3 eye 2,5, 8, 9,10,17 ,18,22,28 artifact
cfg.component=[2 3 5 15 26]; %ii=28; 's' 2,3,5 eye  15, 26 artifact
cfg.component=[1 2 5 6 10 17 24 28 30]; %ii=29; 's' 1 2 5 6 10 eye 17 24 28 30  artifact
cfg.component=[1 2 6 21 28 29 30]; %ii=30; 's' 1 2 6 eye 21 28 29 30  artifact
cfg.component=[1 2 12 18 20 24 27 28]; %ii=31; 's' 1 2 eye 12 18 20 artifact
cfg.component=[1 16 18 19 20 22 29]; %ii=32; 's' 1 eye  16 18 19 20 22 29 artifact
raw_all_ica=ft_rejectcomponent(cfg,comp30,raw_all_rej);
load tmp.mat 
cfg=[];
cfg.artfctdef=artfctdef;
cfg.artfctdef.reject='partial';
raw_all_ica_rej=ft_rejectartifact(cfg,raw_all_ica);

cfg=[];cfg.layout='EEG1010.lay';cfg=ft_databrowser(cfg,raw_all_ica_rej);

cfg.artfctdef.reject='partial';
raw_all_ica_rej=ft_rejectartifact(cfg,raw_all_ica_rej);

save(['raw_all_ica_rej_' sub{ii}],'raw_all_ica_rej','-v7.3');
delete raw_all_rej.mat 

%  Trial sorting

concat_stim_files(ii+100,1);
keyboard
close all

% See eeg_legomagic_epoching.m next

