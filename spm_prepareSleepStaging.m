% spm preproc: convert, concat, montage


clear all
if ispc
    sdir='D:\sleep_EEG_example\';
    edir='D:\audtac\eeg_data\';
%     ddir='D:\audtac\legomagic\diaries\';
else
    sdir='/mnt/hgfs/D/sleep_EEG_example/';
    edir='/mnt/hgfs/D/audtac/eeg_data/';
%     ddir='/mnt/hgfs/D/audtac/legomagic/diaries/';
end
cd(sdir)

sub{1}='sub01_sleep1'; %
sub{2}='sub01_sleep2'; %

ii=1; % change for each

%%

% load into SPM format
sfiles=dir([sub{ii} '.eeg']);
for ll=1:length(sfiles) % should only be one
    S.dataset=sfiles(ll).name;
    S.continuous=1;
    D=spm_eeg_convert(S);
end

%  files=dir(['spm8_' sub{ii} '.mat']);

event=ft_read_event([sub{ii} '.vhdr']);

% Do chunking tool in Fasst using event timings

% if 0
% 
% S.D=D;
% S.timewin=[D.time(event(2).sample) D.time(end)];
% S.channels='all';
% D = spm_eeg_crop(S);
% end

% Create *.mat sleep montage using: eeg_sleep_MEEGmontage_AnalyzerJo.m
%
% See also eeg_sleep_montage_AnalyzerJo.m and eeg_sleep_montageTil_ExG.m and eeg_sleep_montageTil_EOGawake.m 

load([edir 'sleepMEEGAnalyzerMontage.mat']);
S=[];
S.D=['chk1_spm8_' sub{ii} '.mat'];
S.montage=montage;
S.keepothers=0;
[D,montage]=spm_eeg_montage(S);

% Do sleep staging in FASST in Windows Matlab
fasst % opens software

S=[];
S.D=['M' D.fname];
S.timewin=[D.time(event(2).sample) D.time(end)];
S.channels='all';
D = spm_eeg_crop(S);







