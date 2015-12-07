clear all
if ispc
  edir='D:\audtac\eeg_data\';
  ddir='D:\audtac\legomagic\diaries\';
else
  edir='/mnt/hgfs/D/audtac/eeg_data/';
  ddir='/mnt/hgfs/D/audtac/legomagic/diaries/';
end
cd(edir)

sub{1}='p01'; % ma.a. 03/04/14
sub{2}='e01'; % ab.m. 21/05/14
sub{3}='e02'; % ma.a. 04/06/14
sub{4}='e03'; % ag.m. 10/06/14
sub{5}='e04'; % re.g. 17/06/14


ii=4;

cd([edir sub{ii} ])

%% Load data

file=dir('e*s.eeg')

spm('defaults', 'eeg');

for ll=1:length(file)
S = [];
S.dataset = ['D:\audtac\eeg_data\' sub{ii} '\' file(ll).name]
S.outfile = ['spm8_' file(ll).name(1:end-4)];
S.channels = 'all';
% S.timewindow = [0.001 520.72];
% S.blocksize = 3276800;
S.checkboundary = 1;
S.usetrials = 1;
S.datatype = 'float32-le';
S.eventpadding = 0;
S.saveorigheader = 0;
S.conditionlabel = {'Undefined'};
S.inputformat = [];
S.continuous = true;
D = spm_eeg_convert(S);
end

%% 
% Using FASST, crc_concatenate, make one file for all awake
