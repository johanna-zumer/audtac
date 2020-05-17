clear all
close all
[~,hostname]=system('hostname');
if ispc
  if ~isempty(strfind(hostname,'12807'))  % COLLES-12807
    edir='C:\Bham_PSYCHL132432_D\audtac\eeg_data\';
    esdir='C:\Bham_PSYCHL132432_D\audtac\source_data\';
    ddir='C:\Bham_PSYCHL132432_D\audtac\legomagic\diaries\';
    bdir='C:\Bham_PSYCHL132432_D\audtac\behav_data\';
    sdir='C:\Bham_PSYCHL132432_D\audtac\spss_stuff\';
    fdir='C:\Bham_PSYCHL132432_D\audtac\figs\';
    mdir='C:\Bham_PSYCHL132432_D\audtac\structural_MRI\';
    pdir='C:\Bham_PSYCHL132432_D\audtac\polhemus\';
  else
    edir='D:\audtac\eeg_data\';
    esdir='D:\audtac\source_data\';
    ddir='D:\audtac\legomagic\diaries\';
    bdir='D:\audtac\behav_data\';
    sdir='D:\audtac\spss_stuff\';
    fdir='D:\audtac\figs\';
    mdir='D:\audtac\structural_MRI\';
    pdir='D:\audtac\polhemus\';
  end
else
  if ~isempty(strfind(hostname,'151401')) % COLLES-151401
    edir='/home/adf/zumerj/audtac/eeg_data/';
    esdir='/home/adf/zumerj/audtac/source_data/';
    ddir='/home/adf/zumerj/audtac/legomagic/diaries/';
    bdir='/home/adf/zumerj/audtac/behav_data/';
    sdir='/home/adf/zumerj/audtac/spss_stuff/';
    fdir='/home/adf/zumerj/audtac/figs/';
    mdir='/home/adf/zumerj/audtac/structural_MRI/';
    pdir='/home/adf/zumerj/audtac/polhemus/';
  elseif ~isempty(strfind(hostname,'les')) | ~isempty(strfind(hostname,'LES')) % LES-LINUX_FS3
    edir='/home/zumerj/audtac/eeg_data/';
    esdir='/home/zumerj/audtac/source_data/';
    ddir='/home/zumerj/audtac/legomagic/diaries/';
    bdir='/home/zumerj/audtac/behav_data/';
    sdir='/home/zumerj/audtac/spss_stuff/';
    fdir='/home/zumerj/audtac/figs/';
    mdir='/home/zumerj/audtac/structural_MRI/';
    pdir='/home/zumerj/audtac/polhemus/';
  else % assume on VM linux of psychl-132432
    edir='/mnt/hgfs/D/audtac/eeg_data/';
    esdir='/mnt/hgfs/D/audtac/source_data/';
    ddir='/mnt/hgfs/D/audtac/legomagic/diaries/';
    bdir='/mnt/hgfs/D/audtac/behav_data/';
    sdir='/mnt/hgfs/D/audtac/spss_stuff/';
    fdir='/mnt/hgfs/D/audtac/figs/';
    mdir='/mnt/hgfs/D/audtac/structural_MRI/';
    pdir='/mnt/hgfs/D/audtac/polhemus/';
  end
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
  warning off
  if ~isempty(strfind(hostname,'12807'))  % COLLES-12807
    rmpath(genpath('C:\D_132432\Matlab\spm12b\external\fieldtrip\'))
    rmpath(genpath('C:\D_132432\Matlab\spm12b\external\mne\'))
%     rmpath(genpath('D:\fieldtrip_svn\'))
    rmpath(genpath('C:\Users\zumerj\fieldtrip_git\'))
    warning on
    addpath('C:\Users\zumerj\fieldtrip_git\')
  else
    rmpath(genpath('D:\matlab\spm8\external\fieldtrip\'))
    rmpath(genpath('D:\fieldtrip_svn\'))
    rmpath(genpath('D:\fieldtrip_git\'))
    warning on
    addpath('D:\fieldtrip_git\')
  end
else
  if ~isempty(strfind(hostname,'LES')) | ~isempty(strfind(hostname,'les'))
    warning off
    rmpath(genpath('~/matlab/spm8/external/fieldtrip/'))
    rmpath(genpath('~/fieldtrip_svn/'))
    rmpath(genpath('~/fieldtrip_git/'))
    warning on
    addpath('~/fieldtrip_git/')
  else
    warning off
    rmpath(genpath('/mnt/hgfs/D/matlab/spm8/external/fieldtrip/'))
    rmpath(genpath('/mnt/hgfs/D/fieldtrip_svn/'))
    rmpath(genpath('/mnt/hgfs/D/fieldtrip_git/'))
    warning on
    addpath('/mnt/hgfs/D/fieldtrip_git/')
  end
end
which ft_defaults.m -ALL
ft_defaults;
load([edir 'iikeep.mat'])
soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];
iiBuse=setdiff(iiBuse,3:7);
