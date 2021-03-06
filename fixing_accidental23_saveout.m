clear all
close all
if ispc
  edir='D:\audtac\eeg_data\';
  ddir='D:\audtac\legomagic\diaries\';
  bdir='D:\audtac\behav_data\';
  sdir='D:\audtac\spss_stuff\';
  fdir='D:\audtac\figs\';
  mdir='D:\audtac\structural_MRI\';
  pdir='D:\audtac\polhemus\';
else
  [~,hostname]=system('hostname');
  if ~isempty(strfind(hostname,'les')) | ~isempty(strfind(hostname,'LES')) % either COLLES-151401 or LES-LINUX_FS3
    edir='/home/zumerj/audtac/eeg_data/';
    ddir='/home/zumerj/audtac/legomagic/diaries/';
    bdir='/home/zumerj/audtac/behav_data/';
    sdir='/home/zumerj/audtac/spss_stuff/';
    fdir='/home/zumerj/audtac/figs/';
    mdir='/home/zumerj/audtac/structural_MRI/';
    pdir='/home/zumerj/audtac/polhemus/';
  else % assume on VM linux of psychl-132432
    edir='/mnt/hgfs/D/audtac/eeg_data/';
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
  rmpath(genpath('D:\matlab\spm8\external\fieldtrip\'))
  rmpath(genpath('D:\fieldtrip_svn\'))
  addpath('D:\fieldtrip_svn\')
else
  if strfind(hostname,'les')
  else
    rmpath(genpath('/mnt/hgfs/D/matlab/spm8/external/fieldtrip/'))
    rmpath(genpath('/mnt/hgfs/D/fieldtrip_svn/'))
    addpath('/mnt/hgfs/D/fieldtrip_svn/')
  end
end
which ft_defaults.m
ft_defaults;
load([edir 'iikeep.mat'])

%%
sleep=0;
tt=3;
tacaud=1;
soalist=[1 3 4 5 6 7 9];

for ii=setdiff(iiSuse,3:11)
  cd([edir sub{ii} ])
  clearvars -except ii sub *dir ii*use sleep tt tacaud soalist
  load(['freq_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_tt' num2str(tt) '.mat'])
  whovars=whos;
  for ww=1:length(whovars)
    if length(whovars(ww).size)>2
      wwvv.(whovars(ww).name)=eval(whovars(ww).name);
    end
  end
  wwfields=fieldnames(wwvv);  
  for ww=1:length(wwfields)
    for ll=soalist
      for ss=10:size(wwvv.(wwfields{ww}),3)
        for ff=1:size(wwvv.(wwfields{ww}),4)
          try
            wwvv.(wwfields{ww}){ll,2,ss,ff}=[];
          catch
            wwvv.(wwfields{ww})(ll,2,ss,ff)=0;
          end
        end % ff
      end % ss
    end % ll
  end
  
  save(['freq_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_tt' num2str(tt) '.mat'],'-struct','wwvv','-v7.3')
  
end