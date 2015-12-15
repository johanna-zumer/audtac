% testing/comparing manual artifact rejection of awake data versus
% automated for muscle etc.


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

%%

for ii=30
  for sleep=[0 1]
    
    clearvars -except ii sub  edir ddir sleep
    cd([edir sub{ii} ])
    if sleep
      files=dir('cc*b.mat');
    else
      files=dir('cc*s.mat');
    end
    
    cfg=[];
    cfg.dataset=files.name;
    raw_all=ft_preprocessing(cfg);
    cfg=[];
    cfg.demean='yes';
    cfg.channel={'all' '-ECG' '-EMG' '-EOG'};
    cfg.bsfilter='yes';
    cfg.bsfreq=[49 51; 99 101; 149 151];
    raw_all_demean=ft_preprocessing(cfg,raw_all);
    clear raw_all
    
    % apply filter to unepoched data
    cfg=[];
    cfg.hpfilter='yes';
    cfg.hpfiltord=3;
    cfg.hpfreq=0.2;  % Could do higher for awake data but not want higher for sleep.
    %   raw_all_ica_rej=ft_preprocessing(cfg,raw_all_ica_rej);
    tmp2=ft_preprocessing(cfg,raw_all_demean);
    clear raw_all_demean
    
    fltpadding=1;
    trlpadding=1.5;
    % create trl with 30s segments
    for ll=1:ceil( tmp2.sampleinfo(2)/30000)
      if ll==ceil( tmp2.sampleinfo(2)/30000)
        trl(ll,:)=[(30000*(ll-1))+1 tmp2.sampleinfo(2) 0];
      elseif ll==1
        trl(ll,:)=[200 (30000*(ll-1))+30000 0];
      else
        trl(ll,:)=[(30000*(ll-1))+1 (30000*(ll-1))+30000 0];
      end
    end
    if tmp2.fsample*(fltpadding+trlpadding) + trl(end-1,2) > tmp2.sampleinfo(2)
      trl=trl(1:end-1,:);
    end
    
    for zthresh=[2 3 4 5]
      martifact{zthresh}=[];
      cfg=[];
      cfg.trl=trl(2:end-1,:); % 30s segments
      cfg.continuous='yes';
      cfg.artfctdef.muscle.channel={'FT7' 'FT8' 'AF7' 'AF8' 'TP7' 'TP8'};
      cfg.artfctdef.muscle.bpfilter    = 'yes';
      cfg.artfctdef.muscle.bpfreq      = [100 140];
      cfg.artfctdef.muscle.bpfiltord   = 8;
      cfg.artfctdef.muscle.bpfilttype  = 'but';
      cfg.artfctdef.muscle.hilbert     = 'yes';
      cfg.artfctdef.muscle.interactive = 'yes';
      cfg.artfctdef.muscle.trlpadding  = trlpadding;
      cfg.artfctdef.muscle.fltpadding  = fltpadding;
      cfg.artfctdef.muscle.artpadding  = 1;
      cfg.artfctdef.muscle.cutoff = zthresh;
      [mcfg{zthresh},martifact{zthresh}]=ft_artifact_muscle(cfg,tmp2);
      
      percentkept(zthresh)=nansum(martifact{zthresh}(:,2)-martifact{zthresh}(:,1))/tmp2.sampleinfo(2)
      if sleep==0
        try
          load tmp.mat
        catch
          disp('no tmp file')
          keyboard
          continue
        end
        [mx,mnd]=max([size(artfctdef.visual.artifact,1) size(martifact{zthresh},1)]);
        if mnd==2,
          artfctdef.visual.artifact(size(artfctdef.visual.artifact,1)+1:mx,:)=NaN;
        elseif mnd==1,
          martifact{zthresh}(size(martifact{zthresh},1)+1:mx,:)=NaN;
        end
        [artfctdef.visual.artifact martifact{zthresh}]
      end
      keyboard
    end
    
    if sleep
      save('martifact_sleep.mat', 'martifact','mcfg','percentkept')
    else
      save('martifact.mat', 'martifact','mcfg','percentkept')
    end
  end
end

%%  Notes on zthresh, which seems best for each, for awake data

% check percentkept distance from 0.05.  does this match with 'best' as
% described below?  y means yes.

% e03: 2 too low, 3 too low, 4 tradeoff, 5 tradeoff (better)
% e05: 5; y5
% e06: 3 ok, 4 ok, 5 maybe not enough; y4
% e07: 2 low, 3 ok could be higher, 4 ok, 5 ok; y4;
% e08: 2 little too low, 3 ok and good corr, 4 ok and good corr, 5 too high; y4
% e09: 2 too low, 3 too low, 4 great, 5 slightly too high but ok; n5 but 4great
% e10: 2 good, 3 good, 4 ok slightly too high, 5 too high; n4 but 3great
% e11: 2 too low,3 could be higher, 4 ok, 5 ok or too high; n5 4 so-so
% e12: 2 too low, 3 too low, 4 compromise, 5 ok; y5
% e13: 2 great, 3 too high, 4 too high, 5 too high; y2
% e14: 2 too low, 3 better but too low, 4 fairly good, 5 ok too high; n3 but 4ok
% e15: 2 too low, 3 still too low, 4 ok, 5 ok; y5
% e16: 2 too low, 3 tradeoff, 4 tradeoff, 5 too high; n5 4 not great
% e17: 2 too low, 3 ok, 4 ok, 5 slightly too high; n5 4 so-so
% e18: 2 too low, 3 sort of too low but tradeoff, 4 ok, 5 too high; n5 4 not great
% e19: 2 ok/low, 3 fine/tradeoff, 4 fine, 5 fine; y5
% e20: 2 too low, 3 ok/tradeoff, 4 tradeoff matches prevtmp well, 5 too high; n5 4 so-so
% e21: 2 great almost too high, 3 misses some, 4 misses some, 5 misses some; y2
% e22: 2 too low, 3 tradeoff, 4 tradeoff/too high, 5 too high; y3
% e23: 2 ok/tradeoff, 3 too high and prct agrees, 4 too high, 5 too high; y2
% e24: 2 too low, 3 too low, 4 ok/tradeoff, 5 ok/tradeoff; y4
% e25: 2 ok/tradeoff, 3 ok, 4 ok/too high, 5 too high; n2 3ok
% e26: 2 too low, 3 too low/tradeoff, 4 ok/tradeoff, 5 too high/tradeoff; y4
% e27: 2 too low, 3 ok/good, 4 ok/tradeoff, 5 too high; n5 but 4good
% e28: 2 too low, 3 ok/tradeoff, 4 ok, 5 too high; n3 but 4good
% e29: 2 low/tradeoff, 3 ok/tradeoff, 4 ok/too high, 5 too high; n4 3 so-so
% e30: 2 too low, 3 ok/tradeoff, 4 ok/too high, 5 too high; n5 4 so-so 3 not good
% e31: 2 too low, 3 ok/tradeoff, 4 ok/too high, 5 too high; n4 3 so-so
% e32: 2 ok, 3 ok, 4 too high, 5 too high; n4 3 great

% 13 y
% 7 ok+
% 8 so-so-
musthresh([3 5:32])=[5 5 4 4 4 4 3 4 5 2 4 5 4 4 4 5 4 2 3 2 4 3 4 4 4 3 3 3 3 ];

%% Notes on zthresh, for sleep data

% e03:2 too low, 3 ok/could be higher, 4 good, 5 too high
% e05: 2 tradeoff, 3 ok, 4 ok / too high, 5 too high; n2 3ok
% e06: 2 too low, 3 perfect, 4 too high, 5 too high; n2 3ok
% e07: 2 great, 3 too high, 4 too high, 5 too high; y2
% e08: 2 great, 3 too high, 4 too high, 5 too high; y2
% e09: 2 great, 3 ok/high, 4 too high, 5 too high; n3 2ok
% e10: 2 ok, 3 ok (better?), 4 too high, 5 too high; n2 3ok
% e11: 2 too low, 2 ok / could be higher, 4 good/tradeoff, 5 too high; n5 4so-so
% e12: 2 too low, 3 ok / could be higher, 4 good, 5 too high; n2 4so-so
% e13: ok/tradeoff, 3 ok, 4 too high, 5 too high; n2 3great
% e14: 2 too low, 3 ok / could be higher, 4 ok/ too high, 5 too high;
% e15: 2 too low, 3 ok / could be higher, 4 ok, 5 ok
% e16: 2 ok/could be higher, 3 ok/tradeoff, 4 ok/too high, 5 too high
% e17: 2 good/tradeoff/could be higher, 3 good, 4 ok/tradeoff, 5 too high
% e18: 2 could be higher, 3 ok/tradeoff, 4 ok/tradeoff, 5 ok/tradeoff
% e19: 2 too low, 3 great, 4 ok, 5 ok
% e20: 2 tradeoff, 3 great, 4 great, 5 fine
% e21: 2 too low, 3 tradeoff, 4 good, 5 ok
% e22: 2 great, 3 ok, 4 too high, 5 too high
% e23: 2 perfect/no higher, 3 too high, 4 too high, 5 too high
% e24: 2 too low, 3 great, 4 tradeoff, 5 too high
% e25: 2 great / could be lower!, 3 too high, 4 too high, 5 too high,
% e26: 2 too low/tradeoff, 3 tradeoff/too low, 4 ok/tradeoff, 5 too high
% e27: 2 good tradeoff, 3 too high, 4 too high, 5 too high
% e28: 2 good tradeoff, 3 ok/slightly too high, 4 too high, 5 too high
% e29: 2 good tradeoff, 3 ok/too high, 4 too high, 5 too high
% e30: 2 good tradeoff, 3 good tradeoff, 4 ok/slightly too high, 5 too high
% e31: 2 great / could go lower!, 3 too high, 4 too high, 5 too high
% e32: 2 tradeoff, 3 tradeoff/slightly too high, 4 too high, 5 too high

% based on just above, select best for each pt
persleep =[0.0294    0.0225    0.0275    0.0328    0.0328    0.0520    0.0316    0.0749    0.0144    0.0363    0.0256    0.0452    0.0332
  0.0244    0.0418    0.0607    0.0774    0.0405    0.0428    0.0286    0.0965    0.0702    0.0454    0.1088    0.0415
  0.0326    0.0500    0.0455    0.0315];
% mean(persleep)=    0.0452

musthresh_sleep([3 5:32] )=[4 3 3 2 2  2 3 4 4  3 3 4 3  3 3 3 3  4 2 2 3  2 4 2 2  2 2 2 2];





