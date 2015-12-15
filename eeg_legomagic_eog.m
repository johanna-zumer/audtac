% computing eog artifact segments (saving them to be rejected later)


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
load iikeep.mat

%% wake or sleep

for ii=intersect(iiSuse,iiBuse)
% for ii=30
  for sleep=[0 1]
    clearvars -except ii sub  edir ddir sleep
%     pack
    cd([edir sub{ii} ])
    if sleep
      files=dir('cc*b.mat');
    else
      files=dir('cc*s.mat');
    end
    
    cfg=[];
    cfg.dataset=files.name;
    cfg.demean='yes';
    cfg.channel={'VEOG' 'HEOG' 'F7' 'F8' 'Fp2' 'FT9' 'Fp1' 'FT10'};
    cfg.bsfilter='yes';
    cfg.bsfreq=[49 51; 99 101; 149 151];
    cfg.hpfilter='yes';
    cfg.hpfiltord=3;
    cfg.hpfreq=0.2;  % Could do higher for awake data but not want higher for sleep.
    raw_hpf=ft_preprocessing(cfg);
    
    
    fltpadding=0.3;
    trlpadding=1;
    trllen=20000;
    % create trl with 30s segments
    for ll=1:ceil( raw_hpf.sampleinfo(2)/trllen)
      if ll==ceil( raw_hpf.sampleinfo(2)/trllen)
        trl(ll,:)=[(trllen*(ll-1))+1 raw_hpf.sampleinfo(2) 0];
      elseif ll==1
        trl(ll,:)=[200 (trllen*(ll-1))+trllen 0];
      else
        trl(ll,:)=[(trllen*(ll-1))+1 (trllen*(ll-1))+trllen 0];
      end
    end
    if raw_hpf.fsample*(fltpadding+trlpadding) + trl(end-1,2) > raw_hpf.sampleinfo(2)
      trl=trl(1:end-1,:);
    end
    
    load([edir 'sleepmontageT.mat']);
    
    montage1.labelorg=montage.labelorg(match_str(montage.labelorg,raw_hpf.label));
    montage1.labelnew=montage.labelnew(match_str(montage.labelnew,{ 'RVEOG Fp2-FT9'    'LVEOG Fp1-FT10'  'HEOG F7-F8' 'VEOG' 'HEOG'}));
    montage1.tra=montage.tra(match_str(montage.labelnew,{ 'RVEOG Fp2-FT9'    'LVEOG Fp1-FT10'  'HEOG F7-F8' 'VEOG' 'HEOG'}),match_str(montage.labelorg,raw_hpf.label));
    
    
    cfg=[];
    cfg.montage=montage1;
    raw_reref=ft_preprocessing(cfg,raw_hpf);

    %     raw_all=ft_appenddata([],raw_hpf,raw_reref);
    if ii<11
      channel{1}={'HEOG F7-F8' 'RVEOG Fp2-FT9'  'LVEOG Fp1-FT10'};
    else
      cfg=[];
      cfg.channel={'VEOG' 'HEOG'};
      raw_eog=ft_selectdata(cfg,raw_hpf);
      raw_reref=ft_appenddata([],raw_reref,raw_eog);
      channel{1}={'VEOG' 'HEOG F7-F8' 'RVEOG Fp2-FT9'  'LVEOG Fp1-FT10'}; % HEOG contaminated by tapper on cheek
      channel{2}={'VEOG' 'HEOG' 'HEOG F7-F8' 'RVEOG Fp2-FT9'  'LVEOG Fp1-FT10'}; % try HEOG anyway
    end
    
    for cc=1:length(channel)
      for zthresh=[3 4 5 6]
        eogartifact{zthresh}=[];
        cfg=[];
        cfg.trl=trl(2:end-1,:); % 30s segments
        cfg.continuous='yes';
        cfg.artfctdef.eog.channel=channel{cc};
        cfg.artfctdef.eog.bpfilter    = 'yes';
        cfg.artfctdef.eog.bpfreq      = [1 16]; % this spans Kc and spindles
        cfg.artfctdef.eog.bpfiltord   = 3;
        cfg.artfctdef.eog.bpfilttype  = 'but';
        cfg.artfctdef.eog.hilbert     = 'yes';
        %       cfg.artfctdef.eog.interactive = 'yes';
        cfg.artfctdef.eog.trlpadding  = trlpadding;
        cfg.artfctdef.eog.fltpadding  = fltpadding;
        cfg.artfctdef.eog.artpadding  = 0.2;
        cfg.artfctdef.eog.cutoff = zthresh;
        %     if ii<11
        %       [eogcfg{zthresh},eogartifact{zthresh}]=ft_artifact_eog(cfg,raw_reref);
        %     else
        [eogcfg{zthresh,cc},eogartifact{zthresh,cc}]=ft_artifact_eog(cfg,raw_reref);
        %     end
        
        percentkept(zthresh,cc)=nansum(eogartifact{zthresh,cc}(:,2)-eogartifact{zthresh,cc}(:,1))/raw_hpf.sampleinfo(2)
        %       keyboard
      end
    end
    
    if sleep
      save('eogartifact_sleep.mat', 'eogartifact','eogcfg','percentkept')
    else
      save('eogartifact.mat', 'eogartifact','eogcfg','percentkept')
    end
  end
  
end

%%  Notes on zthresh, which seems best for each, for awake data

%% Notes on zthresh, for sleep data




