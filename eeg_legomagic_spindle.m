% computing spindle segments


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

if strfind(which('butter'),'fieldtrip')
  rmpath(fileparts(which('butter')));
end

load iikeep.mat

%% wake or sleep

for ii=union(iiBuse(iiBuse>7),iiSuse)
% for ii=16:17
  
  
  for sleep=[1 0] % this needs to be done in order of Bed first, since it sets thresholds for detection based on N2
%   for sleep=[0]
    
    clearvars -except ii sub  edir ddir sleep ii*use
    %     pack
    cd([edir sub{ii} ])
    rundetect=1; % run spindle detection by default unless find no stage 2 in sleep1 case
%     runslow=1;  % run slow spindle detection by default unless find no stage 2 in sleep1 case
%     runfast=1;  % run fast spindle detection by default unless find no stage 2 in sleep1 case
    
    if sleep
      files=dir('cc*b.mat');
    else
      files=dir('cc*s.mat');
    end
    
    if exist(['M' files.name])
      mkdir backupM
      movefile(['M' files.name(1:end-3) '*'],'backupM')
    end
    
    S=[];
    S.D=files.name;
    S.refchan={'TP9' 'TP10'};
    [D,S]=spm_eeg_reref_eeg(S);
    
    if sleep==1 && ~any(D.CRC.score{1}==2)
      rundetect=0;
    end
    
    if sleep==0 % use threshold of spindle stage 2 detection from asleep data
      load spindles_sleep.mat;
      try % if sp_* isfield
        slow_threshold=[artfctdef.sp_slow_Fz.thresh artfctdef.sp_slow_Cz.thresh artfctdef.sp_slow_Pz.thresh];
        fast_threshold=[artfctdef.sp_fast_Fz.thresh artfctdef.sp_fast_Cz.thresh artfctdef.sp_fast_Pz.thresh];
      catch
        rundetect=0; %sleep1 case didn't have any stage2 to create threshold so don't run on sleep0 at all (anecdoteally, no one reached N2 in sitting that didn't also reach N2 in bed)
      end
      clear artfctdef crc_def % clearing what was just loaded in this if-statement
    end
    
    if rundetect
      
      
      global crc_def
      crc_defaults;
      % crc_def.sp.highfc          = 8;
      % crc_def.sp.lowfc           = 20;
      % crc_def.sp.elecoi          = {'Fz' 'Cz' 'Pz'};
      % crc_def.sp.stagesp         = [2 3 4];% stages to extract from original scored file
      % crc_def.sp.prcthresh       = 95; %percentile used to perform detection of spindles
      % crc_def.sp.stagethresh     = 2;% stage to compute threshold of detection
      % crc_def.sp.threshold       = 0;% user chosing the threshold of detection
      % crc_def.sp.lengthsp        = 0.4;% spindles of 400 ms duration minimum (in s) [% note, this is ok as it says it misses out the first and last wave so actually about 600ms]
      % crc_def.sp.succsp          = 1;% 1000ms between 2 successive sp (in s)
      % crc_def.sp.type            = 'EEG';% type of channels to detect spindles
      % crc_def.sp.usetheor        = 0;% do not use theoretical positions if localisation file available (1 for always use theoretical positions).
      crc_def.sp.highfc          = 11; % after schabus2007
      crc_def.sp.lowfc           = 13;
      crc_def.sp.stagesp  = [1 2 3 4 5];% stages to extract from original scored file
      crc_def.sp.stagethresh     = [2];% stage to compute threshold of detection
      crc_def.sp.usetheor=1;
      crc_def.sp.succsp          = 0.1;% 100ms between 2 successive sp (in s); % changing from 1000ms to 100ms reduces number spindles by half!
      crc_def.sp.lengthsp        = 0.5;% spindles of 400 ms duration minimum (in s) [% note, this is ok as it says it misses out the first and last wave so actually about 600ms]
      %     crc_def.sp.prcthresh       = 99; %percentile used to perform detection of spindles
      if sleep==0
        crc_def.sp.threshold = slow_threshold;
      end
      
      cfg=[];
      cfg.highfc=[];
      cfg.lowfc=[];
      cfg.review=0;
      cfg.fname=[pwd '\' D.fname];
      cfg.reref=1;
      cfg.scorer=1;
      cfg.analyse=2;
      cfg.Begpts=0;
      cfg.Endpts=1;
      cfg.wav=0;
      cfg.save=0 ;% JZ created this option; rather output D than save to disk
      %     profile on
      D=crc_SP_detect_jz(cfg);
      %     profile viewer
      
      %     load(['M' files.name])
      
      clear artfctdef
      
      % D.CRC.spindles will always exist since the threshold info is saved
      % there even if no spindles detected
      artfctdef.sp_slow_Fz.thresh=D.CRC.spindles.spinthres(match_str(D.CRC.spindles.threshchan,'Fz'));
      artfctdef.sp_slow_Cz.thresh=D.CRC.spindles.spinthres(match_str(D.CRC.spindles.threshchan,'Cz'));
      artfctdef.sp_slow_Pz.thresh=D.CRC.spindles.spinthres(match_str(D.CRC.spindles.threshchan,'Pz'));
      
      %       if isfield(D,'CRC') && isfield(D.CRC,'spindles')
      if cfg.wav==0
        fcnt=1;ccnt=1;pcnt=1;
        if isfield(D.CRC.spindles,'maxelectrode')
          for sp=1:size(D.CRC.spindles.maxelectrode,1)
            if strcmp(D.CRC.spindles.maxelectrode(sp,:),'Fz')
              artfctdef.sp_slow_Fz.artifact(fcnt,:)=D.CRC.spindles.bounds(sp,:);
              artfctdef.sp_slow_Fz.amplitude(fcnt)=D.CRC.spindles.amplitude(sp);
              artfctdef.sp_slow_Fz.duration(fcnt)=D.CRC.spindles.duration(sp);
              fcnt=fcnt+1;
            end
            if strcmp(D.CRC.spindles.maxelectrode(sp,:),'Cz')
              artfctdef.sp_slow_Cz.artifact(ccnt,:)=D.CRC.spindles.bounds(sp,:);
              artfctdef.sp_slow_Cz.amplitude(ccnt)=D.CRC.spindles.amplitude(sp);
              artfctdef.sp_slow_Cz.duration(ccnt)=D.CRC.spindles.duration(sp);
              ccnt=ccnt+1;
            end
            if strcmp(D.CRC.spindles.maxelectrode(sp,:),'Pz')
              artfctdef.sp_slow_Pz.artifact(pcnt,:)=D.CRC.spindles.bounds(sp,:);
              artfctdef.sp_slow_Pz.amplitude(pcnt)=D.CRC.spindles.amplitude(sp);
              artfctdef.sp_slow_Pz.duration(pcnt)=D.CRC.spindles.duration(sp);
              pcnt=pcnt+1;
            end
          end
        end
      elseif cfg.wav==1 %  crc_SP_detect_jz does some extra wavelet analysis of the spindles
        if ~isempty(D.CRC.spindles.index_antsp)
          artfctdef.sp_slow_ant.artifact=D.CRC.spindles.bounds(D.CRC.spindles.index_antsp,:);
        end
        if ~isempty(D.CRC.spindles.index_postsp)
          artfctdef.sp_slow_post.artifact=D.CRC.spindles.bounds(D.CRC.spindles.index_postsp,:);
        end
        if ~isempty(D.CRC.spindles.index_undefsp)
          artfctdef.sp_slow_undef.artifact=D.CRC.spindles.bounds(D.CRC.spindles.index_undefsp,:);
        end
      end
      %       end
      
      crc_def.sp.highfc          = 13; % after schabus2007
      crc_def.sp.lowfc           = 15;
      if sleep==0
        crc_def.sp.threshold = fast_threshold;
      end
      
      
      cfg=[];
      cfg.highfc=[];
      cfg.lowfc=[];
      cfg.review=0;
      cfg.fname=[pwd '\' D.fname];
      cfg.reref=1;
      cfg.scorer=1;
      cfg.analyse=2;
      cfg.Begpts=0;
      cfg.Endpts=1;
      cfg.wav=0;
      cfg.save=0;
      D=crc_SP_detect_jz(cfg);
      
      %     load(['M' files.name])
      
      artfctdef.sp_fast_Fz.thresh=D.CRC.spindles.spinthres(match_str(D.CRC.spindles.threshchan,'Fz'));
      artfctdef.sp_fast_Cz.thresh=D.CRC.spindles.spinthres(match_str(D.CRC.spindles.threshchan,'Cz'));
      artfctdef.sp_fast_Pz.thresh=D.CRC.spindles.spinthres(match_str(D.CRC.spindles.threshchan,'Pz'));
      
      %       if isfield(D,'CRC') && isfield(D.CRC,'spindles')
      if cfg.wav==0
        fcnt=1;ccnt=1;pcnt=1;
        if isfield(D.CRC.spindles,'maxelectrode')
          for sp=1:size(D.CRC.spindles.maxelectrode,1)
            if strcmp(D.CRC.spindles.maxelectrode(sp,:),'Fz')
              artfctdef.sp_fast_Fz.artifact(fcnt,:)=D.CRC.spindles.bounds(sp,:);
              artfctdef.sp_fast_Fz.amplitude(fcnt)=D.CRC.spindles.amplitude(sp);
              artfctdef.sp_fast_Fz.duration(fcnt)=D.CRC.spindles.duration(sp);
              fcnt=fcnt+1;
            end
            if strcmp(D.CRC.spindles.maxelectrode(sp,:),'Cz')
              artfctdef.sp_fast_Cz.artifact(ccnt,:)=D.CRC.spindles.bounds(sp,:);
              artfctdef.sp_fast_Cz.amplitude(ccnt)=D.CRC.spindles.amplitude(sp);
              artfctdef.sp_fast_Cz.duration(ccnt)=D.CRC.spindles.duration(sp);
              ccnt=ccnt+1;
            end
            if strcmp(D.CRC.spindles.maxelectrode(sp,:),'Pz')
              artfctdef.sp_fast_Pz.artifact(pcnt,:)=D.CRC.spindles.bounds(sp,:);
              artfctdef.sp_fast_Pz.amplitude(pcnt)=D.CRC.spindles.amplitude(sp);
              artfctdef.sp_fast_Pz.duration(pcnt)=D.CRC.spindles.duration(sp);
              pcnt=pcnt+1;
            end
          end
        end
      elseif cfg.wav==1
        if ~isempty(D.CRC.spindles.index_antsp)
          artfctdef.sp_fast_ant.artifact=D.CRC.spindles.bounds(D.CRC.spindles.index_antsp,:);
        end
        if ~isempty(D.CRC.spindles.index_postsp)
          artfctdef.sp_fast_post.artifact=D.CRC.spindles.bounds(D.CRC.spindles.index_postsp,:);
        end
        if ~isempty(D.CRC.spindles.index_undefsp)
          artfctdef.sp_fast_undef.artifact=D.CRC.spindles.bounds(D.CRC.spindles.index_undefsp,:);
        end
      end
      %       end
      
      if 0 % visualise results, useful for when testing code, but for routine running, this will be visualised later with other 'artifacts' in eeg_legomagic_viewFeatures.m
        channel={'C3' 'C4' 'Fp2' 'Fp1' 'Fz' 'Oz' 'Cz' 'Pz'};
        if ~exist('raw_reref','var')
          cfg=[];
          cfg.dataset=files.name;
          cfg.demean='yes';
          cfg.channel=[channel {'TP9' 'TP10'}];
          cfg.bsfilter='yes';
          cfg.bsfreq=[49 51; 99 101; 149 151];
          cfg.hpfilter='yes';
          cfg.hpfiltord=3;
          cfg.hpfreq=0.2;  % Could do higher for awake data but not want higher for sleep.
          raw_hpf=ft_preprocessing(cfg);
          
          cfg=[];
          cfg.reref='yes';
          cfg.refchannel={'TP9' 'TP10'};
          raw_reref=ft_preprocessing(cfg,raw_hpf);
        end
        
        cfg=[];
        cfg.layout='EEG1010.lay';
        cfg.artfctdef=artfctdef;
        cfg.viewmode='vertical';
        cfg.blocksize=30;
        %       cfg.preproc.bpfilter='yes';
        %       cfg.preproc.bpfreq=[10 16];
        cfg.preproc.lpfilter='yes';
        cfg.preproc.lpfreq=[20];
        cfg.channel=channel;
        cfg=ft_databrowser(cfg,raw_reref);
      end
      
    else
          artfctdef.sp_slow_Cz=[];
          artfctdef.sp_slow_Fz=[];
          artfctdef.sp_slow_Pz=[];
          artfctdef.sp_fast_Cz=[];
          artfctdef.sp_fast_Fz=[];
          artfctdef.sp_fast_Pz=[];
          crc_def=[];
    end
    
    %     if ~exist('artfctdef','var')
    %       artfctdef=[];
    %     end
    
    delete(['M' files.name(1:end-3) '*'])
    
    if sleep
      save('spindles_sleep.mat', 'artfctdef','crc_def')
    else
      save('spindles.mat', 'artfctdef','crc_def')
    end
  end
  
end




