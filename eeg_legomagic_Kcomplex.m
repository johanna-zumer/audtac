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
load([edir 'iikeep.mat'])

%% wake or sleep

% for ii=union(iiBuse(iiBuse>7),iiSuse)
% for ii=setdiff(union(iiBuse(iiBuse>7),iiSuse),[8:16])
for ii=16:17
  
  
  for sleep=[0 1]
%   for sleep=[1]
    
    clearvars -except ii sub  edir ddir sleep ii*use
%     pack
    cd([edir sub{ii} ])
    if sleep
      files=dir('cc*b.mat');
    else
      files=dir('cc*s.mat');
    end
    
    if exist(['M' files.name]) && ~isdir('backupM')
        mkdir backupM
        movefile(['M' files.name(1:end-3) '*'],'backupM')
    elseif exist(['M' files.name])
      delete(['M' files.name(1:end-3) '*']);
    end
    
    S=[];
    S.D=files.name;
    S.refchan={'TP9' 'TP10'};
    [D,S]=spm_eeg_reref_eeg(S);
    
    if any(isnan(D.CRC.score{1}))
      error('fix scores')
    end
    
    global crc_def
    crc_defaults;
% % % crc_def.swsd.dispscale     = 200; % microV
% % % crc_def.swsd.bpfc          = 200; % lowpass cutoff of bandpass filter
% % % crc_def.swsd.ROI_centers   = ... % position of ROI centers on the 2D unit disk
% % %     [0.5    0.6341;  %Fz
% % %     0.3536 0.5013;  %C3
% % %     0.6464 0.5013;  %C4
% % %     0.5    0.3575]; %Pz
% % % crc_def.swsd.SWlength      = [250 1250 1500 200];
% % % % minimum and maximum time durations between down zero crossing and up zero crossing
% % % crc_def.swsd.SWmAmpl       = [-40 -80 75 140];
% % % % Massimini criteria of magnitude for SWS (-80 140) and delta waves (-40 75)
% % % crc_def.swsd.highfc        = 0.2;  % highpass frequency cutoff
% % % crc_def.swsd.lowfc         = 4;    % lowpass frequency cutoff
% % % crc_def.swsd.stagesw       = [3 4];% stages to extract from original scored file
% % % crc_def.swsd.butterorder   = 4;    % order of butterworth filter
% % % crc_def.swsd.prcentile     = 90;    % 
% % % crc_def.swsd.type          = 'EEG';    % type of channels to detect SWs
% % % crc_def.swsd.usetheor        = 0;% do not use theoretical positions if localisation file available (1 for always use theoretical positions).

    % crc_def.swsd.SWmAmpl       = [-20 -40 37.5 70];
    crc_def.swsd.SWmAmpl       = [-30 -40 50 75 -80 140];
    crc_def.swsd.stagesw       = [0 1 2 3 4 5];% stages to extract from original scored file
    % crc_def.swsd.SWlength      = [250 1250 1500 200]; % assuming in ms
    crc_def.swsd.SWlength      = [250 1500 1500 200]; % assuming in ms
    % minimum and maximum time durations between down zero crossing and up zero crossing
    % SWlength(1): duration of UZC-DZC (c) is greater than this
    % SWlength(2): duration of UZC-DZC (c) is less than this
    % SWlength(3): duration beyond UZC that is relevant for features (e.g. posmax)
    %            : this value divided by 5 for distance between successive detected waves
    % SWlength(4): window either side of posmin to search for min power (reset to SW.negmax)
    crc_def.swsd.prcentile       = 70; % this is to go with JZ way of thresholding, so lower number makes more sense
    crc_def.swsd.usetheor=1;  % JZ has found this is critical to get channel mapping correct
    
    cfg=[];
    cfg.analyse=2;
    cfg.highfc=0.2;
    cfg.lowfc=4;
    cfg.roisel=1;
    cfg.review=0;
    cfg.sensauto=1;
    cfg.sensfname=[];
    cfg.maps=1;
    cfg.fmri=0;
    cfg.Begpts=0;
    cfg.Endpts=1;
    cfg.fname=[pwd '\' D.fname];
    cfg.marker=[];
    cfg.TR=[];
    cfg.reref=1;
    cfg.scorer=1;
    D=crc_SWS_detect_jz(cfg);
    
    clear artfctdef
    kcnt=1;scnt=1;dcnt=1;
    if isfield(D.CRC,'SW')
      for ss=1:length(D.CRC.SW.SW)
%         if floor(D.CRC.SW.SW(ss).code/100)==4
%           artfctdef.delta.artifact(dcnt,1)=D.CRC.SW.SW(ss).down;
%           artfctdef.delta.artifact(dcnt,2)=D.CRC.SW.SW(ss).upend;
%           artfctdef.delta.negmax(dcnt)=D.CRC.SW.SW(ss).negmax;
%           artfctdef.delta.amplitude(dcnt)=D.CRC.SW.SW(ss).amplitude;
%           artfctdef.delta.neg_slope(dcnt)=D.CRC.SW.SW(ss).neg_slope;
%           artfctdef.delta.maxslope(dcnt)=D.CRC.SW.SW(ss).maxslope;
%           dcnt=dcnt+1;
%         elseif floor(D.CRC.SW.SW(ss).code/100)==3
%           artfctdef.kc.artifact(kcnt,1)=D.CRC.SW.SW(ss).down;
%           artfctdef.kc.artifact(kcnt,2)=D.CRC.SW.SW(ss).upend;
%           artfctdef.kc.negmax(kcnt)=D.CRC.SW.SW(ss).negmax;
%           artfctdef.kc.amplitude(kcnt)=D.CRC.SW.SW(ss).amplitude;
%           artfctdef.kc.neg_slope(kcnt)=D.CRC.SW.SW(ss).neg_slope;
%           artfctdef.kc.maxslope(kcnt)=D.CRC.SW.SW(ss).maxslope;
%           kcnt=kcnt+1;
%         elseif floor(D.CRC.SW.SW(ss).code/100)==2
%           artfctdef.sw.artifact(scnt,1)=D.CRC.SW.SW(ss).down;
%           artfctdef.sw.artifact(scnt,2)=D.CRC.SW.SW(ss).upend;
%           artfctdef.sw.negmax(scnt)=D.CRC.SW.SW(ss).negmax;
%           artfctdef.sw.amplitude(scnt)=D.CRC.SW.SW(ss).amplitude;
%           artfctdef.sw.neg_slope(scnt)=D.CRC.SW.SW(ss).neg_slope;
%           artfctdef.sw.maxslope(scnt)=D.CRC.SW.SW(ss).maxslope;
%           scnt=scnt+1;
%         end
        if floor(D.CRC.SW.SW(ss).code/100)==4
          artfctdef.delta.artifact(dcnt,1)=D.CRC.SW.SW(ss).down;
          artfctdef.delta.artifact(dcnt,2)=D.CRC.SW.SW(ss).upend;
          artfctdef.delta.negmax(dcnt)=D.CRC.SW.SW(ss).negmax;
          artfctdef.delta.amplitude(dcnt)=D.CRC.SW.SW(ss).amplitude;
          artfctdef.delta.struct(dcnt)=D.CRC.SW.SW(ss);
          dcnt=dcnt+1;
        elseif floor(D.CRC.SW.SW(ss).code/100)==3
          artfctdef.kc.artifact(kcnt,1)=D.CRC.SW.SW(ss).down;
          artfctdef.kc.artifact(kcnt,2)=D.CRC.SW.SW(ss).upend;
          artfctdef.kc.negmax(kcnt)=D.CRC.SW.SW(ss).negmax;
          artfctdef.kc.amplitude(kcnt)=D.CRC.SW.SW(ss).amplitude;
          artfctdef.kc.struct(kcnt)=D.CRC.SW.SW(ss);
          kcnt=kcnt+1;
        elseif floor(D.CRC.SW.SW(ss).code/100)==2
          artfctdef.sw.artifact(scnt,1)=D.CRC.SW.SW(ss).down;
          artfctdef.sw.artifact(scnt,2)=D.CRC.SW.SW(ss).upend;
          artfctdef.sw.negmax(scnt)=D.CRC.SW.SW(ss).negmax;
          artfctdef.sw.amplitude(scnt)=D.CRC.SW.SW(ss).amplitude;
          artfctdef.sw.struct(scnt)=D.CRC.SW.SW(ss);
          scnt=scnt+1;
        end
      end
    else
      artfctdef=[];
    end
    
    if 0 % visualise results
      if ~exist('raw_reref','var')
        cfg=[];
        cfg.dataset=files.name;
        cfg.demean='yes';
        cfg.channel={'all' '-ECG'};
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
      cfg.preproc.lpfilter='yes';
      cfg.preproc.lpfreq=4;
      %     cfg.channel={'C3' 'C4' 'Fp2' 'Fp1' 'Fz' 'Oz' 'Pz'};
      cfg=ft_databrowser(cfg,raw_reref);
    end
    
    delete(['M' files.name(1:end-3) '*'])
    
    if sleep
      save('Kc_sleep.mat', 'artfctdef','crc_def')
    else
      save('Kc.mat', 'artfctdef','crc_def')
    end
  end
  
end

return;

%% load/view saved results

for ii=[3 5:32]
  for sleep=[0 1]
    clearvars -except ii sub  edir ddir sleep
    cd([edir sub{ii} ])
    if sleep
      load('Kc_sleep.mat')
      files=dir('cc*b.mat');
      load('eogartifact_sleep.mat')
    else
      load('Kc.mat');
      files=dir('cc*s.mat');
      load('eogartifact.mat')
    end
    
    artfctdef.eog.artifact=eogartifact{4};    
    
    cfg=[];
    cfg.dataset=files.name;
    cfg.demean='yes';
    cfg.channel={'all' '-ECG'};
    cfg.bsfilter='yes';
    cfg.bsfreq=[49 51; 99 101; 149 151];
    cfg.hpfilter='yes';
    cfg.hpfiltord=3;
    cfg.hpfreq=0.2;  % Could do higher for awake data but not want higher for sleep.
    raw_hpf=ft_preprocessing(cfg);
    
    cfg=[];cfg.channel={'HEOG' 'VEOG' 'EMG'};
    exg_hpf=ft_selectdata(cfg,raw_hpf);
    
    cfg=[];
    cfg.reref='yes';
    cfg.refchannel={'TP9' 'TP10'};
    cfg.channel='EEG';
    raw_reref=ft_preprocessing(cfg,raw_hpf);
    
    raw_view=ft_appenddata([],raw_reref,exg_hpf);
    
    
    if isfield('raw_hpf.hdr.orig.other','CRC')
    ww=1;n1=1;n2=1;n3=1;rr=1;
    for ss=1:length(raw_hpf.hdr.orig.other.CRC.score{1})
       if raw_hpf.hdr.orig.other.CRC.score{1}(ss)==0
         artfctdef.wake.artifact(ww,:)=[(ss-1)*30*raw_hpf.fsample+100 (ss-1)*30*raw_hpf.fsample+110];
         ww=ww+1;
       elseif raw_hpf.hdr.orig.other.CRC.score{1}(ss)==1
         artfctdef.n1.artifact(n1,:)=[(ss-1)*30*raw_hpf.fsample+100 (ss-1)*30*raw_hpf.fsample+110];
         n1=n1+1;
       elseif raw_hpf.hdr.orig.other.CRC.score{1}(ss)==2
         artfctdef.n2.artifact(n2,:)=[(ss-1)*30*raw_hpf.fsample+100 (ss-1)*30*raw_hpf.fsample+110];
         n2=n2+1;
       elseif raw_hpf.hdr.orig.other.CRC.score{1}(ss)==3
         artfctdef.n3.artifact(n3,:)=[(ss-1)*30*raw_hpf.fsample+100 (ss-1)*30*raw_hpf.fsample+110];
         n3=n3+1;
       elseif raw_hpf.hdr.orig.other.CRC.score{1}(ss)==5
         artfctdef.rem.artifact(rr,:)=[(ss-1)*30*raw_hpf.fsample+100 (ss-1)*30*raw_hpf.fsample+110];
         rr=rr+1;
       end
    end
    end

    % reject waves if they are eyeblinks, but how to choose threshold?
    if 0
    for ff=1:size(artfctdef.delta.artifact,1),
      tmp=corr(raw_view.trial{1}([match_str(raw_view.label,'Fp2') match_str(raw_view.label,'VEOG')],artfctdef.delta.artifact(ff,1):artfctdef.delta.artifact(ff,2))');
      corrff(ff)=tmp(1,2);
    end
    end
    
    cfg=[];
    cfg.layout='EEG1010.lay';
    cfg.artfctdef=artfctdef;
    cfg.viewmode='vertical';
    cfg.blocksize=30;
        cfg.preproc.lpfilter='yes';
        cfg.preproc.lpfreq=4;
    if 0
      cfg.eegscale=3;
    else
      cfg.eegscale=0.25;
      cfg.channel={'C3' 'C4' 'Fp2' 'Fp1' 'Fz' 'Oz' 'Pz' 'VEOG' 'HEOG' 'EMG'};
    end
    cfg=ft_databrowser(cfg,raw_view);

    
    artfctdef=cfg.artfctdef;
    if sleep
      save('Kc_visedit_sleep.mat', 'artfctdef')
    else
      save('Kc_visedit.mat', 'artfctdef')
    end
    
  end  
end




