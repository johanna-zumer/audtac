% load/view saved feature detection

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

% from eeg_legomagic_muscle_artifact.m
musthresh([3 5:32])=[5  5 4 4 4  4 3 4 5  2 4 5 4  4 4 5 4  2 3 2 4  3 4 4 4  3 3 3 3 ];
musthresh_sleep([3 5:32] )=[4  3 3 2 2  2 3 4 4  3 3 4 3  3 3 3 3  4 2 2 3  2 4 2 2  2 2 2 2];
load([edir 'iikeep.mat'])


%%

% for ii=setdiff(union(iiBuse,iiSuse),[3:7 8:19])
% for ii=union(iiBuse(iiBuse>7),iiSuse)
for ii=[ 16 17]
  for sleep=[0 1]
%   for sleep=[0]
    clearvars -except ii sub  edir ddir sleep musthresh* ii*use
    cd([edir sub{ii} ])
    if sleep
      files=dir('cc*b.mat');
      load spindles_sleep.mat
      try
        atmp=load('Kc_visedit_sleep.mat');
      catch
        atmp=load('Kc_sleep.mat');
      end
      load('eogartifact_sleep.mat')
      load martifact_sleep.mat
      
      artfctdef.muscle.artifact=martifact{musthresh_sleep(ii)};
    else
      files=dir('cc*s.mat');
      load spindles.mat
      try
        atmp=load('Kc_visedit.mat');
      catch
        atmp=load('Kc.mat');
      end
      load('eogartifact.mat')
      load martifact.mat
      artfctdef.muscle.artifact=martifact{musthresh(ii)};
    end
    if isfield(atmp.artfctdef,'delta'),artfctdef.delta=atmp.artfctdef.delta;end
    if isfield(atmp.artfctdef,'kc'),artfctdef.kc=atmp.artfctdef.kc;end
    if isfield(atmp.artfctdef,'sw'),artfctdef.sw=atmp.artfctdef.sw;end
    artfctdef.eog.artifact=eogartifact{4};
    
    % artfctdef.sp_fast_Cz should exist with at least .thresh in it even if not any artifacts detected
    if ~isfield(artfctdef,'sp_fast_Cz')
      error('rerun eeg_legomagic_spindle')
    end
    if isfield(artfctdef.sp_fast_Cz','artifact')
      artfctdef.sp_fast.artifact=[artfctdef.sp_fast_Cz.artifact];
      artfctdef.sp_fast.amplitude=[artfctdef.sp_fast_Cz.amplitude];
      artfctdef.sp_fast.duration=[artfctdef.sp_fast_Cz.duration];
    else
      artfctdef.sp_fast.artifact=[];
      artfctdef.sp_fast.amplitude=[];
      artfctdef.sp_fast.duration=[];
    end
    artfctdef=rmfield(artfctdef,'sp_fast_Cz');
    if isfield(artfctdef.sp_fast_Fz,'artifact')
      artfctdef.sp_fast.artifact=[artfctdef.sp_fast.artifact; artfctdef.sp_fast_Fz.artifact];
      artfctdef.sp_fast.amplitude=[artfctdef.sp_fast.amplitude, artfctdef.sp_fast_Fz.amplitude];
      artfctdef.sp_fast.duration=[artfctdef.sp_fast.duration, artfctdef.sp_fast_Fz.duration];
    end
    artfctdef=rmfield(artfctdef,'sp_fast_Fz');
    if isfield(artfctdef.sp_fast_Pz,'artifact')
      artfctdef.sp_fast.artifact=[artfctdef.sp_fast.artifact; artfctdef.sp_fast_Pz.artifact];
      artfctdef.sp_fast.amplitude=[artfctdef.sp_fast.amplitude, artfctdef.sp_fast_Pz.amplitude];
      artfctdef.sp_fast.duration=[artfctdef.sp_fast.duration, artfctdef.sp_fast_Pz.duration];
    end
    artfctdef=rmfield(artfctdef,'sp_fast_Pz');
    if isempty(artfctdef.sp_fast)
      artfctdef=rmfield(artfctdef,'sp_fast');
    end
    
    if isfield(artfctdef.sp_slow_Cz,'artifact')
      artfctdef.sp_slow.artifact=[artfctdef.sp_slow_Cz.artifact];
      artfctdef.sp_slow.amplitude=[artfctdef.sp_slow_Cz.amplitude];
      artfctdef.sp_slow.duration=[artfctdef.sp_slow_Cz.duration];
    else
      artfctdef.sp_slow.artifact=[];
      artfctdef.sp_slow.amplitude=[];
      artfctdef.sp_slow.duration=[];
    end
    artfctdef=rmfield(artfctdef,'sp_slow_Cz');
    if isfield(artfctdef.sp_slow_Fz,'artifact')
      artfctdef.sp_slow.artifact=[artfctdef.sp_slow.artifact; artfctdef.sp_slow_Fz.artifact];
      artfctdef.sp_slow.amplitude=[artfctdef.sp_slow.amplitude, artfctdef.sp_slow_Fz.amplitude];
      artfctdef.sp_slow.duration=[artfctdef.sp_slow.duration, artfctdef.sp_slow_Fz.duration];
    end
    artfctdef=rmfield(artfctdef,'sp_slow_Fz');
    if isfield(artfctdef.sp_slow_Pz,'artifact')
      artfctdef.sp_slow.artifact=[artfctdef.sp_slow.artifact; artfctdef.sp_slow_Pz.artifact];
      artfctdef.sp_slow.amplitude=[artfctdef.sp_slow.amplitude, artfctdef.sp_slow_Pz.amplitude];
      artfctdef.sp_slow.duration=[artfctdef.sp_slow.duration, artfctdef.sp_slow_Pz.duration];
    end
    artfctdef=rmfield(artfctdef,'sp_slow_Pz');
    if isempty(artfctdef.sp_slow)
      artfctdef=rmfield(artfctdef,'sp_slow');
    end
    
    if ii<11
      load ../sleepMontage_combo_noExG.mat
    else
      load ../sleepMontage_combo_ExG.mat
    end
    cfg=[];
    cfg.dataset=files.name;
    cfg.demean='yes';
    cfg.bsfilter='yes';
    cfg.bsfreq=[49 51; 99 101; 149 151];
    cfg.hpfilter='yes';
    cfg.hpfiltord=3;
    cfg.hpfreq=0.2;  % Could do higher for awake data but not want higher for sleep.
    cfg.montage=montage;
    raw_view=ft_preprocessing(cfg);
    
    %     cfg=[];cfg.channel={'HEOG' 'VEOG' 'EMG'};
    %     exg_hpf=ft_selectdata(cfg,raw_hpf);
    %
    %     cfg=[];
    %     cfg.reref='yes';
    %     cfg.refchannel={'TP9' 'TP10'};
    %     cfg.channel='EEG';
    %     raw_reref=ft_preprocessing(cfg,raw_hpf);
    %     raw_view=ft_appenddata([],raw_reref,exg_hpf);
    
    
    
    if isfield(raw_view.hdr.orig.other,'CRC')
      ww=1;nn1=1;nn2=1;nn3=1;rr=1;
      for ss=1:length(raw_view.hdr.orig.other.CRC.score{1})
        if raw_view.hdr.orig.other.CRC.score{1}(ss)==0
          artfctdef.wake.artifact(ww,:)=[(ss-1)*30*raw_view.fsample+100 (ss-1)*30*raw_view.fsample+110];
          %           wake.artifact(ww,:)=[(ss-1)*30*raw_view.fsample+100 (ss-1)*30*raw_view.fsample+110];
          ww=ww+1;
        elseif raw_view.hdr.orig.other.CRC.score{1}(ss)==1
          artfctdef.n1.artifact(nn1,:)=[(ss-1)*30*raw_view.fsample+100 (ss-1)*30*raw_view.fsample+110];
          %           n1.artifact(nn1,:)=[(ss-1)*30*raw_view.fsample+100 (ss-1)*30*raw_view.fsample+110];
          nn1=nn1+1;
        elseif raw_view.hdr.orig.other.CRC.score{1}(ss)==2
          artfctdef.n2.artifact(nn2,:)=[(ss-1)*30*raw_view.fsample+100 (ss-1)*30*raw_view.fsample+110];
          nn2=nn2+1;
        elseif raw_view.hdr.orig.other.CRC.score{1}(ss)==3
          artfctdef.n3.artifact(nn3,:)=[(ss-1)*30*raw_view.fsample+100 (ss-1)*30*raw_view.fsample+110];
          %           n3.artifact(nn3,:)=[(ss-1)*30*raw_view.fsample+100 (ss-1)*30*raw_view.fsample+110];
          nn3=nn3+1;
        elseif raw_view.hdr.orig.other.CRC.score{1}(ss)==5
          artfctdef.REM.artifact(rr,:)=[(ss-1)*30*raw_view.fsample+100 (ss-1)*30*raw_view.fsample+110];
          %           REM.artifact(rr,:)=[(ss-1)*30*raw_view.fsample+100 (ss-1)*30*raw_view.fsample+110];
          rr=rr+1;
        end
      end
    elseif sleep && exist('CRC_sleep.mat')
      keyboard % why got here?
      load CRC_sleep
    elseif ~sleep && exist('CRC_wake.mat')
      keyboard % why got here?
      load CRC_wake
    end
    
    %     % reject waves if they are eyeblinks, but how to choose threshold?
    %     if 0
    %       for ff=1:size(artfctdef.delta.artifact,1),
    %         tmp=corr(raw_view.trial{1}([match_str(raw_view.label,'Fp2') match_str(raw_view.label,'VEOG')],artfctdef.delta.artifact(ff,1):artfctdef.delta.artifact(ff,2))');
    %         corrff(ff)=tmp(1,2);
    %       end
    %     end
    
    %     cfg=[];
    %     cfg.layout='EEG1010.lay';
    %     cfg.artfctdef=artfctdef;
    %     cfg.viewmode='vertical';
    %     cfg.blocksize=30;
    % %         cfg.preproc.lpfilter='yes';
    % %         cfg.preproc.lpfreq=20;
    %     cfg=ft_databrowser(cfg,raw_view);
    
    %%
    
    % automatically correct confusion over EOG vs Kc/SW
    
    %     artfctdef.wake=wake;
    %     artfctdef.n1=n1;
    %     artfctdef.n3=n3;
    %     if exist('r')
    %       artfctdef.r=r;
    %     end
    
    artnames=fieldnames(artfctdef);
    for aa=1:length(artnames)
      if isfield(artfctdef.(artnames{aa}),'artifact')
        % first, some/rare situation they might be NaN
        if ~isempty(artfctdef.(artnames{aa}).artifact)
          aakeep=all(~isnan(artfctdef.(artnames{aa}).artifact)');
          artfctdef.(artnames{aa}).artifact=artfctdef.(artnames{aa}).artifact(aakeep,:);
          if strcmp(artnames{aa},'kc') || strcmp(artnames{aa},'delta') || strcmp(artnames{aa},'sw')
            artfctdef.(artnames{aa}).negmax=artfctdef.(artnames{aa}).negmax(aakeep);
            artfctdef.(artnames{aa}).amplitude=artfctdef.(artnames{aa}).amplitude(aakeep);
            artfctdef.(artnames{aa}).struct=artfctdef.(artnames{aa}).struct(aakeep);
          elseif strcmp(artnames{aa},'sp_fast') || strcmp(artnames{aa},'sp_slow')
            artfctdef.(artnames{aa}).duration=artfctdef.(artnames{aa}).duration(aakeep);
            artfctdef.(artnames{aa}).amplitude=artfctdef.(artnames{aa}).amplitude(aakeep);
          end
        end
        bin.(strcat(artnames{aa},'_bin'))=binarise_artifact_begendpoints(artfctdef.(artnames{aa}).artifact,raw_view.sampleinfo(2));
      end
    end
    stage_bin=zeros(size(bin.eog_bin));
    if isfield(bin,'wake_bin')
      stage_bin=stage_bin +bin.wake_bin*-1;
    end
    if isfield(bin,'n1_bin')
      stage_bin=stage_bin +bin.n1_bin*10;
    end
    if isfield(bin,'n2_bin')
      stage_bin=stage_bin +bin.n2_bin*20;
    end
    if isfield(bin,'n3_bin')
      stage_bin=stage_bin +bin.n3_bin*30;
    end
    if isfield(bin,'REM_bin')
      stage_bin=stage_bin +bin.REM_bin*50;
    end
    
    % QUESTION:
    % automatically remove all EOG from N2 or N3  ??? (then don't need
    % more explicit rejection later, below then would only be for rejecting
    % kc/delta/sw from W or REM if overlaps with EOG.
    % ANSWER:  NO!  I have seen examples of eye movements at end of N2 to
    % indicate awakening
    % BUT: have created new category called 'sleepblink' since sometimes an
    % ambiguous more-eye-than-frontal wave occurs
    
    % they are not in order as the were found over 4 different sensor clusters
    if isfield(artfctdef,'delta')
      [sss,iii]=sort(artfctdef.delta.artifact(:,1));
      artfctdef.delta.artifact=artfctdef.delta.artifact(iii,:);
      artfctdef.delta.negmax=artfctdef.delta.negmax(iii);
      artfctdef.delta.amplitude=artfctdef.delta.amplitude(iii);
      artfctdef.delta.struct=artfctdef.delta.struct(iii);
    end
    if isfield(artfctdef,'kc')
      [sss,iii]=sort(artfctdef.kc.artifact(:,1));
      artfctdef.kc.artifact=artfctdef.kc.artifact(iii,:);
      artfctdef.kc.negmax=artfctdef.kc.negmax(iii);
      artfctdef.kc.amplitude=artfctdef.kc.amplitude(iii);
      artfctdef.kc.struct=artfctdef.kc.struct(iii);
    end
    if isfield(artfctdef,'sw')
      [sss,iii]=sort(artfctdef.sw.artifact(:,1));
      artfctdef.sw.artifact=artfctdef.sw.artifact(iii,:);
      artfctdef.sw.negmax=artfctdef.sw.negmax(iii);
      artfctdef.sw.amplitude=artfctdef.sw.amplitude(iii);
      artfctdef.sw.struct=artfctdef.sw.struct(iii);
    end
    
    
    % reject for EOG & Kc
    if isfield(bin,'kc_bin') && isfield(bin,'eog_bin')
      eogkc=find(bin.kc_bin & bin.eog_bin);
      eogkcstage=zeros(1,length(eogkc));
      kcreject=[];
      eogreject=[];
      kcnotreject=[];
      eognotreject=[];
      for kk=1:length(eogkc)
        if ~isempty(setdiff(unique(stage_bin(eogkc(kk):eogkc(kk)+100)),0))
          eogkcstage(kk)=setdiff(unique(stage_bin(eogkc(kk):eogkc(kk)+100)),0);
        else
          eogkcstage(kk)=setdiff(unique(stage_bin(eogkc(kk):-1:max(eogkc(kk)-30000+100,1))),0);
        end
        [mn,mnd]=min(eogkc(kk)-artfctdef.eog.artifact(artfctdef.eog.artifact(:,1)<=eogkc(kk),1));
        if eogkcstage(kk)==20 || eogkcstage(kk)==30 % keep Kc and remove EOG
          eogreject=[eogreject mnd];
        else
          eognotreject=[eognotreject mnd];
        end
        
        [mn,mnd]=min(eogkc(kk)-artfctdef.kc.artifact(artfctdef.kc.artifact(:,1)<=eogkc(kk),1));
        if eogkcstage(kk)==-1 || eogkcstage(kk)==50 % keep EOG and remove Kc
          kcreject=[kcreject mnd];
        else
          kcnotreject=[kcnotreject mnd];
        end
      end
      
      rej=intersect(unique(eognotreject),unique(eogreject));
      clear rejstage
      if ~isempty(rej) % can happen if an artifact straddles a 30s seg break
        disp('nonempty intersection')
        for jj=1:length(rej)
          rejstage(jj)=unique(stage_bin(artfctdef.eog.artifact(rej(jj),1)-1+find(stage_bin(artfctdef.eog.artifact(rej(jj),1):artfctdef.eog.artifact(rej(jj),2)+100))));
          if rejstage(jj)==20 || rejstage(jj)==30
            eognotreject=setdiff(eognotreject,rej(jj)); % rej(jj) shoudl stay in reject
          else
            eogreject=setdiff(eogreject,rej(jj)); % rej(jj) should stay in notreject
          end
        end
      end
      rej=intersect(unique(kcnotreject),unique(kcreject));
      clear rejstage
      if ~isempty(rej) % can happen if an artifact straddles a 30s seg break
        disp('nonempty intersection')
        for jj=1:length(rej)
          rejstage(jj)=unique(stage_bin(artfctdef.kc.artifact(rej(jj),1)-1+find(stage_bin(artfctdef.kc.artifact(rej(jj),1):artfctdef.kc.artifact(rej(jj),2)+100))));
          if rejstage(jj)==-1 || rejstage(jj)==50
            kcnotreject=setdiff(kcnotreject,rej(jj)); % rej(jj) shoudl stay in reject
          else
            kcreject=setdiff(kcreject,rej(jj)); % rej(jj) should stay in notreject
          end
        end
      end
      
      artfctdef.eog.artifact(unique(eogreject),:)=[];
      artfctdef.kc.artifact(unique(kcreject),:)=[];
      artfctdef.kc.negmax(unique(kcreject))=[];
      artfctdef.kc.amplitude(unique(kcreject))=[];
      artfctdef.kc.struct(unique(kcreject))=[];
      if size(artfctdef.kc.artifact,1)==0
        artfctdef=rmfield(artfctdef,'kc');
      end
    end
    bin.eog_bin=binarise_artifact_begendpoints(artfctdef.eog.artifact,raw_view.sampleinfo(2));
    
    % reject for EOG & sw
    if isfield(bin,'sw_bin') && isfield(bin,'eog_bin')
      eogsw=find(bin.sw_bin & bin.eog_bin);
      eogswstage=zeros(1,length(eogsw));
      swreject=[];
      eogreject=[];
      swnotreject=[];
      eognotreject=[];
      for kk=1:length(eogsw)
        if ~isempty(setdiff(unique(stage_bin(eogsw(kk):eogsw(kk)+100)),0))
          eogswstage(kk)=setdiff(unique(stage_bin(eogsw(kk):eogsw(kk)+100)),0);
        else
          eogswstage(kk)=setdiff(unique(stage_bin(eogsw(kk):-1:max(eogsw(kk)-30000+100,1))),0);
        end
        [mn,mnd]=min(eogsw(kk)-artfctdef.eog.artifact(artfctdef.eog.artifact(:,1)<=eogsw(kk),1));
        if eogswstage(kk)==20 || eogswstage(kk)==30 % keep sw and remove EOG
          eogreject=[eogreject mnd];
        else
          eognotreject=[eognotreject mnd];
        end
        
        [mn,mnd]=min(eogsw(kk)-artfctdef.sw.artifact(artfctdef.sw.artifact(:,1)<=eogsw(kk),1));
        if eogswstage(kk)==-1 || eogswstage(kk)==50 % keep EOG and remove sw
          swreject=[swreject mnd];
        else
          swnotreject=[swnotreject mnd];
        end
      end
      rej=intersect(unique(eognotreject),unique(eogreject));
      clear rejstage
      if ~isempty(rej) % can happen if an artifact straddles a 30s seg break
        disp('nonempty intersection')
        for jj=1:length(rej)
          rejstage(jj)=unique(stage_bin(artfctdef.eog.artifact(rej(jj),1)-1+find(stage_bin(artfctdef.eog.artifact(rej(jj),1):artfctdef.eog.artifact(rej(jj),2)+100))));
          if rejstage(jj)==20 || rejstage(jj)==30
            eognotreject=setdiff(eognotreject,rej(jj)); % rej(jj) shoudl stay in reject
          else
            eogreject=setdiff(eogreject,rej(jj)); % rej(jj) should stay in notreject
          end
        end
      end
      rej=intersect(unique(swnotreject),unique(swreject));
      clear rejstage
      if ~isempty(rej) % can happen if an artifact straddles a 30s seg break
        disp('nonempty intersection')
        for jj=1:length(rej)
          rejstage(jj)=unique(stage_bin(artfctdef.sw.artifact(rej(jj),1)-1+find(stage_bin(artfctdef.sw.artifact(rej(jj),1):artfctdef.sw.artifact(rej(jj),2)+100))));
          if rejstage(jj)==-1 || rejstage(jj)==50
            swnotreject=setdiff(swnotreject,rej(jj)); % rej(jj) shoudl stay in reject
          else
            swreject=setdiff(swreject,rej(jj)); % rej(jj) should stay in notreject
          end
        end
      end
      
      artfctdef.eog.artifact(unique(eogreject),:)=[];
      artfctdef.sw.artifact(unique(swreject),:)=[];
      artfctdef.sw.negmax(unique(swreject))=[];
      artfctdef.sw.amplitude(unique(swreject))=[];
      artfctdef.sw.struct(unique(swreject))=[];
      if size(artfctdef.sw.artifact,1)==0
        artfctdef=rmfield(artfctdef,'sw');
      end
    end
    bin.eog_bin=binarise_artifact_begendpoints(artfctdef.eog.artifact,raw_view.sampleinfo(2));
    
    % reject for EOG & delta
    if isfield(bin,'delta_bin') && isfield(bin,'eog_bin')
      eogdelta=find(bin.delta_bin & bin.eog_bin);
      eogdeltastage=zeros(1,length(eogdelta));
      deltareject=[];
      deltanotreject=[];
      eogreject=[];
      eognotreject=[];
      for kk=1:length(eogdelta)
        if ~isempty(setdiff(unique(stage_bin(eogdelta(kk):eogdelta(kk)+100)),0))
          eogdeltastage(kk)=setdiff(unique(stage_bin(eogdelta(kk):eogdelta(kk)+100)),0);
        else
          eogdeltastage(kk)=setdiff(unique(stage_bin(eogdelta(kk):-1:max(eogdelta(kk)-30000+100,1))),0);
        end
        [mn,mnd]=min(eogdelta(kk)-artfctdef.eog.artifact(artfctdef.eog.artifact(:,1)<=eogdelta(kk),1));
        if eogdeltastage(kk)==20 || eogdeltastage(kk)==30 % keep delta and remove EOG
          eogreject=[eogreject mnd];
        else
          eognotreject=[eognotreject mnd];
        end
        
        [mn,mnd]=min(eogdelta(kk)-artfctdef.delta.artifact(artfctdef.delta.artifact(:,1)<=eogdelta(kk),1));
        if eogdeltastage(kk)==-1 || eogdeltastage(kk)==50 % keep EOG and remove delta
          deltareject=[deltareject mnd];
        else
          deltanotreject=[deltanotreject mnd];
        end
      end
      
      rej=intersect(unique(eognotreject),unique(eogreject));
      clear rejstage
      if ~isempty(rej) % can happen if an artifact straddles a 30s seg break
        disp('nonempty intersection')
        for jj=1:length(rej)
          rejstage(jj)=unique(stage_bin(artfctdef.eog.artifact(rej(jj),1)-1+find(stage_bin(artfctdef.eog.artifact(rej(jj),1):artfctdef.eog.artifact(rej(jj),2)+100))));
          if rejstage(jj)==20 || rejstage(jj)==30
            eognotreject=setdiff(eognotreject,rej(jj)); % rej(jj) shoudl stay in reject
          else
            eogreject=setdiff(eogreject,rej(jj)); % rej(jj) should stay in notreject
          end
        end
      end
      rej=intersect(unique(deltanotreject),unique(deltareject));
      clear rejstage
      if ~isempty(rej) % can happen if an artifact straddles a 30s seg break
        disp('nonempty intersection')
        for jj=1:length(rej)
          rejstage(jj)=unique(stage_bin(artfctdef.delta.artifact(rej(jj),1)-1+find(stage_bin(artfctdef.delta.artifact(rej(jj),1):artfctdef.delta.artifact(rej(jj),2)+100))));
          if rejstage(jj)==-1 || rejstage(jj)==50
            deltanotreject=setdiff(deltanotreject,rej(jj)); % rej(jj) shoudl stay in reject
          else
            deltareject=setdiff(deltareject,rej(jj)); % rej(jj) should stay in notreject
          end
        end
      end
      
      artfctdef.eog.artifact(unique(eogreject),:)=[];
      artfctdef.delta.artifact(unique(deltareject),:)=[];
      artfctdef.delta.negmax(unique(deltareject))=[];
      artfctdef.delta.amplitude(unique(deltareject))=[];
      artfctdef.delta.struct(unique(deltareject))=[];
      if size(artfctdef.delta.artifact,1)==0
        artfctdef=rmfield(artfctdef,'delta');
      end
    end
    
    % changing EOG during N2 or N3 to a new category 'sleepblink' (which
    % might be later kept or rejected)
    if size(artfctdef.eog.artifact,1)==0
      artfctdef=rmfield(artfctdef,'eog');
    else
      bin.eog_bin=binarise_artifact_begendpoints(artfctdef.eog.artifact,raw_view.sampleinfo(2));
      
      eogfind=find(bin.eog_bin);
      eogreject=[];
      eognotreject=[];
      for kk=1:length(eogfind)
        if ~isempty(setdiff(unique(stage_bin(eogfind(kk):eogfind(kk)+100)),0))
          eogfindstage(kk)=setdiff(unique(stage_bin(eogfind(kk):eogfind(kk)+100)),0);
        else
          eogfindstage(kk)=setdiff(unique(stage_bin(eogfind(kk):-1:max(eogfind(kk)-30000+100,1))),0);
        end
        [mn,mnd]=min(eogfind(kk)-artfctdef.eog.artifact(artfctdef.eog.artifact(:,1)<=eogfind(kk),1));
        if eogfindstage(kk)==20 || eogfindstage(kk)==30 % change EOG to sleepblink
          eogreject=[eogreject mnd];
        else
          eognotreject=[eognotreject mnd];
        end
      end
      
      rej=intersect(unique(eognotreject),unique(eogreject));
      clear rejstage
      if ~isempty(rej) % can happen if an artifact straddles a 30s seg break
        disp('nonempty intersection')
        for jj=1:length(rej)
          rejstage(jj)=unique(stage_bin(artfctdef.eog.artifact(rej(jj),1)-1+find(stage_bin(artfctdef.eog.artifact(rej(jj),1):artfctdef.eog.artifact(rej(jj),2)+100))));
          if rejstage(jj)==20 || rejstage(jj)==30 % change EOG to sleepblink
            eognotreject=setdiff(eognotreject,rej(jj)); % rej(jj) shoudl stay in reject
          else
            eogreject=setdiff(eogreject,rej(jj)); % rej(jj) should stay in notreject
          end
        end
      end
      
      if ~isempty(eogreject)
        artfctdef.sleepblink.artifact=artfctdef.eog.artifact(unique(eogreject),:);
        artfctdef.eog.artifact(unique(eogreject),:)=[];
        if size(artfctdef.eog.artifact,1)==0
          artfctdef=rmfield(artfctdef,'eog');
        end
      end
    end
    
    % hope this isn't too strict but if the session never enters N2 *at all*
    % then there aren't really any true spindles
    if ~isfield(artfctdef,'n2') && isfield(artfctdef,'sp_slow')
      artfctdef=rmfield(artfctdef,'sp_slow');
    end
    if ~isfield(artfctdef,'n2') && isfield(artfctdef,'sp_fast')
      artfctdef=rmfield(artfctdef,'sp_fast');
    end
    
    % databrowser only wants up to 9 artifacts to display at once
    if length(fieldnames(artfctdef))>9
      if isfield(artfctdef,'REM')
        artfctdef=rmfield(artfctdef,'REM');
      end
    end
    if length(fieldnames(artfctdef))>9
      if isfield(artfctdef,'n3')
        artfctdef=rmfield(artfctdef,'n3');
      end
    end
    if length(fieldnames(artfctdef))>9
      if isfield(artfctdef,'wake')
        artfctdef=rmfield(artfctdef,'wake');
      end
    end
    if length(fieldnames(artfctdef))>9
      if isfield(artfctdef,'n1')
        artfctdef=rmfield(artfctdef,'n1');
      end
    end
    
    %     % keep correct feature stats for the BigWaves
    %     if isfield(artfctdef,'delta')
    %       featkeep=[];
    %       for dd=1:size(artfctdef.delta.artifact,1)
    %         featkeep=[featkeep find([artfctdef.delta.negmax-artfctdef.delta.artifact(dd,1)]<1000 & [artfctdef.delta.negmax-artfctdef.delta.artifact(dd,1)]>0) ];
    %       end
    %       artfctdef.delta.negmax=artfctdef.delta.negmax(featkeep);
    %       artfctdef.delta.amplitude=artfctdef.delta.amplitude(featkeep);
    %     end
    %     if isfield(artfctdef,'kc')
    %       featkeep=[];
    %       for dd=1:size(artfctdef.kc.artifact,1)
    %         featkeep=[featkeep find([artfctdef.kc.negmax-artfctdef.kc.artifact(dd,1)]<2000 & [artfctdef.kc.negmax-artfctdef.kc.artifact(dd,1)]>0) ];
    %       end
    %       artfctdef.kc.negmax=artfctdef.kc.negmax(featkeep);
    %       artfctdef.kc.amplitude=artfctdef.kc.amplitude(featkeep);
    %     end
    %     if isfield(artfctdef,'sw')
    %       featkeep=[];
    %       for dd=1:size(artfctdef.sw.artifact,1)
    %         featkeep=[featkeep find([artfctdef.sw.negmax-artfctdef.sw.artifact(dd,1)]<1000 & [artfctdef.sw.negmax-artfctdef.sw.artifact(dd,1)]>0) ];
    %       end
    %       artfctdef.sw.negmax=artfctdef.sw.negmax(featkeep);
    %       artfctdef.sw.amplitude=artfctdef.sw.amplitude(featkeep);
    %     end
    
    % automatic detection
    if sleep
      save('artfctdef_bed_auto.mat','artfctdef')
    else
      save('artfctdef_sit_auto.mat','artfctdef')
    end
    
    
  end
end

return
%%

for sleep=[0 1]
  if sleep
    iiuse=iiBuse;
  else
    iiuse=iiSuse;
  end
  
  for ii=setdiff(iiuse,8)
    %   for ii=iiuse
    clearvars -except ii sub  edir  sleep ii*use
    cd([edir sub{ii} ])
    
    if sleep
      files=dir('cc*b.mat');
    else
      files=dir('cc*s.mat');
    end
    
    if sleep
      load('artfctdef_bed_auto.mat')
    else
      load('artfctdef_sit_auto.mat')
    end
    
    if ii<11
      load ../sleepMontage_combo_noExG.mat
    else
      load ../sleepMontage_combo_ExG.mat
    end
    cfg=[];
    cfg.dataset=files.name;
    cfg.demean='yes';
    cfg.bsfilter='yes';
    cfg.bsfreq=[49 51; 99 101; 149 151];
    cfg.hpfilter='yes';
    cfg.hpfiltord=3;
    cfg.hpfreq=0.2;  % Could do higher for awake data but not want higher for sleep.
    cfg.montage=montage;
    raw_view=ft_preprocessing(cfg);
    
    % trust muscle, kc, sw, delta, sp_slow, sp_fast
    
    while size(fieldnames(artfctdef),1)>8
      artfctdef_orig=artfctdef;
      try
        artfctdef=rmfield(artfctdef,'n2');
        continue
      catch
      end
    end
    
    % must still moderate eog and sleepblink
    cfg=[];
    cfg.layout='EEG1010.lay';
    cfg.artfctdef=artfctdef;
    cfg.viewmode='vertical';
    cfg.selectfeature='eog';
    cfg.blocksize=30;
    %         cfg.preproc.lpfilter='yes';
    %         cfg.preproc.lpfreq=20;
    cfg=ft_databrowser(cfg,raw_view);
    
    % this is now manual tweak of EOG and 'sleepblink'
    artfctdef_new=cfg.artfctdef;
    if sleep
      save('artfctdef_bed_automan.mat','artfctdef_new')
    else
      save('artfctdef_sit_automan.mat','artfctdef_new')
    end
    
  end % ii
  
end % sleep


