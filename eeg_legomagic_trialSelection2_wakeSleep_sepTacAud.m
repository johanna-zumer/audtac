function [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection2_wakeSleep_sepTacAud(ii,sleep,ttin,tacaud,varargin)
% function [tlock_tac,tlock_aud,tlock_nul]=eeg_legomagic_trialSelection(ii)
% Trial selection: equal trial numbers per condition and removing NaN at ends
%
% ii       subject index
% sleep    sleep=0 if sitting, sleep=1 if bed
% ttin     index for threshold of which tactile latency jitter to allow
% tacaud   flag for whether doing TacPlusAud (1) or AudPlusTac (0) which
%          has to do with which stimulus is at latency time=0
tomflag=0; % set to 1 if you are Tom (0 for Johanna)
if ispc
  if tomflag
    edir='G:\jz_sleep_data\data\';
    ddir='E:\JohannaProject\diaries\';
  else
    edir='D:\audtac\eeg_data\';
    ddir='D:\audtac\legomagic\diaries\';
  end
  bdir='D:\audtac\behav_data\';
  sdir='D:\audtac\spss_stuff\';
else
  [~,hostname]=system('hostname');
  if ~isempty(strfind(hostname,'les')) | ~isempty(strfind(hostname,'LES')) % either COLLES-151401 or LES-LINUX_FS3
    edir='/home/zumerj/audtac/eeg_data/';
    ddir='/home/zumerj/audtac/legomagic/diaries/';
    bdir='/home/zumerj/audtac/behav_data/';
    sdir='/home/zumerj/audtac/spss_stuff/';
  else % assume on VM linux of psychl-132432
    edir='/mnt/hgfs/D/audtac/eeg_data/';
    ddir='/mnt/hgfs/D/audtac/legomagic/diaries/';
    bdir='/mnt/hgfs/D/audtac/behav_data/';
    sdir='/mnt/hgfs/D/audtac/spss_stuff/';
  end
end

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

cd([edir sub{ii} ])
%%

if nargin>4
  data_tac_ref=varargin{1};
  data_aud_ref=varargin{2};
  data_nul_ref=varargin{3};
  iter=varargin{4};
  usetr=varargin{5};
  trialkc=varargin{6};  % -1 for ignore delta/kc/sw value (e.g. awake or collapse over all); 0 for non-Kc trials and 1 for Kc-trials
else
  % load(['raw_each_rej_' sub{ii}],'raw_*'); % see eeg_legomagic_epoching.m for difference of these files
  % load(['raw1_each_rej_' sub{ii}],'raw_*'); % data has 0.2Hz highpass filter applied before epoching
  load(['raw2_each_rej_' sub{ii} '_sleep' num2str(sleep)],'raw_*'); % data has 0.2Hz highpass filter applied before epoching
  
  if strcmp(raw_tac.cfg.previous.previous.reref,'no')
    cfg=[];
    cfg.channel= {'all', '-ECG'};
    cfg.reref='yes';
    cfg.refchannel={'all', '-ECG'};
    data_tac_ref=ft_preprocessing(cfg,raw_tac);
    data_aud_ref=ft_preprocessing(cfg,raw_aud);
    data_nul_ref=ft_preprocessing(cfg,raw_nul);
    clear raw_*
  else
    data_tac_ref=raw_tac;
    data_aud_ref=raw_aud;
    data_nul_ref=raw_nul;
    clear raw_*
  end
end
clear varargin

if ii<8
  soalist=[3 4 5 6 7];
else
  soalist=[1 3 4 5 6 7 9];
end


if 1
  %   stageuse=[0:3 5];
  if sleep
    stageuse=[0:2];
  else
    stageuse=[0:1];
  end
else
  if sleep
    stageuse=[2 3];
  else
    stageuse=[0 ];
  end
end


% only loop through stages actually present in the given dataset
for ss=stageuse
  if ~any(data_tac_ref.trialinfo(:,2)==ss) || ~any(data_aud_ref.trialinfo(:,2)==ss) || ~any(data_nul_ref.trialinfo(:,2)==ss)
    stageuse=setdiff(stageuse,ss);
  end
end
tr.stageuse=stageuse;
tr.trialkc=trialkc;

if usetr
  if iter>30
    error('all the usetr flag settings ok for iter<=30. what to do here?')
  end
  try
    load(['trialkept_tt' num2str(ttin) '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_usetr' num2str(0) '_trialkc' num2str(trialkc) '.mat'],'tr');
  catch ME
    disp(ME.message)
    load(['trialkept_tt' num2str(ttin) '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '.mat'],'tr');
  end
  stageuse=tr.stageuse;
else
  for ss=0:3
    if ~any(ismember(stageuse,ss)) % needs to exist for later code in eeg_legomagic_erp_stats2_sepTacAud
      for ll=soalist
        if tacaud
          tr.t10trialkept{ll,ttin,ss+10}=0;
          tr.nllttrialkept{ll,ttin,ss+10}=0;
          tr.tlltrialkept{ll,ttin,ss+10}=0;
          tr.all40trialkept{ll,ttin,ss+10}=0;
        else
          tr.a10trialkept{ll,ttin,ss+10}=0;
          tr.nllatrialkept{ll,ttin,ss+10}=0;
          tr.alltrialkept{ll,ttin,ss+10}=0;
          tr.tll40trialkept{ll,ttin,ss+10}=0;
        end
      end
    end
  end
end
if usetr
  tr_back=tr;
  auwrong=nan(9,max(ttin),max(stageuse)+10);
end

if ~isempty(stageuse)
  % Do tac, then aud, then nul for memory reasons
  for stage=stageuse
    for ll=1:length(soalist)
      %     for tt=3:6 % different thresholds of inclusion of SOA variance around desired SOA (see concat_stim_files.m)
      % in ideal world, still keep/test various thresholds for inclusion of asynchrony (indexed by tt) but realistically only keep/use one
      %     for tt=4:5 % 4 means +/- 10ms (thus full 20ms variation); 5 means +/- 5ms (thus full 10ms variation)
      tt=ttin+2;
      cfg=[];
      cfg.vartrllength=2;
      cfg.keeptrials='yes';
      cfg.trials=[data_tac_ref.trialinfo(:,tt)==soalist(ll) & data_tac_ref.trialinfo(:,2)==stage & ~isnan(data_tac_ref.trialinfo(:,10))];
      cfg.trials=cfg.trials & data_tac_ref.trialinfo(:,13)==0; % exclude EOG.
      cfgtrialstacorig=cfg.trials;
      if 1
        length(find(cfg.trials))-length(find(cfg.trials & data_tac_ref.trialinfo(:,23)==0))
        cfg.trials=cfg.trials & data_tac_ref.trialinfo(:,23)==0; % exclude EOG.
      end
      origtrialindex=find(cfgtrialstacorig);
      if ~isempty(setdiff(find(cfgtrialstacorig),find(cfg.trials)))
        tacllorig_excluded{soalist(ll),tt-2,stage+10}=dsearchn(origtrialindex,setdiff(find(cfgtrialstacorig),find(cfg.trials)));
        origindex=setdiff(1:length(find(cfgtrialstacorig)),tacllorig_excluded{soalist(ll),tt-2,stage+10});
        if tacaud==0
          error('fix tr.tm here also')
        end
        tr.tlltrialkept_new{soalist(ll),tt-2,stage+10}=dsearchn(origindex',tr.tlltrialkept{soalist(ll),tt-2,stage+10}')';
      else
        tacllorig_excluded{soalist(ll),tt-2,stage+10}=[];
      end
      if trialkc>=0
        if trialkc==0
          cfg.trials=cfg.trials & sum(data_tac_ref.trialinfo(:,[15 14 16 25 24 26]),2)<.1; % exclude delta/Kc/SW trials
        elseif trialkc==1 % include evoked Kc trials
          cfg.trials=cfg.trials & sum(data_tac_ref.trialinfo(:,[25 24 26]),2)>sum(data_tac_ref.trialinfo(:,[15 14 16]),2) & sum(data_tac_ref.trialinfo(:,[25 24 26]),2)>.25; % trials where more Kc in poststim than prestim and in at least more than 25% of poststim time period
        else
          error('wrong value for trialkc')
        end
      end
      if ~isempty(find(cfg.trials))
        tlock_tac{soalist(ll),tt-2,stage+10}=ft_timelockanalysis(cfg,data_tac_ref);
        if strfind(lastwarn,'could not determine dimord') % keep information for features if present
          tlock_tac{soalist(ll),tt-2,stage+10}.delta={data_tac_ref.delta{find(cfg.trials)}};
          tlock_tac{soalist(ll),tt-2,stage+10}.kc={data_tac_ref.kc{find(cfg.trials)}};
          tlock_tac{soalist(ll),tt-2,stage+10}.sw={data_tac_ref.sw{find(cfg.trials)}};
          tlock_tac{soalist(ll),tt-2,stage+10}.sp_fast={data_tac_ref.sp_fast{find(cfg.trials)}};
          tlock_tac{soalist(ll),tt-2,stage+10}.sp_slow={data_tac_ref.sp_slow{find(cfg.trials)}};
        end
        lastwarn('')
      end
      %   cfg.trials=data_tac_ref.trialinfo(:,3)==soalist(ll);
      %   tlock_tac_sall{soalist(ll)}=ft_timelockanalysis(cfg,data_tac_ref); % all trials, irrepective if fell asleep
      
      %     end
    end
    cfg=[]; % tac alone 38 (-2)
    cfg.vartrllength=2;
    cfg.keeptrials='yes';
    cfg.trials=[data_tac_ref.trialinfo(:,3)==-2 & data_tac_ref.trialinfo(:,2)==stage];
    cfg.trials=cfg.trials & data_tac_ref.trialinfo(:,13)==0; % exclude EOG.
    cfgtrialstac10orig=cfg.trials;
    if 1
      length(find(cfg.trials))-length(find(cfg.trials & data_tac_ref.trialinfo(:,23)==0))
      cfg.trials=cfg.trials & data_tac_ref.trialinfo(:,23)==0; % exclude EOG.
    end
    origtrialindex=find(cfgtrialstac10orig);
    if ~isempty(setdiff(find(cfgtrialstac10orig),find(cfg.trials)))
      tac10orig_excluded{tt-2,stage+10}=dsearchn(origtrialindex,setdiff(find(cfgtrialstac10orig),find(cfg.trials)));
    else
      tac10orig_excluded{tt-2,stage+10}=[];
    end
    origindex=setdiff(1:length(find(cfgtrialstac10orig)),tac10orig_excluded{tt-2,stage+10});
    for ll=soalist
      if ~isempty(origindex) && ~isempty(tr.t10trialkept{ll,tt-2,stage+10})
        tr.t10trialkept{ll,tt-2,stage+10}=dsearchn(origindex',tr.t10trialkept{ll,tt-2,stage+10}')';
      end
    end
    
    if trialkc>=0
      if trialkc==0
        cfg.trials=cfg.trials & sum(data_tac_ref.trialinfo(:,[15 14 16 25 24 26]),2)<.1; % exclude delta/Kc/SW trials
      elseif trialkc==1 % include evoked Kc trials
        cfg.trials=cfg.trials & sum(data_tac_ref.trialinfo(:,[25 24 26]),2)>sum(data_tac_ref.trialinfo(:,[15 14 16]),2) & sum(data_tac_ref.trialinfo(:,[25 24 26]),2)>.25; % trials where more Kc in poststim than prestim and in at least more than 25% of poststim time period
      else
        error('wrong value for trialkc')
      end
    end
    if ~isempty(find(cfg.trials))
      tlock_tac{10,tt-2,stage+10}=ft_timelockanalysis(cfg,data_tac_ref); %technically this doesn't depend on 'tt' but to make it consistent with other lines of code and this new way of calling 'ttin' as input
      if strfind(lastwarn,'could not determine dimord') % keep information for features if present
        tlock_tac{10,tt-2,stage+10}.delta={data_tac_ref.delta{find(cfg.trials)}};
        tlock_tac{10,tt-2,stage+10}.kc={data_tac_ref.kc{find(cfg.trials)}};
        tlock_tac{10,tt-2,stage+10}.sw={data_tac_ref.sw{find(cfg.trials)}};
        tlock_tac{10,tt-2,stage+10}.sp_fast={data_tac_ref.sp_fast{find(cfg.trials)}};
        tlock_tac{10,tt-2,stage+10}.sp_slow={data_tac_ref.sp_slow{find(cfg.trials)}};
      end
      lastwarn('')
    end
    % average over all Tac, irrespective of SOA (but for each sleep stage)
    cfg=[]; % tac alone 38 (-2)
    cfg.vartrllength=2;
    %         cfg.keeptrials='yes';
    cfg.trials=[data_tac_ref.trialinfo(:,2)==stage];
    cfg.trials=cfg.trials & data_tac_ref.trialinfo(:,13)==0; % exclude EOG.
    if 1
      length(find(cfg.trials))-length(find(cfg.trials & data_tac_ref.trialinfo(:,23)==0))
      cfg.trials=cfg.trials & data_tac_ref.trialinfo(:,23)==0; % exclude EOG.
    end
    % FIXME auwrong: do as tac10orig_excluded above
    if trialkc>=0
      if trialkc==0
        cfg.trials=cfg.trials & sum(data_tac_ref.trialinfo(:,[15 14 16 25 24 26]),2)<.1; % exclude delta/Kc/SW trials
      elseif trialkc==1 % include evoked Kc trials
        cfg.trials=cfg.trials & sum(data_tac_ref.trialinfo(:,[25 24 26]),2)>sum(data_tac_ref.trialinfo(:,[15 14 16]),2) & sum(data_tac_ref.trialinfo(:,[25 24 26]),2)>.25; % trials where more Kc in poststim than prestim and in at least more than 25% of poststim time period
      else
        error('wrong value for trialkc')
      end
    end
    if ~isempty(find(cfg.trials))
      tlock_tac{100,tt-2,stage+10}=ft_timelockanalysis(cfg,data_tac_ref); %technically this doesn't depend on 'tt' but to make it consistent with other lines of code and this new way of calling 'ttin' as input
      if strfind(lastwarn,'could not determine dimord') % keep information for features if present
        tlock_tac{100,tt-2,stage+10}.delta={data_tac_ref.delta{find(cfg.trials)}};
        tlock_tac{100,tt-2,stage+10}.kc={data_tac_ref.kc{find(cfg.trials)}};
        tlock_tac{100,tt-2,stage+10}.sw={data_tac_ref.sw{find(cfg.trials)}};
        tlock_tac{100,tt-2,stage+10}.sp_fast={data_tac_ref.sp_fast{find(cfg.trials)}};
        tlock_tac{100,tt-2,stage+10}.sp_slow={data_tac_ref.sp_slow{find(cfg.trials)}};
      end
      lastwarn('')
    end
    % average over TacAlone and SOA+/-500ms (ignoring more close SOAs to avoid auditory)
    cfg=[]; % tac alone 38 (-2)
    cfg.vartrllength=2;
    %         cfg.keeptrials='yes';
    cfg.trials=[[data_tac_ref.trialinfo(:,3)==-2 | data_tac_ref.trialinfo(:,3)==1 | data_tac_ref.trialinfo(:,3)==9] & data_tac_ref.trialinfo(:,2)==stage];
    cfg.trials=cfg.trials & data_tac_ref.trialinfo(:,13)==0; % exclude EOG.
    if 1
      length(find(cfg.trials))-length(find(cfg.trials & data_tac_ref.trialinfo(:,23)==0))
      cfg.trials=cfg.trials & data_tac_ref.trialinfo(:,23)==0; % exclude EOG.
    end
    % FIXME auwrong
    if trialkc>=0
      if trialkc==0
        cfg.trials=cfg.trials & sum(data_tac_ref.trialinfo(:,[15 14 16 25 24 26]),2)<.1; % exclude delta/Kc/SW trials
      elseif trialkc==1 % include evoked Kc trials
        cfg.trials=cfg.trials & sum(data_tac_ref.trialinfo(:,[25 24 26]),2)>sum(data_tac_ref.trialinfo(:,[15 14 16]),2) & sum(data_tac_ref.trialinfo(:,[25 24 26]),2)>.25; % trials where more Kc in poststim than prestim and in at least more than 25% of poststim time period
      else
        error('wrong value for trialkc')
      end
    end
    if ~isempty(find(cfg.trials))
      tlock_tac{101,tt-2,stage+10}=ft_timelockanalysis(cfg,data_tac_ref); %technically this doesn't depend on 'tt' but to make it consistent with other lines of code and this new way of calling 'ttin' as input
      if strfind(lastwarn,'could not determine dimord') % keep information for features if present
        tlock_tac{101,tt-2,stage+10}.delta={data_tac_ref.delta{find(cfg.trials)}};
        tlock_tac{101,tt-2,stage+10}.kc={data_tac_ref.kc{find(cfg.trials)}};
        tlock_tac{101,tt-2,stage+10}.sw={data_tac_ref.sw{find(cfg.trials)}};
        tlock_tac{101,tt-2,stage+10}.sp_fast={data_tac_ref.sp_fast{find(cfg.trials)}};
        tlock_tac{101,tt-2,stage+10}.sp_slow={data_tac_ref.sp_slow{find(cfg.trials)}};
      end
      lastwarn('')
    end
    % average over TacAlone and T-lead-A by 500ms (to avoid auditory
    % contamination during fro 500ms but doubles available tactile trials
    cfg=[]; % tac alone and T-lead-A-500
    cfg.vartrllength=2;
    cfg.keeptrials='yes';
    cfg.trials=[[data_tac_ref.trialinfo(:,3)==-2 | data_tac_ref.trialinfo(:,3)==9] & data_tac_ref.trialinfo(:,2)==stage];
    cfg.trials=cfg.trials & data_tac_ref.trialinfo(:,13)==0; % exclude EOG.
    if 1
      length(find(cfg.trials))-length(find(cfg.trials & data_tac_ref.trialinfo(:,23)==0))
      cfg.trials=cfg.trials & data_tac_ref.trialinfo(:,23)==0; % exclude EOG.
    end
    % FIXME auwrong
    if trialkc>=0
      if trialkc==0
        cfg.trials=cfg.trials & sum(data_tac_ref.trialinfo(:,[15 14 16 25 24 26]),2)<.1; % exclude delta/Kc/SW trials
      elseif trialkc==1 % include evoked Kc trials
        cfg.trials=cfg.trials & sum(data_tac_ref.trialinfo(:,[25 24 26]),2)>sum(data_tac_ref.trialinfo(:,[15 14 16]),2) & sum(data_tac_ref.trialinfo(:,[25 24 26]),2)>.25; % trials where more Kc in poststim than prestim and in at least more than 25% of poststim time period
      else
        error('wrong value for trialkc')
      end
    end
    if ~isempty(find(cfg.trials))
      tlock_tac{102,tt-2,stage+10}=ft_timelockanalysis(cfg,data_tac_ref); %technically this doesn't depend on 'tt' but to make it consistent with other lines of code and this new way of calling 'ttin' as input
      if strfind(lastwarn,'could not determine dimord') % keep information for features if present
        tlock_tac{102,tt-2,stage+10}.delta={data_tac_ref.delta{find(cfg.trials)}};
        tlock_tac{102,tt-2,stage+10}.kc={data_tac_ref.kc{find(cfg.trials)}};
        tlock_tac{102,tt-2,stage+10}.sw={data_tac_ref.sw{find(cfg.trials)}};
        tlock_tac{102,tt-2,stage+10}.sp_fast={data_tac_ref.sp_fast{find(cfg.trials)}};
        tlock_tac{102,tt-2,stage+10}.sp_slow={data_tac_ref.sp_slow{find(cfg.trials)}};
      end
      lastwarn('')
    end
  end
  taceog_exclude=unique([find(data_tac_ref.trialinfo(:,13)~=0); find(data_tac_ref.trialinfo(:,23)~=0)]);
  clear data_tac_ref
  
  % Do tac, then aud, then nul for memory reasons
  for stage=stageuse
    for ll=1:length(soalist)
      %     for tt=3:6 % different thresholds of inclusion of SOA variance around desired SOA (see concat_stim_files.m)
      tt=ttin+2;
      cfg=[];
      cfg.vartrllength=2;
      cfg.keeptrials='yes';
      cfg.trials=[data_aud_ref.trialinfo(:,tt)==soalist(ll) & data_aud_ref.trialinfo(:,2)==stage & ~isnan(data_aud_ref.trialinfo(:,10))];
      cfg.trials=cfg.trials & data_aud_ref.trialinfo(:,13)==0; % exclude EOG.
      if 1
        length(find(cfg.trials))-length(find(cfg.trials & data_aud_ref.trialinfo(:,23)==0))
        cfg.trials=cfg.trials & data_aud_ref.trialinfo(:,23)==0; % exclude EOG.
      end
      if trialkc>=0
        if trialkc==0
          cfg.trials=cfg.trials & sum(data_aud_ref.trialinfo(:,[15 14 16 25 24 26]),2)<.1; % exclude delta/Kc/SW trials
        elseif trialkc==1 % include evoked Kc trials
          cfg.trials=cfg.trials & sum(data_aud_ref.trialinfo(:,[25 24 26]),2)>sum(data_aud_ref.trialinfo(:,[15 14 16]),2) & sum(data_aud_ref.trialinfo(:,[25 24 26]),2)>.25; % trials where more Kc in poststim than prestim and in at least more than 25% of poststim time period
        else
          error('wrong value for trialkc')
        end
      end
      if ~isempty(find(cfg.trials))
        tlock_aud{soalist(ll),tt-2,stage+10}=ft_timelockanalysis(cfg,data_aud_ref);
        if strfind(lastwarn,'could not determine dimord') % keep information for features if present
          tlock_aud{soalist(ll),tt-2,stage+10}.delta={data_aud_ref.delta{find(cfg.trials)}};
          tlock_aud{soalist(ll),tt-2,stage+10}.kc={data_aud_ref.kc{find(cfg.trials)}};
          tlock_aud{soalist(ll),tt-2,stage+10}.sw={data_aud_ref.sw{find(cfg.trials)}};
          tlock_aud{soalist(ll),tt-2,stage+10}.sp_fast={data_aud_ref.sp_fast{find(cfg.trials)}};
          tlock_aud{soalist(ll),tt-2,stage+10}.sp_slow={data_aud_ref.sp_slow{find(cfg.trials)}};
        end
        lastwarn('')
      end
      %   cfg.trials=data_aud_ref.trialinfo(:,3)==soalist(ll);
      %   tlock_aud_sall{soalist(ll)}=ft_timelockanalysis(cfg,data_aud_ref);
      %     end
    end % ll
    cfg=[]; % aud alone 39 (-1)
    cfg.vartrllength=2;
    cfg.keeptrials='yes';
    cfg.trials=[data_aud_ref.trialinfo(:,3)==-1 & data_aud_ref.trialinfo(:,2)==stage];
    cfg.trials=cfg.trials & data_aud_ref.trialinfo(:,13)==0; % exclude EOG.
    trialstmp=cfg.trials;
    if 1
      length(find(cfg.trials))-length(find(cfg.trials & data_aud_ref.trialinfo(:,23)==0))
      cfg.trials=cfg.trials & data_aud_ref.trialinfo(:,23)==0; % exclude EOG.
    end
    %     origcfgtrials=cfg.trials;
    if usetr
      origtrialindex=find(trialstmp);
      for ll=soalist,
        try
          if isempty(max(tr.autrialkept{ll,tt-2,stage+10})>length(find(cfg.trials)))
            auwrong(ll,tt-2,stage+10)=0;
          else
            auwrong(ll,tt-2,stage+10)=max(tr.autrialkept{ll,tt-2,stage+10})>length(find(cfg.trials));
          end;
        catch
          auwrong(ll,tt-2,stage+10)=0;
        end
      end
      if any(auwrong(soalist,tt-2,stage+10))  % horrible hack to make up for not saving out tr.autrialkept correctly awhile ago.
        poststimeogexcluded=dsearchn(origtrialindex,setdiff(origtrialindex,find(cfg.trials)));
        origindex=setdiff(1:length(origtrialindex),poststimeogexcluded);
        for ll=soalist,  % one of the ll needs to be 10 instead?
          %           cfg.trials=origcfgtrials;
%           prevtotnumused=length(tr_back.t10trialkept{ll,3,10});
          traukeep_oldind=setdiff(tr.autrialkept{ll,tt-2,stage+10},poststimeogexcluded);
          all40excluded{ll,tt-2,stage+10}=dsearchn(tr.autrialkept{ll,tt-2,stage+10}',setdiff(tr.autrialkept{ll,tt-2,stage+10},traukeep_oldind)');
          tr.autrialkept{ll,tt-2,stage+10}=dsearchn(origindex',traukeep_oldind')';
        end
        
        
        %         for ll=soalist,  % one of the ll needs to be 10 instead?
        %           cfg.trials(origtrialindex(tr.autrialkept{ll,tt-2,stage+10}(find(diff(dsearchn(origindex',tr.autrialkept{ll,tt-2,stage+10}')')<1)+1)))=1;
        %         end
        %         origindex2=setdiff(1:length(origtrialindex),dsearchn(origtrialindex,setdiff(origtrialindex,find(cfg.trials))));
        %         for ll=soalist
        %
        %           tr.autrialkept{ll,tt-2,stage+10}=dsearchn(origindex2',tr.autrialkept{ll,tt-2,stage+10}');
        %         end
      end
      %       origtrialindex=find(trialstmp);
      %       for ll=soalist,auwrong(ll)=max(tr.autrialkept{ll,tt-2,stage+10})>length(find(cfg.trials));end
      %       if any(auwrong(soalist))  % horrible hack to make up for not saving out tr.autrialkept correctly awhile ago.
      %         origindex=setdiff(1:length(find(trialstmp)),dsearchn(find(trialstmp),setdiff(find(trialstmp),find(cfg.trials))));
      %         cfg.trials(origtrialindex(tr.autrialkept{ll,tt-2,stage+10}(find(diff(dsearchn(origindex',tr.autrialkept{ll,tt-2,stage+10}')')<1)+1)))=1;
      %         origindex=setdiff(1:length(find(trialstmp)),dsearchn(find(trialstmp),setdiff(find(trialstmp),find(cfg.trials))));
      %         for ll=soalist,tr.autrialkept{ll,tt-2,stage+10}=dsearchn(origindex',tr.autrialkept{ll,tt-2,stage+10}');end
      %       end
    end
    if trialkc>=0
      if trialkc==0
        cfg.trials=cfg.trials & sum(data_aud_ref.trialinfo(:,[15 14 16 25 24 26]),2)<.1; % exclude delta/Kc/SW trials
      elseif trialkc==1 % include evoked Kc trials
        cfg.trials=cfg.trials & sum(data_aud_ref.trialinfo(:,[25 24 26]),2)>sum(data_aud_ref.trialinfo(:,[15 14 16]),2) & sum(data_aud_ref.trialinfo(:,[25 24 26]),2)>.25; % trials where more Kc in poststim than prestim and in at least more than 25% of poststim time period
      else
        error('wrong value for trialkc')
      end
    end
    cfgtrials_aud10=find(cfg.trials);
    if ~isempty(find(cfg.trials))
      tlock_aud{10,tt-2,stage+10}=ft_timelockanalysis(cfg,data_aud_ref);
      if strfind(lastwarn,'could not determine dimord') % keep information for features if present
        tlock_aud{10,tt-2,stage+10}.delta={data_aud_ref.delta{find(cfg.trials)}};
        tlock_aud{10,tt-2,stage+10}.kc={data_aud_ref.kc{find(cfg.trials)}};
        tlock_aud{10,tt-2,stage+10}.sw={data_aud_ref.sw{find(cfg.trials)}};
        tlock_aud{10,tt-2,stage+10}.sp_fast={data_aud_ref.sp_fast{find(cfg.trials)}};
        tlock_aud{10,tt-2,stage+10}.sp_slow={data_aud_ref.sp_slow{find(cfg.trials)}};
      end
      lastwarn('')
    end
    % average over all Aud, irrespective of SOA (but for each sleep stage)
    cfg=[]; % aud alone 39 (-1)
    cfg.vartrllength=2;
    %         cfg.keeptrials='yes';
    cfg.trials=[data_aud_ref.trialinfo(:,2)==stage];
    cfg.trials=cfg.trials & data_aud_ref.trialinfo(:,13)==0; % exclude EOG.
    if 1
      length(find(cfg.trials))-length(find(cfg.trials & data_aud_ref.trialinfo(:,23)==0))
      cfg.trials=cfg.trials & data_aud_ref.trialinfo(:,23)==0; % exclude EOG.
    end
    % FIXME: auwrong
    if trialkc>=0
      if trialkc==0
        cfg.trials=cfg.trials & sum(data_aud_ref.trialinfo(:,[15 14 16 25 24 26]),2)<.1; % exclude delta/Kc/SW trials
      elseif trialkc==1 % include evoked Kc trials
        cfg.trials=cfg.trials & sum(data_aud_ref.trialinfo(:,[25 24 26]),2)>sum(data_aud_ref.trialinfo(:,[15 14 16]),2) & sum(data_aud_ref.trialinfo(:,[25 24 26]),2)>.25; % trials where more Kc in poststim than prestim and in at least more than 25% of poststim time period
      else
        error('wrong value for trialkc')
      end
    end
    if ~isempty(find(cfg.trials))
      tlock_aud{100,tt-2,stage+10}=ft_timelockanalysis(cfg,data_aud_ref);
      if strfind(lastwarn,'could not determine dimord') % keep information for features if present
        tlock_aud{100,tt-2,stage+10}.delta={data_aud_ref.delta{find(cfg.trials)}};
        tlock_aud{100,tt-2,stage+10}.kc={data_aud_ref.kc{find(cfg.trials)}};
        tlock_aud{100,tt-2,stage+10}.sw={data_aud_ref.sw{find(cfg.trials)}};
        tlock_aud{100,tt-2,stage+10}.sp_fast={data_aud_ref.sp_fast{find(cfg.trials)}};
        tlock_aud{100,tt-2,stage+10}.sp_slow={data_aud_ref.sp_slow{find(cfg.trials)}};
      end
      lastwarn('')
    end
    % average over AudAlone and SOA+/-500ms (ignoring more close SOAs to avoid tactile)
    cfg=[]; % aud alone 39 (-1)
    cfg.vartrllength=2;
    %         cfg.keeptrials='yes';
    cfg.trials=[[data_aud_ref.trialinfo(:,3)==-1 | data_aud_ref.trialinfo(:,3)==1 | data_aud_ref.trialinfo(:,3)==9] & data_aud_ref.trialinfo(:,2)==stage];
    cfg.trials=cfg.trials & data_aud_ref.trialinfo(:,13)==0; % exclude EOG.
    if 1
      length(find(cfg.trials))-length(find(cfg.trials & data_aud_ref.trialinfo(:,23)==0))
      cfg.trials=cfg.trials & data_aud_ref.trialinfo(:,23)==0; % exclude EOG.
    end
    % FIXME: auwrong
    if trialkc>=0
      if trialkc==0
        cfg.trials=cfg.trials & sum(data_aud_ref.trialinfo(:,[15 14 16 25 24 26]),2)<.1; % exclude delta/Kc/SW trials
      elseif trialkc==1 % include evoked Kc trials
        cfg.trials=cfg.trials & sum(data_aud_ref.trialinfo(:,[25 24 26]),2)>sum(data_aud_ref.trialinfo(:,[15 14 16]),2) & sum(data_aud_ref.trialinfo(:,[25 24 26]),2)>.25; % trials where more Kc in poststim than prestim and in at least more than 25% of poststim time period
      else
        error('wrong value for trialkc')
      end
    end
    if ~isempty(find(cfg.trials))
      tlock_aud{101,tt-2,stage+10}=ft_timelockanalysis(cfg,data_aud_ref);
      if strfind(lastwarn,'could not determine dimord') % keep information for features if present
        tlock_aud{101,tt-2,stage+10}.delta={data_aud_ref.delta{find(cfg.trials)}};
        tlock_aud{101,tt-2,stage+10}.kc={data_aud_ref.kc{find(cfg.trials)}};
        tlock_aud{101,tt-2,stage+10}.sw={data_aud_ref.sw{find(cfg.trials)}};
        tlock_aud{101,tt-2,stage+10}.sp_fast={data_aud_ref.sp_fast{find(cfg.trials)}};
        tlock_aud{101,tt-2,stage+10}.sp_slow={data_aud_ref.sp_slow{find(cfg.trials)}};
      end
      lastwarn('')
    end
    % average over AudAlone and A-lead-T by 500ms
    cfg=[]; %
    cfg.vartrllength=2;
    cfg.keeptrials='yes';
    cfg.trials=[[data_aud_ref.trialinfo(:,3)==-1 | data_aud_ref.trialinfo(:,3)==1] & data_aud_ref.trialinfo(:,2)==stage];
    cfg.trials=cfg.trials & data_aud_ref.trialinfo(:,13)==0; % exclude EOG.
    if 1
      length(find(cfg.trials))-length(find(cfg.trials & data_aud_ref.trialinfo(:,23)==0))
      cfg.trials=cfg.trials & data_aud_ref.trialinfo(:,23)==0; % exclude EOG.
    end
    % FIXME: auwrong
    if trialkc>=0
      if trialkc==0
        cfg.trials=cfg.trials & sum(data_aud_ref.trialinfo(:,[15 14 16 25 24 26]),2)<.1; % exclude delta/Kc/SW trials
      elseif trialkc==1 % include evoked Kc trials
        cfg.trials=cfg.trials & sum(data_aud_ref.trialinfo(:,[25 24 26]),2)>sum(data_aud_ref.trialinfo(:,[15 14 16]),2) & sum(data_aud_ref.trialinfo(:,[25 24 26]),2)>.25; % trials where more Kc in poststim than prestim and in at least more than 25% of poststim time period
      else
        error('wrong value for trialkc')
      end
    end
    if ~isempty(find(cfg.trials))
      tlock_aud{102,tt-2,stage+10}=ft_timelockanalysis(cfg,data_aud_ref);
      if strfind(lastwarn,'could not determine dimord') % keep information for features if present
        tlock_aud{102,tt-2,stage+10}.delta={data_aud_ref.delta{find(cfg.trials)}};
        tlock_aud{102,tt-2,stage+10}.kc={data_aud_ref.kc{find(cfg.trials)}};
        tlock_aud{102,tt-2,stage+10}.sw={data_aud_ref.sw{find(cfg.trials)}};
        tlock_aud{102,tt-2,stage+10}.sp_fast={data_aud_ref.sp_fast{find(cfg.trials)}};
        tlock_aud{102,tt-2,stage+10}.sp_slow={data_aud_ref.sp_slow{find(cfg.trials)}};
      end
      lastwarn('')
    end
  end
  audeog_exclude=unique([find(data_aud_ref.trialinfo(:,13)~=0); find(data_aud_ref.trialinfo(:,23)~=0)]);
  clear data_aud_ref
  
  % Do tac, then aud, then nul for memory reasons
  for stage=stageuse
    cfg=[];
    cfg.vartrllength=2;
    cfg.keeptrials='yes';
    cfg.trials=data_nul_ref.trialinfo(:,2)==stage;
    cfg.trials=cfg.trials & data_nul_ref.trialinfo(:,13)==0; % exclude EOG.
    cfgtrialsnul10orig=cfg.trials;
    if 1
      length(find(cfg.trials))-length(find(cfg.trials & data_nul_ref.trialinfo(:,23)==0))
      cfg.trials=cfg.trials & data_nul_ref.trialinfo(:,23)==0; % exclude EOG.
    end
    origtrialindex=find(cfgtrialsnul10orig);
    if ~isempty(origtrialindex) && ~isempty(setdiff(find(cfgtrialsnul10orig),find(cfg.trials)))
      nul10orig_excluded{tt-2,stage+10}=dsearchn(origtrialindex,setdiff(find(cfgtrialsnul10orig),find(cfg.trials)));
      origindex=setdiff(1:length(find(cfgtrialsnul10orig)),nul10orig_excluded{tt-2,stage+10});
      for ll=soalist
        tr.nllttrialkept{ll,tt-2,stage+10}=dsearchn(origindex',tr.nllttrialkept{ll,tt-2,stage+10}')';
      end
    else
      nul10orig_excluded{tt-2,stage+10}=[];
    end
    
    if trialkc>=0
      if trialkc==0
        cfg.trials=cfg.trials & sum(data_nul_ref.trialinfo(:,[15 14 16 25 24 26]),2)<.1; % exclude delta/Kc/SW trials
      elseif trialkc==1 % include evoked Kc trials
        cfg.trials=cfg.trials & sum(data_nul_ref.trialinfo(:,[25 24 26]),2)>sum(data_nul_ref.trialinfo(:,[15 14 16]),2) & sum(data_nul_ref.trialinfo(:,[25 24 26]),2)>.25; % trials where more Kc in poststim than prestim and in at least more than 25% of poststim time period
      else
        error('wrong value for trialkc')
      end
    end
    if ~isempty(find(cfg.trials))
      tlock_nul{10,tt-2,stage+10}=ft_timelockanalysis(cfg,data_nul_ref);
      if strfind(lastwarn,'could not determine dimord') % keep information for features if present
        tlock_nul{10,tt-2,stage+10}.delta={data_nul_ref.delta{find(cfg.trials)}};
        tlock_nul{10,tt-2,stage+10}.kc={data_nul_ref.kc{find(cfg.trials)}};
        tlock_nul{10,tt-2,stage+10}.sw={data_nul_ref.sw{find(cfg.trials)}};
        tlock_nul{10,tt-2,stage+10}.sp_fast={data_nul_ref.sp_fast{find(cfg.trials)}};
        tlock_nul{10,tt-2,stage+10}.sp_slow={data_nul_ref.sp_slow{find(cfg.trials)}};
      end
      lastwarn('')
    end
    % cfg.trials='all';
    % tlock_nul_sall{10}=ft_timelockanalysis(cfg,data_nul_ref_ref);
  end
  nuleog_exclude=unique([find(data_nul_ref.trialinfo(:,13)~=0); find(data_nul_ref.trialinfo(:,23)~=0)]);
  clear data_nul_ref
  
  
  
  
  %%
  
  
  soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];
  ttimwin=[min(-.2,soades-.2); max(0.5,soades+.5)]';
  ttimwin(10,:)=[min(ttimwin(1:9,1)) max(ttimwin(1:9,2))];
  atimwin=[min(-.2,fliplr(soades)-.2); max(0.5,fliplr(soades)+.5)]';
  atimwin(10,:)=ttimwin(10,:);
  ntimwin(10,:)=ttimwin(10,:);
  
  
  for ss=stageuse+10
    if exist('tlock_tac','var') && exist('tlock_aud','var') && exist('tlock_nul','var') && ss<=size(tlock_tac,3) && ss<=size(tlock_aud,3) && ss<=size(tlock_nul,3) && ~isempty(tlock_tac{4,ttin,ss}) && ~isempty(tlock_aud{5,ttin,ss}) && ~isempty(tlock_nul{10,ttin,ss})
      try
        chanuse=match_str(tlock_tac{3,ttin,ss}.label,'Fz');
      catch
        chanuse=match_str(tlock_nul{10,ttin,ss}.label,'Fz');
      end
    else
      stageuse=setdiff(stageuse,ss-10);
      tr.stageuse=stageuse;
      continue
    end
    
    for ll=[soalist 10] % this loop is to remove trials with NaN in (e.g. if artifact rejection cut through middle of a trial)
      %     for tt=1:4
      tt=ttin;
      
      if ~isempty(tlock_tac{ll,tt,ss})
        cfg=[];
        %       cfg.keeptrials='yes';
        %         try
        cfg.trials=find(~any(isnan((tlock_tac{ll,tt,ss}.trial(:,chanuse,dsearchn(tlock_tac{ll,tt,ss}.time',ttimwin(ll,1)):dsearchn(tlock_tac{ll,tt,ss}.time',ttimwin(ll,2)) ))) ,3));
        %       tlock_tac{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
        if size(tlock_tac{ll,tt,ss}.trial,1)==1 && length(cfg.trials)<2
        else
          tmp=tlock_tac{ll,tt,ss};
          tlock_tac{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{ll,tt,ss});
          if strfind(lastwarn,'could not determine dimord') % keep information for features if present
            tlock_tac{ll,tt,ss}.delta={tmp.delta{(cfg.trials)}};
            tlock_tac{ll,tt,ss}.kc={tmp.kc{(cfg.trials)}};
            tlock_tac{ll,tt,ss}.sw={tmp.sw{(cfg.trials)}};
            tlock_tac{ll,tt,ss}.sp_fast={tmp.sp_fast{(cfg.trials)}};
            tlock_tac{ll,tt,ss}.sp_slow={tmp.sp_slow{(cfg.trials)}};
          end
          lastwarn('')
        end
        tmp=tlock_tac{ll,tt,ss};
        tlock_tac{ll,tt,ss}=trim_nans(tlock_tac{ll,tt,ss});
        if strfind(lastwarn,'could not determine dimord') % keep information for features if present
          tlock_tac{ll,tt,ss}.delta=tmp.delta;
          tlock_tac{ll,tt,ss}.kc=tmp.kc;
          tlock_tac{ll,tt,ss}.sw=tmp.sw;
          tlock_tac{ll,tt,ss}.sp_fast={tmp.sp_fast{(cfg.trials)}};
          tlock_tac{ll,tt,ss}.sp_slow={tmp.sp_slow{(cfg.trials)}};
        end
        lastwarn('')
      else
      end
      
      %         catch
      %           disp(['no tlock_tac for ' num2str(ll) num2str(tt) num2str(ss)])
      %         end
      
      
      if ~isempty(tlock_aud{ll,tt,ss})
        cfg=[];
        %       cfg.keeptrials='yes';
        %       cfg.trials=find(sum(squeeze(isnan(tlock_aud{ll,tt,ss}.trial(:,10,:)))')<.05*size(tlock_aud{ll,tt,ss}.trial,3));
        %         try
        cfg.trials=find(~any(isnan((tlock_aud{ll,tt,ss}.trial(:,chanuse,dsearchn(tlock_aud{ll,tt,ss}.time',atimwin(ll,1)):dsearchn(tlock_aud{ll,tt,ss}.time',atimwin(ll,2)) ))),3));
        %       tlock_aud{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_aud{ll,tt,ss});
        if size(tlock_aud{ll,tt,ss}.trial,1)==1 && length(cfg.trials)<2
        else
          tmp=tlock_aud{ll,tt,ss};
          tlock_aud{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{ll,tt,ss});
          if strfind(lastwarn,'could not determine dimord') % keep information for features if present
            tlock_aud{ll,tt,ss}.delta={tmp.delta{(cfg.trials)}};
            tlock_aud{ll,tt,ss}.kc={tmp.kc{(cfg.trials)}};
            tlock_aud{ll,tt,ss}.sw={tmp.sw{(cfg.trials)}};
            tlock_aud{ll,tt,ss}.sp_fast={tmp.sp_fast{(cfg.trials)}};
            tlock_aud{ll,tt,ss}.sp_slow={tmp.sp_slow{(cfg.trials)}};
          end
          lastwarn('')
        end
        %       tmp=squeeze(mean(tlock_aud{ll,tt,ss}.trial,1));
        %       if sum(isnan(tmp(1,:)))>.07*length(tmp(1,:))
        %         error('this averaging didnt work');
        %       end
        tmp=tlock_aud{ll,tt,ss};
        tlock_aud{ll,tt,ss}=trim_nans(tlock_aud{ll,tt,ss});
        if strfind(lastwarn,'could not determine dimord') % keep information for features if present
          tlock_aud{ll,tt,ss}.delta=tmp.delta;
          tlock_aud{ll,tt,ss}.kc=tmp.kc;
          tlock_aud{ll,tt,ss}.sw=tmp.sw;
          tlock_aud{ll,tt,ss}.sp_fast={tmp.sp_fast{(cfg.trials)}};
          tlock_aud{ll,tt,ss}.sp_slow={tmp.sp_slow{(cfg.trials)}};
        end
        lastwarn('')
      else
      end
      %         catch
      %           disp(['no tlock_aud for ' num2str(ll) num2str(tt) num2str(ss)])
      %         end
      %       end
      
      if ll==10
        %       if ll==10 && tt==1
        cfg=[];
        %       cfg.keeptrials='yes';
        %       cfg.trials=find(sum(squeeze(isnan(tlock_nul{ll}.trial(:,10,:)))')<.05*size(tlock_nul{ll}.trial,3));
        %         try
        cfg.trials=find(~any(isnan((tlock_nul{ll,tt,ss}.trial(:,chanuse,dsearchn(tlock_nul{ll,tt,ss}.time',ntimwin(ll,1)):dsearchn(tlock_nul{ll,tt,ss}.time',ntimwin(ll,2)) ))),3));
        %       tlock_nul{ll}=ft_timelockanalysis(cfg,tlock_nul{ll});
        if size(tlock_nul{ll,tt,ss}.trial,1)==1 && length(cfg.trials)<2
        else
          tmp=tlock_nul{ll,tt,ss};
          tlock_nul{ll,tt,ss}=ft_selectdata(cfg,tlock_nul{ll,tt,ss});
          if strfind(lastwarn,'could not determine dimord') % keep information for features if present
            tlock_nul{ll,tt,ss}.delta={tmp.delta{(cfg.trials)}};
            tlock_nul{ll,tt,ss}.kc={tmp.kc{(cfg.trials)}};
            tlock_nul{ll,tt,ss}.sw={tmp.sw{(cfg.trials)}};
            tlock_nul{ll,tt,ss}.sp_fast={tmp.sp_fast{(cfg.trials)}};
            tlock_nul{ll,tt,ss}.sp_slow={tmp.sp_slow{(cfg.trials)}};
          end
          lastwarn('')
        end
        %       tmp=squeeze(mean(tlock_nul{ll,tt,ss}.trial,1));
        %       if sum(isnan(tmp(1,:)))>.07*length(tmp(1,:))
        %         error('this averaging didnt work');
        %       end
        tmp=tlock_nul{ll,tt,ss};
        tlock_nul{ll,tt,ss}=trim_nans(tlock_nul{ll,tt,ss});
        if strfind(lastwarn,'could not determine dimord') % keep information for features if present
          tlock_nul{ll,tt,ss}.delta=tmp.delta;
          tlock_nul{ll,tt,ss}.kc=tmp.kc;
          tlock_nul{ll,tt,ss}.sw=tmp.sw;
          tlock_nul{ll,tt,ss}.sp_fast={tmp.sp_fast{(cfg.trials)}};
          tlock_nul{ll,tt,ss}.sp_slow={tmp.sp_slow{(cfg.trials)}};
        end
        lastwarn('')
        %         catch
        %           disp(['no tlock_nul for ' num2str(ll) num2str(tt) num2str(ss)])
        %         end
        
      end
      if ll<10
        % make sure the multisensory auditory-locked and ms tactile-locked have same trials; sometimes not exactly!
        %         try
        if ~isempty(tlock_tac{ll,tt,ss}) && ~isempty(tlock_aud{ll,tt,ss})
          trcom=intersect(tlock_tac{ll,tt,ss}.trialinfo(:,11),tlock_aud{ll,tt,ss}.trialinfo(:,11));
          [~,tbb,~]=intersect(tlock_tac{ll,tt,ss}.trialinfo(:,11),trcom);
          [~,abb,~]=intersect(tlock_aud{ll,tt,ss}.trialinfo(:,11),trcom);
          cfg=[];
          cfg.trials=tbb;
          if size(tlock_tac{ll,tt,ss}.trial,1)==1 && length(cfg.trials)<2
          else
            tmp=tlock_tac{ll,tt,ss};
            tlock_tac{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{ll,tt,ss});
            if strfind(lastwarn,'could not determine dimord') % keep information for features if present
              tlock_tac{ll,tt,ss}.delta={tmp.delta{(cfg.trials)}};
              tlock_tac{ll,tt,ss}.kc={tmp.kc{(cfg.trials)}};
              tlock_tac{ll,tt,ss}.sw={tmp.sw{(cfg.trials)}};
              tlock_tac{ll,tt,ss}.sp_fast={tmp.sp_fast{(cfg.trials)}};
              tlock_tac{ll,tt,ss}.sp_slow={tmp.sp_slow{(cfg.trials)}};
            end
            lastwarn('')
          end
          cfg=[];
          cfg.trials=abb;
          if size(tlock_aud{ll,tt,ss}.trial,1)==1 && length(cfg.trials)<2
          else
            tmp=tlock_aud{ll,tt,ss};
            tlock_aud{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{ll,tt,ss});
            if strfind(lastwarn,'could not determine dimord') % keep information for features if present
              tlock_aud{ll,tt,ss}.delta={tmp.delta{(cfg.trials)}};
              tlock_aud{ll,tt,ss}.kc={tmp.kc{(cfg.trials)}};
              tlock_aud{ll,tt,ss}.sw={tmp.sw{(cfg.trials)}};
              tlock_aud{ll,tt,ss}.sp_fast={tmp.sp_fast{(cfg.trials)}};
              tlock_aud{ll,tt,ss}.sp_slow={tmp.sp_slow{(cfg.trials)}};
            end
            lastwarn('')
          end
        end
        %         catch
        %           disp(['no tlock_tack or tlock_aud for ' num2str(ll) num2str(tt) num2str(ss)])
        %         end
      end
      %     end
    end
  end
  
  if 0  %why had we been loading this?
    load([ddir sub{ii} '_audtac.mat']);
  end
  
  tnumtr=zeros(max(soalist),ttin,ss); % use ss as max(stageuse)
  anumtr=zeros(max(soalist),ttin,ss);
  
  % % fix not checking ss if not exist
  for ss=stageuse+10
    try
      fsample=1/diff(tlock_tac{10,tt,ss}.time(1:2)); % sampling rate
    catch
      try
        fsample=1/diff(tlock_aud{10,tt,ss}.time(1:2)); % sampling rate
      catch
        fsample=1/diff(tlock_nul{10,tt,ss}.time(1:2)); % sampling rate
      end
    end
    
    for ll=soalist
      %     for tt=1:4 % limits of plus/minus:  [no-limit .01 .005 .002]
      tt=ttin;
      try
        tnumtr(ll,tt,ss)=min([size(tlock_tac{ll,tt,ss}.trial,1) size(tlock_tac{10,tt,ss}.trial,1)]);
      catch
        tnumtr(ll,tt,ss)=0;
      end
      try
        anumtr(ll,tt,ss)=min([size(tlock_aud{ll,tt,ss}.trial,1) size(tlock_aud{10,tt,ss}.trial,1)]);
      catch
        anumtr(ll,tt,ss)=0;
      end
      %     %     end
    end
    
    % this is for Uta's idea of comparing (TA70 + AT70_shifted70ms) vs (TA0 + TA0_shifted70ms)
    %     try
    %       tnumtrMSshift(1,tt,ss)=min([size(tlock_tac{1,tt,ss}.trial,1) size(tlock_tac{5,tt,ss}.trial,1) size(tlock_tac{9,tt,ss}.trial,1)]);
    %     catch
    %       tnumtrMSshift(1,tt,ss)=0;
    %     end
    %     try
    %       tnumtrMSshift(3,tt,ss)=min([size(tlock_tac{3,tt,ss}.trial,1) size(tlock_tac{5,tt,ss}.trial,1) size(tlock_tac{7,tt,ss}.trial,1)]);
    %     catch
    %       tnumtrMSshift(3,tt,ss)=0;
    %     end
    %     try
    %       tnumtrMSshift(4,tt,ss)=min([size(tlock_tac{4,tt,ss}.trial,1) size(tlock_tac{5,tt,ss}.trial,1) size(tlock_tac{6,tt,ss}.trial,1)]);
    %     catch
    %       tnumtrMSshift(4,tt,ss)=0;
    %     end
    
    % MUst do shifting first before freqanalysis.
    
    % shifting unisensory to match for later sum of unisensory to match a multisensory condition
    for ll=soalist
      %     for tt=1:4
      tt=ttin;
      
      %       if ~isnan(tnumtr(ll,ss,tt))
      if ~tnumtr(ll,tt,ss) || ~anumtr(ll,tt,ss)
        tlock_aud{ll+40,tt,ss}=[];
        tlock_tac{ll+40,tt,ss}=[];
        continue
      else
        if ~usetr
          tmrandtr=Shuffle(1:size(tlock_tac{ll,tt,ss}.trial,1));
          turandtr=Shuffle(1:size(tlock_tac{10,tt,ss}.trial,1));
          
          tr.tmtrialkept{ll,tt,ss}=sort(tmrandtr(1:tnumtr(ll,tt,ss)));
          tr.tutrialkept{ll,tt,ss}=sort(turandtr(1:tnumtr(ll,tt,ss)));
        end
        
        %         if ll<5
        %           trandtrMSshift1=Shuffle(1:size(tlock_tac{ll,tt,ss}.trial,1));
        %           trandtrMSshift2=Shuffle(1:size(tlock_tac{10-ll,tt,ss}.trial,1));
        %           trandtrMSshift3=Shuffle(1:size(tlock_tac{5,tt,ss}.trial,1));
        %           trandtrMSshift4=Shuffle(1:size(tlock_tac{5,tt,ss}.trial,1)); % do this twice to get 2 different random subsamples.
        %
        %           ttrialkeptMSshift1{ll,tt,ss}=sort(trandtrMSshift1(1:tnumtrMSshift(ll,tt,ss)));
        %           ttrialkeptMSshift2{ll,tt,ss}=sort(trandtrMSshift2(1:tnumtrMSshift(ll,tt,ss)));
        %           ttrialkeptMSshift3{ll,tt,ss}=sort(trandtrMSshift3(1:tnumtrMSshift(ll,tt,ss)));
        %           ttrialkeptMSshift4{ll,tt,ss}=sort(trandtrMSshift4(1:tnumtrMSshift(ll,tt,ss)));
        %         end
        
      end
      
      if tacaud==0
        error('fix me with tutrialkept')
        cfg=[];
        % I have thought it through, and no minus sign needed here.
        cfg.offset=round(fsample*(tlock_tac{ll,tt,ss}.trialinfo(tr.tmtrialkept{ll,tt,ss},10)-soades(ll) -soades(ll))); % offset is in samples not seconds; -2*soades required for tactile not auditory
        cfg.trials=tr.tutrialkept{ll,tt,ss};
        if length(cfg.trials)==1
          tlock_tac{ll+40,tt,ss}=[];
        else
          tmp=tlock_tac{10,tt,ss};
          warning off
          tlock_tac{ll+40,tt,ss}=ft_redefinetrial(cfg,tlock_tac{10,tt,ss}); % Talone with shift in time matching SOA jitter
          warning on
          if isfield(tmp,'delta')
            tlock_tac{ll+40,tt,ss}.delta={tmp.delta{(cfg.trials)}};
            tlock_tac{ll+40,tt,ss}.kc={tmp.kc{(cfg.trials)}};
            tlock_tac{ll+40,tt,ss}.sw={tmp.sw{(cfg.trials)}};
            tlock_tac{ll+40,tt,ss}.sp_fast={tmp.sp_fast{(cfg.trials)}};
            tlock_tac{ll+40,tt,ss}.sp_slow={tmp.sp_slow{(cfg.trials)}};
          end
          lastwarn('')
          if size(tlock_tac{ll+40,tt,ss}.trialinfo,2)>240 % as it will if doing 'features'
            % need to also shift the time of the saved out negmax of Kc
            keyboard % shift anything for 'struct' as well?
            for nn=[50:59 70:79 90:99 110:119 130:139 150:159]
              tmpchg=find(tlock_tac{ll+40,tt,ss}.trialinfo(:,nn));
              tlock_tac{ll+40,tt,ss}.trialinfo(tmpchg,nn+10) = cfg.offset(tmpchg)/1000 + tlock_tac{ll+40,tt,ss}.trialinfo(tmpchg,nn+10);
            end
          end
        end
      end
      
      
      % Note to self: even though this will unnecesarily cut a few trials
      % from multisensory condition (due to possible fewer matching
      % unisensory) it's easier for RAM reasons to just go with it and use
      % this subset of tlock_tac{ll,tt,ss} for the Multisensory-Shifted
      % comparison.
      
      cfg=[];
      cfg.trials=tr.tmtrialkept{ll,tt,ss}; % only keep ms trials that have a paired shifted corresponding unisensory trial
      if length(cfg.trials)==1
        tlock_tac{ll,tt,ss}=[];
      else
        tmp=tlock_tac{ll,tt,ss};
        %         tlock_tac{ll,tt,ss}=ft_redefinetrial(cfg,tlock_tac{ll,tt,ss});
        tlock_tac{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{ll,tt,ss});
        if strfind(lastwarn,'trial')
          keyboard
        end
        if isfield(tmp,'delta')
          tlock_tac{ll,tt,ss}.delta={tmp.delta{(cfg.trials)}};
          tlock_tac{ll,tt,ss}.kc={tmp.kc{(cfg.trials)}};
          tlock_tac{ll,tt,ss}.sw={tmp.sw{(cfg.trials)}};
          tlock_tac{ll,tt,ss}.sp_fast={tmp.sp_fast{(cfg.trials)}};
          tlock_tac{ll,tt,ss}.sp_slow={tmp.sp_slow{(cfg.trials)}};
        end
        lastwarn('')
      end
      %       end
      
      %       if ~isnan(anumtr(ll,ss,tt))
      if ~usetr
        amrandtr=Shuffle(1:size(tlock_aud{ll,tt,ss}.trial,1));
        aurandtr=Shuffle(1:size(tlock_aud{10,tt,ss}.trial,1));
        
        tr.amtrialkept{ll,tt,ss}=sort(amrandtr(1:anumtr(ll,tt,ss)));
        tr.autrialkept{ll,tt,ss}=sort(aurandtr(1:anumtr(ll,tt,ss)));
      end
      
      if tacaud==1
        cfg=[];
        cfg.trials=tr.autrialkept{ll,tt,ss};
        if usetr && any(auwrong(soalist,tt,ss))
          tr.amtrialkept{ll,tt,ss}=tr.amtrialkept{ll,tt,ss}(setdiff(1:length(tr.amtrialkept{ll,tt,ss}),all40excluded{ll,tt,ss}));
        end
        % I have thought it through, and no minus sign needed here.
        cfg.offset=round(fsample*(tlock_aud{ll,tt,ss}.trialinfo(tr.amtrialkept{ll,tt,ss},10)-soades(ll) +soades(ll))); % offset is in samples not seconds; -2*soades required for tactile not auditory
        if length(cfg.trials)==1
          tlock_aud{ll+40,tt,ss}=[];
        else
          tmp=tlock_aud{10,tt,ss};
          warning off
          tlock_aud{ll+40,tt,ss}=ft_redefinetrial(cfg,tlock_aud{10,tt,ss}); % Talone with shift in time matching SOA jitter
          warning on
          if isfield(tmp,'delta')
            tlock_aud{ll+40,tt,ss}.delta={tmp.delta{(cfg.trials)}};
            tlock_aud{ll+40,tt,ss}.kc={tmp.kc{(cfg.trials)}};
            tlock_aud{ll+40,tt,ss}.sw={tmp.sw{(cfg.trials)}};
            tlock_aud{ll+40,tt,ss}.sp_fast={tmp.sp_fast{(cfg.trials)}};
            tlock_aud{ll+40,tt,ss}.sp_slow={tmp.sp_slow{(cfg.trials)}};
          end
          if size(tlock_aud{ll+40,tt,ss}.trialinfo,2)>240 % as it will if doing 'features'
            % need to also shift the time of the saved out negmax of Kc
            %   keyboard % shift anything for 'struct' as well?
            for nn=[50:59 70:79 90:99 110:119 130:139 150:159]
              tmpchg=find(tlock_aud{ll+40,tt,ss}.trialinfo(:,nn));
              tlock_aud{ll+40,tt,ss}.trialinfo(tmpchg,nn+10) = cfg.offset(tmpchg)/1000 + tlock_aud{ll+40,tt,ss}.trialinfo(tmpchg,nn+10);
            end
            %             for nn=1:size(tlock_aud{ll+40,tt,ss}.trialinfo,1)
            %               if ~isempty(tlock_aud{ll+40,tt,ss}.delta{nn})
            %                 [aaa,bbb]=sort([tlock_aud{ll+40,tt,ss}.delta{nn}.down])
            %               end
            %
            %             end
            
          end
        end
      end
      
      cfg=[];
      cfg.trials=tr.amtrialkept{ll,tt,ss};
      if length(cfg.trials)==1
        tlock_aud{ll,tt,ss}=[];
      else
        tmp=tlock_aud{ll,tt,ss};
        %         tlock_aud{ll,tt,ss}=ft_redefinetrial(cfg,tlock_aud{ll,tt,ss});
        tlock_aud{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{ll,tt,ss});
        if isfield(tmp,'delta')
          tlock_aud{ll,tt,ss}.delta={tmp.delta{(cfg.trials)}};
          tlock_aud{ll,tt,ss}.kc={tmp.kc{(cfg.trials)}};
          tlock_aud{ll,tt,ss}.sw={tmp.sw{(cfg.trials)}};
          tlock_aud{ll,tt,ss}.sp_fast={tmp.sp_fast{(cfg.trials)}};
          tlock_aud{ll,tt,ss}.sp_slow={tmp.sp_slow{(cfg.trials)}};
        end
        lastwarn('')
      end
      %       end
      
      %     end
    end
    
    % timwin40=[timwin(1)+0.5 timwin(2)-0.5];
    
    for ll=[soalist+40]
      %     for tt=1:4
      tt=ttin;
      if tacaud==0
        % need to use atimwin to get right range
        if ~isempty(tlock_tac{ll,tt,ss})
          cfg=[];
          cfg.trials=find(~any(isnan((tlock_tac{ll,tt,ss}.trial(:,chanuse,dsearchn(tlock_tac{ll,tt,ss}.time',atimwin(ll-40,1)):dsearchn(tlock_tac{ll,tt,ss}.time',atimwin(ll-40,2)) ))),3));
          if size(tlock_tac{ll,tt,ss}.trial,1)==1 && length(cfg.trials)<2
          else
            tmp=tlock_tac{ll,tt,ss};
            tlock_tac{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{ll,tt,ss});
            if strfind(lastwarn,'could not determine dimord') % keep information for features if present
              tlock_tac{ll,tt,ss}.delta={tmp.delta{(cfg.trials)}};
              tlock_tac{ll,tt,ss}.kc={tmp.kc{(cfg.trials)}};
              tlock_tac{ll,tt,ss}.sw={tmp.sw{(cfg.trials)}};
              tlock_tac{ll,tt,ss}.sp_fast={tmp.sp_fast{(cfg.trials)}};
              tlock_tac{ll,tt,ss}.sp_slow={tmp.sp_slow{(cfg.trials)}};
            end
            lastwarn('')
          end
        else
          tlock_tac{ll,tt,ss}=[];
        end
      end
      %     tlock_tac{ll,tt,ss}=trim_nans(tlock_tac{ll,tt,ss});
      
      if tacaud==1
        % need to use ttimwin to get right range
        if ~isempty(tlock_aud{ll,tt,ss})
          cfg=[];
          cfg.trials=find(~any(isnan((tlock_aud{ll,tt,ss}.trial(:,chanuse,dsearchn(tlock_aud{ll,tt,ss}.time',ttimwin(ll-40,1)):dsearchn(tlock_aud{ll,tt,ss}.time',ttimwin(ll-40,2)) ))),3));
          if size(tlock_aud{ll,tt,ss}.trial,1)==1 && length(cfg.trials)<2
          else
            tmp=tlock_aud{ll,tt,ss};
            tlock_aud{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{ll,tt,ss});
            if strfind(lastwarn,'could not determine dimord') % keep information for features if present
              tlock_aud{ll,tt,ss}.delta={tmp.delta{(cfg.trials)}};
              tlock_aud{ll,tt,ss}.kc={tmp.kc{(cfg.trials)}};
              tlock_aud{ll,tt,ss}.sw={tmp.sw{(cfg.trials)}};
              tlock_aud{ll,tt,ss}.sp_fast={tmp.sp_fast{(cfg.trials)}};
              tlock_aud{ll,tt,ss}.sp_slow={tmp.sp_slow{(cfg.trials)}};
            end
            lastwarn('')
          end
        else
          tlock_aud{ll,tt,ss}=[];
        end
      end
      %     tlock_aud{ll,tt,ss}=trim_nans(tlock_aud{ll,tt,ss});
      
      %     end
    end
    
    for ll=soalist
      %     for tt=1:4
      tt=ttin;
      if [tacaud==1 && [isempty(tlock_aud{ll+40,tt,ss}) || isempty(tlock_tac{ll,tt,ss}) || isempty(tlock_tac{10,tt,ss}) || isempty(tlock_nul{10,tt,ss})]  ] || [tacaud==0 && [isempty(tlock_tac{ll+40,tt,ss}) || isempty(tlock_aud{ll,tt,ss}) || isempty(tlock_aud{10,tt,ss}) || isempty(tlock_nul{10,tt,ss})]  ]
        if tacaud==1
          tr.t10trialkept{ll,tt,ss}=0;
          tr.all40trialkept{ll,tt,ss}=0;
          tr.tlltrialkept{ll,tt,ss}=0;
          tr.nllttrialkept{ll,tt,ss}=0;
        elseif tacaud==0
          tr.a10trialkept{ll,tt,ss}=0;
          tr.tll40trialkept{ll,tt,ss}=0;
          tr.alltrialkept{ll,tt,ss}=0;
          tr.nllatrialkept{ll,tt,ss}=0;
        end
        continue
      else
        if tacaud==1
          numcondtfinal(ll,tt,ss)=min([size(tlock_tac{10,tt,ss}.trial,1) size(tlock_aud{ll+40,tt,ss}.trial,1) size(tlock_tac{ll,tt,ss}.trial,1) size(tlock_nul{10,tt,ss}.trial,1)]);
        elseif tacaud==0
          numcondafinal(ll,tt,ss)=min([size(tlock_aud{10,tt,ss}.trial,1) size(tlock_tac{ll+40,tt,ss}.trial,1) size(tlock_aud{ll,tt,ss}.trial,1) size(tlock_nul{10,tt,ss}.trial,1)]);
        end
        
        if tacaud==1
          if ~usetr
            t10randtr=Shuffle(1:size(tlock_tac{10,tt,ss}.trial,1));
            all40randtr=Shuffle(1:size(tlock_aud{ll+40,tt,ss}.trial,1));
            tllrandtr=Shuffle(1:size(tlock_tac{ll,tt,ss}.trial,1));
            nllrandtr=Shuffle(1:size(tlock_nul{10,tt,ss}.trial,1));
            
            tr.t10trialkept{ll,tt,ss}=sort(t10randtr(1:numcondtfinal(ll,tt,ss)));
            tr.all40trialkept{ll,tt,ss}=sort(all40randtr(1:numcondtfinal(ll,tt,ss)));
            tr.tlltrialkept{ll,tt,ss}=sort(tllrandtr(1:numcondtfinal(ll,tt,ss)));
            tr.nllttrialkept{ll,tt,ss}=sort(nllrandtr(1:numcondtfinal(ll,tt,ss)));
          elseif usetr && any(auwrong(soalist,tt,ss))
            if numcondtfinal(ll,tt,ss)==size(tlock_aud{ll+40,tt,ss}.trial,1)
              tr.t10trialkept{ll,tt,ss}=tr.t10trialkept{ll,tt,ss}(setdiff(1:length(tr.t10trialkept{ll,tt,ss}),all40excluded{ll,tt,ss}));
              tr.nllttrialkept{ll,tt,ss}=tr.nllttrialkept{ll,tt,ss}(setdiff(1:length(tr.nllttrialkept{ll,tt,ss}),all40excluded{ll,tt,ss}));
              tr.all40trialkept{ll,tt,ss}=1:numcondtfinal(ll,tt,ss);
              tr.tlltrialkept{ll,tt,ss}=tr.tlltrialkept{ll,tt,ss}(setdiff(1:length(tr.tlltrialkept{ll,tt,ss}),all40excluded{ll,tt,ss}));
            elseif numcondtfinal(ll,tt,ss)==size(tlock_tac{10,tt,ss}.trial,1)
              tr.t10trialkept{ll,tt,ss}=1:numcondtfinal(ll,tt,ss);
              tr.nllttrialkept{ll,tt,ss}=tr.nllttrialkept{ll,tt,ss}(setdiff(1:length(tr.nllttrialkept{ll,tt,ss}),tac10orig_excluded{tt,ss}));              
              tr.all40trialkept{ll,tt,ss}=tr.all40trialkept{ll,tt,ss}(setdiff(1:length(tr.all40trialkept{ll,tt,ss}),tac10orig_excluded{tt,ss}));
              tr.tlltrialkept{ll,tt,ss}=tr.tlltrialkept{ll,tt,ss}(setdiff(1:length(tr.tlltrialkept{ll,tt,ss}),tac10orig_excluded{tt,ss}));
            elseif numcondtfinal(ll,tt,ss)==size(tlock_tac{ll,tt,ss}.trial,1)
              tr.t10trialkept{ll,tt,ss}=tr.t10trialkept{ll,tt,ss}(setdiff(1:length(tr.t10trialkept{ll,tt,ss}),tacllorig_excluded{ll,tt,ss}));
              tr.nllttrialkept{ll,tt,ss}=tr.nllttrialkept{ll,tt,ss}(setdiff(1:length(tr.nllttrialkept{ll,tt,ss}),tacllorig_excluded{ll,tt,ss}));              
              tr.all40trialkept{ll,tt,ss}=tr.all40trialkept{ll,tt,ss}(setdiff(1:length(tr.all40trialkept{ll,tt,ss}),tacllorig_excluded{ll,tt,ss}));
              tr.tlltrialkept{ll,tt,ss}=1:numcondtfinal(ll,tt,ss);
            elseif numcondtfinal(ll,tt,ss)==size(tlock_nul{10,tt,ss}.trial,1)
              tr.t10trialkept{ll,tt,ss}=tr.t10trialkept{ll,tt,ss}(setdiff(1:length(tr.t10trialkept{ll,tt,ss}),nul10orig_excluded{tt,ss}));
              tr.nllttrialkept{ll,tt,ss}=1:numcondtfinal(ll,tt,ss);
              tr.all40trialkept{ll,tt,ss}=tr.all40trialkept{ll,tt,ss}(setdiff(1:length(tr.all40trialkept{ll,tt,ss}),nul10orig_excluded{tt,ss}));
              tr.tlltrialkept{ll,tt,ss}=tr.tlltrialkept{ll,tt,ss}(setdiff(1:length(tr.tlltrialkept{ll,tt,ss}),nul10orig_excluded{tt,ss}));
            else
              error('now what?')
            end
            if max(tr.all40trialkept{ll,tt,ss})>size(tlock_aud{ll+40,tt,ss}.trial,1)
              tr.all40trialkept{ll,tt,ss}=dsearchn(tr.autrialkept{ll,tt,ss}',tr.all40trialkept{ll,tt,ss}')';
            end

            tr.t10trialkept{ll,tt,ss}=unique(tr.t10trialkept{ll,tt,ss});
            tr.nllttrialkept{ll,tt,ss}=unique(tr.nllttrialkept{ll,tt,ss});
            tr.all40trialkept{ll,tt,ss}=unique(tr.all40trialkept{ll,tt,ss});
            tr.tlltrialkept{ll,tt,ss}=unique(tr.tlltrialkept{ll,tt,ss});
            newmin=min([length(tr.t10trialkept{ll,tt,ss}) length(tr.nllttrialkept{ll,tt,ss}) length(tr.all40trialkept{ll,tt,ss}) length(tr.tlltrialkept{ll,tt,ss})]);
            if numcondtfinal(ll,tt,ss)>newmin
              numcondtfinal(ll,tt,ss)=newmin;
              tr.t10trialkept{ll,tt,ss}=tr.t10trialkept{ll,tt,ss}(1:newmin);
              tr.nllttrialkept{ll,tt,ss}=tr.nllttrialkept{ll,tt,ss}(1:newmin);
              tr.all40trialkept{ll,tt,ss}=tr.all40trialkept{ll,tt,ss}(1:newmin);
              tr.tlltrialkept{ll,tt,ss}=tr.tlltrialkept{ll,tt,ss}(1:newmin);
            end
          end
        end
        
        if tacaud==0
          error('fix me for auwrong')
          if ~usetr
            a10randtr=Shuffle(1:size(tlock_aud{10,tt,ss}.trial,1));
            tll40randtr=Shuffle(1:size(tlock_tac{ll+40,tt,ss}.trial,1));
            allrandtr=Shuffle(1:size(tlock_aud{ll,tt,ss}.trial,1));
            nllrandtr=Shuffle(1:size(tlock_nul{10,tt,ss}.trial,1));
            
            tr.a10trialkept{ll,tt,ss}=sort(a10randtr(1:numcondafinal(ll,tt,ss)));
            tr.tll40trialkept{ll,tt,ss}=sort(tll40randtr(1:numcondafinal(ll,tt,ss)));
            tr.alltrialkept{ll,tt,ss}=sort(allrandtr(1:numcondafinal(ll,tt,ss)));
            tr.nllatrialkept{ll,tt,ss}=sort(nllrandtr(1:numcondafinal(ll,tt,ss)));
          end
        end
        
        %     cfgtrialstac10orig=cfg.trials;
        %     if 1
        %       length(find(cfg.trials))-length(find(cfg.trials & data_tac_ref.trialinfo(:,23)==0))
        %       cfg.trials=cfg.trials & data_tac_ref.trialinfo(:,23)==0; % exclude EOG.
        %     end
        %     origtrialindex=find(cfgtrialstac10orig);
        %     tac10orig_excluded=dsearchn(origtrialindex,setdiff(find(cfgtrialstac10orig),find(cfg.trials)));
        %     origindex=setdiff(1:length(find(cfgtrialstac10orig)),tac10orig_excluded);
        %     for ll=soalist
        %       tr.t10trialkept{ll,tt-2,stage+10}=dsearchn(origindex',tr.t10trialkept{ll,tt-2,stage+10}')';
        %     end
        
        %       origtrialindex=find(trialstmp);
        %       for ll=soalist,
        %         if isempty(max(tr.autrialkept{ll,tt-2,stage+10})>length(find(cfg.trials)))
        %           auwrong(ll)=0;
        %         else
        %           auwrong(ll)=max(tr.autrialkept{ll,tt-2,stage+10})>length(find(cfg.trials));
        %         end;
        %       end
        %       if any(auwrong(soalist))  % horrible hack to make up for not saving out tr.autrialkept correctly awhile ago.
        %         poststimeogexcluded=dsearchn(origtrialindex,setdiff(origtrialindex,find(cfg.trials)));
        %         origindex=setdiff(1:length(origtrialindex),poststimeogexcluded);
        %         origcfgtrials=cfg.trials;
        %         for ll=soalist,  % one of the ll needs to be 10 instead?
        %           %           cfg.trials=origcfgtrials;
        %           cfg.trials(origtrialindex(tr.autrialkept{ll,tt-2,stage+10}(find(diff(dsearchn(origindex',tr.autrialkept{ll,tt-2,stage+10}')')<1)+1)))=1;
        %         end
        %         origindex2=setdiff(1:length(origtrialindex),dsearchn(origtrialindex,setdiff(origtrialindex,find(cfg.trials))));
        %         for ll=soalist
        %           tr.autrialkept{ll,tt-2,stage+10}=dsearchn(origindex2',tr.autrialkept{ll,tt-2,stage+10}');
        %         end
        %       end
        
        
        if tacaud==1
          cfg=[];
          if usetr && max(tr.all40trialkept{ll,tt,ss})>size(tlock_aud{ll+40,tt,ss}.trial,1)
            error('fix me')
%             tr.all40trialkept{ll,tt,ss}=dsearchn(tr.autrialkept{ll,tt,ss}',tr.all40trialkept{ll,tt,ss}')';
          end
          cfg.trials=tr.all40trialkept{ll,tt,ss};
          %           if usetr && any(auwrong(soalist)) && length(cfg.trials)>size(tlock_aud{ll+40,tt,ss}.trial,1)
          %             % continuation of horrible hack above
          %             cfg.trials=cfg.trials(1:size(tlock_aud{ll+40,tt,ss}.trial,1));
          %             tr.all40trialkept{ll,tt,ss}=tr.all40trialkept{ll,tt,ss}(1:size(tlock_aud{ll+40,tt,ss}.trial,1));
          %           end
          %           if usetr && any(auwrong(soalist))
          %             cfgtrialsorig=cfg.trials;
          %           end
          
          if size(tlock_aud{ll,tt,ss}.trial,1)==1 && length(cfg.trials)<2
          else
            tmp=tlock_aud{ll+40,tt,ss};
            tlock_aud{ll+40,tt,ss}=ft_selectdata(cfg,tlock_aud{ll+40,tt,ss});
            if strfind(lastwarn,'could not determine dimord') % keep information for features if present
              tlock_aud{ll+40,tt,ss}.delta={tmp.delta{(cfg.trials)}};
              tlock_aud{ll+40,tt,ss}.kc={tmp.kc{(cfg.trials)}};
              tlock_aud{ll+40,tt,ss}.sw={tmp.sw{(cfg.trials)}};
              tlock_aud{ll+40,tt,ss}.sp_fast={tmp.sp_fast{(cfg.trials)}};
              tlock_aud{ll+40,tt,ss}.sp_slow={tmp.sp_slow{(cfg.trials)}};
            end
            lastwarn('')
          end
          cfg=[];
          cfg.trials=tr.tlltrialkept{ll,tt,ss};
          if size(tlock_aud{ll,tt,ss}.trial,1)==1 && length(cfg.trials)<2
          else
            tmp=tlock_tac{ll,tt,ss};
            tlock_tac{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{ll,tt,ss});
            if strfind(lastwarn,'could not determine dimord') % keep information for features if present
              tlock_tac{ll,tt,ss}.delta={tmp.delta{(cfg.trials)}};
              tlock_tac{ll,tt,ss}.kc={tmp.kc{(cfg.trials)}};
              tlock_tac{ll,tt,ss}.sw={tmp.sw{(cfg.trials)}};
              tlock_tac{ll,tt,ss}.sp_fast={tmp.sp_fast{(cfg.trials)}};
              tlock_tac{ll,tt,ss}.sp_slow={tmp.sp_slow{(cfg.trials)}};
            end
            lastwarn('')
          end
          if 0 % too memory intensive to do here
            cfg=[];
            cfg.trials=tr.t10trialkept{ll,tt,ss};
            tmp=tlock_tac{10,tt,ss};
            tlock_tac{ll+20,tt,ss}=ft_selectdata(cfg,tlock_tac{10,tt,ss}); % create +20 here as 10 will be different for every ll,tt,ss
            if strfind(lastwarn,'could not determine dimord') % keep information for features if present
              tlock_tac{ll+20,tt,ss}.delta={tmp.delta{(cfg.trials)}};
              tlock_tac{ll+20,tt,ss}.kc={tmp.kc{(cfg.trials)}};
              tlock_tac{ll+20,tt,ss}.sw={tmp.sw{(cfg.trials)}};
              tlock_tac{ll+20,tt,ss}.sp_fast={tmp.sp_fast{(cfg.trials)}};
              tlock_tac{ll+20,tt,ss}.sp_slow={tmp.sp_slow{(cfg.trials)}};
            end
            lastwarn('')
            cfg=[];
            cfg.trials=tr.nllttrialkept{ll,tt,ss};
            tmp=tlock_nul{10,tt,ss};
            tlock_nul{ll+50,tt,ss}=ft_selectdata(cfg,tlock_nul{10,tt,ss});
            if strfind(lastwarn,'could not determine dimord') % keep information for features if present
              tlock_nul{ll+50,tt,ss}.delta={tmp.delta{(cfg.trials)}};
              tlock_nul{ll+50,tt,ss}.kc={tmp.kc{(cfg.trials)}};
              tlock_nul{ll+50,tt,ss}.sw={tmp.sw{(cfg.trials)}};
              tlock_nul{ll+50,tt,ss}.sp_fast={tmp.sp_fast{(cfg.trials)}};
              tlock_nul{ll+50,tt,ss}.sp_slow={tmp.sp_slow{(cfg.trials)}};
            end
            lastwarn('')
          end
        end
        
        if tacaud==0
          cfg=[];
          cfg.trials=tr.tll40trialkept{ll,tt,ss};
          if size(tlock_aud{ll,tt,ss}.trial,1)==1 && length(cfg.trials)<2
          else
            tmp=tlock_tac{ll+40,tt,ss};
            tlock_tac{ll+40,tt,ss}=ft_selectdata(cfg,tlock_tac{ll+40,tt,ss});
            if strfind(lastwarn,'could not determine dimord') % keep information for features if present
              tlock_tac{ll+40,tt,ss}.delta={tmp.delta{(cfg.trials)}};
              tlock_tac{ll+40,tt,ss}.kc={tmp.kc{(cfg.trials)}};
              tlock_tac{ll+40,tt,ss}.sw={tmp.sw{(cfg.trials)}};
              tlock_tac{ll+40,tt,ss}.sp_fast={tmp.sp_fast{(cfg.trials)}};
              tlock_tac{ll+40,tt,ss}.sp_slow={tmp.sp_slow{(cfg.trials)}};
            end
            lastwarn('')
          end
          cfg=[];
          cfg.trials=tr.alltrialkept{ll,tt,ss};
          if size(tlock_aud{ll,tt,ss}.trial,1)==1 && length(cfg.trials)<2
          else
            tmp=tlock_aud{ll,tt,ss};
            tlock_aud{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{ll,tt,ss});
            if strfind(lastwarn,'could not determine dimord') % keep information for features if present
              tlock_aud{ll,tt,ss}.delta={tmp.delta{(cfg.trials)}};
              tlock_aud{ll,tt,ss}.kc={tmp.kc{(cfg.trials)}};
              tlock_aud{ll,tt,ss}.sw={tmp.sw{(cfg.trials)}};
              tlock_aud{ll,tt,ss}.sp_fast={tmp.sp_fast{(cfg.trials)}};
              tlock_aud{ll,tt,ss}.sp_slow={tmp.sp_slow{(cfg.trials)}};
            end
            lastwarn('')
          end
          if 0 % too memory intensive to do here
            cfg=[];
            cfg.trials=tr.a10trialkept{ll,tt,ss};
            tmp=tlock_aud{10,tt,ss};
            tlock_aud{ll+20,tt,ss}=ft_selectdata(cfg,tlock_aud{10,tt,ss}); % create +20 here as 10 will be different for every ll,tt,ss
            if strfind(lastwarn,'could not determine dimord') % keep information for features if present
              tlock_aud{ll+20,tt,ss}.delta={tmp.delta{(cfg.trials)}};
              tlock_aud{ll+20,tt,ss}.kc={tmp.kc{(cfg.trials)}};
              tlock_aud{ll+20,tt,ss}.sw={tmp.sw{(cfg.trials)}};
              tlock_aud{ll+20,tt,ss}.sp_fast={tmp.sp_fast{(cfg.trials)}};
              tlock_aud{ll+20,tt,ss}.sp_slow={tmp.sp_slow{(cfg.trials)}};
            end
            lastwarn('')
            cfg=[];
            cfg.trials=tr.nllatrialkept{ll,tt,ss};
            tmp=tlock_nul{10,tt,ss};
            tlock_nul{ll+60,tt,ss}=ft_selectdata(cfg,tlock_nul{10,tt,ss});
            if strfind(lastwarn,'could not determine dimord') % keep information for features if present
              tlock_nul{ll+60,tt,ss}.delta={tmp.delta{(cfg.trials)}};
              tlock_nul{ll+60,tt,ss}.kc={tmp.kc{(cfg.trials)}};
              tlock_nul{ll+60,tt,ss}.sw={tmp.sw{(cfg.trials)}};
              tlock_nul{ll+60,tt,ss}.sp_fast={tmp.sp_fast{(cfg.trials)}};
              tlock_nul{ll+60,tt,ss}.sp_slow={tmp.sp_slow{(cfg.trials)}};
            end
            lastwarn('')
          end
        end
        %     end
      end % if
    end % ll
    
    for ll=[1 3 4] % for multisensory-shifted comparison
      tt=ttin;
      if [tacaud==1 && [isempty(tlock_tac{ll,tt,ss}) || isempty(tlock_tac{10-ll,tt,ss}) || isempty(tlock_tac{5,tt,ss})]] || [tacaud==0 && [isempty(tlock_aud{ll,tt,ss}) || isempty(tlock_aud{10-ll,tt,ss}) || isempty(tlock_aud{5,tt,ss})] ]
        if tacaud==1
          tr.tmstrialkept1{ll,tt,ss}=0;
          tr.tmstrialkept2{ll,tt,ss}=0;
          tr.tmstrialkept3{ll,tt,ss}=0;
          tr.tmstrialkept4{ll,tt,ss}=0;
        elseif tacaud==0
          tr.amstrialkept1{ll,tt,ss}=0;
          tr.amstrialkept2{ll,tt,ss}=0;
          tr.amstrialkept3{ll,tt,ss}=0;
          tr.amstrialkept4{ll,tt,ss}=0;
        end
        continue
      else
        if tacaud==1
          numcondmstfinal(ll,tt,ss)=min([size(tlock_tac{ll,tt,ss}.trial,1) size(tlock_tac{10-ll,tt,ss}.trial,1) size(tlock_tac{5,tt,ss}.trial,1)]);
        elseif tacaud==0
          numcondmsafinal(ll,tt,ss)=min([size(tlock_aud{ll,tt,ss}.trial,1) size(tlock_aud{10-ll,tt,ss}.trial,1) size(tlock_aud{5,tt,ss}.trial,1)]);
        end
        
        if tacaud==1
          if ~usetr
            tms1randtr=Shuffle(1:size(tlock_tac{ll,tt,ss}.trial,1));
            tms2randtr=Shuffle(1:size(tlock_tac{10-ll,tt,ss}.trial,1));
            tms3randtr=Shuffle(1:size(tlock_tac{5,tt,ss}.trial,1));
            tms4randtr=Shuffle(1:size(tlock_tac{5,tt,ss}.trial,1)); % do this twice to get 2nd random subsampling
            tr.tmstrialkept1{ll,tt,ss}=sort(tms1randtr(1:numcondmstfinal(ll,tt,ss)));
            tr.tmstrialkept2{ll,tt,ss}=sort(tms2randtr(1:numcondmstfinal(ll,tt,ss)));
            tr.tmstrialkept3{ll,tt,ss}=sort(tms3randtr(1:numcondmstfinal(ll,tt,ss)));
            tr.tmstrialkept4{ll,tt,ss}=sort(tms4randtr(1:numcondmstfinal(ll,tt,ss)));
          end
        end
        
        if tacaud==0
          if ~usetr
            ams1randtr=Shuffle(1:size(tlock_aud{ll,tt,ss}.trial,1));
            ams2randtr=Shuffle(1:size(tlock_aud{10-ll,tt,ss}.trial,1));
            ams3randtr=Shuffle(1:size(tlock_aud{5,tt,ss}.trial,1));
            ams4randtr=Shuffle(1:size(tlock_aud{5,tt,ss}.trial,1)); % do this twice to get 2nd random subsampling
            tr.amstrialkept1{ll,tt,ss}=sort(ams1randtr(1:numcondmsafinal(ll,tt,ss)));
            tr.amstrialkept2{ll,tt,ss}=sort(ams2randtr(1:numcondmsafinal(ll,tt,ss)));
            tr.amstrialkept3{ll,tt,ss}=sort(ams3randtr(1:numcondmsafinal(ll,tt,ss)));
            tr.amstrialkept4{ll,tt,ss}=sort(ams4randtr(1:numcondmsafinal(ll,tt,ss)));
          end
        end
      end
    end % ll
    
  end % ss
  
else % if stageuse empty
  warning(['No useful stages for this call to TrialSelection' num2str(ii) num2str(ttin) num2str(sleep) num2str(tacaud)])
  tlock_tac=[];
  tlock_aud=[];
  tlock_nul=[];
  numcond=0;
  eog_exclude=0;
end

if isempty(whos('numcond*'))
  numcond=0;
end

if ~exist('tlock_tac','var')
  tlock_tac=[];
end
if ~exist('tlock_aud','var')
  tlock_aud=[];
end
if ~exist('tlock_nul','var')
  tlock_nul=[];
end

tr.tacaud=tacaud;
% save('trialkept.mat','tr','numcond*')

calling=dbstack; % find out which function called this one
if tomflag
  
  if strcmp(calling(end).name,'eeg_freqanalysis_sensor1')
    save(['trialkeptTFR_tt' num2str(ttin) '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '.mat'],'tr','numcond*','*eog_exclude')
  elseif strcmp(calling(end).name,'eeg_legomagic_erp_stats2_sepTacAud')
    save(['trialkept_tt' num2str(ttin) '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '.mat'],'tr','numcond*','*eog_exclude')
  elseif strcmp(calling(end).name,'eeg_legomagic_featurestats_sepTacAud')
    save(['trialkeptFeature_tt' num2str(ttin) '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '.mat'],'tr','numcond*','*eog_exclude')
  elseif strcmp(calling(end).name,'eeg_legomagic_featurestats1_sepTacAud') %% comment out later - TW
    save(['trialkeptFeature_tt' num2str(ttin) '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '.mat'],'tr','numcond*','*eog_exclude')
  elseif strcmp(calling(end).name,'eeg_legomagic_featurestats2_sepTacAud')%% comment out later - TW
    save(['trialkeptFeature_tt' num2str(ttin) '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '.mat'],'tr','numcond*','*eog_exclude')
  else
    warning('which function called this? nothing saved out')
  end
  
else
  if strcmp(calling(end).name,'eeg_freqanalysis_sensor1')
    save(['trialkeptTFR_tt' num2str(ttin) '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_usetr' num2str(usetr) '_trialkc' num2str(trialkc) '.mat'],'tr','numcond*','*eog_exclude')
  elseif strcmp(calling(end).name,'eeg_legomagic_erp_stats2_sepTacAud')
    save(['trialkept_tt' num2str(ttin) '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_usetr' num2str(usetr) '_trialkc' num2str(trialkc) '.mat'],'tr','numcond*','*eog_exclude')
  elseif strcmp(calling(end).name,'eeg_legomagic_featurestats_sepTacAud')
    save(['trialkeptFeature_tt' num2str(ttin) '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_usetr' num2str(usetr) '_trialkc' num2str(trialkc) '.mat'],'tr','numcond*','*eog_exclude')
  else
    warning('which function called this? nothing saved out')
  end
end




