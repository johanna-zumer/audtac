function [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection2_wakeSleep(ii,sleep,ttin)
% function [tlock_tac,tlock_aud,tlock_nul]=eeg_legomagic_trialSelection(ii)
% Trial selection: equal trial numbers per condition and removing NaN at ends

if ispc
  edir='D:\audtac\eeg_data\';
  ddir='D:\audtac\legomagic\diaries\';
else
  edir='/mnt/hgfs/D/audtac/eeg_data/';
  ddir='/mnt/hgfs/D/audtac/legomagic/diaries/';
end
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

cd([edir sub{ii} ])
%%

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

if ii<8
  soalist=[3 4 5 6 7];
else
  soalist=[1 3 4 5 6 7 9];
end


% only loop through stages actually present in the given dataset
if 0 % ideal world with more memory process for all stages
  stageuse=[0:3 5];
else
  if sleep
    stageuse=[2 3];
  else
    stageuse=[0 ];
  end
end
for ss=stageuse
  if ~any(data_tac_ref.trialinfo(:,2)==ss)  % testing tac should be good enough (rare chance if only auditory and nul trials are in a segment of sleep that not includes tactile, then we can't use those anyway, right?)
    stageuse=setdiff(stageuse,ss);
  end
end
tr.stageuse=stageuse;

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
      if ~isempty(find(cfg.trials))
        tlock_tac{soalist(ll),tt-2,stage+10}=ft_timelockanalysis(cfg,data_tac_ref);
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
    if ~isempty(find(cfg.trials))
      tlock_tac{10,tt-2,stage+10}=ft_timelockanalysis(cfg,data_tac_ref); %technically this doesn't depend on 'tt' but to make it consistent with other lines of code and this new way of calling 'ttin' as input
    end
    % cfg.trials=data_tac_ref.trialinfo(:,3)==-2;
    % tlock_tac_sall{10}=ft_timelockanalysis(cfg,data_tac_ref);
  end
  taceog_exclude=find(data_tac_ref.trialinfo(:,13)==0);
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
      if ~isempty(find(cfg.trials))
        tlock_aud{soalist(ll),tt-2,stage+10}=ft_timelockanalysis(cfg,data_aud_ref);
      end
      %   cfg.trials=data_aud_ref.trialinfo(:,3)==soalist(ll);
      %   tlock_aud_sall{soalist(ll)}=ft_timelockanalysis(cfg,data_aud_ref);
      %     end
    end
    cfg=[]; % aud alone 39 (-1)
    cfg.vartrllength=2;
    cfg.keeptrials='yes';
    cfg.trials=[data_aud_ref.trialinfo(:,3)==-1 & data_aud_ref.trialinfo(:,2)==stage];
    cfg.trials=cfg.trials & data_aud_ref.trialinfo(:,13)==0; % exclude EOG.
    if ~isempty(find(cfg.trials))
      tlock_aud{10,tt-2,stage+10}=ft_timelockanalysis(cfg,data_aud_ref);
    end
    % cfg.trials=data_aud_ref.trialinfo(:,3)==-1;
    % tlock_aud_sall{10}=ft_timelockanalysis(cfg,data_aud_ref);
  end
  audeog_exclude=find(data_aud_ref.trialinfo(:,13)==0);
  clear data_aud_ref
  
  % Do tac, then aud, then nul for memory reasons
  for stage=stageuse
    cfg=[];
    cfg.vartrllength=2;
    cfg.keeptrials='yes';
    cfg.trials=data_nul_ref.trialinfo(:,2)==stage;
    cfg.trials=cfg.trials & data_nul_ref.trialinfo(:,13)==0; % exclude EOG.
    if ~isempty(find(cfg.trials))
      tlock_nul{10,tt-2,stage+10}=ft_timelockanalysis(cfg,data_nul_ref);
    end
    % cfg.trials='all';
    % tlock_nul_sall{10}=ft_timelockanalysis(cfg,data_nul_ref_ref);
  end
  nuleog_exclude=find(data_nul_ref.trialinfo(:,13)==0);
  clear data_nul_ref
  
  
  
  
  %%
  
  
  soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];
  ttimwin=[min(-.2,soades-.2); max(0.5,soades+.5)]';
  ttimwin(10,:)=[min(ttimwin(1:9,1)) max(ttimwin(1:9,2))];
  atimwin=[min(-.2,fliplr(soades)-.2); max(0.5,fliplr(soades)+.5)]';
  atimwin(10,:)=ttimwin(10,:);
  ntimwin(10,:)=ttimwin(10,:);
  
  chanuse=match_str(tlock_tac{3,ttin,stage+10}.label,'Fz');
  
  for ss=stageuse+10
    for ll=[soalist 10] % this loop is to remove trials with NaN in (e.g. if artifact rejection cut through middle of a trial)
      %     for tt=1:4
      tt=ttin;
      
      %       if ll==10 && tt>1
      %         continue
      %       else
      if ~isempty(tlock_tac{ll,tt,ss})
        cfg=[];
        %       cfg.keeptrials='yes';
        %         try
        cfg.trials=find(~any(isnan((tlock_tac{ll,tt,ss}.trial(:,chanuse,dsearchn(tlock_tac{ll,tt,ss}.time',ttimwin(ll,1)):dsearchn(tlock_tac{ll,tt,ss}.time',ttimwin(ll,2)) ))) ,3));
        %       tlock_tac{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
        if size(tlock_tac{ll,tt,ss}.trial,1)==1 && cfg.trials==1
        else
          tlock_tac{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{ll,tt,ss});
        end
        
        tlock_tac{ll,tt,ss}=trim_nans(tlock_tac{ll,tt,ss});
      else
        continue
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
        if size(tlock_aud{ll,tt,ss}.trial,1)==1 && cfg.trials==1
        else
          tlock_aud{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{ll,tt,ss});
        end
        %       tmp=squeeze(mean(tlock_aud{ll,tt,ss}.trial,1));
        %       if sum(isnan(tmp(1,:)))>.07*length(tmp(1,:))
        %         error('this averaging didnt work');
        %       end
        tlock_aud{ll,tt,ss}=trim_nans(tlock_aud{ll,tt,ss});
      else
        continue
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
        if size(tlock_nul{ll,tt,ss}.trial,1)==1 && cfg.trials==1
        else
          tlock_nul{ll,tt,ss}=ft_selectdata(cfg,tlock_nul{ll,tt,ss});
        end
        %       tmp=squeeze(mean(tlock_nul{ll,tt,ss}.trial,1));
        %       if sum(isnan(tmp(1,:)))>.07*length(tmp(1,:))
        %         error('this averaging didnt work');
        %       end
        tlock_nul{ll,tt,ss}=trim_nans(tlock_nul{ll,tt,ss});
        %         catch
        %           disp(['no tlock_nul for ' num2str(ll) num2str(tt) num2str(ss)])
        %         end
        
      end
      if ll<10
        % make sure the multisensory auditory-locked and ms tactile-locked have same trials; sometimes not exactly!
        %         try
        trcom=intersect(tlock_tac{ll,tt,ss}.trialinfo(:,11),tlock_aud{ll,tt,ss}.trialinfo(:,11));
        [~,tbb,~]=intersect(tlock_tac{ll,tt,ss}.trialinfo(:,11),trcom);
        [~,abb,~]=intersect(tlock_aud{ll,tt,ss}.trialinfo(:,11),trcom);
        cfg=[];
        cfg.trials=tbb;
        if size(tlock_tac{ll,tt,ss}.trial,1)==1 && cfg.trials==1
        else
          tlock_tac{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{ll,tt,ss});
        end
        cfg=[];
        cfg.trials=abb;
        if size(tlock_aud{ll,tt,ss}.trial,1)==1 && cfg.trials==1
        else
          tlock_aud{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{ll,tt,ss});
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
  
    fsample=1/diff(tlock_tac{10,tt,ss}.time(1:2)); % sampling rate
  
  tnumtr=zeros(max(soalist),ttin,15);
  anumtr=zeros(max(soalist),ttin,15);
  % % fix not checking ss if not exist
  for ss=stageuse+10
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
        tmrandtr=Shuffle(1:size(tlock_tac{ll,tt,ss}.trial,1));
        turandtr=Shuffle(1:size(tlock_tac{10,tt,ss}.trial,1));
        
        tmtrialkept{ll,tt,ss}=sort(tmrandtr(1:tnumtr(ll,tt,ss)));
        tutrialkept{ll,tt,ss}=sort(turandtr(1:tnumtr(ll,tt,ss)));
      end
      
      cfg=[];
      % I have thought it through, and no minus sign needed here.
      cfg.offset=round(fsample*(tlock_tac{ll,tt,ss}.trialinfo(tmtrialkept{ll,tt,ss},10)-soades(ll) -soades(ll))); % offset is in samples not seconds; -2*soades required for tactile not auditory
      cfg.trials=tutrialkept{ll,tt,ss};
      if length(cfg.trials)==1
        tlock_tac{ll+40,tt,ss}=[];
      else
        tlock_tac{ll+40,tt,ss}=ft_redefinetrial(cfg,tlock_tac{10,tt,ss}); % Talone with shift in time matching SOA jitter
      end
      
      cfg=[];
      cfg.trials=tmtrialkept{ll,tt,ss}; % only keep ms trials that have a paired shifted corresponding unisensory trial
      if length(cfg.trials)==1
        tlock_tac{ll,tt,ss}=[];
      else
        tlock_tac{ll,tt,ss}=ft_redefinetrial(cfg,tlock_tac{ll,tt,ss});
      end
      %       end
      
      %       if ~isnan(anumtr(ll,ss,tt))
      amrandtr=Shuffle(1:size(tlock_aud{ll,tt,ss}.trial,1));
      aurandtr=Shuffle(1:size(tlock_aud{10,tt,ss}.trial,1));
      
      amtrialkept{ll,tt,ss}=sort(amrandtr(1:anumtr(ll,tt,ss)));
      autrialkept{ll,tt,ss}=sort(aurandtr(1:anumtr(ll,tt,ss)));
      
      cfg=[];
      % I have thought it through, and no minus sign needed here.
      cfg.offset=round(fsample*(tlock_aud{ll,tt,ss}.trialinfo(amtrialkept{ll,tt,ss},10)-soades(ll) +soades(ll))); % offset is in samples not seconds; -2*soades required for tactile not auditory
      cfg.trials=autrialkept{ll,tt,ss};
      if length(cfg.trials)==1
        tlock_aud{ll+40,tt,ss}=[];
      else
        tlock_aud{ll+40,tt,ss}=ft_redefinetrial(cfg,tlock_aud{10,tt,ss}); % Talone with shift in time matching SOA jitter
      end
      
      cfg=[];
      cfg.trials=amtrialkept{ll,tt,ss};
      if length(cfg.trials)==1
        tlock_aud{ll,tt,ss}=[];
      else
        tlock_aud{ll,tt,ss}=ft_redefinetrial(cfg,tlock_aud{ll,tt,ss});
      end
      %       end
      
      %     end
    end
    
    % timwin40=[timwin(1)+0.5 timwin(2)-0.5];
    
    for ll=[soalist+40]
      %     for tt=1:4
      tt=ttin;
      % need to use atimwin to get right range
      if ~isempty(tlock_tac{ll,tt,ss})
        cfg=[];
        cfg.trials=find(~any(isnan((tlock_tac{ll,tt,ss}.trial(:,chanuse,dsearchn(tlock_tac{ll,tt,ss}.time',atimwin(ll-40,1)):dsearchn(tlock_tac{ll,tt,ss}.time',atimwin(ll-40,2)) ))),3));
        if size(tlock_aud{ll,tt,ss}.trial,1)==1 && cfg.trials==1
        else
          tlock_tac{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{ll,tt,ss});
        end
      else
        tlock_tac{ll,tt,ss}=[];
      end
      %     tlock_tac{ll,tt,ss}=trim_nans(tlock_tac{ll,tt,ss});
      
      % need to use ttimwin to get right range
      if ~isempty(tlock_aud{ll,tt,ss})
        cfg=[];
        cfg.trials=find(~any(isnan((tlock_aud{ll,tt,ss}.trial(:,chanuse,dsearchn(tlock_aud{ll,tt,ss}.time',ttimwin(ll-40,1)):dsearchn(tlock_aud{ll,tt,ss}.time',ttimwin(ll-40,2)) ))),3));
        if size(tlock_aud{ll,tt,ss}.trial,1)==1 && cfg.trials==1
        else
          tlock_aud{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{ll,tt,ss});
        end
      else
        tlock_aud{ll,tt,ss}=[];
      end
      %     tlock_aud{ll,tt,ss}=trim_nans(tlock_aud{ll,tt,ss});
      
      %     end
    end
    
    for ll=soalist
      %     for tt=1:4
      tt=ttin;
      if isempty(tlock_aud{ll+40,tt,ss}) || isempty(tlock_tac{ll+40,tt,ss})
        continue
      else
        numcondtfinal(ll,tt,ss)=min([size(tlock_tac{10,tt,ss}.trial,1) size(tlock_aud{ll+40,tt,ss}.trial,1) size(tlock_tac{ll,tt,ss}.trial,1) size(tlock_nul{10,tt,ss}.trial,1)]);
        numcondafinal(ll,tt,ss)=min([size(tlock_aud{10,tt,ss}.trial,1) size(tlock_tac{ll+40,tt,ss}.trial,1) size(tlock_aud{ll,tt,ss}.trial,1) size(tlock_nul{10,tt,ss}.trial,1)]);
        
        t10randtr=Shuffle(1:size(tlock_tac{10,tt,ss}.trial,1));
        all40randtr=Shuffle(1:size(tlock_aud{ll+40,tt,ss}.trial,1));
        tllrandtr=Shuffle(1:size(tlock_tac{ll,tt,ss}.trial,1));
        nllrandtr=Shuffle(1:size(tlock_nul{10,tt,ss}.trial,1));
        
        tr.t10trialkept{ll,tt,ss}=sort(t10randtr(1:numcondtfinal(ll,tt,ss)));
        tr.all40trialkept{ll,tt,ss}=sort(all40randtr(1:numcondtfinal(ll,tt,ss)));
        tr.tlltrialkept{ll,tt,ss}=sort(tllrandtr(1:numcondtfinal(ll,tt,ss)));
        tr.nllttrialkept{ll,tt,ss}=sort(nllrandtr(1:numcondtfinal(ll,tt,ss)));
        
        a10randtr=Shuffle(1:size(tlock_aud{10,tt,ss}.trial,1));
        tll40randtr=Shuffle(1:size(tlock_tac{ll+40,tt,ss}.trial,1));
        allrandtr=Shuffle(1:size(tlock_aud{ll,tt,ss}.trial,1));
        nllrandtr=Shuffle(1:size(tlock_nul{10,tt,ss}.trial,1));
        
        tr.a10trialkept{ll,tt,ss}=sort(a10randtr(1:numcondafinal(ll,tt,ss)));
        tr.tll40trialkept{ll,tt,ss}=sort(tll40randtr(1:numcondafinal(ll,tt,ss)));
        tr.alltrialkept{ll,tt,ss}=sort(allrandtr(1:numcondafinal(ll,tt,ss)));
        tr.nllatrialkept{ll,tt,ss}=sort(nllrandtr(1:numcondafinal(ll,tt,ss)));
        
        cfg=[];
        cfg.trials=tr.all40trialkept{ll,tt,ss};
        if size(tlock_aud{ll,tt,ss}.trial,1)==1 && cfg.trials==1
        else
          tlock_aud{ll+40,tt,ss}=ft_selectdata(cfg,tlock_aud{ll+40,tt,ss});
        end
        cfg=[];
        cfg.trials=tr.tlltrialkept{ll,tt,ss};
        if size(tlock_aud{ll,tt,ss}.trial,1)==1 && cfg.trials==1
        else
          tlock_tac{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{ll,tt,ss});
        end
        if 0 % too memory intensive to do here
          cfg=[];
          cfg.trials=tr.t10trialkept{ll,tt,ss};
          tlock_tac{ll+20,tt,ss}=ft_selectdata(cfg,tlock_tac{10,tt,ss}); % create +20 here as 10 will be different for every ll,tt,ss
          cfg=[];
          cfg.trials=tr.nllttrialkept{ll,tt,ss};
          tlock_nul{ll+50,tt,ss}=ft_selectdata(cfg,tlock_nul{10,tt,ss});
        end
        
        cfg=[];
        cfg.trials=tr.tll40trialkept{ll,tt,ss};
        if size(tlock_aud{ll,tt,ss}.trial,1)==1 && cfg.trials==1
        else
          tlock_tac{ll+40,tt,ss}=ft_selectdata(cfg,tlock_tac{ll+40,tt,ss});
        end
        cfg=[];
        cfg.trials=tr.alltrialkept{ll,tt,ss};
        if size(tlock_aud{ll,tt,ss}.trial,1)==1 && cfg.trials==1
        else
          tlock_aud{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{ll,tt,ss});
        end
        if 0 % too memory intensive to do here
          cfg=[];
          cfg.trials=tr.a10trialkept{ll,tt,ss};
          tlock_aud{ll+20,tt,ss}=ft_selectdata(cfg,tlock_aud{10,tt,ss}); % create +20 here as 10 will be different for every ll,tt,ss
          cfg=[];
          cfg.trials=tr.nllatrialkept{ll,tt,ss};
          tlock_nul{ll+60,tt,ss}=ft_selectdata(cfg,tlock_nul{10,tt,ss});
        end
        %     end
      end
    end
  end
  
else % if stageuse empty
  tlock_tac=[];
  tlock_aud=[];
  tlock_nul=[];
  numcond=0;
  eog_exclude=0;
end

% save('trialkept.mat','tr','numcond*')
save(['trialkept_tt' num2str(ttin) '_sleep' num2str(sleep) '.mat'],'tr','numcond*','*eog_exclude')





