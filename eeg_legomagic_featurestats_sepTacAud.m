% statistics of EEG awake data

eeg_legomagic_preamble

%%
usetr=1;
sleep=1;
if sleep
  iiuse=setdiff(iiBuse,3:7);
  trialkc=0;  % CHANGE ME: -1 (use all), 0 (only non-Kc), or 1 (only Kc)
  iter=11;
else
  iiuse=iiSuse;
  trialkc=-1; % always -1
  iter=27; % for main analysis
end


for ii=setdiff(iiBuse,[3:7])
% for ii=8
% for ii=iiBuse(7:end)
  clear featurestats*
  %   for sleep=[0 1]
  %   for sleep=[1]
  cd([edir sub{ii} ])
  clearvars -except ii sub edir ddir ii*use sleep featurestats* iter* usetr trialkc*
  %     [raw_tac, raw_aud, raw_nul, artflag]=eeg_legomagic_epoching2_featurefunction(ii,sleep);
  [raw_tac, raw_aud, raw_nul, artflag]=eeg_legomagic_epoching2(ii,sleep,0,0); % featfull=0, saveflag =0;
  
  % For memory reasons, do TacPlusAud separately from AudPlusTac
  %% Do TacPlusAud first
  tacaud=1;
  %   for tt=2 % refers to lim=[no-limit .01 .005 .002];
  tt=3;
  [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection2_wakeSleep_sepTacAud(ii,sleep,tt,tacaud,raw_tac, raw_aud, raw_nul, iter, usetr, trialkc);
  
  featurestats_tacMSpN=zeros(8,9,3,23,8);
  featurestats_tacPaud=zeros(8,9,3,23,8);
  
  if ii<8
    soalist=[3 4 5 6 7];
  else
    soalist=[1 3 4 5 6 7 9];
  end
  
  if sleep
    ssuse=[tr.stageuse+10 23];
  else
    ssuse=tr.stageuse+10;
  end
  
  for ss=ssuse
    for ll=[soalist soalist+20 soalist+40]
      featureKcHist_tac{ll,tt,ss}=[];
      featureDeltaHist_tac{ll,tt,ss}=[];
      featureSwHist_tac{ll,tt,ss}=[];
      featureKcHist_aud{ll,tt,ss}=[];
      featureDeltaHist_aud{ll,tt,ss}=[];
      featureSwHist_aud{ll,tt,ss}=[];
      featureKcHist_nul{ll,tt,ss}=[];
      featureDeltaHist_nul{ll,tt,ss}=[];
      featureSwHist_nul{ll,tt,ss}=[];
      %           try
      if ll<10
        if ss==23 % concatenate over N2 and N3 together
          if length(tr.t10trialkept{ll,tt,12})<2 && length(tr.t10trialkept{ll,tt,13})>1
            cfg=[];
            cfg.trials=tr.t10trialkept{ll,tt,13};
            tlock_tac{ll+20,tt,ss}=ft_selectdata(cfg,tlock_tac{10,tt,13}); % create +20 here as 10 will be different for every ll,tt
          elseif length(tr.t10trialkept{ll,tt,13})<2 && length(tr.t10trialkept{ll,tt,12})>1
            cfg=[];
            cfg.trials=tr.t10trialkept{ll,tt,12};
            tlock_tac{ll+20,tt,ss}=ft_selectdata(cfg,tlock_tac{10,tt,12}); % create +20 here as 10 will be different for every ll,tt
          elseif length(tr.t10trialkept{ll,tt,13})>1 && length(tr.t10trialkept{ll,tt,12})>1
            cfg=[];
            cfg.trials=tr.t10trialkept{ll,tt,12};
            tmp12=ft_selectdata(cfg,tlock_tac{10,tt,12}); % create +20 here as 10 will be different for every ll,tt
            cfg.trials=tr.t10trialkept{ll,tt,13};
            tmp13=ft_selectdata(cfg,tlock_tac{10,tt,13}); % create +20 here as 10 will be different for every ll,tt
            cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
            cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
            tlock_tac{ll+20,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
          else
            tlock_tac{ll+20,tt,ss}=[];
            continue
          end
          if ll==max(soalist)
            tlock_tac{10,tt,12}=[];
            tlock_tac{10,tt,13}=[];
          end
          
          %               cfg=[];
          %               cfg.trials=tr.a10trialkept{ll,tt,12};
          %               tmp12=ft_selectdata(cfg,tlock_aud{10,tt,12}); % create +20 here as 10 will be different for every ll,tt
          %               cfg.trials=tr.a10trialkept{ll,tt,13};
          %               tmp13=ft_selectdata(cfg,tlock_aud{10,tt,13}); % create +20 here as 10 will be different for every ll,tt
          %               cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
          %               cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
          %               tlock_aud{ll+20,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
          if ll==max(soalist)
            tlock_aud{10,tt,12}=[];
            tlock_aud{10,tt,13}=[];
          end
          
          if length(tr.nllttrialkept{ll,tt,12})<2 && length(tr.nllttrialkept{ll,tt,13})>1
            cfg=[];
            cfg.trials=tr.nllttrialkept{ll,tt,13};
            tlock_nul{ll+50,tt,ss}=ft_selectdata(cfg,tlock_nul{10,tt,13}); % create +20 here as 10 will be different for every ll,tt
            tr.nllttrialkept{ll,tt,ss}=[tr.nllttrialkept{ll,tt,13}];
          elseif length(tr.nllttrialkept{ll,tt,13})<2 && length(tr.nllttrialkept{ll,tt,12})>1
            cfg=[];
            cfg.trials=tr.nllttrialkept{ll,tt,12};
            tlock_nul{ll+50,tt,ss}=ft_selectdata(cfg,tlock_nul{10,tt,12}); % create +20 here as 10 will be different for every ll,tt
            tr.nllttrialkept{ll,tt,ss}=[tr.nllttrialkept{ll,tt,12}];
          elseif length(tr.nllttrialkept{ll,tt,13})>1 && length(tr.nllttrialkept{ll,tt,12})>1
            cfg=[];
            cfg.trials=tr.nllttrialkept{ll,tt,12};
            tmp12=ft_selectdata(cfg,tlock_nul{10,tt,12});
            cfg.trials=tr.nllttrialkept{ll,tt,13};
            tmp13=ft_selectdata(cfg,tlock_nul{10,tt,13});
            cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
            cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
            tlock_nul{ll+50,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
            tr.nllttrialkept{ll,tt,ss}=[tr.nllttrialkept{ll,tt,12} tr.nllttrialkept{ll,tt,13}];
          else
            tlock_nul{ll+50,tt,ss}=[];
            tr.nllttrialkept{ll,tt,ss}=0;
            continue
          end
          
          %               cfg=[];
          %               cfg.trials=tr.nllatrialkept{ll,tt,12};
          %               tmp12=ft_selectdata(cfg,tlock_nul{10,tt,12});
          %               cfg.trials=tr.nllatrialkept{ll,tt,13};
          %               tmp13=ft_selectdata(cfg,tlock_nul{10,tt,13});
          %               cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
          %               cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
          %               tlock_nul{ll+60,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
          if ll==max(soalist)
            tlock_nul{10,tt,12}=[];
            tlock_nul{10,tt,13}=[];
          end
          
        else
          if ~tr.t10trialkept{ll,tt,ss}
            numt_trials(ll,tt,ss)=0;
            continue
          else
            cfg=[];
            cfg.trials=tr.t10trialkept{ll,tt,ss};
            if length(cfg.trials)<2
              tlock_tac{ll+20,tt,ss}=[]; % create +20 here as 10 will be different for every ll,tt
            else
              tlock_tac{ll+20,tt,ss}=ft_selectdata(cfg,tlock_tac{10,tt,ss}); % create +20 here as 10 will be different for every ll,tt
            end
            if ll==max(soalist) && ss<12
              tlock_tac{10,tt,ss}=[];
            end
            
            %               cfg=[];
            %               cfg.trials=tr.a10trialkept{ll,tt,ss};
            %               tlock_aud{ll+20,tt,ss}=ft_selectdata(cfg,tlock_aud{10,tt,ss}); % create +20 here as 10 will be different for every ll,tt
            if ll==max(soalist) && ss<12
              tlock_aud{10,tt,ss}=[];
            end
            
            cfg=[];
            cfg.trials=tr.nllttrialkept{ll,tt,ss};
            if length(cfg.trials)<2
              tlock_nul{ll+50,tt,ss}=[];
            else
              tlock_nul{ll+50,tt,ss}=ft_selectdata(cfg,tlock_nul{10,tt,ss});
            end
            %               cfg=[];
            %               cfg.trials=tr.nllatrialkept{ll,tt,ss};
            %               tlock_nul{ll+60,tt,ss}=ft_selectdata(cfg,tlock_nul{10,tt,ss});
            if ll==max(soalist) && ss<12
              tlock_nul{10,tt,ss}=[];
            end
          end
          
        end
        
      end
      
      if ll<10
        if isempty(tlock_tac{ll+20,tt,ss})
          numt_trials(ll,tt,ss)=0;
        else
          numt_trials(ll,tt,ss)=size(tlock_tac{ll+20,tt,ss}.trial,1); % this should be the same for tacll20, audll40, tacll, nulll50
        end
        %             numa_trials(ll,tt,ss)=size(tlock_aud{ll+20,tt,ss}.trial,1); % this should be the same for audll20, tacll40, audll, nulll60
      end
      
      if ss==23 && ll<10
        if length(tr.tlltrialkept{ll,tt,12})>=2 && length(tr.tlltrialkept{ll,tt,13})<2
          tlock_tac{ll,tt,ss}=tlock_tac{ll,tt,12};
        elseif length(tr.tlltrialkept{ll,tt,12})<2 && length(tr.tlltrialkept{ll,tt,13})>=2
          tlock_tac{ll,tt,ss}=tlock_tac{ll,tt,13};
        elseif length(tr.tlltrialkept{ll,tt,12})>=2 && length(tr.tlltrialkept{ll,tt,13})>=2
          tmp12=tlock_tac{ll,tt,12};
          tmp13=tlock_tac{ll,tt,13};
          cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
          cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
          tlock_tac{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
        else
          tlock_tac{ll,tt,ss}=[];
          continue
        end
        tlock_tac{ll,tt,12}=[];
        tlock_tac{ll,tt,13}=[];
        tlock_aud{ll,tt,12}=[];
        tlock_aud{ll,tt,13}=[];
      end
      
      if ss==23 && ll>40
        if length(tr.all40trialkept{ll-40,tt,12})>=2 && length(tr.all40trialkept{ll-40,tt,13})<2
          tlock_aud{ll,tt,ss}=tlock_aud{ll,tt,12};
        elseif length(tr.all40trialkept{ll-40,tt,12})<2 && length(tr.all40trialkept{ll-40,tt,13})>=2
          tlock_aud{ll,tt,ss}=tlock_aud{ll,tt,13};
        elseif length(tr.all40trialkept{ll-40,tt,12})>=2 && length(tr.all40trialkept{ll-40,tt,13})>=2
          tmp12=tlock_aud{ll,tt,12};
          tmp13=tlock_aud{ll,tt,13};
          cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
          cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
          tlock_aud{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
        else
          tlock_aud{ll,tt,ss}=[];
          continue
        end
        
        tlock_aud{ll,tt,12}=[];
        tlock_aud{ll,tt,13}=[];
        tlock_tac{ll,tt,12}=[];
        tlock_tac{ll,tt,13}=[];
      end
      
      if ll>40
        %             if ~isempty(tlock_tac{ll,tt,ss})
        %               tlock_tac{ll,tt,ss}=trim_nans(tlock_tac{ll,tt,ss});
        %             end
        if ~isempty(tlock_aud{ll,tt,ss})
          tlock_aud{ll,tt,ss}=trim_nans(tlock_aud{ll,tt,ss});
        end
      end
      
      %up to here, same as for TFR
      
      %%
      
      % 40Hz lowpass filter (already 0.2Hz highpass on non-epoched data)
      % also do baseline correction here
      if ll<40
        if isfield(tlock_tac{ll,tt,ss},'trial') && size(tlock_tac{ll,tt,ss}.trial,1)
          featstruct_tac{ll,tt,ss}=rmfield(tlock_tac{ll,tt,ss},intersect(fieldnames(tlock_tac{ll,tt,ss}),{'trial' 'avg' 'var' 'dof'}));
          
          cfg=[];
          cfg.lpfilter='yes';
          cfg.lpfreq=40;
          cfg.demean='yes';
          cfg.baselinewindow=[-1.7 -0.6];
          disp('ft_preprocessing tac')
          tlock_tac{ll,tt,ss}=ft_preprocessing(cfg,tlock_tac{ll,tt,ss});
          featurestats_tac(:,ll,tt,ss,1)=mean(tlock_tac{ll,tt,ss}.trialinfo(:,12:19),1); % 5s prestim, mean of means
          featurestats_tac(:,ll,tt,ss,2)=mean(tlock_tac{ll,tt,ss}.trialinfo(:,22:29),1); % 5s poststim, mean of means
          featurestats_tac(:,ll,tt,ss,3)=mean(ceil(tlock_tac{ll,tt,ss}.trialinfo(:,12:19)),1); % 5s prestim, mean of binary
          featurestats_tac(:,ll,tt,ss,4)=mean(ceil(tlock_tac{ll,tt,ss}.trialinfo(:,22:29)),1); % 5s poststim, mean of binary
          featurestats_tac(:,ll,tt,ss,5)=mean(tlock_tac{ll,tt,ss}.trialinfo(:,32:39),1); % 1s prestim, mean of means
          featurestats_tac(:,ll,tt,ss,6)=mean(tlock_tac{ll,tt,ss}.trialinfo(:,42:49),1); % 1s poststim, mean of means
          featurestats_tac(:,ll,tt,ss,7)=mean(ceil(tlock_tac{ll,tt,ss}.trialinfo(:,32:39)),1); % 1s prestim, mean of binary
          featurestats_tac(:,ll,tt,ss,8)=mean(ceil(tlock_tac{ll,tt,ss}.trialinfo(:,42:49)),1); % 1s poststim, mean of binary
          
          %               featureKcPre_tac{ll,tt,ss}    =tlock_tac{ll,tt,ss}.trialinfo(:,50:59);
          %               featureKcPost_tac{ll,tt,ss}   =tlock_tac{ll,tt,ss}.trialinfo(:,60:69);
          %               featureDeltaPre_tac{ll,tt,ss} =tlock_tac{ll,tt,ss}.trialinfo(:,70:79);
          %               featureDeltaPost_tac{ll,tt,ss}=tlock_tac{ll,tt,ss}.trialinfo(:,80:89);
          %               featureSwPre_tac{ll,tt,ss}    =tlock_tac{ll,tt,ss}.trialinfo(:,90:99);
          %               featureSwPost_tac{ll,tt,ss}   =tlock_tac{ll,tt,ss}.trialinfo(:,100:109);
          
          for hh=[50:59 70:79]
            featureKcHist_tac{ll,tt,ss}=[featureKcHist_tac{ll,tt,ss}; tlock_tac{ll,tt,ss}.trialinfo(find(tlock_tac{ll,tt,ss}.trialinfo(:,hh)),hh+10)];
          end
          for hh=[90:99 110:119]
            featureDeltaHist_tac{ll,tt,ss}=[featureDeltaHist_tac{ll,tt,ss}; tlock_tac{ll,tt,ss}.trialinfo(find(tlock_tac{ll,tt,ss}.trialinfo(:,hh)),hh+10)];
          end
          for hh=[130:139 150:159]
            featureSwHist_tac{ll,tt,ss}=[featureSwHist_tac{ll,tt,ss}; tlock_tac{ll,tt,ss}.trialinfo(find(tlock_tac{ll,tt,ss}.trialinfo(:,hh)),hh+10)];
          end
          
          %               keyboard; % what about spindles? save that info out too?
          % Here, the 4th index means:
          % 1: prestim, vector over trials in which feature exists of amplitude of most recent preceding feature.
          % 2: poststim, vector over trials in which feature exists of amplitude of immediately following feature.
          featureSpFAmp_tac{ll,tt,ss,2}=tlock_tac{ll,tt,ss}.trialinfo(find(tlock_tac{ll,tt,ss}.trialinfo(:,190)),190);
          for ff=1:size(tlock_tac{ll,tt,ss}.trialinfo,1)
            for hh=179:-1:170 % start with most recent before stim onset
              if find(tlock_tac{ll,tt,ss}.trialinfo(ff,hh))
                featureSpFAmp_tac{ll,tt,ss,1}=[featureSpFAmp_tac{ll,tt,ss,1} tlock_tac{ll,tt,ss}.trialinfo(ff,hh)];
                break % don't save out amplitude of previous ones
              end
            end
          end
          featureSpSAmp_tac{ll,tt,ss,2}=tlock_tac{ll,tt,ss}.trialinfo(find(tlock_tac{ll,tt,ss}.trialinfo(:,230)),230);
          for ff=1:size(tlock_tac{ll,tt,ss}.trialinfo,1)
            for hh=219:-1:210 % start with most recent before stim onset
              if find(tlock_tac{ll,tt,ss}.trialinfo(ff,hh))
                featureSpSAmp_tac{ll,tt,ss,1}=[featureSpSAmp_tac{ll,tt,ss,1} tlock_tac{ll,tt,ss}.trialinfo(ff,hh)];
                break % don't save out amplitude of previous ones
              end
            end
          end
          featureSpFDur_tac{ll,tt,ss,2}=tlock_tac{ll,tt,ss}.trialinfo(find(tlock_tac{ll,tt,ss}.trialinfo(:,200)),200);
          for ff=1:size(tlock_tac{ll,tt,ss}.trialinfo,1)
            for hh=189:-1:180 % start with most recent before stim onset
              if find(tlock_tac{ll,tt,ss}.trialinfo(ff,hh))
                featureSpFDur_tac{ll,tt,ss,1}=[featureSpFDur_tac{ll,tt,ss,1} tlock_tac{ll,tt,ss}.trialinfo(ff,hh)];
                break % don't save out amplitude of previous ones
              end
            end
          end
          featureSpSDur_tac{ll,tt,ss,2}=tlock_tac{ll,tt,ss}.trialinfo(find(tlock_tac{ll,tt,ss}.trialinfo(:,240)),240);
          for ff=1:size(tlock_tac{ll,tt,ss}.trialinfo,1)
            for hh=229:-1:220 % start with most recent before stim onset
              if find(tlock_tac{ll,tt,ss}.trialinfo(ff,hh))
                featureSpSDur_tac{ll,tt,ss,1}=[featureSpSDur_tac{ll,tt,ss,1} tlock_tac{ll,tt,ss}.trialinfo(ff,hh)];
                break % don't save out amplitude of previous ones
              end
            end
          end
          
          % Here, the 4th index means:
          % 1: prestim, vector over trials in which feature exists of amplitude of most recent preceding feature.
          % 2: poststim, vector over trials in which feature exists of amplitude of immediately following feature.
          % 3: prestim, scalar of percentage trials with that feature (redundant with above?)
          % 4: poststim, scalar of percentage trials with that feature (redundant with above?)
          featureKcAmp_tac{ll,tt,ss,2}=tlock_tac{ll,tt,ss}.trialinfo(find(tlock_tac{ll,tt,ss}.trialinfo(:,70)),70);
          featureKcAmp_tac{ll,tt,ss,4}=length(find(tlock_tac{ll,tt,ss}.trialinfo(:,70)))/size(tlock_tac{ll,tt,ss}.trialinfo,1); % is this redundant with featurestats_tac(.....,4) ?
          for ff=1:size(tlock_tac{ll,tt,ss}.trialinfo,1)
            for hh=59:-1:50 % start with most recent before stim onset
              if find(tlock_tac{ll,tt,ss}.trialinfo(ff,hh))
                featureKcAmp_tac{ll,tt,ss,1}=[featureKcAmp_tac{ll,tt,ss,1} tlock_tac{ll,tt,ss}.trialinfo(ff,hh)];
                break % don't save out amplitude of previous ones
              end
            end
          end
          featureKcAmp_tac{ll,tt,ss,3}=length(featureKcAmp_tac{ll,tt,ss,1})/size(tlock_tac{ll,tt,ss}.trialinfo,1); % is this redundant with featurestats_tac(.....
          
          featureDeltaAmp_tac{ll,tt,ss,2}=tlock_tac{ll,tt,ss}.trialinfo(find(tlock_tac{ll,tt,ss}.trialinfo(:,110)),110);
          featureDeltaAmp_tac{ll,tt,ss,4}=length(find(tlock_tac{ll,tt,ss}.trialinfo(:,110)))/size(tlock_tac{ll,tt,ss}.trialinfo,1); % is this redundant with featurestats_tac(.....,4) ?
          for ff=1:size(tlock_tac{ll,tt,ss}.trialinfo,1)
            for hh=99:-1:90 % start with most recent before stim onset
              if find(tlock_tac{ll,tt,ss}.trialinfo(ff,hh))
                featureDeltaAmp_tac{ll,tt,ss,1}=[featureDeltaAmp_tac{ll,tt,ss,1} tlock_tac{ll,tt,ss}.trialinfo(ff,hh)];
                break % don't save out amplitude of previous ones
              end
            end
          end
          featureDeltaAmp_tac{ll,tt,ss,3}=length(featureDeltaAmp_tac{ll,tt,ss,1})/size(tlock_tac{ll,tt,ss}.trialinfo,1); % is this redundant with featurestats_tac(.....
          
          featureSwAmp_tac{ll,tt,ss,2}=tlock_tac{ll,tt,ss}.trialinfo(find(tlock_tac{ll,tt,ss}.trialinfo(:,150)),150);
          featureSwAmp_tac{ll,tt,ss,4}=length(find(tlock_tac{ll,tt,ss}.trialinfo(:,150)))/size(tlock_tac{ll,tt,ss}.trialinfo,1); % is this redundant with featurestats_tac(.....,4) ?
          for ff=1:size(tlock_tac{ll,tt,ss}.trialinfo,1)
            for hh=139:-1:130 % start with most recent before stim onset
              if find(tlock_tac{ll,tt,ss}.trialinfo(ff,hh))
                featureSwAmp_tac{ll,tt,ss,1}=[featureSwAmp_tac{ll,tt,ss,1} tlock_tac{ll,tt,ss}.trialinfo(ff,hh)];
                break % don't save out amplitude of previous ones
              end
            end
          end
          featureSwAmp_tac{ll,tt,ss,3}=length(featureSwAmp_tac{ll,tt,ss,1})/size(tlock_tac{ll,tt,ss}.trialinfo,1); % is this redundant with featurestats_tac(.....
          
          
          cfg=[];
          tlock_tactlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
        else
          tlock_tactlock{ll,tt,ss}=[];
        end
      else % we for sure do not need ll>40 for tac here; it shouldn't even exist actually
        tlock_tac{ll,tt,ss}=[]; % clearing to help save memory as running
      end
      if ss==12 || ss==13
      else
        tlock_tac{ll,tt,ss}=[]; % clearing to help save memory as running
      end
      
      if ll>40 % can we clear ll<40 earlier?
        if isfield(tlock_aud{ll,tt,ss},'trial') && size(tlock_aud{ll,tt,ss}.trial,1)
          featstruct_aud{ll,tt,ss}=rmfield(tlock_aud{ll,tt,ss},intersect(fieldnames(tlock_aud{ll,tt,ss}),{'trial' 'avg' 'var' 'dof'}));
          
          cfg=[];
          cfg.lpfilter='yes';
          cfg.lpfreq=40;
          %         cfg.bpfilter='yes';
          %         cfg.bpfreq=[1 40];
          cfg.demean='yes';
          cfg.baselinewindow=[-1.7 -0.6];
          disp('ft_preprocessing aud')
          tlock_aud{ll,tt,ss}=ft_preprocessing(cfg,tlock_aud{ll,tt,ss});
          %               featurestats_aud(:,ll,tt,ss,ii)=mean(tlock_aud{ll,tt,ss}.trialinfo(:,12:19));
          featurestats_aud(:,ll,tt,ss,1)=mean(tlock_aud{ll,tt,ss}.trialinfo(:,12:19)); % prestim, mean of means
          featurestats_aud(:,ll,tt,ss,2)=mean(tlock_aud{ll,tt,ss}.trialinfo(:,22:29)); % poststim, mean of means
          featurestats_aud(:,ll,tt,ss,3)=mean(ceil(tlock_aud{ll,tt,ss}.trialinfo(:,12:19))); % prestim, mean of binary
          featurestats_aud(:,ll,tt,ss,4)=mean(ceil(tlock_aud{ll,tt,ss}.trialinfo(:,22:29))); % poststim, mean of binary
          featurestats_aud(:,ll,tt,ss,5)=mean(tlock_aud{ll,tt,ss}.trialinfo(:,32:39)); % prestim, mean of means
          featurestats_aud(:,ll,tt,ss,6)=mean(tlock_aud{ll,tt,ss}.trialinfo(:,42:49)); % poststim, mean of means
          featurestats_aud(:,ll,tt,ss,7)=mean(ceil(tlock_aud{ll,tt,ss}.trialinfo(:,32:39))); % prestim, mean of binary
          featurestats_aud(:,ll,tt,ss,8)=mean(ceil(tlock_aud{ll,tt,ss}.trialinfo(:,42:49))); % poststim, mean of binary
          
          %               featureKcPre_aud{ll,tt,ss}    =tlock_aud{ll,tt,ss}.trialinfo(:,50:59);
          %               featureKcPost_aud{ll,tt,ss}   =tlock_aud{ll,tt,ss}.trialinfo(:,60:69);
          %               featureDeltaPre_aud{ll,tt,ss} =tlock_aud{ll,tt,ss}.trialinfo(:,70:79);
          %               featureDeltaPost_aud{ll,tt,ss}=tlock_aud{ll,tt,ss}.trialinfo(:,80:89);
          %               featureSwPre_aud{ll,tt,ss}    =tlock_aud{ll,tt,ss}.trialinfo(:,90:99);
          %               featureSwPost_aud{ll,tt,ss}   =tlock_aud{ll,tt,ss}.trialinfo(:,100:109);
          
          %               featureBWtime_aud(1:3,ll,tt,ss,:)=tlock_aud{ll,tt,ss}.trialinfo(:,34:36)';
          %               featureBWamp_aud(1:3,ll,tt,ss,:)=tlock_aud{ll,tt,ss}.trialinfo(:,44:46)';
          for hh=[50:59 70:79]
            featureKcHist_aud{ll,tt,ss}=[featureKcHist_aud{ll,tt,ss}; tlock_aud{ll,tt,ss}.trialinfo(find(tlock_aud{ll,tt,ss}.trialinfo(:,hh)),hh+10)];
          end
          for hh=[90:99 110:119]
            featureDeltaHist_aud{ll,tt,ss}=[featureDeltaHist_aud{ll,tt,ss}; tlock_aud{ll,tt,ss}.trialinfo(find(tlock_aud{ll,tt,ss}.trialinfo(:,hh)),hh+10)];
          end
          for hh=[130:139 150:159]
            featureSwHist_aud{ll,tt,ss}=[featureSwHist_aud{ll,tt,ss}; tlock_aud{ll,tt,ss}.trialinfo(find(tlock_aud{ll,tt,ss}.trialinfo(:,hh)),hh+10)];
          end
          
          featureKcAmp_aud{ll,tt,ss,2}=tlock_aud{ll,tt,ss}.trialinfo(find(tlock_aud{ll,tt,ss}.trialinfo(:,70)),70);
          featureKcAmp_aud{ll,tt,ss,4}=length(find(tlock_aud{ll,tt,ss}.trialinfo(:,70)))/size(tlock_aud{ll,tt,ss}.trialinfo,1); % is this redundant with featurestats_aud(.....,4) ?
          for ff=1:size(tlock_aud{ll,tt,ss}.trialinfo,1)
            for hh=59:-1:50 % start with most recent before stim onset
              if find(tlock_aud{ll,tt,ss}.trialinfo(ff,hh))
                featureKcAmp_aud{ll,tt,ss,1}=[featureKcAmp_aud{ll,tt,ss,1} tlock_aud{ll,tt,ss}.trialinfo(ff,hh)];
                break % don't save out amplitude of previous ones
              end
            end
          end
          featureKcAmp_aud{ll,tt,ss,3}=length(featureKcAmp_aud{ll,tt,ss,1})/size(tlock_aud{ll,tt,ss}.trialinfo,1); % is this redundant with featurestats_aud(.....
          
          featureDeltaAmp_aud{ll,tt,ss,2}=tlock_aud{ll,tt,ss}.trialinfo(find(tlock_aud{ll,tt,ss}.trialinfo(:,110)),110);
          featureDeltaAmp_aud{ll,tt,ss,4}=length(find(tlock_aud{ll,tt,ss}.trialinfo(:,110)))/size(tlock_aud{ll,tt,ss}.trialinfo,1); % is this redundant with featurestats_aud(.....,4) ?
          for ff=1:size(tlock_aud{ll,tt,ss}.trialinfo,1)
            for hh=99:-1:90 % start with most recent before stim onset
              if find(tlock_aud{ll,tt,ss}.trialinfo(ff,hh))
                featureDeltaAmp_aud{ll,tt,ss,1}=[featureDeltaAmp_aud{ll,tt,ss,1} tlock_aud{ll,tt,ss}.trialinfo(ff,hh)];
                break % don't save out amplitude of previous ones
              end
            end
          end
          featureDeltaAmp_aud{ll,tt,ss,3}=length(featureDeltaAmp_aud{ll,tt,ss,1})/size(tlock_aud{ll,tt,ss}.trialinfo,1); % is this redundant with featurestats_aud(.....
          
          featureSwAmp_aud{ll,tt,ss,2}=tlock_aud{ll,tt,ss}.trialinfo(find(tlock_aud{ll,tt,ss}.trialinfo(:,150)),150);
          featureSwAmp_aud{ll,tt,ss,4}=length(find(tlock_aud{ll,tt,ss}.trialinfo(:,150)))/size(tlock_aud{ll,tt,ss}.trialinfo,1); % is this redundant with featurestats_aud(.....,4) ?
          for ff=1:size(tlock_aud{ll,tt,ss}.trialinfo,1)
            for hh=139:-1:130 % start with most recent before stim onset
              if find(tlock_aud{ll,tt,ss}.trialinfo(ff,hh))
                featureSwAmp_aud{ll,tt,ss,1}=[featureSwAmp_aud{ll,tt,ss,1} tlock_aud{ll,tt,ss}.trialinfo(ff,hh)];
                break % don't save out amplitude of previous ones
              end
            end
          end
          featureSwAmp_aud{ll,tt,ss,3}=length(featureSwAmp_aud{ll,tt,ss,1})/size(tlock_aud{ll,tt,ss}.trialinfo,1); % is this redundant with featurestats_aud(.....
          
          % spindles
          featureSpFAmp_aud{ll,tt,ss,2}=tlock_aud{ll,tt,ss}.trialinfo(find(tlock_aud{ll,tt,ss}.trialinfo(:,190)),190);
          for ff=1:size(tlock_aud{ll,tt,ss}.trialinfo,1)
            for hh=179:-1:170 % start with most recent before stim onset
              if find(tlock_aud{ll,tt,ss}.trialinfo(ff,hh))
                featureSpFAmp_aud{ll,tt,ss,1}=[featureSpFAmp_aud{ll,tt,ss,1} tlock_aud{ll,tt,ss}.trialinfo(ff,hh)];
                break % don't save out amplitude of previous ones
              end
            end
          end
          featureSpSAmp_aud{ll,tt,ss,2}=tlock_aud{ll,tt,ss}.trialinfo(find(tlock_aud{ll,tt,ss}.trialinfo(:,230)),230);
          for ff=1:size(tlock_aud{ll,tt,ss}.trialinfo,1)
            for hh=219:-1:210 % start with most recent before stim onset
              if find(tlock_aud{ll,tt,ss}.trialinfo(ff,hh))
                featureSpSAmp_aud{ll,tt,ss,1}=[featureSpSAmp_aud{ll,tt,ss,1} tlock_aud{ll,tt,ss}.trialinfo(ff,hh)];
                break % don't save out amplitude of previous ones
              end
            end
          end
          featureSpFDur_aud{ll,tt,ss,2}=tlock_aud{ll,tt,ss}.trialinfo(find(tlock_aud{ll,tt,ss}.trialinfo(:,200)),200);
          for ff=1:size(tlock_aud{ll,tt,ss}.trialinfo,1)
            for hh=189:-1:180 % start with most recent before stim onset
              if find(tlock_aud{ll,tt,ss}.trialinfo(ff,hh))
                featureSpFDur_aud{ll,tt,ss,1}=[featureSpFDur_aud{ll,tt,ss,1} tlock_aud{ll,tt,ss}.trialinfo(ff,hh)];
                break % don't save out amplitude of previous ones
              end
            end
          end
          featureSpSDur_aud{ll,tt,ss,2}=tlock_aud{ll,tt,ss}.trialinfo(find(tlock_aud{ll,tt,ss}.trialinfo(:,240)),240);
          for ff=1:size(tlock_aud{ll,tt,ss}.trialinfo,1)
            for hh=229:-1:220 % start with most recent before stim onset
              if find(tlock_aud{ll,tt,ss}.trialinfo(ff,hh))
                featureSpSDur_aud{ll,tt,ss,1}=[featureSpSDur_aud{ll,tt,ss,1} tlock_aud{ll,tt,ss}.trialinfo(ff,hh)];
                break % don't save out amplitude of previous ones
              end
            end
          end
          
          cfg=[];
          tlock_audtlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_aud{ll,tt,ss});
        else
          tlock_audtlock{ll,tt,ss}=[];
        end
        %           elseif ll<40 && ss~=12 && ss~=13
        %             tlock_aud{ll,tt,ss}=[];
      end
      if ss==12 || ss==13
      else
        tlock_aud{ll,tt,ss}=[];
      end
      %           catch
      %             disp('didnt work for this ll')
      %           end
    end % end ll
    
    
    %         for ll=[soalist+50 soalist+60]
    for ll=[soalist+50 ]
      featureKcHist_nul{ll,tt,ss}=[];
      featureDeltaHist_nul{ll,tt,ss}=[];
      featureSwHist_nul{ll,tt,ss}=[];
      if length(tr.nllttrialkept{ll-50,tt,ss})>=2 && isfield(tlock_nul{ll,tt,ss},'trial') && size(tlock_nul{ll,tt,ss}.trial,1)
        featstruct_nul{ll,tt,ss}=rmfield(tlock_nul{ll,tt,ss},intersect(fieldnames(tlock_nul{ll,tt,ss}),{'trial' 'avg' 'var' 'dof'}));
        
        cfg=[];
        cfg.lpfilter='yes';
        cfg.lpfreq=40;
        cfg.demean='yes';
        cfg.baselinewindow=[-1.7 -0.6];
        %           cfg.bpfilter='yes';
        %           cfg.bpfreq=[1 40];
        disp('ft_preprocessing nul')
        tlock_nul{ll,tt,ss}=ft_preprocessing(cfg,tlock_nul{ll,tt,ss});
        featurestats_nul(:,ll,tt,ss,1)=mean(tlock_nul{ll,tt,ss}.trialinfo(:,12:19)); % prestim, mean of means
        featurestats_nul(:,ll,tt,ss,2)=mean(tlock_nul{ll,tt,ss}.trialinfo(:,22:29)); % poststim, mean of means
        featurestats_nul(:,ll,tt,ss,3)=mean(ceil(tlock_nul{ll,tt,ss}.trialinfo(:,12:19))); % prestim, mean of binary
        featurestats_nul(:,ll,tt,ss,4)=mean(ceil(tlock_nul{ll,tt,ss}.trialinfo(:,22:29))); % poststim, mean of binary
        featurestats_nul(:,ll,tt,ss,5)=mean(tlock_nul{ll,tt,ss}.trialinfo(:,32:39)); % prestim, mean of means
        featurestats_nul(:,ll,tt,ss,6)=mean(tlock_nul{ll,tt,ss}.trialinfo(:,42:49)); % poststim, mean of means
        featurestats_nul(:,ll,tt,ss,7)=mean(ceil(tlock_nul{ll,tt,ss}.trialinfo(:,32:39))); % prestim, mean of binary
        featurestats_nul(:,ll,tt,ss,8)=mean(ceil(tlock_nul{ll,tt,ss}.trialinfo(:,42:49))); % poststim, mean of binary
        
        %               featureKcPre_nul{ll,tt,ss}    =tlock_nul{ll,tt,ss}.trialinfo(:,50:59);
        %               featureKcPost_nul{ll,tt,ss}   =tlock_nul{ll,tt,ss}.trialinfo(:,60:69);
        %               featureDeltaPre_nul{ll,tt,ss} =tlock_nul{ll,tt,ss}.trialinfo(:,70:79);
        %               featureDeltaPost_nul{ll,tt,ss}=tlock_nul{ll,tt,ss}.trialinfo(:,80:89);
        %               featureSwPre_nul{ll,tt,ss}    =tlock_nul{ll,tt,ss}.trialinfo(:,90:99);
        %               featureSwPost_nul{ll,tt,ss}   =tlock_nul{ll,tt,ss}.trialinfo(:,100:109);
        
        %               %             featurestats_nul(:,ll,tt,ss)=mean(tlock_nul{ll,tt,ss}.trialinfo(:,12:19));
        %               featureBWtime_nul(1:3,ll,tt,ss,:)=tlock_nul{ll,tt,ss}.trialinfo(:,34:36)';
        %               featureBWamp_nul(1:3,ll,tt,ss,:)=tlock_nul{ll,tt,ss}.trialinfo(:,44:46)';
        for hh=[50:59 70:79]
          featureKcHist_nul{ll,tt,ss}=[featureKcHist_nul{ll,tt,ss}; tlock_nul{ll,tt,ss}.trialinfo(find(tlock_nul{ll,tt,ss}.trialinfo(:,hh)),hh+10)];
        end
        for hh=[90:99 110:119]
          featureDeltaHist_nul{ll,tt,ss}=[featureDeltaHist_nul{ll,tt,ss}; tlock_nul{ll,tt,ss}.trialinfo(find(tlock_nul{ll,tt,ss}.trialinfo(:,hh)),hh+10)];
        end
        for hh=[130:139 150:159]
          featureSwHist_nul{ll,tt,ss}=[featureSwHist_nul{ll,tt,ss}; tlock_nul{ll,tt,ss}.trialinfo(find(tlock_nul{ll,tt,ss}.trialinfo(:,hh)),hh+10)];
        end
        
        featureKcAmp_nul{ll,tt,ss,2}=tlock_nul{ll,tt,ss}.trialinfo(find(tlock_nul{ll,tt,ss}.trialinfo(:,70)),70);
        featureKcAmp_nul{ll,tt,ss,4}=length(find(tlock_nul{ll,tt,ss}.trialinfo(:,70)))/size(tlock_nul{ll,tt,ss}.trialinfo,1); % is this redundant with featurestats_nul(.....,4) ?
        for ff=1:size(tlock_nul{ll,tt,ss}.trialinfo,1)
          for hh=59:-1:50 % start with most recent before stim onset
            if find(tlock_nul{ll,tt,ss}.trialinfo(ff,hh))
              featureKcAmp_nul{ll,tt,ss,1}=[featureKcAmp_nul{ll,tt,ss,1} tlock_nul{ll,tt,ss}.trialinfo(ff,hh)];
              break % don't save out amplitude of previous ones
            end
          end
        end
        featureKcAmp_nul{ll,tt,ss,3}=length(featureKcAmp_nul{ll,tt,ss,1})/size(tlock_nul{ll,tt,ss}.trialinfo,1); % is this redundant with featurestats_nul(.....
        
        featureDeltaAmp_nul{ll,tt,ss,2}=tlock_nul{ll,tt,ss}.trialinfo(find(tlock_nul{ll,tt,ss}.trialinfo(:,110)),110);
        featureDeltaAmp_nul{ll,tt,ss,4}=length(find(tlock_nul{ll,tt,ss}.trialinfo(:,110)))/size(tlock_nul{ll,tt,ss}.trialinfo,1); % is this redundant with featurestats_nul(.....,4) ?
        for ff=1:size(tlock_nul{ll,tt,ss}.trialinfo,1)
          for hh=99:-1:90 % start with most recent before stim onset
            if find(tlock_nul{ll,tt,ss}.trialinfo(ff,hh))
              featureDeltaAmp_nul{ll,tt,ss,1}=[featureDeltaAmp_nul{ll,tt,ss,1} tlock_nul{ll,tt,ss}.trialinfo(ff,hh)];
              break % don't save out amplitude of previous ones
            end
          end
        end
        featureDeltaAmp_nul{ll,tt,ss,3}=length(featureDeltaAmp_nul{ll,tt,ss,1})/size(tlock_nul{ll,tt,ss}.trialinfo,1); % is this redundant with featurestats_nul(.....
        
        featureSwAmp_nul{ll,tt,ss,2}=tlock_nul{ll,tt,ss}.trialinfo(find(tlock_nul{ll,tt,ss}.trialinfo(:,150)),150);
        featureSwAmp_nul{ll,tt,ss,4}=length(find(tlock_nul{ll,tt,ss}.trialinfo(:,150)))/size(tlock_nul{ll,tt,ss}.trialinfo,1); % is this redundant with featurestats_nul(.....,4) ?
        for ff=1:size(tlock_nul{ll,tt,ss}.trialinfo,1)
          for hh=139:-1:130 % start with most recent before stim onset
            if find(tlock_nul{ll,tt,ss}.trialinfo(ff,hh))
              featureSwAmp_nul{ll,tt,ss,1}=[featureSwAmp_nul{ll,tt,ss,1} tlock_nul{ll,tt,ss}.trialinfo(ff,hh)];
              break % don't save out amplitude of previous ones
            end
          end
        end
        featureSwAmp_nul{ll,tt,ss,3}=length(featureSwAmp_nul{ll,tt,ss,1})/size(tlock_nul{ll,tt,ss}.trialinfo,1); % is this redundant with featurestats_nul(.....
        
        
        % spindles
        featureSpFAmp_nul{ll,tt,ss,2}=tlock_nul{ll,tt,ss}.trialinfo(find(tlock_nul{ll,tt,ss}.trialinfo(:,190)),190);
        for ff=1:size(tlock_nul{ll,tt,ss}.trialinfo,1)
          for hh=179:-1:170 % start with most recent before stim onset
            if find(tlock_nul{ll,tt,ss}.trialinfo(ff,hh))
              featureSpFAmp_nul{ll,tt,ss,1}=[featureSpFAmp_nul{ll,tt,ss,1} tlock_nul{ll,tt,ss}.trialinfo(ff,hh)];
              break % don't save out amplitude of previous ones
            end
          end
        end
        featureSpSAmp_nul{ll,tt,ss,2}=tlock_nul{ll,tt,ss}.trialinfo(find(tlock_nul{ll,tt,ss}.trialinfo(:,230)),230);
        for ff=1:size(tlock_nul{ll,tt,ss}.trialinfo,1)
          for hh=219:-1:210 % start with most recent before stim onset
            if find(tlock_nul{ll,tt,ss}.trialinfo(ff,hh))
              featureSpSAmp_nul{ll,tt,ss,1}=[featureSpSAmp_nul{ll,tt,ss,1} tlock_nul{ll,tt,ss}.trialinfo(ff,hh)];
              break % don't save out amplitude of previous ones
            end
          end
        end
        featureSpFDur_nul{ll,tt,ss,2}=tlock_nul{ll,tt,ss}.trialinfo(find(tlock_nul{ll,tt,ss}.trialinfo(:,200)),200);
        for ff=1:size(tlock_nul{ll,tt,ss}.trialinfo,1)
          for hh=189:-1:180 % start with most recent before stim onset
            if find(tlock_nul{ll,tt,ss}.trialinfo(ff,hh))
              featureSpFDur_nul{ll,tt,ss,1}=[featureSpFDur_nul{ll,tt,ss,1} tlock_nul{ll,tt,ss}.trialinfo(ff,hh)];
              break % don't save out amplitude of previous ones
            end
          end
        end
        featureSpSDur_nul{ll,tt,ss,2}=tlock_nul{ll,tt,ss}.trialinfo(find(tlock_nul{ll,tt,ss}.trialinfo(:,240)),240);
        for ff=1:size(tlock_nul{ll,tt,ss}.trialinfo,1)
          for hh=229:-1:220 % start with most recent before stim onset
            if find(tlock_nul{ll,tt,ss}.trialinfo(ff,hh))
              featureSpSDur_nul{ll,tt,ss,1}=[featureSpSDur_nul{ll,tt,ss,1} tlock_nul{ll,tt,ss}.trialinfo(ff,hh)];
              break % don't save out amplitude of previous ones
            end
          end
        end
        
        
        cfg=[];
        tlock_nultlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_nul{ll,tt,ss});
      else
        tlock_nultlock{ll,tt,ss}=[];
      end
      tlock_nul{ll,tt,ss}=[];
    end
    
    
    for ll=soalist
      % create sum of unisensory conditions
      cfg=[];
      cfg.operation='add';
      cfg.parameter='avg';
      % the call to ft_timelockanalysis is b/c .avg field not present after ft_redefinetrial
      if numt_trials(ll,tt,ss)>=2
        tlock_tacPaud{ll,tt,ss}=ft_math(cfg,tlock_tactlock{ll+20,tt,ss},tlock_audtlock{ll+40,tt,ss}); % keeping tac centred but jitter aud (to match tlock_tac{X,tt})
        if ss<14
          featstruct_tacPaud{ll,tt,ss}=add_sleep_fields(featstruct_tac{ll+20,tt,ss},featstruct_aud{ll+40,tt,ss});
        end
        
        featurestats_tacPaud(:,ll,tt,ss,:)=mean([featurestats_tac(:,ll+20,tt,ss,:), featurestats_aud(:,ll+40,tt,ss,:)],2);
        featureKcHist_tacPaud{ll,tt,ss}=[featureKcHist_tac{ll+20,tt,ss}; featureKcHist_aud{ll+40,tt,ss}];
        featureDeltaHist_tacPaud{ll,tt,ss}=[featureDeltaHist_tac{ll+20,tt,ss}; featureDeltaHist_aud{ll+40,tt,ss}];
        featureSwHist_tacPaud{ll,tt,ss}=[featureSwHist_tac{ll+20,tt,ss}; featureSwHist_aud{ll+40,tt,ss}];
        featureKcAmp_tacPaud{ll,tt,ss,1}=[featureKcAmp_tac{ll+20,tt,ss,1}, featureKcAmp_aud{ll+40,tt,ss,1}]';
        featureKcAmp_tacPaud{ll,tt,ss,2}=[featureKcAmp_tac{ll+20,tt,ss,2}; featureKcAmp_aud{ll+40,tt,ss,2}];
        featureDeltaAmp_tacPaud{ll,tt,ss,1}=[featureDeltaAmp_tac{ll+20,tt,ss,1}, featureDeltaAmp_aud{ll+40,tt,ss,1}]';
        featureDeltaAmp_tacPaud{ll,tt,ss,2}=[featureDeltaAmp_tac{ll+20,tt,ss,2}; featureDeltaAmp_aud{ll+40,tt,ss,2}];
        featureSwAmp_tacPaud{ll,tt,ss,1}=[featureSwAmp_tac{ll+20,tt,ss,1}, featureSwAmp_aud{ll+40,tt,ss,1}]';
        featureSwAmp_tacPaud{ll,tt,ss,2}=[featureSwAmp_tac{ll+20,tt,ss,2}; featureSwAmp_aud{ll+40,tt,ss,2}];
        
        featureSpFAmp_tacPaud{ll,tt,ss,1}=[featureSpFAmp_tac{ll+20,tt,ss,1}, featureSpFAmp_aud{ll+40,tt,ss,1}]';
        featureSpFAmp_tacPaud{ll,tt,ss,2}=[featureSpFAmp_tac{ll+20,tt,ss,2}; featureSpFAmp_aud{ll+40,tt,ss,2}];
        featureSpSAmp_tacPaud{ll,tt,ss,1}=[featureSpSAmp_tac{ll+20,tt,ss,1}, featureSpSAmp_aud{ll+40,tt,ss,1}]';
        featureSpSAmp_tacPaud{ll,tt,ss,2}=[featureSpSAmp_tac{ll+20,tt,ss,2}; featureSpSAmp_aud{ll+40,tt,ss,2}];
        featureSpFDur_tacPaud{ll,tt,ss,1}=[featureSpFDur_tac{ll+20,tt,ss,1}, featureSpFDur_aud{ll+40,tt,ss,1}]';
        featureSpFDur_tacPaud{ll,tt,ss,2}=[featureSpFDur_tac{ll+20,tt,ss,2}; featureSpFDur_aud{ll+40,tt,ss,2}];
        featureSpSDur_tacPaud{ll,tt,ss,1}=[featureSpSDur_tac{ll+20,tt,ss,1}, featureSpSDur_aud{ll+40,tt,ss,1}]';
        featureSpSDur_tacPaud{ll,tt,ss,2}=[featureSpSDur_tac{ll+20,tt,ss,2}; featureSpSDur_aud{ll+40,tt,ss,2}];
        
        tlock_tTacAlone{ll,tt,ss}=tlock_tactlock{ll+20,tt,ss}; % keeping tac centred but jitter aud (to match tlock_tac{X,tt})
        tlock_tAudAlone{ll,tt,ss}=tlock_audtlock{ll+40,tt,ss}; % keeping tac centred but jitter aud (to match tlock_tac{X,tt})
      else
        %             featurestats_tacPaud(:,ll,tt,ss,:)=[];
        tlock_tacPaud{ll,tt,ss}=[];
        tlock_tTacAlone{ll,tt,ss}=[]; % keeping tac centred but jitter aud (to match tlock_tac{X,tt})
        tlock_tAudAlone{ll,tt,ss}=[]; % keeping tac centred but jitter aud (to match tlock_tac{X,tt})
        featstruct_tacPaud{ll,tt,ss}=[];
        
        featureSpFAmp_tacPaud{ll,tt,ss,1}=[]';
        featureSpFAmp_tacPaud{ll,tt,ss,2}=[];
        featureSpSAmp_tacPaud{ll,tt,ss,1}=[]';
        featureSpSAmp_tacPaud{ll,tt,ss,2}=[];
        featureSpFDur_tacPaud{ll,tt,ss,1}=[]';
        featureSpFDur_tacPaud{ll,tt,ss,2}=[];
        featureSpSDur_tacPaud{ll,tt,ss,1}=[]';
        featureSpSDur_tacPaud{ll,tt,ss,2}=[];
        
        featureKcHist_tacPaud{ll,tt,ss}=[];
        featureDeltaHist_tacPaud{ll,tt,ss}=[];
        featureSwHist_tacPaud{ll,tt,ss}=[];
        featureKcAmp_tacPaud{ll,tt,ss,1}=[];
        featureKcAmp_tacPaud{ll,tt,ss,2}=[];
        featureDeltaAmp_tacPaud{ll,tt,ss,1}=[];
        featureDeltaAmp_tacPaud{ll,tt,ss,2}=[];
        featureSwAmp_tacPaud{ll,tt,ss,1}=[];
        featureSwAmp_tacPaud{ll,tt,ss,2}=[];
      end
      %           if numa_trials(ll,tt,ss)
      %             tlock_audPtac{ll,tt,ss}=ft_math(cfg,tlock_audtlock{ll+20,tt,ss},tlock_tactlock{ll+40,tt,ss}); % keeping aud centred but jitter tac (to match tlock_aud{X,tt})
      %             featurestats_audPtac(:,ll,tt,ss,ii)=mean([featurestats_aud(:,ll+20,tt,ss,ii), featurestats_tac(:,ll+40,tt,ss,ii)]');
      %           else
      %             tlock_audPtac{ll,tt,ss}=[];
      %           end
      
      % create sum of multisensory plus nul conditions
      cfg=[];
      cfg.operation='add';
      cfg.parameter='avg';
      if numt_trials(ll,tt,ss)>=2
        tlock_tacMSpN{ll,tt,ss}=ft_math(cfg,tlock_tactlock{ll,tt,ss},tlock_nultlock{ll+50,tt,ss}); % TA+N
        if ss<14
          featstruct_tacMSpN{ll,tt,ss}=add_sleep_fields(featstruct_tac{ll,tt,ss},featstruct_nul{ll+50,tt,ss});
        end
        
        featurestats_tacMSpN(:,ll,tt,ss,:)=mean([featurestats_tac(:,ll,tt,ss,:), featurestats_nul(:,ll+50,tt,ss,:)],2);
        featureKcHist_tacMSpN{ll,tt,ss}=[featureKcHist_tac{ll,tt,ss}; featureKcHist_nul{ll+50,tt,ss}];
        featureDeltaHist_tacMSpN{ll,tt,ss}=[featureDeltaHist_tac{ll,tt,ss}; featureDeltaHist_nul{ll+50,tt,ss}];
        featureSwHist_tacMSpN{ll,tt,ss}=[featureSwHist_tac{ll,tt,ss}; featureSwHist_nul{ll+50,tt,ss}];
        featureKcAmp_tacMSpN{ll,tt,ss,1}=[featureKcAmp_tac{ll,tt,ss,1}, featureKcAmp_nul{ll+50,tt,ss,1}]';
        featureKcAmp_tacMSpN{ll,tt,ss,2}=[featureKcAmp_tac{ll,tt,ss,2}; featureKcAmp_nul{ll+50,tt,ss,2}];
        featureDeltaAmp_tacMSpN{ll,tt,ss,1}=[featureDeltaAmp_tac{ll,tt,ss,1}, featureDeltaAmp_nul{ll+50,tt,ss,1}]';
        featureDeltaAmp_tacMSpN{ll,tt,ss,2}=[featureDeltaAmp_tac{ll,tt,ss,2}; featureDeltaAmp_nul{ll+50,tt,ss,2}];
        featureSwAmp_tacMSpN{ll,tt,ss,1}=[featureSwAmp_tac{ll,tt,ss,1}, featureSwAmp_nul{ll+50,tt,ss,1}]';
        featureSwAmp_tacMSpN{ll,tt,ss,2}=[featureSwAmp_tac{ll,tt,ss,2}; featureSwAmp_nul{ll+50,tt,ss,2}];
        
        featureSpFAmp_tacMSpN{ll,tt,ss,1}=[featureSpFAmp_tac{ll,tt,ss,1}, featureSpFAmp_nul{ll+50,tt,ss,1}]';
        featureSpFAmp_tacMSpN{ll,tt,ss,2}=[featureSpFAmp_tac{ll,tt,ss,2}; featureSpFAmp_nul{ll+50,tt,ss,2}];
        featureSpSAmp_tacMSpN{ll,tt,ss,1}=[featureSpSAmp_tac{ll,tt,ss,1}, featureSpSAmp_nul{ll+50,tt,ss,1}]';
        featureSpSAmp_tacMSpN{ll,tt,ss,2}=[featureSpSAmp_tac{ll,tt,ss,2}; featureSpSAmp_nul{ll+50,tt,ss,2}];
        featureSpFDur_tacMSpN{ll,tt,ss,1}=[featureSpFDur_tac{ll,tt,ss,1}, featureSpFDur_nul{ll+50,tt,ss,1}]';
        featureSpFDur_tacMSpN{ll,tt,ss,2}=[featureSpFDur_tac{ll,tt,ss,2}; featureSpFDur_nul{ll+50,tt,ss,2}];
        featureSpSDur_tacMSpN{ll,tt,ss,1}=[featureSpSDur_tac{ll,tt,ss,1}, featureSpSDur_nul{ll+50,tt,ss,1}]';
        featureSpSDur_tacMSpN{ll,tt,ss,2}=[featureSpSDur_tac{ll,tt,ss,2}; featureSpSDur_nul{ll+50,tt,ss,2}];
        
        tlock_tMSAlone{ll,tt,ss}=tlock_tactlock{ll,tt,ss}; % TA+N
        tlock_tNulAlone{ll,tt,ss}=tlock_nultlock{ll+50,tt,ss}; % TA+N
      else
        tlock_tacMSpN{ll,tt,ss}=[];
        featstruct_tacMSpN{ll,tt,ss}=[];
        %             featurestats_tacMSpN(:,ll,tt,ss,:)=[];
        tlock_tMSAlone{ll,tt,ss}=[]; % TA+N
        tlock_tNulAlone{ll,tt,ss}=[]; % TA+N
        featureKcHist_tacMSpN{ll,tt,ss}=[];
        featureDeltaHist_tacMSpN{ll,tt,ss}=[];
        featureSwHist_tacMSpN{ll,tt,ss}=[];
        featureKcAmp_tacMSpN{ll,tt,ss,1}=[];
        featureKcAmp_tacMSpN{ll,tt,ss,2}=[];
        featureDeltaAmp_tacMSpN{ll,tt,ss,1}=[];
        featureDeltaAmp_tacMSpN{ll,tt,ss,2}=[];
        featureSwAmp_tacMSpN{ll,tt,ss,1}=[];
        featureSwAmp_tacMSpN{ll,tt,ss,2}=[];
        
        featureSpFAmp_tacMSpN{ll,tt,ss,1}=[]';
        featureSpFAmp_tacMSpN{ll,tt,ss,2}=[];
        featureSpSAmp_tacMSpN{ll,tt,ss,1}=[]';
        featureSpSAmp_tacMSpN{ll,tt,ss,2}=[];
        featureSpFDur_tacMSpN{ll,tt,ss,1}=[]';
        featureSpFDur_tacMSpN{ll,tt,ss,2}=[];
        featureSpSDur_tacMSpN{ll,tt,ss,1}=[]';
        featureSpSDur_tacMSpN{ll,tt,ss,2}=[];
      end
      %           if numa_trials(ll,tt,ss)
      %             tlock_audMSpN{ll,tt,ss}=ft_math(cfg,tlock_audtlock{ll,tt,ss},tlock_nultlock{ll+60,tt,ss}); % AT+N
      %             featurestats_audMSpN(:,ll,tt,ss,ii)=mean([featurestats_aud(:,ll,tt,ss,ii), featurestats_nul(:,ll+60,tt,ss,ii)]');
      %           else
      %             tlock_audMSpN{ll,tt,ss}=[];
      %           end
    end % end ll
    
    
    
  end  % end ss
  
  clear tlock_tac tlock_aud tlock_nul tr tlock*tlock
  %   end % end tt
  
  
  %% Do AudPlusTac second
  tacaud=0;
  if 0
    for tt=2:3 % refers to lim=[no-limit .01 .005 .002];
      %   profile on
      [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection2_wakeSleep_sepTacAud(ii,sleep,tt,tacaud);
      %     profile viewer
      
      % % would love to save out here but really can't; it would be over 13GB for
      % % sitting and 20GB for bed, *per* tt
      %     save(['tlock_trialSel_' sub{ii} '.mat'],'tlock*','tr','-v7.3');
      
      
      %   end
      %   if ~exist('tlock_tac_s0','var')
      %     [tlock_tac_s0,tlock_aud_s0,tlock_nul_s0,tr]=eeg_legomagic_trialSelection1(ii);
      %     save(['tlock_trialSel_' sub{ii} '.mat'],'tlock*s0','tr','-v7.3');
      %   end
      %   if ~exist('tlock_aud_s0','var')
      %     [tlock_tac_s0,tlock_aud_s0,tlock_nul_s0,tr]=eeg_legomagic_trialSelection1(ii);
      %     save(['tlock_trialSel_' sub{ii} '.mat'],'tlock*s0','tr','-v7.3');
      %   end
      
      
      if ii<8
        soalist=[3 4 5 6 7];
      else
        soalist=[1 3 4 5 6 7 9];
      end
      
      if sleep
        ssuse=[tr.stageuse+10 23];
      else
        ssuse=tr.stageuse+10;
      end
      
      for ss=ssuse
        %   for tt=1:4
        for ll=[soalist soalist+20 soalist+40]
          %           try
          if ll<10
            if ss==23 % concatenate over N2 and N3 together
              %               cfg=[];
              %               cfg.trials=tr.t10trialkept{ll,tt,12};
              %               tmp12=ft_selectdata(cfg,tlock_tac{10,tt,12}); % create +20 here as 10 will be different for every ll,tt
              %               cfg.trials=tr.t10trialkept{ll,tt,13};
              %               tmp13=ft_selectdata(cfg,tlock_tac{10,tt,13}); % create +20 here as 10 will be different for every ll,tt
              %               cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
              %               cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
              %               tlock_tac{ll+20,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
              if ll==max(soalist)
                tlock_tac{10,tt,12}=[];
                tlock_tac{10,tt,13}=[];
              end
              
              if length(tr.a10trialkept{ll,tt,12})<2 && length(tr.a10trialkept{ll,tt,13})>1
                cfg=[];
                cfg.trials=tr.a10trialkept{ll,tt,13};
                tlock_aud{ll+20,tt,ss}=ft_selectdata(cfg,tlock_aud{10,tt,13}); % create +20 here as 10 will be different for every ll,tt
              elseif length(tr.a10trialkept{ll,tt,13})<2 && length(tr.a10trialkept{ll,tt,12})>1
                cfg=[];
                cfg.trials=tr.a10trialkept{ll,tt,12};
                tlock_aud{ll+20,tt,ss}=ft_selectdata(cfg,tlock_aud{10,tt,12}); % create +20 here as 10 will be different for every ll,tt
              elseif length(tr.a10trialkept{ll,tt,13})>1 && length(tr.a10trialkept{ll,tt,12})>1
                cfg=[];
                cfg.trials=tr.a10trialkept{ll,tt,12};
                tmp12=ft_selectdata(cfg,tlock_aud{10,tt,12}); % create +20 here as 10 will be different for every ll,tt
                cfg.trials=tr.a10trialkept{ll,tt,13};
                tmp13=ft_selectdata(cfg,tlock_aud{10,tt,13}); % create +20 here as 10 will be different for every ll,tt
                cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
                cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
                tlock_aud{ll+20,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
              else
                tlock_aud{ll+20,tt,ss}=[];
              end
              if ll==max(soalist)
                tlock_aud{10,tt,12}=[];
                tlock_aud{10,tt,13}=[];
              end
              
              %               cfg=[];
              %               cfg.trials=tr.nllttrialkept{ll,tt,12};
              %               tmp12=ft_selectdata(cfg,tlock_nul{10,tt,12});
              %               cfg.trials=tr.nllttrialkept{ll,tt,13};
              %               tmp13=ft_selectdata(cfg,tlock_nul{10,tt,13});
              %               cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
              %               cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
              %               tlock_nul{ll+50,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
              
              if length(tr.nllatrialkept{ll,tt,12})<2 && length(tr.nllatrialkept{ll,tt,13})>1
                cfg=[];
                cfg.trials=tr.nllatrialkept{ll,tt,13};
                tlock_nul{ll+60,tt,ss}=ft_selectdata(cfg,tlock_nul{10,tt,13}); % create +20 here as 10 will be different for every ll,tt
                tr.nllatrialkept{ll,tt,ss}=[tr.nllatrialkept{ll,tt,13}];
              elseif length(tr.nllatrialkept{ll,tt,13})<2 && length(tr.nllatrialkept{ll,tt,12})>1
                cfg=[];
                cfg.trials=tr.nllatrialkept{ll,tt,12};
                tlock_nul{ll+60,tt,ss}=ft_selectdata(cfg,tlock_nul{10,tt,12}); % create +20 here as 10 will be different for every ll,tt
                tr.nllatrialkept{ll,tt,ss}=[tr.nllatrialkept{ll,tt,12}];
              elseif length(tr.a10trialkept{ll,tt,13})>1 && length(tr.a10trialkept{ll,tt,12})>1
                cfg=[];
                cfg.trials=tr.nllatrialkept{ll,tt,12};
                tmp12=ft_selectdata(cfg,tlock_nul{10,tt,12});
                cfg.trials=tr.nllatrialkept{ll,tt,13};
                tmp13=ft_selectdata(cfg,tlock_nul{10,tt,13});
                cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
                cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
                tlock_nul{ll+60,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
                tr.nllatrialkept{ll,tt,ss}=[tr.nllatrialkept{ll,tt,12} tr.nllatrialkept{ll,tt,13}];
              else
                tlock_nul{ll+60,tt,ss}=[];
                tr.nllatrialkept{ll,tt,ss}=0;
              end
              
              if ll==max(soalist)
                tlock_nul{10,tt,12}=[];
                tlock_nul{10,tt,13}=[];
              end
              
            else
              if ~tr.a10trialkept{ll,tt,ss}
                numa_trials(ll,tt,ss)=0;
                continue
              else
                %               cfg=[];
                %               cfg.trials=tr.t10trialkept{ll,tt,ss};
                %               tlock_tac{ll+20,tt,ss}=ft_selectdata(cfg,tlock_tac{10,tt,ss}); % create +20 here as 10 will be different for every ll,tt
                if ll==max(soalist) && ss<12
                  tlock_tac{10,tt,ss}=[];
                end
                
                cfg=[];
                cfg.trials=tr.a10trialkept{ll,tt,ss};
                if length(cfg.trials)<2
                  tlock_aud{ll+20,tt,ss}=[]; % create +20 here as 10 will be different for every ll,tt
                else
                  tlock_aud{ll+20,tt,ss}=ft_selectdata(cfg,tlock_aud{10,tt,ss}); % create +20 here as 10 will be different for every ll,tt
                end
                if ll==max(soalist) && ss<12
                  tlock_aud{10,tt,ss}=[];
                end
                
                %               cfg=[];
                %               cfg.trials=tr.nllttrialkept{ll,tt,ss};
                %               tlock_nul{ll+50,tt,ss}=ft_selectdata(cfg,tlock_nul{10,tt,ss});
                cfg=[];
                cfg.trials=tr.nllatrialkept{ll,tt,ss};
                if length(cfg.trials)<2
                  tlock_nul{ll+60,tt,ss}=[];
                else
                  tlock_nul{ll+60,tt,ss}=ft_selectdata(cfg,tlock_nul{10,tt,ss});
                end
                if ll==max(soalist) && ss<12
                  tlock_nul{10,tt,ss}=[];
                end
              end
              
            end
          end
          
          if ll<10
            %             numt_trials(ll,tt,ss)=size(tlock_tac{ll+20,tt,ss}.trial,1); % this should be the same for tacll20, audll40, tacll, nulll50
            if isempty(tlock_aud{ll+20,tt,ss})
              numa_trials(ll,tt,ss)=0;
            else
              numa_trials(ll,tt,ss)=size(tlock_aud{ll+20,tt,ss}.trial,1); % this should be the same for audll20, tacll40, audll, nulll60
            end
          end
          
          if ss==23 && ll<10
            %             tmp12=tlock_tac{ll,tt,12};
            %             tmp13=tlock_tac{ll,tt,13};
            %             cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
            %             cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
            %             tlock_tac{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
            tlock_tac{ll,tt,12}=[];
            tlock_tac{ll,tt,13}=[];
            
            if length(tr.alltrialkept{ll,tt,12})>=2 && length(tr.alltrialkept{ll,tt,13})<2
              tlock_aud{ll,tt,ss}=tlock_aud{ll,tt,12};
            elseif length(tr.alltrialkept{ll,tt,12})<2 && length(tr.alltrialkept{ll,tt,13})>=2
              tlock_aud{ll,tt,ss}=tlock_aud{ll,tt,13};
            elseif length(tr.alltrialkept{ll,tt,12})>=2 && length(tr.alltrialkept{ll,tt,13})>=2
              tmp12=tlock_aud{ll,tt,12};
              tmp13=tlock_aud{ll,tt,13};
              cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
              cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
              tlock_aud{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
            else
              tlock_aud{ll,tt,ss}=[];
              continue
            end
            tlock_aud{ll,tt,12}=[];
            tlock_aud{ll,tt,13}=[];
          end
          if ss==23 && ll>40
            if length(tr.tll40trialkept{ll-40,tt,12})>=2 && length(tr.tll40trialkept{ll-40,tt,13})<2
              tlock_tac{ll,tt,ss}=tlock_tac{ll,tt,12};
            elseif length(tr.tll40trialkept{ll-40,tt,12})<2 && length(tr.tll40trialkept{ll-40,tt,13})>=2
              tlock_tac{ll,tt,ss}=tlock_tac{ll,tt,13};
            elseif length(tr.tll40trialkept{ll-40,tt,12})>=2 && length(tr.tll40trialkept{ll-40,tt,13})>=2
              tmp12=tlock_tac{ll,tt,12};
              tmp13=tlock_tac{ll,tt,13};
              cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
              cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
              tlock_tac{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
            else
              tlock_tac{ll,tt,ss}=[];
              continue
            end
            tlock_tac{ll,tt,12}=[];
            tlock_tac{ll,tt,13}=[];
            
            %             tmp12=tlock_aud{ll,tt,12};
            %             tmp13=tlock_aud{ll,tt,13};
            %             cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
            %             cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
            %             tlock_aud{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
            tlock_aud{ll,tt,12}=[];
            tlock_aud{ll,tt,13}=[];
          end
          
          if ll>40
            if ~isempty(tlock_tac{ll,tt,ss})
              tlock_tac{ll,tt,ss}=trim_nans(tlock_tac{ll,tt,ss});
            end
            %             if ~isempty(tlock_aud{ll,tt,ss})
            %               tlock_aud{ll,tt,ss}=trim_nans(tlock_aud{ll,tt,ss});
            %             end
          end
          
          % up to here same as for TFR
          
          %%
          % 40Hz lowpass filter (already 0.2Hz highpass on non-epoched data)
          % also do baseline correction here
          if ll>40
            if isfield(tlock_tac{ll,tt,ss},'trial') && size(tlock_tac{ll,tt,ss}.trial,1)
              featstruct_tac{ll,tt,ss}=rmfield(tlock_tac{ll,tt,ss},intersect(fieldnames(tlock_tac{ll,tt,ss}),{'trial' 'avg' 'var' 'dof'}));
              
              cfg=[];
              cfg.lpfilter='yes';
              cfg.lpfreq=40;
              cfg.demean='yes';
              cfg.baselinewindow=[-1.7 -0.6];
              disp('ft_preprocessing tac')
              tlock_tac{ll,tt,ss}=ft_preprocessing(cfg,tlock_tac{ll,tt,ss});
              featurestats_tac(:,ll,tt,ss)=mean(tlock_tac{ll,tt,ss}.trialinfo(:,12:19));
              cfg=[];
              tlock_tactlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
            else
              tlock_tactlock{ll,tt,ss}=[];
            end
          end
          if ss==12 || ss==13
          else
            tlock_tac{ll,tt,ss}=[]; % clearing to help save memory as running
          end
          
          if ll<40
            if isfield(tlock_aud{ll,tt,ss},'trial') && size(tlock_aud{ll,tt,ss}.trial,1)
              featstruct_tac{ll,tt,ss}=rmfield(tlock_tac{ll,tt,ss},intersect(fieldnames(tlock_tac{ll,tt,ss}),{'trial' 'avg' 'var' 'dof'}));
              cfg=[];
              cfg.lpfilter='yes';
              cfg.lpfreq=40;
              %         cfg.bpfilter='yes';
              %         cfg.bpfreq=[1 40];
              cfg.demean='yes';
              cfg.baselinewindow=[-1.7 -0.6];
              disp('ft_preprocessing aud')
              tlock_aud{ll,tt,ss}=ft_preprocessing(cfg,tlock_aud{ll,tt,ss});
              featurestats_aud(:,ll,tt,ss)=mean(tlock_aud{ll,tt,ss}.trialinfo(:,12:19));
              cfg=[];
              tlock_audtlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_aud{ll,tt,ss});
            else
              tlock_audtlock{ll,tt,ss}=[];
            end
          end
          if ss==12 || ss==13
          else
            tlock_aud{ll,tt,ss}=[];
          end
          %           catch
          %             disp('didnt work for this ll')
          %           end
        end % end ll
        
        
        %         for ll=[soalist+50 soalist+60]
        for ll=[soalist+60]
          if length(tr.nllatrialkept{ll-60,tt,ss})>=2 && isfield(tlock_nul{ll,tt,ss},'trial') && size(tlock_nul{ll,tt,ss}.trial,1)
            featstruct_tac{ll,tt,ss}=rmfield(tlock_tac{ll,tt,ss},intersect(fieldnames(tlock_tac{ll,tt,ss}),{'trial' 'avg' 'var' 'dof'}));
            cfg=[];
            cfg.lpfilter='yes';
            cfg.lpfreq=40;
            cfg.demean='yes';
            cfg.baselinewindow=[-1.7 -0.6];
            %           cfg.bpfilter='yes';
            %           cfg.bpfreq=[1 40];
            disp('ft_preprocessing nul')
            tlock_nul{ll,tt,ss}=ft_preprocessing(cfg,tlock_nul{ll,tt,ss});
            featurestats_nul(:,ll,tt,ss)=mean(tlock_nul{ll,tt,ss}.trialinfo(:,12:19));
            cfg=[];
            tlock_nultlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_nul{ll,tt,ss});
          else
            tlock_nultlock{ll,tt,ss}=[];
          end
          tlock_nul{ll,tt,ss}=[];
        end
        
        
        for ll=soalist
          keyboard
          % create sum of unisensory conditions
          cfg=[];
          cfg.operation='add';
          cfg.parameter='avg';
          % the call to ft_timelockanalysis is b/c .avg field not present after ft_redefinetrial
          %           if numt_trials(ll,tt,ss)
          %             tlock_tacPaud{ll,tt,ss}=ft_math(cfg,tlock_tactlock{ll+20,tt,ss},tlock_audtlock{ll+40,tt,ss}); % keeping tac centred but jitter aud (to match tlock_tac{X,tt})
          %             featurestats_tacPaud(:,ll,tt,ss,ii)=mean([featurestats_tac(:,ll+20,tt,ss,ii), featurestats_aud(:,ll+40,tt,ss,ii)]');
          %           else
          %             tlock_tacPaud{ll,tt,ss}=[];
          %           end
          if numa_trials(ll,tt,ss)>=2
            tlock_audPtac{ll,tt,ss}=ft_math(cfg,tlock_audtlock{ll+20,tt,ss},tlock_tactlock{ll+40,tt,ss}); % keeping aud centred but jitter tac (to match tlock_aud{X,tt})
            featurestats_audPtac(:,ll,tt,ss)=mean([featurestats_aud(:,ll+20,tt,ss,ii), featurestats_tac(:,ll+40,tt,ss,ii)]');
            tlock_aAudAlone{ll,tt,ss}=tlock_audtlock{ll+20,tt,ss}; % keeping aud centred but jitter tac (to match tlock_aud{X,tt})
            tlock_aTacAlone{ll,tt,ss}=tlock_tactlock{ll+40,tt,ss}; % keeping aud centred but jitter tac (to match tlock_aud{X,tt})
          else
            tlock_audPtac{ll,tt,ss}=[];
            tlock_aAudAlone{ll,tt,ss}=[]; % keeping aud centred but jitter tac (to match tlock_aud{X,tt})
            tlock_aTacAlone{ll,tt,ss}=[]; % keeping aud centred but jitter tac (to match tlock_aud{X,tt})
          end
          
          % create sum of multisensory plus nul conditions
          cfg=[];
          cfg.operation='add';
          cfg.parameter='avg';
          %           if numt_trials(ll,tt,ss)
          %             tlock_tacMSpN{ll,tt,ss}=ft_math(cfg,tlock_tactlock{ll,tt,ss},tlock_nultlock{ll+50,tt,ss}); % TA+N
          %             featurestats_tacMSpN(:,ll,tt,ss,ii)=mean([featurestats_tac(:,ll,tt,ss,ii), featurestats_nul(:,ll+50,tt,ss,ii)]');
          %           else
          %             tlock_tacMSpN{ll,tt,ss}=[];
          %           end
          if numa_trials(ll,tt,ss)>=2
            tlock_audMSpN{ll,tt,ss}=ft_math(cfg,tlock_audtlock{ll,tt,ss},tlock_nultlock{ll+60,tt,ss}); % AT+N
            featurestats_audMSpN(:,ll,tt,ss)=mean([featurestats_aud(:,ll,tt,ss,ii), featurestats_nul(:,ll+60,tt,ss,ii)]');
            tlock_aMSAlone{ll,tt,ss}=tlock_audtlock{ll,tt,ss}; % AT+N
            tlock_aNulAlone{ll,tt,ss}=tlock_nultlock{ll+60,tt,ss}; % AT+N
          else
            tlock_audMSpN{ll,tt,ss}=[];
            tlock_aMSAlone{ll,tt,ss}=[]; % AT+N
            tlock_aNulAlone{ll,tt,ss}=[]; % AT+N
          end
        end % end ll
        
        tlock_tacAll{tt,ss}=tlock_tac{100,tt,ss};
        tlock_audAll{tt,ss}=tlock_aud{100,tt,ss};
        tlock_tac19T{tt,ss}=tlock_tac{101,tt,ss};
        tlock_aud19A{tt,ss}=tlock_aud{101,tt,ss};
        
      end  % end ss
      
      clear tlock_tac tlock_aud tlock_nul tr tlock*tlock
    end % end tt
  end
  
  
  %%
  
  
  %     save(['tlock_diffs_features_' sub{ii} '_sleep' num2str(sleep) '.mat'],'tlock_*P*','tlock_*N*','tlock_t*','tlock_a*','num*trials','featurestats_*')
  save(['tlock_diffs_features_' sub{ii} '_sleep' num2str(sleep) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'num*trials','feat*')
  
  
  %   end % sleep
end

return

