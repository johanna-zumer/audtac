% statistics of EEG awake data

clear all
close all
if ispc
  edir='D:\audtac\eeg_data\';
  ddir='D:\audtac\legomagic\diaries\';
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

%%
for ii=setdiff(iiBuse,[3:7])
%   for ii=8
  % for ii=setdiff(6:32,[iiuse 19]) % for awake-alone
  clear featurestats*
  %   for sleep=[0 1]
  for sleep=[1]
    cd([edir sub{ii} ])
    clearvars -except ii sub edir ddir ii*use sleep featurestats*
    %     [raw_tac, raw_aud, raw_nul, artflag]=eeg_legomagic_epoching2_featurefunction(ii,sleep);
    [raw_tac, raw_aud, raw_nul, artflag]=eeg_legomagic_epoching2(ii,sleep,0,0); % featfull=0, saveflag =0;
    
    % For memory reasons, do TacPlusAud separately from AudPlusTac
    %% Do TacPlusAud first
    tacaud=1;
    featurestats_tacMSpN=zeros(8,9,3,23,8);
    featurestats_tacPaud=zeros(8,9,3,23,8);
    for tt=2 % refers to lim=[no-limit .01 .005 .002];
      [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection2_wakeSleep_sepTacAud(ii,sleep,tt,tacaud,raw_tac, raw_aud, raw_nul);
      
      
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
    end % end tt
    
    
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
    save(['tlock_diffs_features_' sub{ii} '_sleep' num2str(sleep) '.mat'],'num*trials','feature*')
    
    
  end
end

return



%%  Stats on features

printflag=0;
plotflag=0;
dostats=0;
audtacflag=0;
histbins=[-4:.05:4];

soalist=[1 3 4 5 6 7 9];
ttuse=2;

kcts=nan(9,2,23,3,2,19);
kctt=nan(9,2,23,3,2,19);
kcta=nan(9,2,23,3,2,19);
kctn=nan(9,2,23,3,2,19);

for sleep=1 % for Kc/spindles, only makes sense to look at 'bed' data
  
  if sleep
    subuse=setdiff(iiBuse,[3:7]);
  else
    subuse=iiSuse;
  end
  
  submin=subuse(1)-1;
  subuseind=0;
  for ii=subuse
    clearvars -except sub* edir ddir sdir ii* soalist *flag dostats fstats* histbins k* sleep ttuse
    %       for ii=setdiff(subuse,[8 9 10 12 14 15 16 17 18])
    cd([edir sub{ii} ])
    %         load(['tlock_diffs_features_' sub{ii} '_sleep' num2str(sleep) '.mat'])
    load(['tlock_diffs_features_' sub{ii} '_sleep' num2str(sleep) '.mat']);
    if exist('numa_trials','var')
      audtacflag=1;
    end
    
%     tk0=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(0) '.mat']);
%     tk1=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '.mat']);
    
    %         % THis is preliminary...how best to included all stages later on?
    %         if sleep==0
    %           ss=10; % awake
    %         elseif sleep==1
    %           ss=23; % this is concatenation of N2 and N3
    %         end
    if sleep==1
      ssuse=[10 11 12 13 23];
    else
      ssuse=[10 11];
    end
    
    for tt=ttuse
      for ll=soalist
        for ss=ssuse
          numtrt(ll,tt,ss,ii-submin)=numt_trials(ll,tt,ss); % does this need to be per 'sleep01' as well?
          if audtacflag
            numtra(ll,tt,ss,ii-submin)=numa_trials(ll,tt,ss);
          end
        end
      end
    end
    
    if [sleep==1 && any(numtrt(soalist,ttuse,12,ii-submin)<20)] || [sleep==0 && any(numtrt(soalist,ttuse,10,ii-submin)<20)]
      subuse=setdiff(subuse,ii);
      continue
    else
      subuseind=subuseind+1;
    end
    
    for tt=ttuse
      for ll=soalist
        for ss=ssuse
          %         if sleep==0 %load for both at sleep0 then will still be in memory for sleep1
          %           numtrt(ll,tt,10,ii-submin)=numt_trials(ll,tt,10); % does this need to be per 'sleep01' as well?
          %           numtra(ll,tt,10,ii-submin)=numa_trials(ll,tt,10);
          %           load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(1) '.mat'],'num*trials')
          %           load(['tlock_diffs_features_' sub{ii} '_sleep' num2str(1) '.mat']);
          %           numtrt(ll,tt,23,ii-submin)=numt_trials(ll,tt,23); % does this need to be per 'sleep01' as well?
          %           numtra(ll,tt,23,ii-submin)=numa_trials(ll,tt,23);
          %         end
          
          % discard from both if either sleep/wake doesn't have 'enough' (what is enough? 20?)
          %         if numtrt(ll,tt,10,ii-submin)<20 || numtrt(ll,tt,23,ii-submin)<20
          %           tlock_tacPaud_each{subuseind}=tlock_tacPaud{ll,tt,ss};
          %           tlock_audPtac_each{subuseind}=tlock_audPtac{ll,tt,ss};
          %           tlock_tacMSpN_each{subuseind}=tlock_tacMSpN{ll,tt,ss};
          %           tlock_audMSpN_each{subuseind}=tlock_audMSpN{ll,tt,ss};
          %           if sleep==1
          
          % Histogram of time of negmax
          khtm(ll,tt,ss,:,1,subuseind)=histcjz(featureDeltaHist_tacMSpN{ll,tt,ss},histbins); % multisens plus nul
          khtm(ll,tt,ss,:,2,subuseind)=histcjz(featureKcHist_tacMSpN{ll,tt,ss},histbins);
          khtm(ll,tt,ss,:,3,subuseind)=histcjz(featureSwHist_tacMSpN{ll,tt,ss},histbins);
          khtu(ll,tt,ss,:,1,subuseind)=histcjz(featureDeltaHist_tacPaud{ll,tt,ss},histbins); % unisensory sum
          khtu(ll,tt,ss,:,2,subuseind)=histcjz(featureKcHist_tacPaud{ll,tt,ss},histbins);
          khtu(ll,tt,ss,:,3,subuseind)=histcjz(featureSwHist_tacPaud{ll,tt,ss},histbins);
          khts(ll,tt,ss,:,1,subuseind)=histcjz(featureDeltaHist_tac{ll,tt,ss},histbins); % simult mulitsens (not plus nul)
          khts(ll,tt,ss,:,2,subuseind)=histcjz(featureKcHist_tac{ll,tt,ss},histbins);
          khts(ll,tt,ss,:,3,subuseind)=histcjz(featureSwHist_tac{ll,tt,ss},histbins);
          khtt(ll,tt,ss,:,1,subuseind)=histcjz(featureDeltaHist_tac{ll+20,tt,ss},histbins); % tactile alone
          khtt(ll,tt,ss,:,2,subuseind)=histcjz(featureKcHist_tac{ll+20,tt,ss},histbins);
          khtt(ll,tt,ss,:,3,subuseind)=histcjz(featureSwHist_tac{ll+20,tt,ss},histbins);
          khta(ll,tt,ss,:,1,subuseind)=histcjz(featureDeltaHist_aud{ll+40,tt,ss},histbins); % auditory alone
          khta(ll,tt,ss,:,2,subuseind)=histcjz(featureKcHist_aud{ll+40,tt,ss},histbins);
          khta(ll,tt,ss,:,3,subuseind)=histcjz(featureSwHist_aud{ll+40,tt,ss},histbins);
          khtn(ll,tt,ss,:,1,subuseind)=histcjz(featureDeltaHist_nul{ll+50,tt,ss},histbins); % nul alone
          khtn(ll,tt,ss,:,2,subuseind)=histcjz(featureKcHist_nul{ll+50,tt,ss},histbins);
          khtn(ll,tt,ss,:,3,subuseind)=histcjz(featureSwHist_nul{ll+50,tt,ss},histbins);
          
          kahtm(ll,tt,ss,:,subuseind)=histcjz([featureDeltaHist_tacMSpN{ll,tt,ss}; featureKcHist_tacMSpN{ll,tt,ss}; featureSwHist_tacMSpN{ll,tt,ss}],histbins);
          kahtu(ll,tt,ss,:,subuseind)=histcjz([featureDeltaHist_tacPaud{ll,tt,ss}; featureKcHist_tacPaud{ll,tt,ss}; featureSwHist_tacPaud{ll,tt,ss}],histbins);
          kahts(ll,tt,ss,:,subuseind)=histcjz([featureDeltaHist_tac{ll,tt,ss}; featureKcHist_tac{ll,tt,ss}; featureSwHist_tac{ll,tt,ss}],histbins);
          kahtt(ll,tt,ss,:,subuseind)=histcjz([featureDeltaHist_tac{ll+20,tt,ss}; featureKcHist_tac{ll+20,tt,ss}; featureSwHist_tac{ll+20,tt,ss}],histbins);
          kahta(ll,tt,ss,:,subuseind)=histcjz([featureDeltaHist_aud{ll+40,tt,ss}; featureKcHist_aud{ll+40,tt,ss}; featureSwHist_aud{ll+40,tt,ss}],histbins);
          kahtn(ll,tt,ss,:,subuseind)=histcjz([featureDeltaHist_nul{ll+50,tt,ss}; featureKcHist_nul{ll+50,tt,ss}; featureSwHist_nul{ll+50,tt,ss}],histbins);
          
          
          % Amplitude
          kaats(ll,tt,ss,1,1,subuseind)=mean([featureDeltaAmp_tac{ll,tt,ss,1}, featureKcAmp_tac{ll,tt,ss,1}, featureSwAmp_tac{ll,tt,ss,1}]); % pre
          kaats(ll,tt,ss,1,2,subuseind)=mean([featureDeltaAmp_tac{ll,tt,ss,2}', featureKcAmp_tac{ll,tt,ss,2}', featureSwAmp_tac{ll,tt,ss,2}']); % post
          kaatt(ll,tt,ss,1,1,subuseind)=mean([featureDeltaAmp_tac{ll+20,tt,ss,1}, featureKcAmp_tac{ll+20,tt,ss,1}, featureSwAmp_tac{ll+20,tt,ss,1}]); % pre
          kaatt(ll,tt,ss,1,2,subuseind)=mean([featureDeltaAmp_tac{ll+20,tt,ss,2}', featureKcAmp_tac{ll+20,tt,ss,2}', featureSwAmp_tac{ll+20,tt,ss,2}']); % post
          kaata(ll,tt,ss,1,1,subuseind)=mean([featureDeltaAmp_aud{ll+40,tt,ss,1}, featureKcAmp_aud{ll+40,tt,ss,1}, featureSwAmp_aud{ll+40,tt,ss,1}]); % pre
          kaata(ll,tt,ss,1,2,subuseind)=mean([featureDeltaAmp_aud{ll+40,tt,ss,2}', featureKcAmp_aud{ll+40,tt,ss,2}', featureSwAmp_aud{ll+40,tt,ss,2}']); % post
          kaatn(ll,tt,ss,1,1,subuseind)=mean([featureDeltaAmp_nul{ll+50,tt,ss,1}, featureKcAmp_nul{ll+50,tt,ss,1}, featureSwAmp_nul{ll+50,tt,ss,1}]); % pre
          kaatn(ll,tt,ss,1,2,subuseind)=mean([featureDeltaAmp_nul{ll+50,tt,ss,2}', featureKcAmp_nul{ll+50,tt,ss,2}', featureSwAmp_nul{ll+50,tt,ss,2}']); % post
          kaatm(ll,tt,ss,1,1,subuseind)=mean([featureDeltaAmp_tacMSpN{ll,tt,ss,1}', featureKcAmp_tacMSpN{ll,tt,ss,1}', featureSwAmp_tacMSpN{ll,tt,ss,1}']); % pre
          kaatm(ll,tt,ss,1,2,subuseind)=mean([featureDeltaAmp_tacMSpN{ll,tt,ss,2}', featureKcAmp_tacMSpN{ll,tt,ss,2}', featureSwAmp_tacMSpN{ll,tt,ss,2}']); % post
          kaatu(ll,tt,ss,1,1,subuseind)=mean([featureDeltaAmp_tacPaud{ll,tt,ss,1}', featureKcAmp_tacPaud{ll,tt,ss,1}', featureSwAmp_tacPaud{ll,tt,ss,1}']); % pre
          kaatu(ll,tt,ss,1,2,subuseind)=mean([featureDeltaAmp_tacPaud{ll,tt,ss,2}', featureKcAmp_tacPaud{ll,tt,ss,2}', featureSwAmp_tacPaud{ll,tt,ss,2}']); % post
          kaats(ll,tt,ss,2,1,subuseind)=length([featureDeltaAmp_tac{ll,tt,ss,1}, featureKcAmp_tac{ll,tt,ss,1}, featureSwAmp_tac{ll,tt,ss,1}]); % pre
          kaats(ll,tt,ss,2,2,subuseind)=length([featureDeltaAmp_tac{ll,tt,ss,2}', featureKcAmp_tac{ll,tt,ss,2}', featureSwAmp_tac{ll,tt,ss,2}']); % post
          kaatt(ll,tt,ss,2,1,subuseind)=length([featureDeltaAmp_tac{ll+20,tt,ss,1}, featureKcAmp_tac{ll+20,tt,ss,1}, featureSwAmp_tac{ll+20,tt,ss,1}]); % pre
          kaatt(ll,tt,ss,2,2,subuseind)=length([featureDeltaAmp_tac{ll+20,tt,ss,2}', featureKcAmp_tac{ll+20,tt,ss,2}', featureSwAmp_tac{ll+20,tt,ss,2}']); % post
          kaata(ll,tt,ss,2,1,subuseind)=length([featureDeltaAmp_aud{ll+40,tt,ss,1}, featureKcAmp_aud{ll+40,tt,ss,1}, featureSwAmp_aud{ll+40,tt,ss,1}]); % pre
          kaata(ll,tt,ss,2,2,subuseind)=length([featureDeltaAmp_aud{ll+40,tt,ss,2}', featureKcAmp_aud{ll+40,tt,ss,2}', featureSwAmp_aud{ll+40,tt,ss,2}']); % post
          kaatn(ll,tt,ss,2,1,subuseind)=length([featureDeltaAmp_nul{ll+50,tt,ss,1}, featureKcAmp_nul{ll+50,tt,ss,1}, featureSwAmp_nul{ll+50,tt,ss,1}]); % pre
          kaatn(ll,tt,ss,2,2,subuseind)=length([featureDeltaAmp_nul{ll+50,tt,ss,2}', featureKcAmp_nul{ll+50,tt,ss,2}', featureSwAmp_nul{ll+50,tt,ss,2}']); % post
          kaatm(ll,tt,ss,2,1,subuseind)=length([featureDeltaAmp_tacMSpN{ll,tt,ss,1}', featureKcAmp_tacMSpN{ll,tt,ss,1}', featureSwAmp_tacMSpN{ll,tt,ss,1}']); % pre
          kaatm(ll,tt,ss,2,2,subuseind)=length([featureDeltaAmp_tacMSpN{ll,tt,ss,2}', featureKcAmp_tacMSpN{ll,tt,ss,2}', featureSwAmp_tacMSpN{ll,tt,ss,2}']); % post
          kaatu(ll,tt,ss,2,1,subuseind)=length([featureDeltaAmp_tacPaud{ll,tt,ss,1}', featureKcAmp_tacPaud{ll,tt,ss,1}', featureSwAmp_tacPaud{ll,tt,ss,1}']); % pre
          kaatu(ll,tt,ss,2,2,subuseind)=length([featureDeltaAmp_tacPaud{ll,tt,ss,2}', featureKcAmp_tacPaud{ll,tt,ss,2}', featureSwAmp_tacPaud{ll,tt,ss,2}']); % post
          

          katm(ll,tt,ss,1,1,1,subuseind)=mean(featureDeltaAmp_tacMSpN{ll,tt,ss,1}); % pre
          katm(ll,tt,ss,1,2,1,subuseind)=mean(featureDeltaAmp_tacMSpN{ll,tt,ss,2}); % post
          katu(ll,tt,ss,1,1,1,subuseind)=mean(featureDeltaAmp_tacPaud{ll,tt,ss,1});
          katu(ll,tt,ss,1,2,1,subuseind)=mean(featureDeltaAmp_tacPaud{ll,tt,ss,2});
          katm(ll,tt,ss,2,1,1,subuseind)=length(featureDeltaAmp_tacMSpN{ll,tt,ss,1});
          katm(ll,tt,ss,2,2,1,subuseind)=length(featureDeltaAmp_tacMSpN{ll,tt,ss,2});
          katu(ll,tt,ss,2,1,1,subuseind)=length(featureDeltaAmp_tacPaud{ll,tt,ss,1});
          katu(ll,tt,ss,2,2,1,subuseind)=length(featureDeltaAmp_tacPaud{ll,tt,ss,2});
          katm(ll,tt,ss,1,1,2,subuseind)=mean(featureKcAmp_tacMSpN{ll,tt,ss,1}); % pre
          katm(ll,tt,ss,1,2,2,subuseind)=mean(featureKcAmp_tacMSpN{ll,tt,ss,2}); % post
          katu(ll,tt,ss,1,1,2,subuseind)=mean(featureKcAmp_tacPaud{ll,tt,ss,1});
          katu(ll,tt,ss,1,2,2,subuseind)=mean(featureKcAmp_tacPaud{ll,tt,ss,2});
          katm(ll,tt,ss,2,1,2,subuseind)=length(featureKcAmp_tacMSpN{ll,tt,ss,1});
          katm(ll,tt,ss,2,2,2,subuseind)=length(featureKcAmp_tacMSpN{ll,tt,ss,2});
          katu(ll,tt,ss,2,1,2,subuseind)=length(featureKcAmp_tacPaud{ll,tt,ss,1});
          katu(ll,tt,ss,2,2,2,subuseind)=length(featureKcAmp_tacPaud{ll,tt,ss,2});
          katm(ll,tt,ss,1,1,3,subuseind)=mean(featureSwAmp_tacMSpN{ll,tt,ss,1}); % pre
          katm(ll,tt,ss,1,2,3,subuseind)=mean(featureSwAmp_tacMSpN{ll,tt,ss,2}); % post
          katu(ll,tt,ss,1,1,3,subuseind)=mean(featureSwAmp_tacPaud{ll,tt,ss,1});
          katu(ll,tt,ss,1,2,3,subuseind)=mean(featureSwAmp_tacPaud{ll,tt,ss,2});
          katm(ll,tt,ss,2,1,3,subuseind)=length(featureSwAmp_tacMSpN{ll,tt,ss,1});
          katm(ll,tt,ss,2,2,3,subuseind)=length(featureSwAmp_tacMSpN{ll,tt,ss,2});
          katu(ll,tt,ss,2,1,3,subuseind)=length(featureSwAmp_tacPaud{ll,tt,ss,1});
          katu(ll,tt,ss,2,2,3,subuseind)=length(featureSwAmp_tacPaud{ll,tt,ss,2});
          
          
%           keyboard;
          % why is [5:8] of 5th dimension always all zero??  and why subuseind going up higher to 32?
%           sum(fstats_tpa([3 5],:,tt,ss,[3 4 7 8],2),2)


          % percentage of trials with this slowwave having negmax in window
          try kcts(ll,tt,ss,1,:,subuseind)=[featureDeltaAmp_tac{ll,tt,ss,3:4}]; end
          try kcts(ll,tt,ss,2,:,subuseind)=[featureKcAmp_tac{ll,tt,ss,3:4}]; end
          try kcts(ll,tt,ss,3,:,subuseind)=[featureSwAmp_tac{ll,tt,ss,3:4}]; end
          try kctt(ll,tt,ss,1,:,subuseind)=[featureDeltaAmp_tac{ll+20,tt,ss,3:4}]; end
          try kctt(ll,tt,ss,2,:,subuseind)=[featureKcAmp_tac{ll+20,tt,ss,3:4}]; end
          try kctt(ll,tt,ss,3,:,subuseind)=[featureSwAmp_tac{ll+20,tt,ss,3:4}]; end
          try kcta(ll,tt,ss,1,:,subuseind)=[featureDeltaAmp_aud{ll+40,tt,ss,3:4}]; end
          try kcta(ll,tt,ss,2,:,subuseind)=[featureKcAmp_aud{ll+40,tt,ss,3:4}]; end
          try kcta(ll,tt,ss,3,:,subuseind)=[featureSwAmp_aud{ll+40,tt,ss,3:4}]; end
          try kctn(ll,tt,ss,1,:,subuseind)=[featureDeltaAmp_nul{ll+50,tt,ss,3:4}]; end
          try kctn(ll,tt,ss,2,:,subuseind)=[featureKcAmp_nul{ll+50,tt,ss,3:4}]; end
          try kctn(ll,tt,ss,3,:,subuseind)=[featureSwAmp_nul{ll+50,tt,ss,3:4}]; end      
          kctu(ll,tt,ss,:,:,subuseind)=nanmean(cat(1,kctt(ll,tt,ss,:,:,subuseind), kcta(ll,tt,ss,:,:,subuseind)),1);
          kctm(ll,tt,ss,:,:,subuseind)=nanmean(cat(1,kcts(ll,tt,ss,:,:,subuseind), kctn(ll,tt,ss,:,:,subuseind)),1);
            
          fstats_tAudAlone(:,ll,tt,ss,:,subuseind)=featurestats_aud(:,ll+40,tt,ss,:);
          fstats_tMSAlone(:,ll,tt,ss,:,subuseind)=featurestats_tac(:,ll,tt,ss,:);
          fstats_tNulAlone(:,ll,tt,ss,:,subuseind)=featurestats_nul(:,ll+50,tt,ss,:);

          % we can save these out, but slightly more biased/wrong count of
          % Delta/Kc/SW, as it includes if any timepoint, not just max, whereas
          % above includes just if negmax included.
          % But these useful for spindle
          fstats_tTacAlone(:,ll,tt,ss,:,subuseind)=featurestats_tac(:,ll+20,tt,ss,:);
          fstats_tAudAlone(:,ll,tt,ss,:,subuseind)=featurestats_aud(:,ll+40,tt,ss,:);
          fstats_tMSAlone(:,ll,tt,ss,:,subuseind)=featurestats_tac(:,ll,tt,ss,:);
          fstats_tNulAlone(:,ll,tt,ss,:,subuseind)=featurestats_nul(:,ll+50,tt,ss,:);


          fstats_tpa(:,ll,tt,ss,:,subuseind)=featurestats_tacPaud(:,ll,tt,ss,:);
          fstats_tmspn(:,ll,tt,ss,:,subuseind)=featurestats_tacMSpN(:,ll,tt,ss,:);
          if audtacflag
            fstats_apt(:,ll,tt,ss,:,subuseind)=featurestats_audPtac(:,ll,tt,ss,:);
            fstats_amspn(:,ll,tt,ss,:,subuseind)=featurestats_audMSpN(:,ll,tt,ss,:);
          end
          %           end
        end % ss
        %       end % ll
        %       clear *_s0
      end % ll
      
    end % tt
    clear tlock*N tlock*tac tlock*aud
    
  end % ii
  subuseindfinal=subuseind;
  
  % do some processing on fstats* here, over all 'll' but within a 'tt'
  
  fstats_tpa(3,:,tt,:,:)
  
end % sleep
save([edir 'tlockSLEEP01_numtrlltt.mat'],'numtr*','fstats*','k*');

%%

% fstats (percentage of trials with features)
for ss=ssuse
  figure;
  % Nul_alone
  subplot(4,5,1);
  bar([squeeze(mean(fstats_tNulAlone(4,soalist,tt,ss,3:4,:),6)) squeeze(mean( fstats_tNulAlone(4,soalist,tt,ss,4,:)-fstats_tNulAlone(4,soalist,tt,ss,3,:),6))']);title('Delta'); ylabel('Nul')
  hold on; [hh,pp]=ttest(squeeze(fstats_tNulAlone(4,soalist,tt,ss,4,:)-fstats_tNulAlone(4,soalist,tt,ss,3,:))');
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,2);
  bar([squeeze(mean(fstats_tNulAlone(3,soalist,tt,ss,3:4,:),6)) squeeze(mean( fstats_tNulAlone(3,soalist,tt,ss,4,:)-fstats_tNulAlone(3,soalist,tt,ss,3,:),6))']);title('Kc')
  hold on; [hh,pp]=ttest(squeeze(fstats_tNulAlone(3,soalist,tt,ss,4,:)-fstats_tNulAlone(3,soalist,tt,ss,3,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,3);
  bar([squeeze(mean(fstats_tNulAlone(5,soalist,tt,ss,3:4,:),6)) squeeze(mean( fstats_tNulAlone(5,soalist,tt,ss,4,:)-fstats_tNulAlone(5,soalist,tt,ss,3,:),6))']);title('SW')
  hold on; [hh,pp]=ttest(squeeze(fstats_tNulAlone(5,soalist,tt,ss,4,:)-fstats_tNulAlone(5,soalist,tt,ss,3,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,4);
  bar([squeeze(mean(fstats_tNulAlone(6,soalist,tt,ss,3:4,:),6)) squeeze(mean( fstats_tNulAlone(6,soalist,tt,ss,4,:)-fstats_tNulAlone(6,soalist,tt,ss,3,:),6))']);title('Spindle fast')
  hold on; [hh,pp]=ttest(squeeze(fstats_tNulAlone(6,soalist,tt,ss,4,:)-fstats_tNulAlone(6,soalist,tt,ss,3,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,5);
  bar([squeeze(mean(fstats_tNulAlone(7,soalist,tt,ss,3:4,:),6)) squeeze(mean( fstats_tNulAlone(7,soalist,tt,ss,4,:)-fstats_tNulAlone(7,soalist,tt,ss,3,:),6))']);title('Spindle slow')
  hold on; [hh,pp]=ttest(squeeze(fstats_tNulAlone(7,soalist,tt,ss,4,:)-fstats_tNulAlone(7,soalist,tt,ss,3,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  % Aud_alone
  subplot(4,5,1+5);
  bar([squeeze(mean(fstats_tAudAlone(4,soalist,tt,ss,3:4,:),6)) squeeze(mean( fstats_tAudAlone(4,soalist,tt,ss,4,:)-fstats_tAudAlone(4,soalist,tt,ss,3,:),6))']);title('Delta'); ylabel('Aud')
  hold on; [hh,pp]=ttest(squeeze(fstats_tAudAlone(4,soalist,tt,ss,4,:)-fstats_tAudAlone(4,soalist,tt,ss,3,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,2+5);
  bar([squeeze(mean(fstats_tAudAlone(3,soalist,tt,ss,3:4,:),6)) squeeze(mean( fstats_tAudAlone(3,soalist,tt,ss,4,:)-fstats_tAudAlone(3,soalist,tt,ss,3,:),6))']);title('Kc')
  hold on; [hh,pp]=ttest(squeeze(fstats_tAudAlone(3,soalist,tt,ss,4,:)-fstats_tAudAlone(3,soalist,tt,ss,3,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,3+5);
  bar([squeeze(mean(fstats_tAudAlone(5,soalist,tt,ss,3:4,:),6)) squeeze(mean( fstats_tAudAlone(5,soalist,tt,ss,4,:)-fstats_tAudAlone(5,soalist,tt,ss,3,:),6))']);title('SW')
  hold on; [hh,pp]=ttest(squeeze(fstats_tAudAlone(5,soalist,tt,ss,4,:)-fstats_tAudAlone(5,soalist,tt,ss,3,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,4+5);
  bar([squeeze(mean(fstats_tAudAlone(6,soalist,tt,ss,3:4,:),6)) squeeze(mean( fstats_tAudAlone(6,soalist,tt,ss,4,:)-fstats_tAudAlone(6,soalist,tt,ss,3,:),6))']);title('Spindle fast')
  hold on; [hh,pp]=ttest(squeeze(fstats_tAudAlone(6,soalist,tt,ss,4,:)-fstats_tAudAlone(6,soalist,tt,ss,3,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,5+5);
  bar([squeeze(mean(fstats_tAudAlone(7,soalist,tt,ss,3:4,:),6)) squeeze(mean( fstats_tAudAlone(7,soalist,tt,ss,4,:)-fstats_tAudAlone(7,soalist,tt,ss,3,:),6))']);title('Spindle slow')
  hold on; [hh,pp]=ttest(squeeze(fstats_tAudAlone(7,soalist,tt,ss,4,:)-fstats_tAudAlone(7,soalist,tt,ss,3,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  % Tac_alone
  subplot(4,5,1+10);
  bar([squeeze(mean(fstats_tTacAlone(4,soalist,tt,ss,3:4,:),6)) squeeze(mean( fstats_tTacAlone(4,soalist,tt,ss,4,:)-fstats_tTacAlone(4,soalist,tt,ss,3,:),6))']);title('Delta'); ylabel('Tac')
  hold on; [hh,pp]=ttest(squeeze(fstats_tTacAlone(4,soalist,tt,ss,4,:)-fstats_tTacAlone(4,soalist,tt,ss,3,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,2+10);
  bar([squeeze(mean(fstats_tTacAlone(3,soalist,tt,ss,3:4,:),6)) squeeze(mean( fstats_tTacAlone(3,soalist,tt,ss,4,:)-fstats_tTacAlone(3,soalist,tt,ss,3,:),6))']);title('Kc')
  hold on; [hh,pp]=ttest(squeeze(fstats_tTacAlone(3,soalist,tt,ss,4,:)-fstats_tTacAlone(3,soalist,tt,ss,3,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,3+10);
  bar([squeeze(mean(fstats_tTacAlone(5,soalist,tt,ss,3:4,:),6)) squeeze(mean( fstats_tTacAlone(5,soalist,tt,ss,4,:)-fstats_tTacAlone(5,soalist,tt,ss,3,:),6))']);title('SW')
  hold on; [hh,pp]=ttest(squeeze(fstats_tTacAlone(5,soalist,tt,ss,4,:)-fstats_tTacAlone(5,soalist,tt,ss,3,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,4+10);
  bar([squeeze(mean(fstats_tTacAlone(6,soalist,tt,ss,3:4,:),6)) squeeze(mean( fstats_tTacAlone(6,soalist,tt,ss,4,:)-fstats_tTacAlone(6,soalist,tt,ss,3,:),6))']);title('Spindle fast')
  hold on; [hh,pp]=ttest(squeeze(fstats_tTacAlone(6,soalist,tt,ss,4,:)-fstats_tTacAlone(6,soalist,tt,ss,3,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,5+10);
  bar([squeeze(mean(fstats_tTacAlone(7,soalist,tt,ss,3:4,:),6)) squeeze(mean( fstats_tTacAlone(7,soalist,tt,ss,4,:)-fstats_tTacAlone(7,soalist,tt,ss,3,:),6))']);title('Spindle slow')
  hold on; [hh,pp]=ttest(squeeze(fstats_tTacAlone(7,soalist,tt,ss,4,:)-fstats_tTacAlone(7,soalist,tt,ss,3,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  % MS_alone
  subplot(4,5,1+15);
  bar([squeeze(mean(fstats_tMSAlone(4,soalist,tt,ss,3:4,:),6)) squeeze(mean( fstats_tMSAlone(4,soalist,tt,ss,4,:)-fstats_tMSAlone(4,soalist,tt,ss,3,:),6))']);title('Delta'); ylabel('MS')
  hold on; [hh,pp]=ttest(squeeze(fstats_tMSAlone(4,soalist,tt,ss,4,:)-fstats_tMSAlone(4,soalist,tt,ss,3,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,2+15);
  bar([squeeze(mean(fstats_tMSAlone(3,soalist,tt,ss,3:4,:),6)) squeeze(mean( fstats_tMSAlone(3,soalist,tt,ss,4,:)-fstats_tMSAlone(3,soalist,tt,ss,3,:),6))']);title('Kc')
  hold on; [hh,pp]=ttest(squeeze(fstats_tMSAlone(3,soalist,tt,ss,4,:)-fstats_tMSAlone(3,soalist,tt,ss,3,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,3+15);
  bar([squeeze(mean(fstats_tMSAlone(5,soalist,tt,ss,3:4,:),6)) squeeze(mean( fstats_tMSAlone(5,soalist,tt,ss,4,:)-fstats_tMSAlone(5,soalist,tt,ss,3,:),6))']);title('SW')
  hold on; [hh,pp]=ttest(squeeze(fstats_tMSAlone(5,soalist,tt,ss,4,:)-fstats_tMSAlone(5,soalist,tt,ss,3,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,4+15);
  bar([squeeze(mean(fstats_tMSAlone(6,soalist,tt,ss,3:4,:),6)) squeeze(mean( fstats_tMSAlone(6,soalist,tt,ss,4,:)-fstats_tMSAlone(6,soalist,tt,ss,3,:),6))']);title('Spindle fast')
  hold on; [hh,pp]=ttest(squeeze(fstats_tMSAlone(6,soalist,tt,ss,4,:)-fstats_tMSAlone(6,soalist,tt,ss,3,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,5+15);
  bar([squeeze(mean(fstats_tMSAlone(7,soalist,tt,ss,3:4,:),6)) squeeze(mean( fstats_tMSAlone(7,soalist,tt,ss,4,:)-fstats_tMSAlone(7,soalist,tt,ss,3,:),6))']);title('Spindle slow')
  hold on; [hh,pp]=ttest(squeeze(fstats_tMSAlone(7,soalist,tt,ss,4,:)-fstats_tMSAlone(7,soalist,tt,ss,3,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
end

% fstats (mean percentage of time spent in feature)
for ss=ssuse
  figure;
  % Nul_alone
  subplot(4,5,1);
  bar([squeeze(mean(fstats_tNulAlone(4,soalist,tt,ss,1:2,:),6)) squeeze(mean( fstats_tNulAlone(4,soalist,tt,ss,2,:)-fstats_tNulAlone(4,soalist,tt,ss,1,:),6))']);title('Delta'); ylabel('Nul')
  hold on; [hh,pp]=ttest(squeeze(fstats_tNulAlone(4,soalist,tt,ss,2,:)-fstats_tNulAlone(4,soalist,tt,ss,1,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,2);
  bar([squeeze(mean(fstats_tNulAlone(3,soalist,tt,ss,1:2,:),6)) squeeze(mean( fstats_tNulAlone(3,soalist,tt,ss,2,:)-fstats_tNulAlone(3,soalist,tt,ss,1,:),6))']);title('Kc')
  hold on; [hh,pp]=ttest(squeeze(fstats_tNulAlone(3,soalist,tt,ss,2,:)-fstats_tNulAlone(3,soalist,tt,ss,1,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,3);
  bar([squeeze(mean(fstats_tNulAlone(5,soalist,tt,ss,1:2,:),6)) squeeze(mean( fstats_tNulAlone(5,soalist,tt,ss,2,:)-fstats_tNulAlone(5,soalist,tt,ss,1,:),6))']);title('SW')
  hold on; [hh,pp]=ttest(squeeze(fstats_tNulAlone(5,soalist,tt,ss,2,:)-fstats_tNulAlone(5,soalist,tt,ss,1,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,4);
  bar([squeeze(mean(fstats_tNulAlone(6,soalist,tt,ss,1:2,:),6)) squeeze(mean( fstats_tNulAlone(6,soalist,tt,ss,2,:)-fstats_tNulAlone(6,soalist,tt,ss,1,:),6))']);title('Spindle fast')
  hold on; [hh,pp]=ttest(squeeze(fstats_tNulAlone(6,soalist,tt,ss,2,:)-fstats_tNulAlone(6,soalist,tt,ss,1,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,5);
  bar([squeeze(mean(fstats_tNulAlone(7,soalist,tt,ss,1:2,:),6)) squeeze(mean( fstats_tNulAlone(7,soalist,tt,ss,2,:)-fstats_tNulAlone(7,soalist,tt,ss,1,:),6))']);title('Spindle slow')
  hold on; [hh,pp]=ttest(squeeze(fstats_tNulAlone(7,soalist,tt,ss,2,:)-fstats_tNulAlone(7,soalist,tt,ss,1,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  % Aud_alone
  subplot(4,5,1+5);
  bar([squeeze(mean(fstats_tAudAlone(4,soalist,tt,ss,1:2,:),6)) squeeze(mean( fstats_tAudAlone(4,soalist,tt,ss,2,:)-fstats_tAudAlone(4,soalist,tt,ss,1,:),6))']);title('Delta'); ylabel('Aud')
  hold on; [hh,pp]=ttest(squeeze(fstats_tAudAlone(4,soalist,tt,ss,2,:)-fstats_tAudAlone(4,soalist,tt,ss,1,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,2+5);
  bar([squeeze(mean(fstats_tAudAlone(3,soalist,tt,ss,1:2,:),6)) squeeze(mean( fstats_tAudAlone(3,soalist,tt,ss,2,:)-fstats_tAudAlone(3,soalist,tt,ss,1,:),6))']);title('Kc')
  hold on; [hh,pp]=ttest(squeeze(fstats_tAudAlone(3,soalist,tt,ss,2,:)-fstats_tAudAlone(3,soalist,tt,ss,1,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,3+5);
  bar([squeeze(mean(fstats_tAudAlone(5,soalist,tt,ss,1:2,:),6)) squeeze(mean( fstats_tAudAlone(5,soalist,tt,ss,2,:)-fstats_tAudAlone(5,soalist,tt,ss,1,:),6))']);title('SW')
  hold on; [hh,pp]=ttest(squeeze(fstats_tAudAlone(5,soalist,tt,ss,2,:)-fstats_tAudAlone(5,soalist,tt,ss,1,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,4+5);
  bar([squeeze(mean(fstats_tAudAlone(6,soalist,tt,ss,1:2,:),6)) squeeze(mean( fstats_tAudAlone(6,soalist,tt,ss,2,:)-fstats_tAudAlone(6,soalist,tt,ss,1,:),6))']);title('Spindle fast')
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  hold on; [hh,pp]=ttest(squeeze(fstats_tAudAlone(6,soalist,tt,ss,2,:)-fstats_tAudAlone(6,soalist,tt,ss,1,:))'); 
  subplot(4,5,5+5);
  bar([squeeze(mean(fstats_tAudAlone(7,soalist,tt,ss,1:2,:),6)) squeeze(mean( fstats_tAudAlone(7,soalist,tt,ss,2,:)-fstats_tAudAlone(7,soalist,tt,ss,1,:),6))']);title('Spindle slow')
  hold on; [hh,pp]=ttest(squeeze(fstats_tAudAlone(7,soalist,tt,ss,2,:)-fstats_tAudAlone(7,soalist,tt,ss,1,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  % Tac_alone
  subplot(4,5,1+10);
  bar([squeeze(mean(fstats_tTacAlone(4,soalist,tt,ss,1:2,:),6)) squeeze(mean( fstats_tTacAlone(4,soalist,tt,ss,2,:)-fstats_tTacAlone(4,soalist,tt,ss,1,:),6))']);title('Delta'); ylabel('Tac')
  hold on; [hh,pp]=ttest(squeeze(fstats_tTacAlone(4,soalist,tt,ss,2,:)-fstats_tTacAlone(4,soalist,tt,ss,1,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,2+10);
  bar([squeeze(mean(fstats_tTacAlone(3,soalist,tt,ss,1:2,:),6)) squeeze(mean( fstats_tTacAlone(3,soalist,tt,ss,2,:)-fstats_tTacAlone(3,soalist,tt,ss,1,:),6))']);title('Kc')
  hold on; [hh,pp]=ttest(squeeze(fstats_tTacAlone(3,soalist,tt,ss,2,:)-fstats_tTacAlone(3,soalist,tt,ss,1,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,3+10);
  bar([squeeze(mean(fstats_tTacAlone(5,soalist,tt,ss,1:2,:),6)) squeeze(mean( fstats_tTacAlone(5,soalist,tt,ss,2,:)-fstats_tTacAlone(5,soalist,tt,ss,1,:),6))']);title('SW')
  hold on; [hh,pp]=ttest(squeeze(fstats_tTacAlone(5,soalist,tt,ss,2,:)-fstats_tTacAlone(5,soalist,tt,ss,1,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,4+10);
  bar([squeeze(mean(fstats_tTacAlone(6,soalist,tt,ss,1:2,:),6)) squeeze(mean( fstats_tTacAlone(6,soalist,tt,ss,2,:)-fstats_tTacAlone(6,soalist,tt,ss,1,:),6))']);title('Spindle fast')
  hold on; [hh,pp]=ttest(squeeze(fstats_tTacAlone(6,soalist,tt,ss,2,:)-fstats_tTacAlone(6,soalist,tt,ss,1,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,5+10);
  bar([squeeze(mean(fstats_tTacAlone(7,soalist,tt,ss,1:2,:),6)) squeeze(mean( fstats_tTacAlone(7,soalist,tt,ss,2,:)-fstats_tTacAlone(7,soalist,tt,ss,1,:),6))']);title('Spindle slow')
  hold on; [hh,pp]=ttest(squeeze(fstats_tTacAlone(7,soalist,tt,ss,2,:)-fstats_tTacAlone(7,soalist,tt,ss,1,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  % MS_alone
  subplot(4,5,1+15);
  bar([squeeze(mean(fstats_tMSAlone(4,soalist,tt,ss,1:2,:),6)) squeeze(mean( fstats_tMSAlone(4,soalist,tt,ss,2,:)-fstats_tMSAlone(4,soalist,tt,ss,1,:),6))']);title('Delta'); ylabel('MS')
  hold on; [hh,pp]=ttest(squeeze(fstats_tMSAlone(4,soalist,tt,ss,2,:)-fstats_tMSAlone(4,soalist,tt,ss,1,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,2+15);
  bar([squeeze(mean(fstats_tMSAlone(3,soalist,tt,ss,1:2,:),6)) squeeze(mean( fstats_tMSAlone(3,soalist,tt,ss,2,:)-fstats_tMSAlone(3,soalist,tt,ss,1,:),6))']);title('Kc')
  hold on; [hh,pp]=ttest(squeeze(fstats_tMSAlone(3,soalist,tt,ss,2,:)-fstats_tMSAlone(3,soalist,tt,ss,1,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,3+15);
  bar([squeeze(mean(fstats_tMSAlone(5,soalist,tt,ss,1:2,:),6)) squeeze(mean( fstats_tMSAlone(5,soalist,tt,ss,2,:)-fstats_tMSAlone(5,soalist,tt,ss,1,:),6))']);title('SW')
  hold on; [hh,pp]=ttest(squeeze(fstats_tMSAlone(5,soalist,tt,ss,2,:)-fstats_tMSAlone(5,soalist,tt,ss,1,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,4+15);
  bar([squeeze(mean(fstats_tMSAlone(6,soalist,tt,ss,1:2,:),6)) squeeze(mean( fstats_tMSAlone(6,soalist,tt,ss,2,:)-fstats_tMSAlone(6,soalist,tt,ss,1,:),6))']);title('Spindle fast')
  hold on; [hh,pp]=ttest(squeeze(fstats_tMSAlone(6,soalist,tt,ss,2,:)-fstats_tMSAlone(6,soalist,tt,ss,1,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
  subplot(4,5,5+15);
  bar([squeeze(mean(fstats_tMSAlone(7,soalist,tt,ss,1:2,:),6)) squeeze(mean( fstats_tMSAlone(7,soalist,tt,ss,2,:)-fstats_tMSAlone(7,soalist,tt,ss,1,:),6))']);title('Spindle slow')
  hold on; [hh,pp]=ttest(squeeze(fstats_tMSAlone(7,soalist,tt,ss,2,:)-fstats_tMSAlone(7,soalist,tt,ss,1,:))'); 
  plot(1:7,.1*(pp<.05).*log10(pp),'r*');
end


% Histograms
histbins=histbins(6:end-6);
khtm=khtm(:,:,:,6:end-6,:,:);
khtu=khtu(:,:,:,6:end-6,:,:);
khtt=khtt(:,:,:,6:end-6,:,:);
khta=khta(:,:,:,6:end-6,:,:);
khtn=khtn(:,:,:,6:end-6,:,:);
khts=khts(:,:,:,6:end-6,:,:);
kahtm=kahtm(:,:,:,6:end-6,:);
kahtu=kahtu(:,:,:,6:end-6,:);
kahtt=kahtt(:,:,:,6:end-6,:);
kahta=kahta(:,:,:,6:end-6,:);
kahtn=kahtn(:,:,:,6:end-6,:);
kahts=kahts(:,:,:,6:end-6,:);


% explore which Stage to use
figure(10);
ll=4; % simultaneous
for ss=ssuse(2:end)
  if ss==23
    subplot(4,4,1+3);
  else
    subplot(4,4,1+ss-11);
  end
  imagesc(squeeze(kahtn(soalist(ll),tt,ss,:,:))');colorbar;colormap('gray')
  if ss==11,title('N1');end
  if ss==12,title('N2');end
  if ss==13,title('N3');end
  if ss==23,title('N2+N3');end
  if ss==11,ylabel('Null');end
  if ss==23
    subplot(4,4,5+3);
  else
    subplot(4,4,5+ss-11);
  end
  imagesc(squeeze(kahtt(soalist(ll),tt,ss,:,:))');colorbar;
  if ss==11,ylabel('Tac');end
  if ss==23
    subplot(4,4,9+3);
  else
    subplot(4,4,9+ss-11);
  end
  imagesc(squeeze(kahta(soalist(ll),tt,ss,:,:))');colorbar;
  if ss==11,ylabel('Aud');end
  if ss==23
    subplot(4,4,13+3);
  else
    subplot(4,4,13+ss-11);
  end
  imagesc(squeeze(kahts(soalist(ll),tt,ss,:,:))');colorbar;
  if ss==11,ylabel('MultSens');end
end
% conclude that N2 alone (not N2+N3) most clear results (low baseline,
% stronge evoked)

alphasubthresh=0.05;
alphathresh=0.005;
for ss=ssuse
  for ll=soalist,ttest(squeeze(khtm(ll,tt,ss,:,1,:))'-squeeze(khtu(ll,tt,ss,:,1,:))',[],'alpha',alphathresh),end
  for ll=soalist,ttest(squeeze(khtm(ll,tt,ss,:,2,:))'-squeeze(khtu(ll,tt,ss,:,2,:))',[],'alpha',alphathresh),end
  for ll=soalist,ttest(squeeze(khtm(ll,tt,ss,:,3,:))'-squeeze(khtu(ll,tt,ss,:,3,:))',[],'alpha',alphathresh),end
  
  for ll=soalist,ttest(squeeze(khta(ll,tt,ss,:,1,:))'-squeeze(khtn(ll,tt,ss,:,1,:))',[],'alpha',alphathresh),end
  
  % multisens+nul versus sum-of-unisens, for each Delta, Kc, SW
  figure(1);for ll=1:length(soalist),subplot(1,7,ll);plot(histbins(:),mean(squeeze(khtm(soalist(ll),tt,ss,:,1,:))'-squeeze(khtu(soalist(ll),tt,ss,:,1,:))'));axis([-inf inf -3 3]);
    testtime=find(ttest(squeeze(khta(soalist(ll),tt,ss,:,1,:))'-squeeze(khtn(soalist(ll),tt,ss,:,1,:))',[],'alpha',alphasubthresh));
    hold on;plot(histbins(:),ttest(squeeze(khtm(soalist(ll),tt,ss,:,1,:))'-squeeze(khtu(soalist(ll),tt,ss,:,1,:))',[],'alpha',alphathresh),'*');
    hold on;plot(histbins(testtime),ttest(squeeze(khtm(soalist(ll),tt,ss,testtime,1,:))'-squeeze(khtu(soalist(ll),tt,ss,testtime,1,:))',[],'alpha',alphasubthresh),'r*');end
  figure(2);for ll=1:length(soalist),subplot(1,7,ll);plot(histbins(:),mean(squeeze(khtm(soalist(ll),tt,ss,:,2,:))'-squeeze(khtu(soalist(ll),tt,ss,:,2,:))'));axis([-inf inf -3 3]);
    testtime=find(ttest(squeeze(khta(soalist(ll),tt,ss,:,2,:))'-squeeze(khtn(soalist(ll),tt,ss,:,2,:))',[],'alpha',alphasubthresh));
    hold on;plot(histbins(:),ttest(squeeze(khtm(soalist(ll),tt,ss,:,2,:))'-squeeze(khtu(soalist(ll),tt,ss,:,2,:))',[],'alpha',alphathresh),'*');
    hold on;plot(histbins(testtime),ttest(squeeze(khtm(soalist(ll),tt,ss,testtime,2,:))'-squeeze(khtu(soalist(ll),tt,ss,testtime,2,:))',[],'alpha',alphasubthresh),'r*');end
  figure(3);for ll=1:length(soalist),subplot(1,7,ll);plot(histbins(:),mean(squeeze(khtm(soalist(ll),tt,ss,:,3,:))'-squeeze(khtu(soalist(ll),tt,ss,:,3,:))'));axis([-inf inf -3 3]);
    testtime=find(ttest(squeeze(khta(soalist(ll),tt,ss,:,3,:))'-squeeze(khtn(soalist(ll),tt,ss,:,3,:))',[],'alpha',alphasubthresh));
    hold on;plot(histbins(:),ttest(squeeze(khtm(soalist(ll),tt,ss,:,3,:))'-squeeze(khtu(soalist(ll),tt,ss,:,3,:))',[],'alpha',alphathresh),'*');
    hold on;plot(histbins(testtime),ttest(squeeze(khtm(soalist(ll),tt,ss,testtime,3,:))'-squeeze(khtu(soalist(ll),tt,ss,testtime,3,:))',[],'alpha',alphasubthresh),'r*');end
  
  % multisens+nul versus sum-of-unisens (combined over Delta, Kc, SW)
  figure(4);for ll=1:length(soalist),subplot(1,7,ll);plot(histbins(:),mean(squeeze(kahtm(soalist(ll),tt,ss,:,:))'-squeeze(kahtu(soalist(ll),tt,ss,:,:))'));axis([-inf inf -3 3]);
    testtime=find(ttest(squeeze(kahta(soalist(ll),tt,ss,:,:))'-squeeze(kahtn(soalist(ll),tt,ss,:,:))',[],'alpha',alphasubthresh));
    hold on;plot(histbins(:),ttest(squeeze(kahtm(soalist(ll),tt,ss,:,:))'-squeeze(kahtu(soalist(ll),tt,ss,:,:))',[],'alpha',alphathresh),'*');
    hold on;plot(histbins(testtime),ttest(squeeze(kahtm(soalist(ll),tt,ss,testtime,:))'-squeeze(kahtu(soalist(ll),tt,ss,testtime,:))',[],'alpha',alphasubthresh),'r*');end
  
  % multisens versus aud-alone
  figure(5);for ll=1:length(soalist),subplot(1,7,ll);plot(histbins(:),mean(squeeze(kahts(soalist(ll),tt,ss,:,:))'-squeeze(kahta(soalist(ll),tt,ss,:,:))'));axis([-inf inf -3 3]);
    testtime=find(ttest(squeeze(kahtt(soalist(ll),tt,ss,:,:))'-squeeze(kahtn(soalist(ll),tt,ss,:,:))',[],'alpha',alphasubthresh));
    hold on;plot(histbins(:),ttest(squeeze(kahts(soalist(ll),tt,ss,:,:))'-squeeze(kahta(soalist(ll),tt,ss,:,:))',[],'alpha',alphathresh),'*');
    hold on;plot(histbins(testtime),ttest(squeeze(kahts(soalist(ll),tt,ss,testtime,:))'-squeeze(kahta(soalist(ll),tt,ss,testtime,:))',[],'alpha',alphasubthresh),'r*');end

  % multisens versus tac-alone
  figure(6);for ll=1:length(soalist),subplot(1,7,ll);plot(histbins(:),mean(squeeze(kahts(soalist(ll),tt,ss,:,:))'-squeeze(kahtt(soalist(ll),tt,ss,:,:))'));axis([-inf inf -3 3]);
    testtime=find(ttest(squeeze(kahta(soalist(ll),tt,ss,:,:))'-squeeze(kahtn(soalist(ll),tt,ss,:,:))',[],'alpha',alphasubthresh));
    hold on;plot(histbins(:),ttest(squeeze(kahts(soalist(ll),tt,ss,:,:))'-squeeze(kahtt(soalist(ll),tt,ss,:,:))',[],'alpha',alphathresh),'*');
    hold on;plot(histbins(testtime),ttest(squeeze(kahts(soalist(ll),tt,ss,testtime,:))'-squeeze(kahtt(soalist(ll),tt,ss,testtime,:))',[],'alpha',alphasubthresh),'r*');end

  % each condition/unisensory on its own, separate for Delta, Kc, SW
  figure(1);for ll=soalist,subplot(1,9,ll);plot(histbins(1:end-1),mean(squeeze(khtt(ll,tt,ss,1:end-1,1,:))'));axis([-inf inf -inf inf]);end
  figure(2);for ll=soalist,subplot(1,9,ll);plot(histbins(1:end-1),mean(squeeze(khtt(ll,tt,ss,1:end-1,2,:))'));axis([-inf inf -inf inf]);end
  figure(3);for ll=soalist,subplot(1,9,ll);plot(histbins(1:end-1),mean(squeeze(khtt(ll,tt,ss,1:end-1,3,:))'));axis([-inf inf -inf inf]);end

  figure(1);for ll=soalist,subplot(1,9,ll);plot(histbins(1:end-1),mean(squeeze(khta(ll,tt,ss,1:end-1,1,:))'));axis([-inf inf -inf inf]);end
  figure(2);for ll=soalist,subplot(1,9,ll);plot(histbins(1:end-1),mean(squeeze(khta(ll,tt,ss,1:end-1,2,:))'));axis([-inf inf -inf inf]);end
  figure(3);for ll=soalist,subplot(1,9,ll);plot(histbins(1:end-1),mean(squeeze(khta(ll,tt,ss,1:end-1,3,:))'));axis([-inf inf -inf inf]);end

  figure(1);for ll=soalist,subplot(1,9,ll);plot(histbins(1:end-1),mean(squeeze(khtn(ll,tt,ss,1:end-1,1,:))'));axis([-inf inf -inf inf]);end
  figure(2);for ll=soalist,subplot(1,9,ll);plot(histbins(1:end-1),mean(squeeze(khtn(ll,tt,ss,1:end-1,2,:))'));axis([-inf inf -inf inf]);end
  figure(3);for ll=soalist,subplot(1,9,ll);plot(histbins(1:end-1),mean(squeeze(khtn(ll,tt,ss,1:end-1,3,:))'));axis([-inf inf -inf inf]);end

  figure(1);for ll=soalist,subplot(1,9,ll);plot(histbins(1:end-1),mean(squeeze(khtm(ll,tt,ss,1:end-1,1,:))'));axis([-inf inf -inf inf]);end
  figure(2);for ll=soalist,subplot(1,9,ll);plot(histbins(1:end-1),mean(squeeze(khtm(ll,tt,ss,1:end-1,2,:))'));axis([-inf inf -inf inf]);end
  figure(3);for ll=soalist,subplot(1,9,ll);plot(histbins(1:end-1),mean(squeeze(khtm(ll,tt,ss,1:end-1,3,:))'));axis([-inf inf -inf inf]);end
  
  % first, collapsed over all 'large waves'
  stimtac=zeros(1,length(histbins));
  stimtac(dsearchn(histbins',0))=2;
  soatime=[-.5 -.07 -.02 0 .02 .07 .5];
  figure(4);
  for ll=1:length(soalist)
    subplot(1,7,ll)
    hold on;
    plot(histbins(:),mean(squeeze(kahtn(soalist(ll),tt,ss,:,:))'),'k');
    plot(histbins(:),mean(squeeze(kahtt(soalist(ll),tt,ss,:,:))'),'r');
    plot(histbins(:),mean(squeeze(kahta(soalist(ll),tt,ss,:,:))'),'g');
    plot(histbins(:),mean(squeeze(kahts(soalist(ll),tt,ss,:,:))'),'b');
      axis([-inf inf 0 4]);
    legend('Null','Tactile','Auditory','Multisensory')
    stimaud=zeros(1,length(histbins));
    stimaud(dsearchn(histbins',soatime(ll)))=2;
    plot(histbins(:),stimtac(:),'r');
    plot(histbins(:),stimaud(:),'g');
  end
%   maxx=[6 8 6]; % ss=23
  maxx=[5 6 4]; % ss=12
  stimtac=zeros(1,length(histbins));
  stimtac(dsearchn(histbins',0))=0.5;
  for kk=1:3
    figure(kk);
    for ll=1:length(soalist)
      subplot(1,7,ll)
      hold on;
      plot(histbins(:),mean(squeeze(khtn(soalist(ll),tt,ss,:,kk,:))'),'k');
      plot(histbins(:),mean(squeeze(khtt(soalist(ll),tt,ss,:,kk,:))'),'r');
      plot(histbins(:),mean(squeeze(khta(soalist(ll),tt,ss,:,kk,:))'),'g');
      plot(histbins(:),mean(squeeze(khts(soalist(ll),tt,ss,:,kk,:))'),'b');
      axis([-inf inf 0 maxx(kk)]);
      legend('Null','Tactile','Auditory','Multisensory')
      stimaud=zeros(1,length(histbins));
      stimaud(dsearchn(histbins',soatime(ll)))=0.5;
      plot(histbins(:),stimtac(:),'r');
      plot(histbins(:),stimaud(:),'g');
    end
  end
end

% amplitude and number total

% 1=mean (ave amplitude)
% 2=length (ave number)
for ss=ssuse
  figure;
  subplot(2,3,1);
  bar([squeeze(mean(kaatn(soalist,tt,ss,1,:,:),6)) squeeze(mean(kaatn(soalist,tt,ss,1,2,:)-kaatn(soalist,tt,ss,1,1,:),6))]);ylabel('Nul')
  hold on; [hh,pp]=ttest(squeeze(kaatn(soalist,tt,ss,1,2,:)-kaatn(soalist,tt,ss,1,1,:))');
  plot(1:7,10*(pp<.05).*log10(pp),'r*');axis([-inf inf -60 110])
  subplot(2,3,2);
  bar([squeeze(mean(kaata(soalist,tt,ss,1,:,:),6)) squeeze(mean(kaata(soalist,tt,ss,1,2,:)-kaata(soalist,tt,ss,1,1,:),6))]);ylabel('Aud')
  hold on; [hh,pp]=ttest(squeeze(kaata(soalist,tt,ss,1,2,:)-kaata(soalist,tt,ss,1,1,:))');
  plot(1:7,10*(pp<.05).*log10(pp),'r*');axis([-inf inf -60 110])
  subplot(2,3,3);
  bar([squeeze(mean(kaatt(soalist,tt,ss,1,:,:),6)) squeeze(mean(kaatt(soalist,tt,ss,1,2,:)-kaatt(soalist,tt,ss,1,1,:),6))]);ylabel('Tac')
  hold on; [hh,pp]=ttest(squeeze(kaatt(soalist,tt,ss,1,2,:)-kaatt(soalist,tt,ss,1,1,:))');
  plot(1:7,10*(pp<.05).*log10(pp),'r*');axis([-inf inf -60 110])
  subplot(2,3,4);
  bar([squeeze(mean(kaats(soalist,tt,ss,1,:,:),6)) squeeze(mean(kaats(soalist,tt,ss,1,2,:)-kaats(soalist,tt,ss,1,1,:),6))]);ylabel('MS')
  hold on; [hh,pp]=ttest(squeeze(kaats(soalist,tt,ss,1,2,:)-kaats(soalist,tt,ss,1,1,:))');
  plot(1:7,10*(pp<.05).*log10(pp),'r*');axis([-inf inf -60 110])
  subplot(2,3,5);
  bar([squeeze(mean(kaatm(soalist,tt,ss,1,:,:),6)) squeeze(mean(kaatm(soalist,tt,ss,1,2,:)-kaatm(soalist,tt,ss,1,1,:),6))]);ylabel('MS+Nul')
  hold on; [hh,pp]=ttest(squeeze(kaatm(soalist,tt,ss,1,2,:)-kaatm(soalist,tt,ss,1,1,:))');
  plot(1:7,10*(pp<.05).*log10(pp),'r*');axis([-inf inf -60 110])
  subplot(2,3,6);
  bar([squeeze(mean(kaatu(soalist,tt,ss,1,:,:),6)) squeeze(mean(kaatu(soalist,tt,ss,1,2,:)-kaatu(soalist,tt,ss,1,1,:),6))]);ylabel('Tac+Aud')
  hold on; [hh,pp]=ttest(squeeze(kaatu(soalist,tt,ss,1,2,:)-kaatu(soalist,tt,ss,1,1,:))');
  plot(1:7,10*(pp<.05).*log10(pp),'r*');axis([-inf inf -60 110])

  % this way may be flawed (see above where created, and below.
  figure; 
  subplot(2,3,1);
  bar([squeeze(mean(kaatn(soalist,tt,ss,2,:,:),6)) squeeze(mean(kaatn(soalist,tt,ss,2,2,:)-kaatn(soalist,tt,ss,2,1,:),6))]);ylabel('Nul')
  hold on; [hh,pp]=ttest(squeeze(kaatn(soalist,tt,ss,2,2,:)-kaatn(soalist,tt,ss,2,1,:))');
  plot(1:7,10*(pp<.05).*log10(pp),'r*');axis([-inf inf -70 150])
  subplot(2,3,2);
  bar([squeeze(mean(kaata(soalist,tt,ss,2,:,:),6)) squeeze(mean(kaata(soalist,tt,ss,2,2,:)-kaata(soalist,tt,ss,2,1,:),6))]);ylabel('Aud')
  hold on; [hh,pp]=ttest(squeeze(kaata(soalist,tt,ss,2,2,:)-kaata(soalist,tt,ss,2,1,:))');
  plot(1:7,10*(pp<.05).*log10(pp),'r*');axis([-inf inf -70 150])
  subplot(2,3,3);
  bar([squeeze(mean(kaatt(soalist,tt,ss,2,:,:),6)) squeeze(mean(kaatt(soalist,tt,ss,2,2,:)-kaatt(soalist,tt,ss,2,1,:),6))]);ylabel('Tac')
  hold on; [hh,pp]=ttest(squeeze(kaatt(soalist,tt,ss,2,2,:)-kaatt(soalist,tt,ss,2,1,:))');
  plot(1:7,10*(pp<.05).*log10(pp),'r*');axis([-inf inf -70 150])
  subplot(2,3,4);
  bar([squeeze(mean(kaats(soalist,tt,ss,2,:,:),6)) squeeze(mean(kaats(soalist,tt,ss,2,2,:)-kaats(soalist,tt,ss,2,1,:),6))]);ylabel('MS')
  hold on; [hh,pp]=ttest(squeeze(kaats(soalist,tt,ss,2,2,:)-kaats(soalist,tt,ss,2,1,:))');
  plot(1:7,10*(pp<.05).*log10(pp),'r*');axis([-inf inf -70 150])
  subplot(2,3,5);
  bar([squeeze(mean(kaatm(soalist,tt,ss,2,:,:),6)) squeeze(mean(kaatm(soalist,tt,ss,2,2,:)-kaatm(soalist,tt,ss,2,1,:),6))]);ylabel('MS+Nul')
  hold on; [hh,pp]=ttest(squeeze(kaatm(soalist,tt,ss,2,2,:)-kaatm(soalist,tt,ss,2,1,:))');
  plot(1:7,10*(pp<.05).*log10(pp),'r*');axis([-inf inf -70 150])
  subplot(2,3,6);
  bar([squeeze(mean(kaatu(soalist,tt,ss,2,:,:),6)) squeeze(mean(kaatu(soalist,tt,ss,2,2,:)-kaatu(soalist,tt,ss,2,1,:),6))]);ylabel('Tac+Aud')
  hold on; [hh,pp]=ttest(squeeze(kaatu(soalist,tt,ss,2,2,:)-kaatu(soalist,tt,ss,2,1,:))');
  plot(1:7,10*(pp<.05).*log10(pp),'r*');axis([-inf inf -70 150])
  
  % better count of percentage of trials
  figure;  % ave over delta/kc/sw (but hides shift from small to big)
  subplot(2,3,1);
  bar([squeeze(mean(mean(kctn(soalist,tt,ss,:,:,:),4),6))  squeeze(mean( mean(kctn(soalist,tt,ss,:,2,:),4)-mean(kctn(soalist,tt,ss,:,1,:),4) ,6))]);ylabel('Nul');
  hold on; [hh,pp]=ttest(squeeze(mean(kctn(soalist,tt,ss,:,2,:),4)-mean(kctn(soalist,tt,ss,:,1,:),4) )');
  plot(1:7,0.1*(pp<.05).*log10(pp),'r*');axis([-inf inf -.3 .3])
  subplot(2,3,2);
  bar([squeeze(mean(mean(kcta(soalist,tt,ss,:,:,:),4),6))  squeeze(mean( mean(kcta(soalist,tt,ss,:,2,:),4)-mean(kcta(soalist,tt,ss,:,1,:),4) ,6))]);ylabel('Aud');
  hold on; [hh,pp]=ttest(squeeze(mean(kcta(soalist,tt,ss,:,2,:),4)-mean(kcta(soalist,tt,ss,:,1,:),4) )');
  plot(1:7,0.1*(pp<.05).*log10(pp),'r*');axis([-inf inf -.3 .3])
  subplot(2,3,3);
  bar([squeeze(mean(mean(kctt(soalist,tt,ss,:,:,:),4),6))  squeeze(mean( mean(kctt(soalist,tt,ss,:,2,:),4)-mean(kctt(soalist,tt,ss,:,1,:),4) ,6))]);ylabel('Tac');
  hold on; [hh,pp]=ttest(squeeze(mean(kctt(soalist,tt,ss,:,2,:),4)-mean(kctt(soalist,tt,ss,:,1,:),4) )');
  plot(1:7,0.1*(pp<.05).*log10(pp),'r*');axis([-inf inf -.3 .3])
  subplot(2,3,4);
  bar([squeeze(mean(mean(kcts(soalist,tt,ss,:,:,:),4),6))  squeeze(mean( mean(kcts(soalist,tt,ss,:,2,:),4)-mean(kcts(soalist,tt,ss,:,1,:),4) ,6))]);ylabel('MS');
  hold on; [hh,pp]=ttest(squeeze(mean(kcts(soalist,tt,ss,:,2,:),4)-mean(kcts(soalist,tt,ss,:,1,:),4) )');
  plot(1:7,0.1*(pp<.05).*log10(pp),'r*');axis([-inf inf -.3 .3])
  subplot(2,3,5);
  bar([squeeze(mean(mean(kctm(soalist,tt,ss,:,:,:),4),6))  squeeze(mean( mean(kctm(soalist,tt,ss,:,2,:),4)-mean(kctm(soalist,tt,ss,:,1,:),4) ,6))]);ylabel('MS+Nul');
  hold on; [hh,pp]=ttest(squeeze(mean(kctm(soalist,tt,ss,:,2,:),4)-mean(kctm(soalist,tt,ss,:,1,:),4) )');
  plot(1:7,0.1*(pp<.05).*log10(pp),'r*');axis([-inf inf -.3 .3])
  subplot(2,3,6);
  bar([squeeze(mean(mean(kctu(soalist,tt,ss,:,:,:),4),6))  squeeze(mean( mean(kctu(soalist,tt,ss,:,2,:),4)-mean(kctu(soalist,tt,ss,:,1,:),4) ,6))]);ylabel('Tac+Aud');
  hold on; [hh,pp]=ttest(squeeze(mean(kctu(soalist,tt,ss,:,2,:),4)-mean(kctu(soalist,tt,ss,:,1,:),4) )');
  plot(1:7,0.1*(pp<.05).*log10(pp),'r*');axis([-inf inf -.3 .3])
  
  for kk=1:3
    figure;  % just delta
    subplot(2,3,1);
    bar([squeeze(mean(mean(kctn(soalist,tt,ss,kk,:,:),4),6))  squeeze(mean( mean(kctn(soalist,tt,ss,kk,2,:),4)-mean(kctn(soalist,tt,ss,kk,1,:),4) ,6))]);ylabel('Nul');
    hold on; [hh,pp]=ttest(squeeze(mean(kctn(soalist,tt,ss,kk,2,:),4)-mean(kctn(soalist,tt,ss,kk,1,:),4) )');
    plot(1:7,0.1*(pp<.05).*log10(pp),'r*');axis([-inf inf -.3 .4])
    subplot(2,3,2);
    bar([squeeze(mean(mean(kcta(soalist,tt,ss,kk,:,:),4),6))  squeeze(mean( mean(kcta(soalist,tt,ss,kk,2,:),4)-mean(kcta(soalist,tt,ss,kk,1,:),4) ,6))]);ylabel('Aud');
    hold on; [hh,pp]=ttest(squeeze(mean(kcta(soalist,tt,ss,kk,2,:),4)-mean(kcta(soalist,tt,ss,kk,1,:),4) )');
    plot(1:7,0.1*(pp<.05).*log10(pp),'r*');axis([-inf inf -.3 .4])
    subplot(2,3,3);
    bar([squeeze(mean(mean(kctt(soalist,tt,ss,kk,:,:),4),6))  squeeze(mean( mean(kctt(soalist,tt,ss,kk,2,:),4)-mean(kctt(soalist,tt,ss,kk,1,:),4) ,6))]);ylabel('Tac');
    hold on; [hh,pp]=ttest(squeeze(mean(kctt(soalist,tt,ss,kk,2,:),4)-mean(kctt(soalist,tt,ss,kk,1,:),4) )');
    plot(1:7,0.1*(pp<.05).*log10(pp),'r*');axis([-inf inf -.3 .4])
    subplot(2,3,4);
    bar([squeeze(mean(mean(kcts(soalist,tt,ss,kk,:,:),4),6))  squeeze(mean( mean(kcts(soalist,tt,ss,kk,2,:),4)-mean(kcts(soalist,tt,ss,kk,1,:),4) ,6))]);ylabel('MS');
    hold on; [hh,pp]=ttest(squeeze(mean(kcts(soalist,tt,ss,kk,2,:),4)-mean(kcts(soalist,tt,ss,kk,1,:),4) )');
    plot(1:7,0.1*(pp<.05).*log10(pp),'r*');axis([-inf inf -.3 .4])
    subplot(2,3,5);
    bar([squeeze(mean(mean(kctm(soalist,tt,ss,kk,:,:),4),6))  squeeze(mean( mean(kctm(soalist,tt,ss,kk,2,:),4)-mean(kctm(soalist,tt,ss,kk,1,:),4) ,6))]);ylabel('MS+Nul');
    hold on; [hh,pp]=ttest(squeeze(mean(kctm(soalist,tt,ss,kk,2,:),4)-mean(kctm(soalist,tt,ss,kk,1,:),4) )');
    plot(1:7,0.1*(pp<.05).*log10(pp),'r*');axis([-inf inf -.3 .4])
    subplot(2,3,6);
    bar([squeeze(mean(mean(kctu(soalist,tt,ss,kk,:,:),4),6))  squeeze(mean( mean(kctu(soalist,tt,ss,kk,2,:),4)-mean(kctu(soalist,tt,ss,kk,1,:),4) ,6))]);ylabel('Tac+Aud');
    hold on; [hh,pp]=ttest(squeeze(mean(kctu(soalist,tt,ss,kk,2,:),4)-mean(kctu(soalist,tt,ss,kk,1,:),4) )');
    plot(1:7,0.1*(pp<.05).*log10(pp),'r*');axis([-inf inf -.3 .4])
  end
  
  
end


for ss=ssuse
  %
  % testing for sum-uni difference of [delta, kc, sw] for post-pre/post, amplitude
  try [hh,pp]=ttest(squeeze([katu(soalist,tt,ss,1,2,1,:)-katu(soalist,tt,ss,1,1,1,:)]./katu(soalist,tt,ss,1,1,1,:))',[],'alpha',.05/7),end
  try [hh,pp]=ttest(squeeze([katu(soalist,tt,ss,1,2,2,:)-katu(soalist,tt,ss,1,1,2,:)]./katu(soalist,tt,ss,1,1,2,:))',[],'alpha',.05/7),end
  try [hh,pp]=ttest(squeeze([katu(soalist,tt,ss,1,2,3,:)-katu(soalist,tt,ss,1,1,3,:)]./katu(soalist,tt,ss,1,1,3,:))',[],'alpha',.05/7),end
  
  % testing for sum-uni difference of [delta, kc, sw] for post-pre/post, number
  try [hh,pp]=ttest(squeeze([katu(soalist,tt,ss,2,2,1,:)-katu(soalist,tt,ss,2,1,1,:)]./katu(soalist,tt,ss,2,1,1,:))',[],'alpha',.05/7),end
  try [hh,pp]=ttest(squeeze([katu(soalist,tt,ss,2,2,2,:)-katu(soalist,tt,ss,2,1,2,:)]./katu(soalist,tt,ss,2,1,2,:))',[],'alpha',.05/7),end
  try [hh,pp]=ttest(squeeze([katu(soalist,tt,ss,2,2,3,:)-katu(soalist,tt,ss,2,1,3,:)]./katu(soalist,tt,ss,2,1,3,:))',[],'alpha',.05/7),end
  
  % testing for multisens difference of [delta, kc, sw] for post-pre/post, amplitude
  try [hh,pp]=ttest(squeeze([katm(soalist,tt,ss,1,2,1,:)-katm(soalist,tt,ss,1,1,1,:)]./katm(soalist,tt,ss,1,1,1,:))',[],'alpha',.05/7),end
  try [hh,pp]=ttest(squeeze([katm(soalist,tt,ss,1,2,2,:)-katm(soalist,tt,ss,1,1,2,:)]./katm(soalist,tt,ss,1,1,2,:))',[],'alpha',.05/7),end
  try [hh,pp]=ttest(squeeze([katm(soalist,tt,ss,1,2,3,:)-katm(soalist,tt,ss,1,1,3,:)]./katm(soalist,tt,ss,1,1,3,:))',[],'alpha',.05/7),end
  
  % testing for multisens difference of [delta, kc, sw] for post-pre/post, number
  try [hh,pp]=ttest(squeeze([katm(soalist,tt,ss,2,2,1,:)-katm(soalist,tt,ss,2,1,1,:)]./katm(soalist,tt,ss,2,1,1,:))',[],'alpha',.05/7),end
  try [hh,pp]=ttest(squeeze([katm(soalist,tt,ss,2,2,2,:)-katm(soalist,tt,ss,2,1,2,:)]./katm(soalist,tt,ss,2,1,2,:))',[],'alpha',.05/7),end
  try [hh,pp]=ttest(squeeze([katm(soalist,tt,ss,2,2,3,:)-katm(soalist,tt,ss,2,1,3,:)]./katm(soalist,tt,ss,2,1,3,:))',[],'alpha',.05/7),end
  
  % interaction: of amplitude
  try [hh,pp]=ttest([squeeze([katu(soalist,tt,ss,1,2,1,:)-katu(soalist,tt,ss,1,1,1,:)]./katu(soalist,tt,ss,1,1,1,:))']-[squeeze([katm(soalist,tt,ss,1,2,1,:)-katm(soalist,tt,ss,1,1,1,:)]./katm(soalist,tt,ss,1,1,1,:))'],[],'alpha',.05/7),end
  try [hh,pp]=ttest([squeeze([katu(soalist,tt,ss,1,2,2,:)-katu(soalist,tt,ss,1,1,2,:)]./katu(soalist,tt,ss,1,1,2,:))']-[squeeze([katm(soalist,tt,ss,1,2,2,:)-katm(soalist,tt,ss,1,1,2,:)]./katm(soalist,tt,ss,1,1,2,:))'],[],'alpha',.05/7),end
  try [hh,pp]=ttest([squeeze([katu(soalist,tt,ss,1,2,3,:)-katu(soalist,tt,ss,1,1,3,:)]./katu(soalist,tt,ss,1,1,3,:))']-[squeeze([katm(soalist,tt,ss,1,2,3,:)-katm(soalist,tt,ss,1,1,3,:)]./katm(soalist,tt,ss,1,1,3,:))'],[],'alpha',.05/7),end
  
  % interaction: of number
  try [hh,pp]=ttest([squeeze([katu(soalist,tt,ss,2,2,1,:)-katu(soalist,tt,ss,2,1,1,:)]./katu(soalist,tt,ss,2,1,1,:))']-[squeeze([katm(soalist,tt,ss,2,2,1,:)-katm(soalist,tt,ss,2,1,1,:)]./katm(soalist,tt,ss,2,1,1,:))'],[],'alpha',.05/7),end
  try [hh,pp]=ttest([squeeze([katu(soalist,tt,ss,2,2,2,:)-katu(soalist,tt,ss,2,1,2,:)]./katu(soalist,tt,ss,2,1,2,:))']-[squeeze([katm(soalist,tt,ss,2,2,2,:)-katm(soalist,tt,ss,2,1,2,:)]./katm(soalist,tt,ss,2,1,2,:))'],[],'alpha',.05/7),end
  try [hh,pp]=ttest([squeeze([katu(soalist,tt,ss,2,2,3,:)-katu(soalist,tt,ss,2,1,3,:)]./katu(soalist,tt,ss,2,1,3,:))']-[squeeze([katm(soalist,tt,ss,2,2,3,:)-katm(soalist,tt,ss,2,1,3,:)]./katm(soalist,tt,ss,2,1,3,:))'],[],'alpha',.05/7),end
  
  % just poststim, not subtracted by pre
  % testing for difference of [delta, kc, sw] for multsens-unisum/unisum, amplitude
  try [hh,pp]=ttest(squeeze([katm(soalist,tt,ss,1,2,1,:)-katu(soalist,tt,ss,1,2,1,:)]./katu(soalist,tt,ss,1,2,1,:))',[],'alpha',.05/7),end
  try [hh,pp]=ttest(squeeze([katm(soalist,tt,ss,1,2,2,:)-katu(soalist,tt,ss,1,2,2,:)]./katu(soalist,tt,ss,1,2,2,:))',[],'alpha',.05/7),end
  try [hh,pp]=ttest(squeeze([katm(soalist,tt,ss,1,2,3,:)-katu(soalist,tt,ss,1,2,3,:)]./katu(soalist,tt,ss,1,2,3,:))',[],'alpha',.05/7),end
  
  % testing for difference of [delta, kc, sw] for multsens-unisum/unisum, number
  try [hh,pp]=ttest(squeeze([katm(soalist,tt,ss,2,2,1,:)-katu(soalist,tt,ss,2,2,1,:)]./katu(soalist,tt,ss,2,2,1,:))',[],'alpha',.05/7),end
  try [hh,pp]=ttest(squeeze([katm(soalist,tt,ss,2,2,2,:)-katu(soalist,tt,ss,2,2,2,:)]./katu(soalist,tt,ss,2,2,2,:))',[],'alpha',.05/7),end
  try [hh,pp]=ttest(squeeze([katm(soalist,tt,ss,2,2,3,:)-katu(soalist,tt,ss,2,2,3,:)]./katu(soalist,tt,ss,2,2,3,:))',[],'alpha',.05/7),end
  
end



% export to SPSS
cd(sdir)
spss_exporter(squeeze([fstats_tmspn(3,[1 3 4 5 6 7 9],2,2,:); fstats_tpa(3,[1 3 4 5 6 7 9],2,2,:)]),{'tKc','MultSens','SOA'},0);  % no
spss_exporter(squeeze([fstats_tmspn(4,[1 3 4 5 6 7 9],2,2,:); fstats_tpa(4,[1 3 4 5 6 7 9],2,2,:)]),{'tDelta','MultSens','SOA'},0) % no
spss_exporter(squeeze([fstats_tmspn(5,[1 3 4 5 6 7 9],2,2,:); fstats_tpa(5,[1 3 4 5 6 7 9],2,2,:)]),{'tSW','MultSens','SOA'},0) % no
spss_exporter(squeeze([fstats_tmspn(6,[1 3 4 5 6 7 9],2,2,:); fstats_tpa(6,[1 3 4 5 6 7 9],2,2,:)]),{'tSpindleF','MultSens','SOA'},0) % possible soa and/or interaction
spss_exporter(squeeze([fstats_tmspn(7,[1 3 4 5 6 7 9],2,2,:); fstats_tpa(7,[1 3 4 5 6 7 9],2,2,:)]),{'tSpindleS','MultSens','SOA'},0) % error

spss1factor_exporter(squeeze(fstats_tmspn(3,[1 3 4 5 6 7 9],2,2,:)-fstats_tpa(3,[1 3 4 5 6 7 9],2,2,:)),{'tKcDiff','SOA','MultSens'},0) % no
spss1factor_exporter(squeeze(fstats_tmspn(4,[1 3 4 5 6 7 9],2,2,:)-fstats_tpa(4,[1 3 4 5 6 7 9],2,2,:)),{'tDeltaDiff','SOA','MultSens'},0) % no
spss1factor_exporter(squeeze(fstats_tmspn(5,[1 3 4 5 6 7 9],2,2,:)-fstats_tpa(5,[1 3 4 5 6 7 9],2,2,:)),{'tSWDiff','SOA','MultSens'},0) % no
spss1factor_exporter(squeeze(fstats_tmspn(6,[1 3 4 5 6 7 9],2,2,:)-fstats_tpa(6,[1 3 4 5 6 7 9],2,2,:)),{'tSpindleFDiff','SOA','MultSens'},0) % yes
spss1factor_exporter(squeeze(fstats_tmspn(7,[1 3 4 5 6 7 9],2,2,:)-fstats_tpa(7,[1 3 4 5 6 7 9],2,2,:)),{'tSpindleSDiff','SOA','MultSens'},0) % no

spss_exporter(squeeze([fstats_amspn(3,[1 3 4 5 6 7 9],2,2,:); fstats_apt(3,[1 3 4 5 6 7 9],2,2,:)]),{'aKc','MultSens','SOA'},0)
spss_exporter(squeeze([fstats_amspn(4,[1 3 4 5 6 7 9],2,2,:); fstats_apt(4,[1 3 4 5 6 7 9],2,2,:)]),{'aDelta','MultSens','SOA'},0)
spss_exporter(squeeze([fstats_amspn(5,[1 3 4 5 6 7 9],2,2,:); fstats_apt(5,[1 3 4 5 6 7 9],2,2,:)]),{'aSW','MultSens','SOA'},0)
spss_exporter(squeeze([fstats_amspn(6,[1 3 4 5 6 7 9],2,2,:); fstats_apt(6,[1 3 4 5 6 7 9],2,2,:)]),{'aSpindleF','MultSens','SOA'},0)
spss_exporter(squeeze([fstats_amspn(7,[1 3 4 5 6 7 9],2,2,:); fstats_apt(7,[1 3 4 5 6 7 9],2,2,:)]),{'aSpindleS','MultSens','SOA'},0)

spss1factor_exporter(squeeze(fstats_amspn(3,[1 3 4 5 6 7 9],2,2,:)-fstats_apt(3,[1 3 4 5 6 7 9],2,2,:)),{'aKcDiff','SOA','MultSens'},0)
spss1factor_exporter(squeeze(fstats_amspn(4,[1 3 4 5 6 7 9],2,2,:)-fstats_apt(4,[1 3 4 5 6 7 9],2,2,:)),{'aDeltaDiff','SOA','MultSens'},0)
spss1factor_exporter(squeeze(fstats_amspn(5,[1 3 4 5 6 7 9],2,2,:)-fstats_apt(5,[1 3 4 5 6 7 9],2,2,:)),{'aSWDiff','SOA','MultSens'},0)
spss1factor_exporter(squeeze(fstats_amspn(6,[1 3 4 5 6 7 9],2,2,:)-fstats_apt(6,[1 3 4 5 6 7 9],2,2,:)),{'aSpindleFDiff','SOA','MultSens'},0) % no
spss1factor_exporter(squeeze(fstats_amspn(7,[1 3 4 5 6 7 9],2,2,:)-fstats_apt(7,[1 3 4 5 6 7 9],2,2,:)),{'aSpindleSDiff','SOA','MultSens'},0)

spss_exporter(squeeze([mean(fstats_tmspn(3:5,[1 3 4 5 6 7 9],2,2,:),1); mean(fstats_tpa(3:5,[1 3 4 5 6 7 9],2,2,:),1)]),{'tBigWaves','MultSens','SOA'},0) % no
spss_exporter(squeeze([mean(fstats_tmspn(6:7,[1 3 4 5 6 7 9],2,2,:),1); mean(fstats_tpa(6:7,[1 3 4 5 6 7 9],2,2,:),1)]),{'tSpindleAll','MultSens','SOA'},0) % possible soa

spss1factor_exporter(squeeze(mean(fstats_tmspn(3:5,[1 3 4 5 6 7 9],2,2,:)-fstats_tpa(3:5,[1 3 4 5 6 7 9],2,2,:),1)),{'tBigWavesDiff','SOA','MultSens'},0) % no
spss1factor_exporter(squeeze(mean(fstats_tmspn(6:7,[1 3 4 5 6 7 9],2,2,:)-fstats_tpa(6:7,[1 3 4 5 6 7 9],2,2,:),1)),{'tSpindleAllDiff','SOA','MultSens'},0) % no

spss_exporter(squeeze([mean(fstats_amspn(3:5,[1 3 4 5 6 7 9],2,2,:),1); mean(fstats_apt(3:5,[1 3 4 5 6 7 9],2,2,:),1)]),{'aBigWaves','MultSens','SOA'},0)
spss_exporter(squeeze([mean(fstats_amspn(6:7,[1 3 4 5 6 7 9],2,2,:),1); mean(fstats_apt(6:7,[1 3 4 5 6 7 9],2,2,:),1)]),{'aSpindleAll','MultSens','SOA'},0)

spss1factor_exporter(squeeze(mean(fstats_amspn(3:5,[1 3 4 5 6 7 9],2,2,:)-fstats_apt(3:5,[1 3 4 5 6 7 9],2,2,:),1)),{'aBigWavesDiff','SOA','MultSens'},0)
spss1factor_exporter(squeeze(mean(fstats_amspn(6:7,[1 3 4 5 6 7 9],2,2,:)-fstats_apt(6:7,[1 3 4 5 6 7 9],2,2,:),1)),{'aSpindleAllDiff','SOA','MultSens'},0)

% test correlations of BigWaves and Spindles to behavioural outcomes
load([bdir 'rtgroup_pcb_diffms.mat'])
pcb=pcb(~isnan(pcb));
diffms=diffms(~isnan(diffms));

%collapsing over Bigwaves, diff of SOA456 to SOA19, split by medican of PCB
[hh,pp]=ttest2(squeeze([mean(mean(fstats_tmspn(3:5,[4 5 6],2,2,pcb>median(pcb)),2),1) - mean(mean(fstats_tmspn(3:5,[1 9],2,2,pcb>median(pcb)),2),1)]),squeeze([mean(mean(fstats_tmspn(3:5,[4 5 6],2,2,pcb<median(pcb)),2),1) - mean(mean(fstats_tmspn(3:5,[1 9],2,2,pcb<median(pcb)),2),1)]))
[hh,pp]=ttest2(squeeze([mean(mean(fstats_tmspn(6:7,[4 5 6],2,2,pcb>median(pcb)),2),1) - mean(mean(fstats_tmspn(6:7,[1 9],2,2,pcb>median(pcb)),2),1)]),squeeze([mean(mean(fstats_tmspn(6:7,[4 5 6],2,2,pcb<median(pcb)),2),1) - mean(mean(fstats_tmspn(6:7,[1 9],2,2,pcb<median(pcb)),2),1)]))

%collapsing over Bigwaves, diff of SOA456 to SOA19, split by medican of DiffMs
[hh,pp]=ttest2(squeeze([mean(mean(fstats_tmspn(3:5,[4 5 6],2,2,diffms>median(diffms)),2),1) - mean(mean(fstats_tmspn(3:5,[1 9],2,2,diffms>median(diffms)),2),1)]),squeeze([mean(mean(fstats_tmspn(3:5,[4 5 6],2,2,diffms<median(diffms)),2),1) - mean(mean(fstats_tmspn(3:5,[1 9],2,2,diffms<median(diffms)),2),1)]))
[hh,pp]=ttest2(squeeze([mean(mean(fstats_tmspn(6:7,[4 5 6],2,2,diffms>median(diffms)),2),1) - mean(mean(fstats_tmspn(6:7,[1 9],2,2,diffms>median(diffms)),2),1)]),squeeze([mean(mean(fstats_tmspn(6:7,[4 5 6],2,2,diffms<median(diffms)),2),1) - mean(mean(fstats_tmspn(6:7,[1 9],2,2,diffms<median(diffms)),2),1)]))


