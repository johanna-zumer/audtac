% preprocessing of EEG data from Hills, 64ch MR-cap
clear all
close all
if ispc
  edir='D:\audtac\eeg_data\';
  esdir='D:\audtac\source_data\';
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
    esdir='/home/zumerj/audtac/source_data/';
    ddir='/home/zumerj/audtac/legomagic/diaries/';
    bdir='/home/zumerj/audtac/behav_data/';
    sdir='/home/zumerj/audtac/spss_stuff/';
    fdir='/home/zumerj/audtac/figs/';
    mdir='/home/zumerj/audtac/structural_MRI/';
    pdir='/home/zumerj/audtac/polhemus/';
  else % assume on VM linux of psychl-132432
    edir='/mnt/hgfs/D/audtac/eeg_data/';
    esdir='/mnt/hgfs/D/audtac/source_data/';
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
soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];
iiBuse=setdiff(iiBuse,3:7);


%%

minnumcomb=1; % in effect just the normal combination of trials
minnumcomb=2; % don't let minnumcomb go higher than 1000.
% minnumcomb=100; % don't let minnumcomb go higher than 1000.

timestepfft=0.02;
timestepplv=0.01;

dofftadd=0;
usetr=2; % see inside eeg_legomagic_trialSelection2_wakeSleep_sepTacAud for what this does/means
synchasynch=0;
use23=0;
phaset0=0;

% for sleep=0
sleep=0;
if sleep
  iiuse=iiBuse;
  %     iiuse=[32];
  iteruse=11;
  trialkc=-1;  % vary this from -1, 0, and 1
else
  iiuse=iiSuse;
  %         iiuse=setdiff(iiSuse,1:30);
  iteruse=31;
  trialkc=-1;
end
for ii=iiuse;
% for ii=setdiff(iiuse,1:8)
  cd([edir sub{ii} ])
  clearvars -except ii sub *dir ii*use sleep minnumcomb hostname timestep* soades dofftadd statst use* synchasynch iteruse trialkc phaset0
  [raw_tac, raw_aud, raw_nul, artflag]=eeg_legomagic_epoching2(ii,sleep,1,0); % featfull=1, saveflag =0;
  
  %% Do TacPlusAud first
  tacaud=1;
  %   for tt=[3]
  tt=3;
  for iter=iteruse
    %       [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection2_wakeSleep_sepTacAud(ii,sleep,tt,tacaud);
    [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection2_wakeSleep_sepTacAud(ii,sleep,tt,tacaud,raw_tac, raw_aud, raw_nul, iter, usetr, trialkc);
    
    if sleep==0
      tlock_aud{10,tt,12}=[];
      tlock_aud{10,tt,13}=[];
    end
    
    if usetr==2
      if length(iteruse)>1
        error('only one iteruse at a time for usetr2')
      end
      iter=tr.iter;
    end
    
    %  Frequency analysis
    %         if ii<8
    %           soalist=[3 4 5 6 7];
    %         else
    soalist=[1 3 4 5 6 7 9];
    %         end
    if sleep
      if use23
        ssuse=[tr.stageuse+10 23];
      else
        ssuse=[tr.stageuse+10];
      end
    else
      ssuse=tr.stageuse+10;
      warning('remove me later')
      ssuse=10;
    end
    
    for ss=ssuse
      if ~any([tr.t10trialkept{:,tt,ss}])
        ssuse=setdiff(ssuse,ss);
        continue
      else
        
        try
          fsample=1/diff(tlock_tac{10,tt,ss}.time(1:2)); % sampling rate
        catch
          try
            fsample=1/diff(tlock_aud{10,tt,ss}.time(1:2)); % sampling rate
          catch
            fsample=1/diff(tlock_nul{10,tt,ss}.time(1:2)); % sampling rate
          end
        end
        %   for tt=1:4
        
        if synchasynch
          % Multisensory-shifted comparison
          for ll=[1 3 4] % need to do this before other loops as the fields get deleted as it loops through ll
            if ss==23
              if length(tr.tmstrialkept1{ll,tt,12})>=2 && length(tr.tmstrialkept1{ll,tt,13})<2
                cfg=[];
                cfg.trials=tr.tmstrialkept1{ll,tt,12};
                tlock_tacMSshift1{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{ll,tt,12});
                cfg=[];
                cfg.trials=tr.tmstrialkept2{ll,tt,12};
                tlock_tacMSshift2{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{10-ll,tt,12});
                cfg=[];
                cfg.trials=tr.tmstrialkept3{ll,tt,12};
                tlock_tacMSshift3{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{5,tt,12});
                cfg=[];
                cfg.trials=tr.tmstrialkept4{ll,tt,12};
                tlock_tacMSshift4{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{5,tt,12});
              elseif length(tr.tmstrialkept1{ll,tt,12})<2 && length(tr.tmstrialkept1{ll,tt,13})>=2
                cfg=[];
                cfg.trials=tr.tmstrialkept1{ll,tt,13};
                tlock_tacMSshift1{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{ll,tt,13});
                cfg=[];
                cfg.trials=tr.tmstrialkept2{ll,tt,13};
                tlock_tacMSshift2{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{10-ll,tt,13});
                cfg=[];
                cfg.trials=tr.tmstrialkept3{ll,tt,13};
                tlock_tacMSshift3{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{5,tt,13});
                cfg=[];
                cfg.trials=tr.tmstrialkept4{ll,tt,13};
                tlock_tacMSshift4{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{5,tt,13});
              elseif length(tr.tmstrialkept1{ll,tt,12})>=2 && length(tr.tmstrialkept1{ll,tt,13})>=2
                cfg=[];
                cfg.trials=tr.tmstrialkept1{ll,tt,12};
                tmp12=ft_selectdata(cfg,tlock_tac{ll,tt,12});
                cfg=[];
                cfg.trials=tr.tmstrialkept1{ll,tt,13};
                tmp13=ft_selectdata(cfg,tlock_tac{ll,tt,13});
                cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
                cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
                tlock_tacMSshift1{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
                cfg=[];
                cfg.trials=tr.tmstrialkept2{ll,tt,12};
                tmp12=ft_selectdata(cfg,tlock_tac{10-ll,tt,12});
                cfg=[];
                cfg.trials=tr.tmstrialkept2{ll,tt,13};
                tmp13=ft_selectdata(cfg,tlock_tac{10-ll,tt,13});
                cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
                cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
                tlock_tacMSshift2{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
                cfg=[];
                cfg.trials=tr.tmstrialkept3{ll,tt,12};
                tmp12=ft_selectdata(cfg,tlock_tac{5,tt,12});
                cfg=[];
                cfg.trials=tr.tmstrialkept3{ll,tt,13};
                tmp13=ft_selectdata(cfg,tlock_tac{5,tt,13});
                cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
                cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
                tlock_tacMSshift3{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
                cfg=[];
                cfg.trials=tr.tmstrialkept4{ll,tt,12};
                tmp12=ft_selectdata(cfg,tlock_tac{5,tt,12});
                cfg=[];
                cfg.trials=tr.tmstrialkept4{ll,tt,13};
                tmp13=ft_selectdata(cfg,tlock_tac{5,tt,13});
                cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
                cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
                tlock_tacMSshift4{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
              else
                %               tlock_tac{ll,tt,ss}=[];
                continue
              end
            else
              if length(tr.tmstrialkept1{ll,tt,ss})>=2 && length(tr.tmstrialkept2{ll,tt,ss})>=2 && length(tr.tmstrialkept3{ll,tt,ss})>=2 && length(tr.tmstrialkept4{ll,tt,ss})>=2
                cfg=[];
                cfg.trials=tr.tmstrialkept1{ll,tt,ss};
                tlock_tacMSshift1{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{ll,tt,ss});
                cfg=[];
                cfg.trials=tr.tmstrialkept2{ll,tt,ss};
                tlock_tacMSshift2{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{10-ll,tt,ss});
                cfg=[];
                cfg.trials=tr.tmstrialkept3{ll,tt,ss};
                tlock_tacMSshift3{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{5,tt,ss});
                cfg=[];
                cfg.trials=tr.tmstrialkept4{ll,tt,ss};
                tlock_tacMSshift4{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{5,tt,ss});
              else
                tlock_tacMSshift1{ll,tt,ss}=[];
                tlock_tacMSshift2{ll,tt,ss}=[];
                tlock_tacMSshift3{ll,tt,ss}=[];
                tlock_tacMSshift4{ll,tt,ss}=[];
              end
            end
            % do time-shifting here
            if ~isempty(tlock_tacMSshift1{ll,tt,ss})
              warning off
              cfg=[];
              cfg.offset=round(fsample*(-soades(ll)));
              tlock_tacMSshift1{ll,tt,ss}=ft_redefinetrial(cfg,tlock_tacMSshift1{ll,tt,ss}); % Aud first shifted so that Aud at time 0
              tlock_tacMSshift4{ll,tt,ss}=ft_redefinetrial(cfg,tlock_tacMSshift4{ll,tt,ss}); % Simult shifted so that both at time of second stim of 'll' condition
              warning on
            end
          end %ll
        end
        
        
        for ll=[soalist soalist+20 soalist+40]
          
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
                if ispc
                  S=memory; % memory doesn't exist on linux
                  gbavail=S.MemAvailableAllArrays/(1024^3);
                else
                  [aa,bb]=system('head /proc/meminfo');
                  if ~isempty(strfind(hostname,'les')) % les-linux-fs3
                    gbavail=str2num(bb(strfind(bb,'MemFree')+8:strfind(bb,'MemFree')+23))/(1024^2);
                  elseif ~isempty(strfind(hostname,'LES')) % COLLES-151401
                    gb1=str2num(bb(strfind(bb,'Inactive:')+9:strfind(bb,'Inactive:')+23))/(1024^2);
                    gb2=str2num(bb(strfind(bb,'MemFree')+8:strfind(bb,'MemFree')+23))/(1024^2);
                    gbavail=max(gb1,gb2);
                  end
                end
                if gbavail < 1.5
                  save tmptlock.mat tlock_aud -v7.3
                  clear tlock_aud
                end
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
                numt_trials(ll,tt,ss)=0;
                continue
              end
              if ll==max(soalist)
                tlock_tac{10,tt,12}=[];
                tlock_tac{10,tt,13}=[];
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
                if ispc
                  S=memory; % memory doesn't exist on linux
                  gbavail=S.MemAvailableAllArrays/(1024^3);
                else
                  [aa,bb]=system('head /proc/meminfo');
                  if ~isempty(strfind(hostname,'les')) % les-linux-fs3
                    gbavail=str2num(bb(strfind(bb,'MemFree')+8:strfind(bb,'MemFree')+23))/(1024^2);
                  elseif ~isempty(strfind(hostname,'LES')) % COLLES-151401
                    gb1=str2num(bb(strfind(bb,'Inactive:')+9:strfind(bb,'Inactive:')+23))/(1024^2);
                    gb2=str2num(bb(strfind(bb,'MemFree')+8:strfind(bb,'MemFree')+23))/(1024^2);
                    gbavail=max(gb1,gb2);
                  end
                end
                if gbavail < 1.5
                  save tmptlock.mat tlock_aud -v7.3
                  clear tlock_aud
                end
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
                if ll==max(soalist) && ( (use23 && ss<12) || ~use23) && ~phaset0
                  tlock_tac{10,tt,ss}=[];
                end
                if ll==max(soalist) && ( (use23 && ss<12) || ~use23) && ~phaset0
                  tlock_aud{10,tt,ss}=[];
                end
                
                cfg=[];
                cfg.trials=tr.nllttrialkept{ll,tt,ss};
                if length(cfg.trials)<2
                  tlock_nul{ll+50,tt,ss}=[];
                else
                  tlock_nul{ll+50,tt,ss}=ft_selectdata(cfg,tlock_nul{10,tt,ss});
                end
                if ll==max(soalist) && ( (use23 && ss<12) || ~use23)
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
          end
          if ~exist('tlock_aud','var')
            load tmptlock.mat
            delete tmptlock.mat
          end
          
          if ss==23 && ll<10
            tlock_aud{ll,tt,12}=[];
            tlock_aud{ll,tt,13}=[];
            if length(tr.tlltrialkept{ll,tt,12})>=2 && length(tr.tlltrialkept{ll,tt,13})<2
              tlock_tac{ll,tt,ss}=tlock_tac{ll,tt,12};
            elseif length(tr.tlltrialkept{ll,tt,12})<2 && length(tr.tlltrialkept{ll,tt,13})>=2
              tlock_tac{ll,tt,ss}=tlock_tac{ll,tt,13};
            elseif length(tr.tlltrialkept{ll,tt,12})>=2 && length(tr.tlltrialkept{ll,tt,13})>=2
              if ispc
                S=memory; % memory doesn't exist on linux
                gbavail=S.MemAvailableAllArrays/(1024^3);
              else
                [aa,bb]=system('head /proc/meminfo');
                if ~isempty(strfind(hostname,'les')) % les-linux-fs3
                  gbavail=str2num(bb(strfind(bb,'MemFree')+8:strfind(bb,'MemFree')+23))/(1024^2);
                elseif ~isempty(strfind(hostname,'LES')) % COLLES-151401
                  gb1=str2num(bb(strfind(bb,'Inactive:')+9:strfind(bb,'Inactive:')+23))/(1024^2);
                  gb2=str2num(bb(strfind(bb,'MemFree')+8:strfind(bb,'MemFree')+23))/(1024^2);
                  gbavail=max(gb1,gb2);
                end
              end
              if gbavail < 1.5
                save tmptlock.mat tlock_nul -v7.3
                clear tlock_nul
              end
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
          end
          
          if ss==23 && ll>40
            if length(tr.all40trialkept{ll-40,tt,12})>=2 && length(tr.all40trialkept{ll-40,tt,13})<2
              tlock_aud{ll,tt,ss}=tlock_aud{ll,tt,12};
            elseif length(tr.all40trialkept{ll-40,tt,12})<2 && length(tr.all40trialkept{ll-40,tt,13})>=2
              tlock_aud{ll,tt,ss}=tlock_aud{ll,tt,13};
            elseif length(tr.all40trialkept{ll-40,tt,12})>=2 && length(tr.all40trialkept{ll-40,tt,13})>=2
              if ispc
                S=memory; % memory doesn't exist on linux
                gbavail=S.MemAvailableAllArrays/(1024^3);
              else
                [aa,bb]=system('head /proc/meminfo');
                if ~isempty(strfind(hostname,'les')) % les-linux-fs3
                  gbavail=str2num(bb(strfind(bb,'MemFree')+8:strfind(bb,'MemFree')+23))/(1024^2);
                elseif ~isempty(strfind(hostname,'LES')) % COLLES-151401
                  gb1=str2num(bb(strfind(bb,'Inactive:')+9:strfind(bb,'Inactive:')+23))/(1024^2);
                  gb2=str2num(bb(strfind(bb,'MemFree')+8:strfind(bb,'MemFree')+23))/(1024^2);
                  gbavail=max(gb1,gb2);
                end
              end
              if gbavail < 1.5
                save tmptlock.mat tlock_nul -v7.3
                clear tlock_nul
              end
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
            if ~isempty(tlock_aud{ll,tt,ss})
              tlock_aud{ll,tt,ss}=trim_nans(tlock_aud{ll,tt,ss});
            end
          end
          
          if ~exist('tlock_nul','var')
            load tmptlock.mat
            delete tmptlock.mat
          end
          
          %up to here, same as for ERP
        end % ll
        
        %%
        
        % Do it two ways:
        
        if dofftadd
          % %%%%      % FFT then add    %%%%%%%%%%%%%
          % Way 1:  FFT+norm of each condition, then add conditions
          
          for ll=[soalist soalist+20]
            %           if 0
            if isfield(tlock_tac{ll,tt,ss},'trial') && size(tlock_tac{ll,tt,ss}.trial,1)
              %         cfg=[];
              %         cfg.lpfilter='yes';
              %         cfg.lpfreq=40;
              %         cfg.demean='yes';
              %         cfg.baselinewindow=[-1.7 -0.6];
              
              % Lower frequencies with Hanning taper
              cfg=[];
              cfg.method='mtmconvol';
              cfg.pad=4;
              cfg.foi=4:2:30;
              cfg.taper='hanning';
              cfg.toi=-1.2:timestepfft:1.3;
              cfg.t_ftimwin=4./cfg.foi;
              cfg.output='powandcsd';  % doesn't make sense to look at PLV for each condition separately in 'way 1'
              % cfg.output='fourier'; % fourier automatically sets keeptrials=yes
              %         cfg.keeptrials='yes';  % not necessary at this stage for keeping trials (yes later for source analysis and phase resetting)
              freqlo_tac{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tac{ll,tt,ss});
              % freqlo_tac{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqlo_tac{ll,tt,ss}.fourierspctrm).^2,1));
              
              % Higher frequencies with dpss multitaper
              cfg=[];
              cfg.method='mtmconvol';
              cfg.pad=4;
              cfg.foi=[30:5:45 55:5:80];
              cfg.taper='dpss';
              cfg.tapsmofrq=7*ones(1,length(cfg.foi));
              cfg.toi=-1.2:timestepfft:1.3;
              cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
              cfg.output='powandcsd';
              freqhi_tac{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tac{ll,tt,ss});
            else
              freqlo_tac{ll,tt,ss}=[];
              freqhi_tac{ll,tt,ss}=[];
            end
          end
          %       tlock_tac{ll,tt,ss}=[]; % clearing to help save memory as running
          for ll=soalist+40
            if numt_trials(ll-40,tt,ss) && isfield(tlock_aud{ll,tt,ss},'trial') && size(tlock_aud{ll,tt,ss}.trial,1)
              cfg=[];
              cfg.method='mtmconvol';
              cfg.pad=4;
              cfg.foi=4:2:30;
              cfg.taper='hanning';
              cfg.toi=-1.2:timestepfft:1.3;
              cfg.t_ftimwin=4./cfg.foi;
              cfg.output='powandcsd';
              freqlo_aud{ll,tt,ss}=ft_freqanalysis(cfg,tlock_aud{ll,tt,ss});
              cfg=[];
              cfg.method='mtmconvol';
              cfg.pad=4;
              cfg.foi=[30:5:45 55:5:80];
              cfg.taper='dpss';
              cfg.tapsmofrq=7*ones(1,length(cfg.foi));
              cfg.toi=-1.2:timestepfft:1.3;
              cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
              cfg.output='powandcsd';
              freqhi_aud{ll,tt,ss}=ft_freqanalysis(cfg,tlock_aud{ll,tt,ss});
            else
              freqlo_aud{ll,tt,ss}=[];
              freqhi_aud{ll,tt,ss}=[];
            end
            %       tlock_aud{ll,tt,ss}=[]; % clearing to help save memory as running
            %           end
          end
          for ll=[soalist+50]
            if numt_trials(ll-50,tt,ss) && isfield(tlock_nul{ll,tt,ss},'trial') && size(tlock_nul{ll,tt,ss}.trial,1)
              cfg=[];
              cfg.method='mtmconvol';
              cfg.pad=4;
              cfg.foi=4:2:30;
              cfg.taper='hanning';
              cfg.toi=-1.2:timestepfft:1.3;
              cfg.t_ftimwin=4./cfg.foi;
              cfg.output='powandcsd';
              freqlo_nul{ll,tt,ss}=ft_freqanalysis(cfg,tlock_nul{ll,tt,ss});
              cfg=[];
              cfg.method='mtmconvol';
              cfg.pad=4;
              cfg.foi=[30:5:45 55:5:80];
              cfg.taper='dpss';
              cfg.tapsmofrq=7*ones(1,length(cfg.foi));
              cfg.toi=-1.2:timestepfft:1.3;
              cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
              cfg.output='powandcsd';
              freqhi_nul{ll,tt,ss}=ft_freqanalysis(cfg,tlock_nul{ll,tt,ss});
            else
              freqlo_nul{ll,tt,ss}=[];
              freqhi_nul{ll,tt,ss}=[];
            end
            %       tlock_nul{ll,tt,ss}=[]; % clearing to help save memory as running
          end
          
          for ll=soalist
            
            % create sum of unisensory conditions
            if numt_trials(ll,tt,ss)
              cfg=[];
              cfg.operation='add';
              cfg.parameter='powspctrm';
              %            cfg.parameter='fourierspctrm';
              %       freqlo_tacPaud{ll,tt,ss}=ft_math(cfg,ft_selectdata(cfgavg,freqlo_tac{10}),ft_selectdata(cfgavg,freqlo_aud{ll+40,tt,ss})); % keeping tac centred but jitter aud (to match tlock_tac{X,tt,ss})
              %       freqlo_audPtac{ll,tt,ss}=ft_math(cfg,ft_selectdata(cfgavg,freqlo_aud{10}),ft_selectdata(cfgavg,freqlo_tac{ll+40,tt,ss})); % keeping aud centred but jitter tac (to match tlock_aud{X,tt,ss})
              %       freqhi_tacPaud{ll,tt,ss}=ft_math(cfg,ft_selectdata(cfgavg,freqhi_tac{10}),ft_selectdata(cfgavg,freqhi_aud{ll+40,tt,ss})); % keeping tac centred but jitter aud (to match tlock_tac{X,tt,ss})
              %       freqhi_audPtac{ll,tt,ss}=ft_math(cfg,ft_selectdata(cfgavg,freqhi_aud{10}),ft_selectdata(cfgavg,freqhi_tac{ll+40,tt,ss})); % keeping aud centred but jitter tac (to match tlock_aud{X,tt,ss})
              freqlo_tacPaud_fftadd{ll,tt,ss}=ft_math(cfg,freqlo_tac{ll+20,tt,ss},freqlo_aud{ll+40,tt,ss}); % keeping tac centred but jitter aud (to match tlock_tac{X,tt,ss})
              freqhi_tacPaud_fftadd{ll,tt,ss}=ft_math(cfg,freqhi_tac{ll+20,tt,ss},freqhi_aud{ll+40,tt,ss}); % keeping tac centred but jitter aud (to match tlock_tac{X,tt,ss})
              freqlo_tTacAlone_fftadd{ll,tt,ss}=freqlo_tac{ll+20,tt,ss};
              freqlo_tAudAlone_fftadd{ll,tt,ss}=freqlo_aud{ll+40,tt,ss};
              freqhi_tTacAlone_fftadd{ll,tt,ss}=freqhi_tac{ll+20,tt,ss};
              freqhi_tAudAlone_fftadd{ll,tt,ss}=freqhi_aud{ll+40,tt,ss};
              cfg.parameter='crsspctrm';
              tmp=ft_math(cfg,freqlo_tac{ll+20,tt,ss},freqlo_aud{ll+40,tt,ss}); % keeping tac centred but jitter aud (to match tlock_tac{X,tt,ss})
              freqlo_tacPaud_fftadd{ll,tt,ss}.crsspctrm=tmp.crsspctrm;
              freqlo_tTacAlone_fftadd{ll,tt,ss}.crsspctrm=freqlo_tac{ll+20,tt,ss}.crsspctrm;
              freqlo_tAudAlone_fftadd{ll,tt,ss}.crsspctrm=freqlo_aud{ll+40,tt,ss}.crsspctrm;
              tmp=ft_math(cfg,freqhi_tac{ll+20,tt,ss},freqhi_aud{ll+40,tt,ss}); % keeping tac centred but jitter aud (to match tlock_tac{X,tt,ss})
              freqhi_tacPaud_fftadd{ll,tt,ss}.crsspctrm=tmp.crsspctrm;
              freqhi_tTacAlone_fftadd{ll,tt,ss}.crsspctrm=freqhi_tac{ll+20,tt,ss}.crsspctrm;
              freqhi_tAudAlone_fftadd{ll,tt,ss}.crsspctrm=freqhi_aud{ll+40,tt,ss}.crsspctrm;
            else
              freqlo_tacPaud_fftadd{ll,tt,ss}=[];
              freqhi_tacPaud_fftadd{ll,tt,ss}=[];
              freqlo_tTacAlone_fftadd{ll,tt,ss}=[];
              freqhi_tTacAlone_fftadd{ll,tt,ss}=[];
              freqlo_tAudAlone_fftadd{ll,tt,ss}=[];
              freqhi_tAudAlone_fftadd{ll,tt,ss}=[];
            end
            
            
            % create sum of mulitsensory and null conditions
            %       freqlo_tacMSpN{ll,tt,ss}=ft_math(cfg,ft_selectdata(cfgavg,freqlo_tac{ll,tt,ss}),ft_selectdata(cfgavg,freqlo_nul{ll+40,tt,ss})); % keeping tac centred but jitter aud (to match tlock_tac{X,tt,ss})
            %       freqlo_audMSpN{ll,tt,ss}=ft_math(cfg,ft_selectdata(cfgavg,freqlo_aud{ll,tt,ss}),ft_selectdata(cfgavg,freqlo_nul{ll+40,tt,ss})); % keeping aud centred but jitter tac (to match tlock_aud{X,tt,ss})
            %       freqhi_tacMSpN{ll,tt,ss}=ft_math(cfg,ft_selectdata(cfgavg,freqhi_tac{ll,tt,ss}),ft_selectdata(cfgavg,freqhi_nul{ll+40,tt,ss})); % keeping tac centred but jitter aud (to match tlock_tac{X,tt,ss})
            %       freqhi_audMSpN{ll,tt,ss}=ft_math(cfg,ft_selectdata(cfgavg,freqhi_aud{ll,tt,ss}),ft_selectdata(cfgavg,freqhi_nul{ll+40,tt,ss})); % keeping aud centred but jitter tac (to match tlock_aud{X,tt,ss})
            if numt_trials(ll,tt,ss)
              cfg=[];
              cfg.operation='add';
              cfg.parameter='powspctrm';
              freqlo_tacMSpN_fftadd{ll,tt,ss}=ft_math(cfg,freqlo_tac{ll,tt,ss},freqlo_nul{ll+50,tt,ss}); % keeping tac centred but jitter aud (to match tlock_tac{X,tt,ss})
              freqhi_tacMSpN_fftadd{ll,tt,ss}=ft_math(cfg,freqhi_tac{ll,tt,ss},freqhi_nul{ll+50,tt,ss}); % keeping tac centred but jitter aud (to match tlock_tac{X,tt,ss})
              freqlo_tMSAlone_fftadd{ll,tt,ss}=freqlo_tac{ll,tt,ss};
              freqhi_tMSAlone_fftadd{ll,tt,ss}=freqhi_tac{ll,tt,ss};
              freqlo_tNulAlone_fftadd{ll,tt,ss}=freqlo_nul{ll+50,tt,ss};
              freqhi_tNulAlone_fftadd{ll,tt,ss}=freqhi_nul{ll+50,tt,ss};
              cfg.parameter='crsspctrm';
              tmp=ft_math(cfg,freqlo_tac{ll,tt,ss},freqlo_nul{ll+50,tt,ss}); % keeping tac centred but jitter aud (to match tlock_tac{X,tt,ss})
              freqlo_tacMSpN_fftadd{ll,tt,ss}.crsspctrm=tmp.crsspctrm;
              freqlo_tMSAlone_fftadd{ll,tt,ss}.crsspctrm=freqlo_tac{ll,tt,ss}.crsspctrm;
              freqlo_tNulAlone_fftadd{ll,tt,ss}.crsspctrm=freqlo_nul{ll+50,tt,ss}.crsspctrm;
              tmp=ft_math(cfg,freqhi_tac{ll,tt,ss},freqhi_nul{ll+50,tt,ss}); % keeping tac centred but jitter aud (to match tlock_tac{X,tt,ss})
              freqhi_tacMSpN_fftadd{ll,tt,ss}.crsspctrm=tmp.crsspctrm;
              freqhi_tMSAlone_fftadd{ll,tt,ss}.crsspctrm=freqhi_tac{ll,tt,ss}.crsspctrm;
              freqhi_tNulAlone_fftadd{ll,tt,ss}.crsspctrm=freqhi_nul{ll+50,tt,ss}.crsspctrm;
            else
              freqlo_tacMSpN_fftadd{ll,tt,ss}=[];
              freqhi_tacMSpN_fftadd{ll,tt,ss}=[];
              freqlo_tMSAlone_fftadd{ll,tt,ss}=[];
              freqhi_tMSAlone_fftadd{ll,tt,ss}=[];
              freqlo_tNulAlone_fftadd{ll,tt,ss}=[];
              freqhi_tNulAlone_fftadd{ll,tt,ss}=[];
            end
          end % ll
          
          if synchasynch
            for ll=[1 3 4]
              if ~isempty(tlock_tacMSshift1{ll,tt,ss}) && size(tlock_tacMSshift1{ll,tt,ss}.trial,1)
                cfg=[];
                cfg.method='mtmconvol';
                cfg.pad=4;
                cfg.foi=4:2:30;
                cfg.taper='hanning';
                cfg.toi=-1.2:timestepfft:1.3;
                cfg.t_ftimwin=4./cfg.foi;
                cfg.output='powandcsd';  % doesn't make sense to look at PLV for each condition separately in 'way 1'
                freqlo_tacMSshift1{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tacMSshift1{ll,tt,ss});
                freqlo_tacMSshift2{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tacMSshift2{ll,tt,ss});
                freqlo_tacMSshift3{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tacMSshift3{ll,tt,ss});
                freqlo_tacMSshift4{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tacMSshift4{ll,tt,ss});
                
                % Higher frequencies with dpss multitaper
                cfg=[];
                cfg.method='mtmconvol';
                cfg.pad=4;
                cfg.foi=[30:5:45 55:5:80];
                cfg.taper='dpss';
                cfg.tapsmofrq=7*ones(1,length(cfg.foi));
                cfg.toi=-1.2:timestepfft:1.3;
                cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
                cfg.output='powandcsd';
                freqhi_tacMSshift1{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tacMSshift1{ll,tt,ss});
                freqhi_tacMSshift2{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tacMSshift2{ll,tt,ss});
                freqhi_tacMSshift3{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tacMSshift3{ll,tt,ss});
                freqhi_tacMSshift4{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tacMSshift4{ll,tt,ss});
              else
                freqlo_tacMSshift1{ll,tt,ss}=[];
                freqlo_tacMSshift2{ll,tt,ss}=[];
                freqlo_tacMSshift3{ll,tt,ss}=[];
                freqlo_tacMSshift4{ll,tt,ss}=[];
                freqhi_tacMSshift1{ll,tt,ss}=[];
                freqhi_tacMSshift2{ll,tt,ss}=[];
                freqhi_tacMSshift3{ll,tt,ss}=[];
                freqhi_tacMSshift4{ll,tt,ss}=[];
              end
              
              if ~isempty(freqlo_tacMSshift1{ll,tt,ss})
                % create sum of asynchronous multisensory conditions (e.g. AT70_shifted plus TA70)
                cfg=[];
                cfg.operation='add';
                cfg.parameter='powspctrm';
                freqlo_tMSasynch_fftadd{ll,tt,ss}=ft_math(cfg,freqlo_tacMSshift1{ll,tt,ss},freqlo_tacMSshift2{ll,tt,ss});
                freqhi_tMSasynch_fftadd{ll,tt,ss}=ft_math(cfg,freqhi_tacMSshift1{ll,tt,ss},freqhi_tacMSshift2{ll,tt,ss});
                cfg.parameter='crsspctrm';
                tmp=ft_math(cfg,freqlo_tacMSshift1{ll,tt,ss},freqlo_tacMSshift2{ll,tt,ss});
                freqlo_tMSasynch_fftadd{ll,tt,ss}.crsspctrm=tmp.crsspctrm;
                tmp=ft_math(cfg,freqhi_tacMSshift1{ll,tt,ss},freqhi_tacMSshift2{ll,tt,ss});
                freqhi_tMSasynch_fftadd{ll,tt,ss}.crsspctrm=tmp.crsspctrm;
                
                % create sum of synchronous multisensory conditions (e.g. TA0_shifted plus TA0)
                cfg=[];
                cfg.operation='add';
                cfg.parameter='powspctrm';
                freqlo_tMSsynch_fftadd{ll,tt,ss}=ft_math(cfg,freqlo_tacMSshift3{ll,tt,ss},freqlo_tacMSshift4{ll,tt,ss});
                freqhi_tMSsynch_fftadd{ll,tt,ss}=ft_math(cfg,freqhi_tacMSshift3{ll,tt,ss},freqhi_tacMSshift4{ll,tt,ss});
                cfg.parameter='crsspctrm';
                tmp=ft_math(cfg,freqlo_tacMSshift3{ll,tt,ss},freqlo_tacMSshift4{ll,tt,ss});
                freqlo_tMSsynch_fftadd{ll,tt,ss}.crsspctrm=tmp.crsspctrm;
                tmp=ft_math(cfg,freqhi_tacMSshift3{ll,tt,ss},freqhi_tacMSshift4{ll,tt,ss});
                freqhi_tMSsynch_fftadd{ll,tt,ss}.crsspctrm=tmp.crsspctrm;
              else
                freqlo_tMSasynch_fftadd{ll,tt,ss}=[];
                freqlo_tMSsynch_fftadd{ll,tt,ss}=[];
                freqhi_tMSasynch_fftadd{ll,tt,ss}=[];
                freqhi_tMSsynch_fftadd{ll,tt,ss}=[];
              end
            end % ll
          end
          
        end %dofftadd
        
        %%
        % Way 2: Add first then FFT
        
        % Done two ways:
        %
        % 2a) freqlo_tacPaud_comb{ll,tt,ss,1} is the result from adding
        % the timelock unaveraged data first, then abs(FFT), with trial
        % pairings as they are originally ordered.
        %
        % 2b) And freqlo_tacPaud_comb{ll,tt,ss,2} is also created from adding
        % the timelock unaveraged data first, then abs(FFT) but with
        % trial pairings to be added randomised and the full distribution
        % computed, followed at last by averaging over all the
        % realisations of this distribution.
        %
        
        if phaset0 && sleep==0
          % For PCA of data, hopefully picking up channels greatly contribution to evoked response
          % This is not dependent on 'll'
          cfg=[];
          cfg.offset=500; % samples
          tmp_ms1=ft_redefinetrial(cfg,tlock_tac{1,tt,ss});
          cfg.offset=70; % samples
          tmp_ms3=ft_redefinetrial(cfg,tlock_tac{3,tt,ss});
          cfg.offset=20; % samples
          tmp_ms4=ft_redefinetrial(cfg,tlock_tac{4,tt,ss});
          cfg=[];
          cfg.latency=[-.1 .8];
          appseldata=ft_selectdata(cfg,ft_appenddata([],tlock_tac{10,tt,ss},tlock_aud{10,tt,ss},tmp_ms1,tmp_ms3,tmp_ms4,tlock_tac{5,tt,ss},tlock_tac{6,tt,ss},tlock_tac{7,tt,ss},tlock_tac{9,tt,ss}));
          cfg=[];
          cfg.method='pca';
          pcatacaud=ft_componentanalysis(cfg,appseldata);
          pcatacaud_tlock=ft_timelockanalysis([],pcatacaud);
          appseldata_tlock=ft_timelockanalysis([],appseldata);
          
          for cc=1:63,corr_erppca(cc)=corr(appseldata_tlock.avg(18,:)',pcatacaud_tlock.avg(cc,:)');end
          [mxc,mind]=max(abs(corr_erppca));
          if mxc<.7
            %                   if     [ii==11 && ll==9 && ss==10 && sleep==0]
            %                   elseif [ii==26 && ll==9 && ss==10 && sleep==0]
            %                   elseif [ii==31 && ll==9 && ss==10 && sleep==0]
            %                     mind=2;
            %                   elseif [ii==8 && ll==3 && ss==10 && sleep==0]
            %                     mind=2;
            %                   else
            keyboard
            %                   end
          end
          montage     = [];
          montage.tra = pcatacaud_tlock.topo(:, mind)';
          montage.labelorg = pcatacaud.topolabel;
          montage.labelnew = pcatacaud.label(mind);
          clear tmp_ms* appseldata pcatacaud
        end
        
        if ~all(numt_trials(soalist,tt,ss))
          ssuse=setdiff(ssuse,ss);
          continue
        end
        
        for ll=soalist
          
          if phaset0  % For phase at time0
            
            if sleep==0
              tlock_tac_pca=ft_apply_montage(tlock_tac{ll+20,tt,ss},montage);
              tlock_tac_pca.label{1}='pca_maxcorrerp';
              tlock_nul_pca=ft_apply_montage(tlock_nul{ll+50,tt,ss},montage);
              tlock_nul_pca.label{1}='pca_maxcorrerp';
              tlock_ms_pca=ft_apply_montage(tlock_tac{ll,tt,ss},montage);
              tlock_ms_pca.label{1}='pca_maxcorrerp';
              tlock_aud_pca=ft_apply_montage(tlock_aud{ll+40,tt,ss},montage);
              tlock_aud_pca.label{1}='pca_maxcorrerp';
              if sign(corr_erppca(mind))==-1
                cfg=[];
                cfg.operation='-1*x1';
                cfg.parameter='trial';
                tlock_tac_pca=ft_math(cfg,tlock_tac_pca);
                tlock_nul_pca=ft_math(cfg,tlock_nul_pca);
                tlock_ms_pca =ft_math(cfg,tlock_ms_pca);
                tlock_aud_pca=ft_math(cfg,tlock_aud_pca);
              end
              cfg=[];
              cfg.keeptrials='yes';
              tlock_tac_pca=ft_timelockanalysis(cfg,tlock_tac_pca);
              tlock_nul_pca=ft_timelockanalysis(cfg,tlock_nul_pca);
              tlock_ms_pca =ft_timelockanalysis(cfg,tlock_ms_pca);
              tlock_aud_pca=ft_timelockanalysis(cfg,tlock_aud_pca);
            end
            
            %                 % add for multisensory interactions
            %                 combuse=Shuffle(1:numt_trials(ll,tt,ss));
            %                 for cc=1
            %                   for at=1:numt_trials(ll,tt,ss)
            %                     if cc==1  % do it as 'normal'
            %                       tind=at;aind=at;
            %                     else
            %                       %                   [tind,aind]=find(combindex==combuse(at));
            %                       tind=at;
            %                       aind=combuse(at);
            %                       %                   [tind,aind]=find(combindex==combuse(at));
            %                     end
            %                     cfg=[];
            %                     cfg.trials=tind;
            %                     tmpt=ft_selectdata(cfg,tlock_tac{ll+20,tt,ss});
            %                     tmpms=ft_selectdata(cfg,tlock_tac{ll,tt,ss});
            %                     cfg.trials=aind;
            %                     tmpa=ft_selectdata(cfg,tlock_aud{ll+40,tt,ss});
            %                     tmpn=ft_selectdata(cfg,tlock_nul{ll+50,tt,ss});
            %                     cfg=[];
            %                     cfg.operation='add';
            %                     cfg.parameter='trial';
            %                     tmpsumU=ft_math(cfg,tmpt,tmpa);
            %                     tmpsumM=ft_math(cfg,tmpms,tmpn);
            %                     cfg=[];
            %                     cfg.operation='subtract';
            %                     cfg.parameter='trial';
            %                     tmpconMT=ft_math(cfg,tmpms,tmpt); % these are different contrasts, since add/subtact is done before FFT (nonlinear) step
            %                     tmpconMA=ft_math(cfg,tmpms,tmpa);
            %                     tmpconAN=ft_math(cfg,tmpa,tmpn);
            %                     tmpconTN=ft_math(cfg,tmpt,tmpn);
            %                     if at==1
            %                       tlock_fakeU=tmpsumU;
            %                       tlock_fakeM=tmpsumM;
            %                       tlock_fakeMT=tmpconMT;
            %                       tlock_fakeMA=tmpconMA;
            %                       tlock_fakeTN=tmpconTN;
            %                       tlock_fakeAN=tmpconAN;
            %                     end
            %                     tlock_fakeU.trial(at,:,:)=tmpsumU.trial(1,:,:);
            %                     tlock_fakeM.trial(at,:,:)=tmpsumM.trial(1,:,:);
            %                     tlock_fakeMT.trial(at,:,:)=tmpconMT.trial(1,:,:);
            %                     tlock_fakeMA.trial(at,:,:)=tmpconMA.trial(1,:,:);
            %                     tlock_fakeTN.trial(at,:,:)=tmpconTN.trial(1,:,:);
            %                     tlock_fakeAN.trial(at,:,:)=tmpconAN.trial(1,:,:);
            %                   end
            %                 end
            
            
            scfg=[];  % for channel subselection to get phase for channel grouping
            scfg.channel={'Fz' 'F1' 'F2' 'FC1' 'FC2' 'C1' 'Cz' 'C2'};
            scfg.avgoverchan='yes';
            
            
            cfg=[];
            cfg.method='mtmconvol';
            cfg.pad=4;
            cfg.foi=1:13;
            cfg.taper='hanning';
            cfg.toi=[-0.5:.05:0]+0;  % use -100ms for awake 10 Hz and -500ms for delta N2
            cfg.t_ftimwin=2./cfg.foi;
            cfg.output='fourier';
            freqlo_tTacAlone_time0{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tac{ll+20,tt,ss});
            freqlo_tNulAlone_time0{ll,tt,ss}=ft_freqanalysis(cfg,tlock_nul{ll+50,tt,ss});
            freqlo_tTacAlone_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_tac{ll+20,tt,ss}));
            freqlo_tNulAlone_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_nul{ll+50,tt,ss}));
            if sleep==0
              freqlo_tTacAlone_chanPCA_time0{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tac_pca);
              freqlo_tNulAlone_chanPCA_time0{ll,tt,ss}=ft_freqanalysis(cfg,tlock_nul_pca);
            end
            if ll==1
              cfg.toi=[-0.5:.05:0]-.5;
            elseif ll==3
              cfg.toi=[-0.5:.05:0]-.07;
            elseif ll==4
              cfg.toi=[-0.5:.05:0]-0.02;
            elseif ll==5 || ll==6 || ll==7 || ll==9
              cfg.toi=[-0.5:.05:0]+0;
            end
            freqlo_tMSAlone_first_time0{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tac{ll,tt,ss});
            freqlo_tMSAlone_first_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_tac{ll,tt,ss}));
            if sleep==0
              freqlo_tMSAlone_first_chanPCA_time0{ll,tt,ss}=ft_freqanalysis(cfg,tlock_ms_pca);
            end
            if ll==1 || ll==3 || ll==4 || ll==5
              cfg.toi=[-0.5:.05:0]+0;
            elseif ll==6
              cfg.toi=[-0.5:.05:0]+0.02;
            elseif ll==7
              cfg.toi=[-0.5:.05:0]+.07;
            elseif ll==9
              cfg.toi=[-0.5:.05:0]+.5;
            end
            freqlo_tMSAlone_second_time0{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tac{ll,tt,ss});
            freqlo_tMSAlone_second_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_tac{ll,tt,ss}));
            if sleep==0
              freqlo_tMSAlone_second_chanPCA_time0{ll,tt,ss}=ft_freqanalysis(cfg,tlock_ms_pca);
            end
            if ll==1
              cfg.toi=[-0.5:.05:0]-.5;
            elseif ll==3
              cfg.toi=[-0.5:.05:0]-.07;
            elseif ll==4
              cfg.toi=[-0.5:.05:0]-0.02;
            elseif ll==5
              cfg.toi=[-0.5:.05:0]+0;
            elseif ll==6
              cfg.toi=[-0.5:.05:0]+0.02;
            elseif ll==7
              cfg.toi=[-0.5:.05:0]+.07;
            elseif ll==9
              cfg.toi=[-0.5:.05:0]+.5;
            end
            freqlo_tAudAlone_time0{ll,tt,ss}=ft_freqanalysis(cfg,tlock_aud{ll+40,tt,ss});
            freqlo_tAudAlone_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_aud{ll+40,tt,ss}));
            if sleep==0
              freqlo_tAudAlone_chanPCA_time0{ll,tt,ss}=ft_freqanalysis(cfg,tlock_aud_pca);
            end
            
            
            if sleep
              if ll==1
                cfg.toi=[-0.5000   -0.2500   -0.1667   -0.1250]-.5;
              elseif ll==3
                cfg.toi=[-0.5000   -0.2500   -0.1667   -0.1250]-.07;
              elseif ll==4
                cfg.toi=[-0.5000   -0.2500   -0.1667   -0.1250]-0.02;
              elseif ll==5 || ll==6 || ll==7 || ll==9
                cfg.toi=[-0.5000   -0.2500   -0.1667   -0.1250]+0;
              end
            else
              if ll==1
                cfg.toi=[-0.1]-.5;
              elseif ll==3
                cfg.toi=[-0.1]-.07;
              elseif ll==4
                cfg.toi=[-0.1]-0.02;
              elseif ll==5 || ll==6 || ll==7 || ll==9
                cfg.toi=[-0.1]+0;
              end
            end
            freqlo_tac4mscon1_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_tac{ll+20,tt,ss}));
            freqlo_nul4mscon1_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_nul{ll+50,tt,ss}));
            freqlo_ms4mscon1_first_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_tac{ll,tt,ss}));
            freqlo_aud4mscon1_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_aud{ll+40,tt,ss}));
            %                 freqlo_tacPaud_first_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_fakeU));
            %                 freqlo_tacMSpN_first_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_fakeM));
            %                 freqlo_tacMSmT_first_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_fakeMT));
            %                 freqlo_tacAmN_first_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_fakeAN));
            %                 freqlo_tacMSmA_first_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_fakeMA));
            %                 freqlo_tacTmN_first_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_fakeTN));
            if sleep
              if ll==1 || ll==3 || ll==4 || ll==5
                cfg.toi=[-0.5000   -0.2500   -0.1667   -0.1250]+0;
              elseif ll==6
                cfg.toi=[-0.5000   -0.2500   -0.1667   -0.1250]+0.02;
              elseif ll==7
                cfg.toi=[-0.5000   -0.2500   -0.1667   -0.1250]+.07;
              elseif ll==9
                cfg.toi=[-0.5000   -0.2500   -0.1667   -0.1250]+.5;
              end
            else
              if ll==1 || ll==3 || ll==4 || ll==5
                cfg.toi=[-0.1]+0;
              elseif ll==6
                cfg.toi=[-0.1]+0.02;
              elseif ll==7
                cfg.toi=[-0.1]+.07;
              elseif ll==9
                cfg.toi=[-0.1]+.5;
              end
            end
            freqlo_tac4mscon2_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_tac{ll+20,tt,ss}));
            freqlo_nul4mscon2_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_nul{ll+50,tt,ss}));
            freqlo_ms4mscon2_first_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_tac{ll,tt,ss}));
            freqlo_aud4mscon2_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_aud{ll+40,tt,ss}));
            %                 freqlo_tacPaud_second_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_fakeU));
            %                 freqlo_tacMSpN_second_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_fakeM));
            %                 freqlo_tacMSmT_second_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_fakeMT));
            %                 freqlo_tacAmN_second_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_fakeAN));
            %                 freqlo_tacMSmA_second_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_fakeMA));
            %                 freqlo_tacTmN_second_chanFC_time0{ll,tt,ss}=ft_freqanalysis(cfg,ft_selectdata(scfg,tlock_fakeTN));
            
            clear tlock_fake*
            
            save(['freqtime0_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_tt' num2str(tt) '_trialkc' num2str(trialkc) '.mat'],'freq*time0','tr','-v7.3')
            continue
          end
          
          tic
          % individual conditions first
          cfg=[];
          cfg.method='mtmconvol';
          cfg.pad=4;
          cfg.foi=4:2:30;
          cfg.taper='hanning';
          cfg.toi=-1.2:timestepplv:1.3;
          cfg.t_ftimwin=4./cfg.foi;
          %         cfg.keeptrials='yes';  % not necessary at this stage for keeping trials (yes later for source analysis and phase resetting)
          cfg.output='fourier';
          freqlo_tTacAlone_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tac{ll+20,tt,ss});
          freqlo_tAudAlone_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_aud{ll+40,tt,ss});
          freqlo_tMSAlone_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tac{ll,tt,ss});
          freqlo_tNulAlone_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_nul{ll+50,tt,ss});
          
          freqlo_tTacAlone_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqlo_tTacAlone_comb{ll,tt,ss}.fourierspctrm).^2,1));
          freqlo_tTacAlone_comb{ll,tt,ss}.plvspctrm=getplv(freqlo_tTacAlone_comb{ll,tt,ss}.fourierspctrm);
          freqlo_tTacAlone_comb{ll,tt,ss}=rmfield(freqlo_tTacAlone_comb{ll,tt,ss},'fourierspctrm');
          freqlo_tAudAlone_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqlo_tAudAlone_comb{ll,tt,ss}.fourierspctrm).^2,1));
          freqlo_tAudAlone_comb{ll,tt,ss}.plvspctrm=getplv(freqlo_tAudAlone_comb{ll,tt,ss}.fourierspctrm);
          freqlo_tAudAlone_comb{ll,tt,ss}=rmfield(freqlo_tAudAlone_comb{ll,tt,ss},'fourierspctrm');
          freqlo_tMSAlone_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqlo_tMSAlone_comb{ll,tt,ss}.fourierspctrm).^2,1));
          freqlo_tMSAlone_comb{ll,tt,ss}.plvspctrm=getplv(freqlo_tMSAlone_comb{ll,tt,ss}.fourierspctrm);
          freqlo_tMSAlone_comb{ll,tt,ss}=rmfield(freqlo_tMSAlone_comb{ll,tt,ss},'fourierspctrm');
          freqlo_tNulAlone_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqlo_tNulAlone_comb{ll,tt,ss}.fourierspctrm).^2,1));
          freqlo_tNulAlone_comb{ll,tt,ss}.plvspctrm=getplv(freqlo_tNulAlone_comb{ll,tt,ss}.fourierspctrm);
          freqlo_tNulAlone_comb{ll,tt,ss}=rmfield(freqlo_tNulAlone_comb{ll,tt,ss},'fourierspctrm');
          
          if synchasynch && ll<5
            freqlo_tacMSshift1_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tacMSshift1{ll,tt,ss});
            freqlo_tacMSshift2_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tacMSshift2{ll,tt,ss});
            freqlo_tacMSshift3_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tacMSshift3{ll,tt,ss});
            freqlo_tacMSshift4_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tacMSshift4{ll,tt,ss});
            freqlo_tacMSshift1_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqlo_tacMSshift1_comb{ll,tt,ss}.fourierspctrm).^2,1));
            freqlo_tacMSshift1_comb{ll,tt,ss}.plvspctrm=getplv(freqlo_tacMSshift1_comb{ll,tt,ss}.fourierspctrm);
            freqlo_tacMSshift1_comb{ll,tt,ss}=rmfield(freqlo_tacMSshift1_comb{ll,tt,ss},'fourierspctrm');
            freqlo_tacMSshift2_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqlo_tacMSshift2_comb{ll,tt,ss}.fourierspctrm).^2,1));
            freqlo_tacMSshift2_comb{ll,tt,ss}.plvspctrm=getplv(freqlo_tacMSshift2_comb{ll,tt,ss}.fourierspctrm);
            freqlo_tacMSshift2_comb{ll,tt,ss}=rmfield(freqlo_tacMSshift2_comb{ll,tt,ss},'fourierspctrm');
            freqlo_tacMSshift3_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqlo_tacMSshift3_comb{ll,tt,ss}.fourierspctrm).^2,1));
            freqlo_tacMSshift3_comb{ll,tt,ss}.plvspctrm=getplv(freqlo_tacMSshift3_comb{ll,tt,ss}.fourierspctrm);
            freqlo_tacMSshift3_comb{ll,tt,ss}=rmfield(freqlo_tacMSshift3_comb{ll,tt,ss},'fourierspctrm');
            freqlo_tacMSshift4_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqlo_tacMSshift4_comb{ll,tt,ss}.fourierspctrm).^2,1));
            freqlo_tacMSshift4_comb{ll,tt,ss}.plvspctrm=getplv(freqlo_tacMSshift4_comb{ll,tt,ss}.fourierspctrm);
            freqlo_tacMSshift4_comb{ll,tt,ss}=rmfield(freqlo_tacMSshift4_comb{ll,tt,ss},'fourierspctrm');
          end
          
          cfg=[];
          cfg.method='mtmconvol';
          cfg.pad=4;
          cfg.foi=[30:5:45 55:5:80];
          cfg.taper='dpss';
          cfg.tapsmofrq=7*ones(1,length(cfg.foi));
          cfg.toi=-1.2:timestepplv:1.3;
          cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
          cfg.output='fourier';
          freqhi_tTacAlone_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tac{ll+20,tt,ss});
          freqhi_tAudAlone_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_aud{ll+40,tt,ss});
          freqhi_tMSAlone_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tac{ll,tt,ss});
          freqhi_tNulAlone_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_nul{ll+50,tt,ss});
          
          freqhi_tTacAlone_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqhi_tTacAlone_comb{ll,tt,ss}.fourierspctrm).^2,1));
          freqhi_tTacAlone_comb{ll,tt,ss}.plvspctrm=getplv(freqhi_tTacAlone_comb{ll,tt,ss}.fourierspctrm);
          freqhi_tTacAlone_comb{ll,tt,ss}=rmfield(freqhi_tTacAlone_comb{ll,tt,ss},'fourierspctrm');
          freqhi_tAudAlone_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqhi_tAudAlone_comb{ll,tt,ss}.fourierspctrm).^2,1));
          freqhi_tAudAlone_comb{ll,tt,ss}.plvspctrm=getplv(freqhi_tAudAlone_comb{ll,tt,ss}.fourierspctrm);
          freqhi_tAudAlone_comb{ll,tt,ss}=rmfield(freqhi_tAudAlone_comb{ll,tt,ss},'fourierspctrm');
          freqhi_tMSAlone_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqhi_tMSAlone_comb{ll,tt,ss}.fourierspctrm).^2,1));
          freqhi_tMSAlone_comb{ll,tt,ss}.plvspctrm=getplv(freqhi_tMSAlone_comb{ll,tt,ss}.fourierspctrm);
          freqhi_tMSAlone_comb{ll,tt,ss}=rmfield(freqhi_tMSAlone_comb{ll,tt,ss},'fourierspctrm');
          freqhi_tNulAlone_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqhi_tNulAlone_comb{ll,tt,ss}.fourierspctrm).^2,1));
          freqhi_tNulAlone_comb{ll,tt,ss}.plvspctrm=getplv(freqhi_tNulAlone_comb{ll,tt,ss}.fourierspctrm);
          freqhi_tNulAlone_comb{ll,tt,ss}=rmfield(freqhi_tNulAlone_comb{ll,tt,ss},'fourierspctrm');
          
          if synchasynch && ll<5
            freqhi_tacMSshift1_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tacMSshift1{ll,tt,ss});
            freqhi_tacMSshift2_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tacMSshift2{ll,tt,ss});
            freqhi_tacMSshift3_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tacMSshift3{ll,tt,ss});
            freqhi_tacMSshift4_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tacMSshift4{ll,tt,ss});
            freqhi_tacMSshift1_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqhi_tacMSshift1_comb{ll,tt,ss}.fourierspctrm).^2,1));
            freqhi_tacMSshift1_comb{ll,tt,ss}.plvspctrm=getplv(freqhi_tacMSshift1_comb{ll,tt,ss}.fourierspctrm);
            freqhi_tacMSshift1_comb{ll,tt,ss}=rmfield(freqhi_tacMSshift1_comb{ll,tt,ss},'fourierspctrm');
            freqhi_tacMSshift2_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqhi_tacMSshift2_comb{ll,tt,ss}.fourierspctrm).^2,1));
            freqhi_tacMSshift2_comb{ll,tt,ss}.plvspctrm=getplv(freqhi_tacMSshift2_comb{ll,tt,ss}.fourierspctrm);
            freqhi_tacMSshift2_comb{ll,tt,ss}=rmfield(freqhi_tacMSshift2_comb{ll,tt,ss},'fourierspctrm');
            freqhi_tacMSshift3_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqhi_tacMSshift3_comb{ll,tt,ss}.fourierspctrm).^2,1));
            freqhi_tacMSshift3_comb{ll,tt,ss}.plvspctrm=getplv(freqhi_tacMSshift3_comb{ll,tt,ss}.fourierspctrm);
            freqhi_tacMSshift3_comb{ll,tt,ss}=rmfield(freqhi_tacMSshift3_comb{ll,tt,ss},'fourierspctrm');
            freqhi_tacMSshift4_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqhi_tacMSshift4_comb{ll,tt,ss}.fourierspctrm).^2,1));
            freqhi_tacMSshift4_comb{ll,tt,ss}.plvspctrm=getplv(freqhi_tacMSshift4_comb{ll,tt,ss}.fourierspctrm);
            freqhi_tacMSshift4_comb{ll,tt,ss}=rmfield(freqhi_tacMSshift4_comb{ll,tt,ss},'fourierspctrm');
          end
          
          %           if ll<10
          % tacPaud
          numcomb=numt_trials(ll,tt,ss)^2;
          numtests=min(numcomb,minnumcomb);
          %             combindex=reshape(1:numcomb,numt_trials(ll,tt,ss),numt_trials(ll,tt,ss));
          freqlo_tacPaud_tmp{1}=[];
          freqhi_tacPaud_tmp{1}=[];
          for cc=1:numtests
            %               tmp=Shuffle(combindex(:));
            %               combuse=tmp(1:numt_trials(ll,tt,ss));
            combuse=Shuffle(1:numt_trials(ll,tt,ss));
            for at=1:numt_trials(ll,tt,ss)
              if cc==1  % do it as 'normal'
                tind=at;aind=at;
              else
                %                   [tind,aind]=find(combindex==combuse(at));
                tind=at;
                aind=combuse(at);
                %                   [tind,aind]=find(combindex==combuse(at));
              end
              cfg=[];
              cfg.trials=tind;
              tmpt=ft_selectdata(cfg,tlock_tac{ll+20,tt,ss});
              cfg.trials=aind;
              tmpa=ft_selectdata(cfg,tlock_aud{ll+40,tt,ss});
              cfg=[];
              cfg.operation='add';
              cfg.parameter='trial';
              tmpsum=ft_math(cfg,tmpt,tmpa);
              if at==1
                tlock_fake=tmpsum;
              end
              tlock_fake.trial(at,:,:)=tmpsum.trial(1,:,:);
            end
            
            cfg=[];
            cfg.method='mtmconvol';
            cfg.pad=4;
            cfg.foi=4:2:30;
            cfg.taper='hanning';
            cfg.toi=-1.2:timestepplv:1.3;
            cfg.t_ftimwin=4./cfg.foi;
            %         cfg.keeptrials='yes';  % not necessary at this stage for keeping trials (yes later for source analysis and phase resetting)
            cfg.output='fourier';
            freqlo_tacPaud_tmp{cc}=ft_freqanalysis(cfg,tlock_fake);
            freqlo_tacPaud_tmp{cc}.powspctrm=squeeze(mean(abs(freqlo_tacPaud_tmp{cc}.fourierspctrm).^2,1));
            freqlo_tacPaud_tmp{cc}.plvspctrm=getplv(freqlo_tacPaud_tmp{cc}.fourierspctrm);
            %             freqlo_tacPaud_tmp{cc}.plvmag=abs(getplv(freqlo_tacPaud_tmp{cc}.fourierspctrm));
            %             freqlo_tacPaud_tmp{cc}.plvang=angle(getplv(freqlo_tacPaud_tmp{cc}.fourierspctrm));
            freqlo_tacPaud_tmp{cc}=rmfield(freqlo_tacPaud_tmp{cc},'fourierspctrm');
            cfg=[];
            cfg.method='mtmconvol';
            cfg.pad=4;
            cfg.foi=[30:5:45 55:5:80];
            cfg.taper='dpss';
            cfg.tapsmofrq=7*ones(1,length(cfg.foi));
            cfg.toi=-1.2:timestepplv:1.3;
            cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
            cfg.output='fourier';
            freqhi_tacPaud_tmp{cc}=ft_freqanalysis(cfg,tlock_fake);
            freqhi_tacPaud_tmp{cc}.powspctrm=squeeze(mean(abs(freqhi_tacPaud_tmp{cc}.fourierspctrm).^2,1));
            freqhi_tacPaud_tmp{cc}.plvspctrm=getplv(freqhi_tacPaud_tmp{cc}.fourierspctrm);
            freqhi_tacPaud_tmp{cc}=rmfield(freqhi_tacPaud_tmp{cc},'fourierspctrm');
            disp(['tacPaud cc ' num2str(cc)]),toc
          end % cc
          freqlo_tacPaud_comb{ll,tt,ss,1}=freqlo_tacPaud_tmp{1};
          clear tmp tmp2 tmp3
          if numtests>1
            for cc=1:numtests,tmp(:,:,:,cc)=freqlo_tacPaud_tmp{cc}.powspctrm;end
            %for cc=1:numtests,tmp2(:,:,:,cc)=freqlo_tacPaud_tmp{cc}.crsspctrm;end
            for cc=1:numtests,tmp3(:,:,:,cc)=freqlo_tacPaud_tmp{cc}.plvspctrm;end
            if exist('tmp','var')
              freqlo_tacPaud_comb{ll,tt,ss,2}=freqlo_tacPaud_tmp{1};
              freqlo_tacPaud_comb{ll,tt,ss,2}.powspctrm=mean(tmp,4);
              freqlo_tacPaud_comb{ll,tt,ss,2}.stdpow=std(tmp,[],4);
              %            freqlo_tacPaud_comb{ll,tt,ss,2}.crsspctrm=mean(tmp2,4);
              %            freqlo_tacPaud_comb{ll,tt,ss,2}.stdcrs=std(tmp2,[],4);
              freqlo_tacPaud_comb{ll,tt,ss,2}.plvspctrm=mean(tmp3,4);
              freqlo_tacPaud_comb{ll,tt,ss,2}.stdplv=std(tmp3,[],4);
              if numtests>=10
                halftest=round(numtests/2);
                freqlo_tacPaud_comb{ll,tt,ss,3}=freqlo_tacPaud_tmp{1};
                freqlo_tacPaud_comb{ll,tt,ss,4}=freqlo_tacPaud_tmp{1};
                freqlo_tacPaud_comb{ll,tt,ss,3}.powspctrm=mean(tmp(:,:,:,1:halftest),4);
                freqlo_tacPaud_comb{ll,tt,ss,3}.stdpow=std(tmp(:,:,:,1:halftest),[],4);
                freqlo_tacPaud_comb{ll,tt,ss,3}.plvspctrm=mean(tmp3(:,:,:,1:halftest),4);
                freqlo_tacPaud_comb{ll,tt,ss,3}.stdplv=std(tmp3(:,:,:,1:halftest),[],4);
                freqlo_tacPaud_comb{ll,tt,ss,4}.powspctrm=mean(tmp(:,:,:,halftest+1:end),4);
                freqlo_tacPaud_comb{ll,tt,ss,4}.stdpow=std(tmp(:,:,:,halftest+1:end),[],4);
                freqlo_tacPaud_comb{ll,tt,ss,4}.plvspctrm=mean(tmp3(:,:,:,halftest+1:end),4);
                freqlo_tacPaud_comb{ll,tt,ss,4}.stdplv=std(tmp3(:,:,:,halftest+1:end),[],4);
              end
            end
          else
            freqlo_tacPaud_comb{ll,tt,ss,2}=[];
          end
          tmplo_tacPaud_comb2=tmp;
          tmplo_tacPaud_comb2plv=tmp3;
          freqhi_tacPaud_comb{ll,tt,ss,1}=freqhi_tacPaud_tmp{1};
          clear tmp tmp2 tmp3
          if numtests>1
            for cc=1:numtests,tmp(:,:,:,cc)=freqhi_tacPaud_tmp{cc}.powspctrm;end
            %          for cc=1:numtests,tmp2(:,:,:,cc)=freqhi_tacPaud_tmp{cc}.crsspctrm;end
            for cc=1:numtests,tmp3(:,:,:,cc)=freqhi_tacPaud_tmp{cc}.plvspctrm;end
            if exist('tmp','var')
              freqhi_tacPaud_comb{ll,tt,ss,2}=freqhi_tacPaud_tmp{1};
              freqhi_tacPaud_comb{ll,tt,ss,2}.powspctrm=mean(tmp,4);
              freqhi_tacPaud_comb{ll,tt,ss,2}.stdpow=std(tmp,[],4);
              %            freqhi_tacPaud_comb{ll,tt,ss,2}.crsspctrm=mean(tmp2,4);
              %            freqhi_tacPaud_comb{ll,tt,ss,2}.stdcsd=std(tmp2,[],4);
              freqhi_tacPaud_comb{ll,tt,ss,2}.plvspctrm=mean(tmp3,4);
              freqhi_tacPaud_comb{ll,tt,ss,2}.stdplv=std(tmp3,[],4);
              if numtests>=10
                halftest=round(numtests/2);
                freqhi_tacPaud_comb{ll,tt,ss,3}=freqhi_tacPaud_tmp{1};
                freqhi_tacPaud_comb{ll,tt,ss,4}=freqhi_tacPaud_tmp{1};
                freqhi_tacPaud_comb{ll,tt,ss,3}.powspctrm=mean(tmp(:,:,:,1:halftest),4);
                freqhi_tacPaud_comb{ll,tt,ss,3}.stdpow=std(tmp(:,:,:,1:halftest),[],4);
                freqhi_tacPaud_comb{ll,tt,ss,3}.plvspctrm=mean(tmp3(:,:,:,1:halftest),4);
                freqhi_tacPaud_comb{ll,tt,ss,3}.stdplv=std(tmp3(:,:,:,1:halftest),[],4);
                freqhi_tacPaud_comb{ll,tt,ss,4}.powspctrm=mean(tmp(:,:,:,halftest+1:end),4);
                freqhi_tacPaud_comb{ll,tt,ss,4}.stdpow=std(tmp(:,:,:,halftest+1:end),[],4);
                freqhi_tacPaud_comb{ll,tt,ss,4}.plvspctrm=mean(tmp3(:,:,:,halftest+1:end),4);
                freqhi_tacPaud_comb{ll,tt,ss,4}.stdplv=std(tmp3(:,:,:,halftest+1:end),[],4);
              end
            end
          else
            freqhi_tacPaud_comb{ll,tt,ss,2}=[];
          end
          tmphi_tacPaud_comb2=tmp;
          tmphi_tacPaud_comb2plv=tmp3;
          clear freqlo_tacPaud_tmp freqhi_tacPaud_tmp
          clear tmp tmp2 tmp3
          
          if (use23 && (ss==12 || ss==13))
          else
            tlock_tac{ll+20,tt,ss}=[]; % clearing to help save memory as running
            tlock_aud{ll+40,tt,ss}=[]; % clearing to help save memory as running
          end
          toc
          disp('finished tacPaud cc'),toc
          
          
          % tacMSpN
          numcomb=numt_trials(ll,tt,ss)^2;
          numtests=min(numcomb,minnumcomb);
          %             combindex=reshape(1:numcomb,numt_trials(ll,tt,ss),numt_trials(ll,tt,ss));
          freqlo_tacMSpN_tmp{1}=[];
          freqhi_tacMSpN_tmp{1}=[];
          for cc=1:numtests
            %               tmp=Shuffle(combindex(:));
            %               combuse=tmp(1:numt_trials(ll,tt,ss));
            combuse=Shuffle(1:numt_trials(ll,tt,ss));
            
            tlock_fake=tlock_tac{ll,tt,ss};
            for at=1:numt_trials(ll,tt,ss)
              if cc==1  % do it as 'normal'
                tind=at;aind=at;
              else
                %                   [tind,aind]=find(combindex==combuse(at));
                tind=at;
                aind=combuse(at);
              end
              cfg=[];
              cfg.trials=tind;
              tmpt=ft_selectdata(cfg,tlock_tac{ll,tt,ss});
              cfg.trials=aind;
              tmpa=ft_selectdata(cfg,tlock_nul{ll+50,tt,ss});
              cfg=[];
              cfg.operation='add';
              cfg.parameter='trial';
              tmpsum=ft_math(cfg,tmpt,tmpa);
              if at==1
                tlock_fake=tmpsum;
              end
              tlock_fake.trial(at,:,:)=tmpsum.trial(1,:,:);
              %           tlock_fake.trial(at,:,:)=tlock_tac{ll,tt,ss}.trial(tind,:,:)+tlock_nul{ll+50,tt,ss}.trial(aind,:,:);
            end
            cfg=[];
            cfg.method='mtmconvol';
            cfg.pad=4;
            cfg.foi=4:2:30;
            cfg.taper='hanning';
            cfg.toi=-1.2:timestepplv:1.3;
            cfg.t_ftimwin=4./cfg.foi;
            %         cfg.keeptrials='yes';  % not necessary at this stage for keeping trials (yes later for source analysis and phase resetting)
            cfg.output='fourier';
            freqlo_tacMSpN_tmp{cc}=ft_freqanalysis(cfg,tlock_fake);
            freqlo_tacMSpN_tmp{cc}.powspctrm=squeeze(mean(abs(freqlo_tacMSpN_tmp{cc}.fourierspctrm).^2,1));
            freqlo_tacMSpN_tmp{cc}.plvspctrm=getplv(freqlo_tacMSpN_tmp{cc}.fourierspctrm);
            freqlo_tacMSpN_tmp{cc}=rmfield(freqlo_tacMSpN_tmp{cc},'fourierspctrm');
            cfg=[];
            cfg.method='mtmconvol';
            cfg.pad=4;
            cfg.foi=[30:5:45 55:5:80];
            cfg.taper='dpss';
            cfg.tapsmofrq=7*ones(1,length(cfg.foi));
            cfg.toi=-1.2:timestepplv:1.3;
            cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
            cfg.output='fourier';
            freqhi_tacMSpN_tmp{cc}=ft_freqanalysis(cfg,tlock_fake);
            freqhi_tacMSpN_tmp{cc}.powspctrm=squeeze(mean(abs(freqhi_tacMSpN_tmp{cc}.fourierspctrm).^2,1));
            freqhi_tacMSpN_tmp{cc}.plvspctrm=getplv(freqhi_tacMSpN_tmp{cc}.fourierspctrm);
            freqhi_tacMSpN_tmp{cc}=rmfield(freqhi_tacMSpN_tmp{cc},'fourierspctrm');
          end
          freqlo_tacMSpN_comb{ll,tt,ss,1}=freqlo_tacMSpN_tmp{1};
          clear tmp tmp2 tmp3
          if numtests>1
            for cc=1:numtests,tmp(:,:,:,cc)=freqlo_tacMSpN_tmp{cc}.powspctrm;end
            %          for cc=1:numtests,tmp2(:,:,:,cc)=freqlo_tacMSpN_tmp{cc}.crsspctrm;end
            for cc=1:numtests,tmp3(:,:,:,cc)=freqlo_tacMSpN_tmp{cc}.plvspctrm;end
            if exist('tmp','var')
              freqlo_tacMSpN_comb{ll,tt,ss,2}=freqlo_tacMSpN_tmp{1};
              freqlo_tacMSpN_comb{ll,tt,ss,2}.powspctrm=mean(tmp,4);
              freqlo_tacMSpN_comb{ll,tt,ss,2}.stdpow=std(tmp,[],4);
              %            freqlo_tacMSpN_comb{ll,tt,ss,2}.crsspctrm=mean(tmp2,4);
              %            freqlo_tacMSpN_comb{ll,tt,ss,2}.stdcrs=std(tmp2,[],4);
              freqlo_tacMSpN_comb{ll,tt,ss,2}.plvspctrm=mean(tmp3,4);
              freqlo_tacMSpN_comb{ll,tt,ss,2}.stdplv=std(tmp3,[],4);
              if numtests>=10
                halftest=round(numtests/2);
                freqlo_tacMSpN_comb{ll,tt,ss,3}=freqlo_tacMSpN_tmp{1};
                freqlo_tacMSpN_comb{ll,tt,ss,4}=freqlo_tacMSpN_tmp{1};
                freqlo_tacMSpN_comb{ll,tt,ss,3}.powspctrm=mean(tmp(:,:,:,1:halftest),4);
                freqlo_tacMSpN_comb{ll,tt,ss,3}.stdpow=std(tmp(:,:,:,1:halftest),[],4);
                freqlo_tacMSpN_comb{ll,tt,ss,3}.plvspctrm=mean(tmp3(:,:,:,1:halftest),4);
                freqlo_tacMSpN_comb{ll,tt,ss,3}.stdplv=std(tmp3(:,:,:,1:halftest),[],4);
                freqlo_tacMSpN_comb{ll,tt,ss,4}.powspctrm=mean(tmp(:,:,:,halftest+1:end),4);
                freqlo_tacMSpN_comb{ll,tt,ss,4}.stdpow=std(tmp(:,:,:,halftest+1:end),[],4);
                freqlo_tacMSpN_comb{ll,tt,ss,4}.plvspctrm=mean(tmp3(:,:,:,halftest+1:end),4);
                freqlo_tacMSpN_comb{ll,tt,ss,4}.stdplv=std(tmp3(:,:,:,halftest+1:end),[],4);
              end
            end
          else
            freqlo_tacMSpN_comb{ll,tt,ss,2}=[];
          end
          tmplo_tacMSpN_comb2=tmp;
          tmplo_tacMSpN_comb2plv=tmp3;
          freqhi_tacMSpN_comb{ll,tt,ss,1}=freqhi_tacMSpN_tmp{1};
          clear tmp tmp2 tmp3
          if numtests>1
            for cc=1:numtests,tmp(:,:,:,cc)=freqhi_tacMSpN_tmp{cc}.powspctrm;end
            %          for cc=1:numtests,tmp2(:,:,:,cc)=freqhi_tacMSpN_tmp{cc}.crsspctrm;end
            for cc=1:numtests,tmp3(:,:,:,cc)=freqhi_tacMSpN_tmp{cc}.plvspctrm;end
            if exist('tmp','var')
              freqhi_tacMSpN_comb{ll,tt,ss,2}=freqhi_tacMSpN_tmp{1};
              freqhi_tacMSpN_comb{ll,tt,ss,2}.powspctrm=mean(tmp,4);
              freqhi_tacMSpN_comb{ll,tt,ss,2}.stdpow=std(tmp,[],4);
              %            freqhi_tacMSpN_comb{ll,tt,ss,2}.crsspctrm=mean(tmp2,4);
              %            freqhi_tacMSpN_comb{ll,tt,ss,2}.stdcrs=std(tmp2,[],4);
              freqhi_tacMSpN_comb{ll,tt,ss,2}.plvspctrm=mean(tmp3,4);
              freqhi_tacMSpN_comb{ll,tt,ss,2}.stdplv=std(tmp3,[],4);
              if numtests>=10
                halftest=round(numtests/2);
                freqhi_tacMSpN_comb{ll,tt,ss,3}=freqhi_tacMSpN_tmp{1};
                freqhi_tacMSpN_comb{ll,tt,ss,4}=freqhi_tacMSpN_tmp{1};
                freqhi_tacMSpN_comb{ll,tt,ss,3}.powspctrm=mean(tmp(:,:,:,1:halftest),4);
                freqhi_tacMSpN_comb{ll,tt,ss,3}.stdpow=std(tmp(:,:,:,1:halftest),[],4);
                freqhi_tacMSpN_comb{ll,tt,ss,3}.plvspctrm=mean(tmp3(:,:,:,1:halftest),4);
                freqhi_tacMSpN_comb{ll,tt,ss,3}.stdplv=std(tmp3(:,:,:,1:halftest),[],4);
                freqhi_tacMSpN_comb{ll,tt,ss,4}.powspctrm=mean(tmp(:,:,:,halftest+1:end),4);
                freqhi_tacMSpN_comb{ll,tt,ss,4}.stdpow=std(tmp(:,:,:,halftest+1:end),[],4);
                freqhi_tacMSpN_comb{ll,tt,ss,4}.plvspctrm=mean(tmp3(:,:,:,halftest+1:end),4);
                freqhi_tacMSpN_comb{ll,tt,ss,4}.stdplv=std(tmp3(:,:,:,halftest+1:end),[],4);
              end
            end
          else
            freqhi_tacMSpN_comb{ll,tt,ss,2}=[];
          end
          tmphi_tacMSpN_comb2=tmp;
          tmphi_tacMSpN_comb2plv=tmp3;
          clear freqlo_tacMSpN_tmp freqhi_tacMSpN_tmp
          clear tmp tmp2 tmp3
          
          if use23 && (ss==12 || ss==13)
          else
            tlock_tac{ll,tt,ss}=[]; % clearing to help save memory as running
            tlock_nul{ll+50,tt,ss}=[]; % clearing to help save memory as running
          end
          disp('finished tacMSpN cc'),toc
          
          if numtests>2
            alpha=.001;
            for cc=1:size(tmplo_tacPaud_comb2,1)
              for ff=1:size(tmplo_tacPaud_comb2,2)
                for ttime=1:size(tmplo_tacPaud_comb2,3)
                  if ~isnan(mean(squeeze(tmplo_tacPaud_comb2(cc,ff,ttime,:)))) && ~isnan(mean(squeeze(tmplo_tacMSpN_comb2(cc,ff,ttime,:))))
                    [statst{ll,tt,ss}.ttestlo_h(cc,ff,ttime),statst{ll,tt,ss}.ttestlo_p(cc,ff,ttime),statst{ll,tt,ss}.ttestlo_ci(cc,ff,ttime,1,:),statst{ll,tt,ss}.ttestlo_stats{cc,ff,ttime,1}]=ttest2(squeeze(tmplo_tacPaud_comb2(cc,ff,ttime,:)),squeeze(tmplo_tacMSpN_comb2(cc,ff,ttime,:)),'alpha',alpha);
                    [statst{ll,tt,ss}.kstestlo_h(cc,ff,ttime),statst{ll,tt,ss}.kstestlo_p(cc,ff,ttime),statst{ll,tt,ss}.kstestlo_stats{cc,ff,ttime,1}]=kstest2(squeeze(tmplo_tacPaud_comb2(cc,ff,ttime,:)),squeeze(tmplo_tacMSpN_comb2(cc,ff,ttime,:)));
                    if numtests>=10
                      [statst{ll,tt,ss}.ttestlo_h(cc,ff,ttime,2),statst{ll,tt,ss}.ttestlo_p(cc,ff,ttime,2),statst{ll,tt,ss}.ttestlo_ci(cc,ff,ttime,2,:),statst{ll,tt,ss}.ttestlo_stats{cc,ff,ttime,2}]=  ttest2(squeeze(tmplo_tacPaud_comb2(cc,ff,ttime,1:halftest)),squeeze(tmplo_tacMSpN_comb2(cc,ff,ttime,1:halftest)),'alpha',alpha);
                      [statst{ll,tt,ss}.kstestlo_h(cc,ff,ttime,2),statst{ll,tt,ss}.kstestlo_p(cc,ff,ttime,2),statst{ll,tt,ss}.kstestlo_stats{cc,ff,ttime,2}]=  kstest2(squeeze(tmplo_tacPaud_comb2(cc,ff,ttime,1:halftest)),squeeze(tmplo_tacMSpN_comb2(cc,ff,ttime,1:halftest)));
                      [statst{ll,tt,ss}.ttestlo_h(cc,ff,ttime,3),statst{ll,tt,ss}.ttestlo_p(cc,ff,ttime,3),statst{ll,tt,ss}.ttestlo_ci(cc,ff,ttime,3,:),statst{ll,tt,ss}.ttestlo_stats{cc,ff,ttime,3}]=  ttest2(squeeze(tmplo_tacPaud_comb2(cc,ff,ttime,halftest+1:end)),squeeze(tmplo_tacMSpN_comb2(cc,ff,ttime,halftest+1:end)),'alpha',alpha);
                      [statst{ll,tt,ss}.kstestlo_h(cc,ff,ttime,3),statst{ll,tt,ss}.kstestlo_p(cc,ff,ttime,3),statst{ll,tt,ss}.kstestlo_stats{cc,ff,ttime,3}]=  kstest2(squeeze(tmplo_tacPaud_comb2(cc,ff,ttime,halftest+1:end)),squeeze(tmplo_tacMSpN_comb2(cc,ff,ttime,halftest+1:end)));
                    end
                    [statstplv{ll,tt,ss}.ttestlo_h(cc,ff,ttime),statstplv{ll,tt,ss}.ttestlo_p(cc,ff,ttime),statstplv{ll,tt,ss}.ttestlo_ci(cc,ff,ttime,1,:),statstplv{ll,tt,ss}.ttestlo_stats{cc,ff,ttime,1}]=ttest2(squeeze(tmplo_tacPaud_comb2plv(cc,ff,ttime,:)),squeeze(tmplo_tacMSpN_comb2plv(cc,ff,ttime,:)),'alpha',alpha);
                    [statstplv{ll,tt,ss}.kstestlo_h(cc,ff,ttime),statstplv{ll,tt,ss}.kstestlo_p(cc,ff,ttime),statstplv{ll,tt,ss}.kstestlo_stats{cc,ff,ttime,1}]=kstest2(abs(squeeze(tmplo_tacPaud_comb2plv(cc,ff,ttime,:))),abs(squeeze(tmplo_tacMSpN_comb2plv(cc,ff,ttime,:))));
                    if numtests>=10
                      [statstplv{ll,tt,ss}.ttestlo_h(cc,ff,ttime,2),statstplv{ll,tt,ss}.ttestlo_p(cc,ff,ttime,2),statstplv{ll,tt,ss}.ttestlo_ci(cc,ff,ttime,2,:),statstplv{ll,tt,ss}.ttestlo_stats{cc,ff,ttime,2}]=  ttest2(squeeze(tmplo_tacPaud_comb2plv(cc,ff,ttime,1:halftest)),squeeze(tmplo_tacMSpN_comb2plv(cc,ff,ttime,1:halftest)),'alpha',alpha);
                      [statstplv{ll,tt,ss}.kstestlo_h(cc,ff,ttime,2),statstplv{ll,tt,ss}.kstestlo_p(cc,ff,ttime,2),statstplv{ll,tt,ss}.kstestlo_stats{cc,ff,ttime,2}]=  kstest2(abs(squeeze(tmplo_tacPaud_comb2plv(cc,ff,ttime,1:halftest))),abs(squeeze(tmplo_tacMSpN_comb2plv(cc,ff,ttime,1:halftest))));
                      [statstplv{ll,tt,ss}.ttestlo_h(cc,ff,ttime,3),statstplv{ll,tt,ss}.ttestlo_p(cc,ff,ttime,3),statstplv{ll,tt,ss}.ttestlo_ci(cc,ff,ttime,3,:),statstplv{ll,tt,ss}.ttestlo_stats{cc,ff,ttime,3}]=  ttest2(squeeze(tmplo_tacPaud_comb2plv(cc,ff,ttime,halftest+1:end)),squeeze(tmplo_tacMSpN_comb2plv(cc,ff,ttime,halftest+1:end)),'alpha',alpha);
                      [statstplv{ll,tt,ss}.kstestlo_h(cc,ff,ttime,3),statstplv{ll,tt,ss}.kstestlo_p(cc,ff,ttime,3),statstplv{ll,tt,ss}.kstestlo_stats{cc,ff,ttime,3}]=  kstest2(abs(squeeze(tmplo_tacPaud_comb2plv(cc,ff,ttime,halftest+1:end))),abs(squeeze(tmplo_tacMSpN_comb2plv(cc,ff,ttime,halftest+1:end))));
                    end
                  end
                end
              end
              for ff=1:size(tmphi_tacPaud_comb2,2)
                for ttime=1:size(tmphi_tacPaud_comb2,3)
                  if ~isnan(mean(squeeze(tmphi_tacPaud_comb2(cc,ff,ttime,:)))) && ~isnan(mean(squeeze(tmphi_tacMSpN_comb2(cc,ff,ttime,:))))
                    [statst{ll,tt,ss}.ttesthi_h(cc,ff,ttime),statst{ll,tt,ss}.ttesthi_p(cc,ff,ttime),statst{ll,tt,ss}.ttesthi_ci(cc,ff,ttime,1,:),statst{ll,tt,ss}.ttesthi_stats{cc,ff,ttime,1}]=ttest2(squeeze(tmphi_tacPaud_comb2(cc,ff,ttime,:)),squeeze(tmphi_tacMSpN_comb2(cc,ff,ttime,:)),'alpha',alpha);
                    [statst{ll,tt,ss}.kstesthi_h(cc,ff,ttime),statst{ll,tt,ss}.kstesthi_p(cc,ff,ttime),statst{ll,tt,ss}.kstesthi_stats{cc,ff,ttime,1}]=kstest2(squeeze(tmphi_tacPaud_comb2(cc,ff,ttime,:)),squeeze(tmphi_tacMSpN_comb2(cc,ff,ttime,:)));
                    if numtests>=10
                      [statst{ll,tt,ss}.ttesthi_h(cc,ff,ttime,2),statst{ll,tt,ss}.ttesthi_p(cc,ff,ttime,2),statst{ll,tt,ss}.ttesthi_ci(cc,ff,ttime,2,:),statst{ll,tt,ss}.ttesthi_stats{cc,ff,ttime,2}]=  ttest2(squeeze(tmphi_tacPaud_comb2(cc,ff,ttime,1:halftest)),squeeze(tmphi_tacMSpN_comb2(cc,ff,ttime,1:halftest)),'alpha',alpha);
                      [statst{ll,tt,ss}.kstesthi_h(cc,ff,ttime,2),statst{ll,tt,ss}.kstesthi_p(cc,ff,ttime,2),statst{ll,tt,ss}.kstesthi_stats{cc,ff,ttime,2}]=  kstest2(squeeze(tmphi_tacPaud_comb2(cc,ff,ttime,1:halftest)),squeeze(tmphi_tacMSpN_comb2(cc,ff,ttime,1:halftest)));
                      [statst{ll,tt,ss}.ttesthi_h(cc,ff,ttime,3),statst{ll,tt,ss}.ttesthi_p(cc,ff,ttime,3),statst{ll,tt,ss}.ttesthi_ci(cc,ff,ttime,3,:),statst{ll,tt,ss}.ttesthi_stats{cc,ff,ttime,3}]=  ttest2(squeeze(tmphi_tacPaud_comb2(cc,ff,ttime,halftest+1:end)),squeeze(tmphi_tacMSpN_comb2(cc,ff,ttime,halftest+1:end)),'alpha',alpha);
                      [statst{ll,tt,ss}.kstesthi_h(cc,ff,ttime,3),statst{ll,tt,ss}.kstesthi_p(cc,ff,ttime,3),statst{ll,tt,ss}.kstesthi_stats{cc,ff,ttime,3}]=  kstest2(squeeze(tmphi_tacPaud_comb2(cc,ff,ttime,halftest+1:end)),squeeze(tmphi_tacMSpN_comb2(cc,ff,ttime,halftest+1:end)));
                    end
                    [statstplv{ll,tt,ss}.ttesthi_h(cc,ff,ttime),statstplv{ll,tt,ss}.ttesthi_p(cc,ff,ttime),statstplv{ll,tt,ss}.ttesthi_ci(cc,ff,ttime,1,:),statstplv{ll,tt,ss}.ttesthi_stats{cc,ff,ttime,1}]=ttest2(squeeze(tmphi_tacPaud_comb2plv(cc,ff,ttime,:)),squeeze(tmphi_tacMSpN_comb2plv(cc,ff,ttime,:)),'alpha',alpha);
                    [statstplv{ll,tt,ss}.kstesthi_h(cc,ff,ttime),statstplv{ll,tt,ss}.kstesthi_p(cc,ff,ttime),statstplv{ll,tt,ss}.kstesthi_stats{cc,ff,ttime,1}]=kstest2(abs(squeeze(tmphi_tacPaud_comb2plv(cc,ff,ttime,:))),abs(squeeze(tmphi_tacMSpN_comb2plv(cc,ff,ttime,:))));
                    if numtests>=10
                      [statstplv{ll,tt,ss}.ttesthi_h(cc,ff,ttime,2),statstplv{ll,tt,ss}.ttesthi_p(cc,ff,ttime,2),statstplv{ll,tt,ss}.ttesthi_ci(cc,ff,ttime,2,:),statstplv{ll,tt,ss}.ttesthi_stats{cc,ff,ttime,2}]=  ttest2(squeeze(tmphi_tacPaud_comb2plv(cc,ff,ttime,1:halftest)),squeeze(tmphi_tacMSpN_comb2plv(cc,ff,ttime,1:halftest)),'alpha',alpha);
                      [statstplv{ll,tt,ss}.kstesthi_h(cc,ff,ttime,2),statstplv{ll,tt,ss}.kstesthi_p(cc,ff,ttime,2),statstplv{ll,tt,ss}.kstesthi_stats{cc,ff,ttime,2}]=  kstest2(abs(squeeze(tmphi_tacPaud_comb2plv(cc,ff,ttime,1:halftest))),abs(squeeze(tmphi_tacMSpN_comb2plv(cc,ff,ttime,1:halftest))));
                      [statstplv{ll,tt,ss}.ttesthi_h(cc,ff,ttime,3),statstplv{ll,tt,ss}.ttesthi_p(cc,ff,ttime,3),statstplv{ll,tt,ss}.ttesthi_ci(cc,ff,ttime,3,:),statstplv{ll,tt,ss}.ttesthi_stats{cc,ff,ttime,3}]=  ttest2(squeeze(tmphi_tacPaud_comb2plv(cc,ff,ttime,halftest+1:end)),squeeze(tmphi_tacMSpN_comb2plv(cc,ff,ttime,halftest+1:end)),'alpha',alpha);
                      [statstplv{ll,tt,ss}.kstesthi_h(cc,ff,ttime,3),statstplv{ll,tt,ss}.kstesthi_p(cc,ff,ttime,3),statstplv{ll,tt,ss}.kstesthi_stats{cc,ff,ttime,3}]=  kstest2(abs(squeeze(tmphi_tacPaud_comb2plv(cc,ff,ttime,halftest+1:end))),abs(squeeze(tmphi_tacMSpN_comb2plv(cc,ff,ttime,halftest+1:end))));
                    end
                  end
                end
              end
              disp(['finished stat cc' num2str(cc)]),toc
            end % cc
          end
          clear tmplo_* tmphi_*
          
          if synchasynch && ll<5
            numtsynch_trials(ll,tt,ss)=size(tlock_tacMSshift1{ll,tt,ss}.trialinfo,1);
            % (AT70 + TA70)
            numcomb=numtsynch_trials(ll,tt,ss)^2;
            numtests=min(numcomb,minnumcomb);
            %             combindex=reshape(1:numcomb,numt_trials(ll,tt,ss),numt_trials(ll,tt,ss));
            freqlo_tMSasynch_tmp{1}=[];
            freqhi_tMSasynch_tmp{1}=[];
            for cc=1:numtests
              %               tmp=Shuffle(combindex(:));
              %               combuse=tmp(1:numt_trials(ll,tt,ss));
              combuse=Shuffle(1:numtsynch_trials(ll,tt,ss));
              
              tlock_fake=tlock_tacMSshift1{ll,tt,ss};
              for at=1:numtsynch_trials(ll,tt,ss)
                if cc==1  % do it as 'normal'
                  tind=at;aind=at;
                else
                  %                   [tind,aind]=find(combindex==combuse(at));
                  tind=at;
                  aind=combuse(at);
                end
                cfg=[];
                cfg.trials=tind;
                tmpt=ft_selectdata(cfg,tlock_tacMSshift1{ll,tt,ss});
                cfg.trials=aind;
                tmpa=ft_selectdata(cfg,tlock_tacMSshift2{ll,tt,ss});
                cfg=[];
                cfg.operation='add';
                cfg.parameter='trial';
                tmpsum=ft_math(cfg,tmpt,tmpa);
                if at==1
                  tlock_fake=tmpsum;
                end
                tlock_fake.trial(at,:,:)=tmpsum.trial(1,:,:);
              end
              cfg=[];
              cfg.method='mtmconvol';
              cfg.pad=4;
              cfg.foi=4:2:30;
              cfg.taper='hanning';
              cfg.toi=-1.2:timestepplv:1.3;
              cfg.t_ftimwin=4./cfg.foi;
              cfg.output='fourier';
              freqlo_tMSasynch_tmp{cc}=ft_freqanalysis(cfg,tlock_fake);
              freqlo_tMSasynch_tmp{cc}.powspctrm=squeeze(mean(abs(freqlo_tMSasynch_tmp{cc}.fourierspctrm).^2,1));
              freqlo_tMSasynch_tmp{cc}.plvspctrm=getplv(freqlo_tMSasynch_tmp{cc}.fourierspctrm);
              freqlo_tMSasynch_tmp{cc}=rmfield(freqlo_tMSasynch_tmp{cc},'fourierspctrm');
              cfg=[];
              cfg.method='mtmconvol';
              cfg.pad=4;
              cfg.foi=[30:5:45 55:5:80];
              cfg.taper='dpss';
              cfg.tapsmofrq=7*ones(1,length(cfg.foi));
              cfg.toi=-1.2:timestepplv:1.3;
              cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
              cfg.output='fourier';
              freqhi_tMSasynch_tmp{cc}=ft_freqanalysis(cfg,tlock_fake);
              freqhi_tMSasynch_tmp{cc}.powspctrm=squeeze(mean(abs(freqhi_tMSasynch_tmp{cc}.fourierspctrm).^2,1));
              freqhi_tMSasynch_tmp{cc}.plvspctrm=getplv(freqhi_tMSasynch_tmp{cc}.fourierspctrm);
              freqhi_tMSasynch_tmp{cc}=rmfield(freqhi_tMSasynch_tmp{cc},'fourierspctrm');
            end
            freqlo_tMSasynch_comb{ll,tt,ss,1}=freqlo_tMSasynch_tmp{1};
            clear tmp tmp2 tmp3
            if numtests>1
              for cc=1:numtests,tmp(:,:,:,cc)=freqlo_tMSasynch_tmp{cc}.powspctrm;end
              for cc=1:numtests,tmp3(:,:,:,cc)=freqlo_tMSasynch_tmp{cc}.plvspctrm;end
              if exist('tmp','var')
                freqlo_tMSasynch_comb{ll,tt,ss,2}=freqlo_tMSasynch_tmp{1};
                freqlo_tMSasynch_comb{ll,tt,ss,2}.powspctrm=mean(tmp,4);
                freqlo_tMSasynch_comb{ll,tt,ss,2}.stdpow=std(tmp,[],4);
                freqlo_tMSasynch_comb{ll,tt,ss,2}.plvspctrm=mean(tmp3,4);
                freqlo_tMSasynch_comb{ll,tt,ss,2}.stdplv=std(tmp3,[],4);
                if numtests>=10
                  halftest=round(numtests/2);
                  freqlo_tMSasynch_comb{ll,tt,ss,3}=freqlo_tMSasynch_tmp{1};
                  freqlo_tMSasynch_comb{ll,tt,ss,4}=freqlo_tMSasynch_tmp{1};
                  freqlo_tMSasynch_comb{ll,tt,ss,3}.powspctrm=mean(tmp(:,:,:,1:halftest),4);
                  freqlo_tMSasynch_comb{ll,tt,ss,3}.stdpow=std(tmp(:,:,:,1:halftest),[],4);
                  freqlo_tMSasynch_comb{ll,tt,ss,3}.plvspctrm=mean(tmp3(:,:,:,1:halftest),4);
                  freqlo_tMSasynch_comb{ll,tt,ss,3}.stdplv=std(tmp3(:,:,:,1:halftest),[],4);
                  freqlo_tMSasynch_comb{ll,tt,ss,4}.powspctrm=mean(tmp(:,:,:,halftest+1:end),4);
                  freqlo_tMSasynch_comb{ll,tt,ss,4}.stdpow=std(tmp(:,:,:,halftest+1:end),[],4);
                  freqlo_tMSasynch_comb{ll,tt,ss,4}.plvspctrm=mean(tmp3(:,:,:,halftest+1:end),4);
                  freqlo_tMSasynch_comb{ll,tt,ss,4}.stdplv=std(tmp3(:,:,:,halftest+1:end),[],4);
                end
              end
            else
              freqlo_tMSasynch_comb{ll,tt,ss,2}=[];
            end
            tmplo_tMSasynch_comb2=tmp;
            tmplo_tMSasynch_comb2plv=tmp3;
            freqhi_tMSasynch_comb{ll,tt,ss,1}=freqhi_tMSasynch_tmp{1};
            clear tmp tmp2 tmp3
            if numtests>1
              for cc=1:numtests,tmp(:,:,:,cc)=freqhi_tMSasynch_tmp{cc}.powspctrm;end
              for cc=1:numtests,tmp3(:,:,:,cc)=freqhi_tMSasynch_tmp{cc}.plvspctrm;end
              if exist('tmp','var')
                freqhi_tMSasynch_comb{ll,tt,ss,2}=freqhi_tMSasynch_tmp{1};
                freqhi_tMSasynch_comb{ll,tt,ss,2}.powspctrm=mean(tmp,4);
                freqhi_tMSasynch_comb{ll,tt,ss,2}.stdpow=std(tmp,[],4);
                freqhi_tMSasynch_comb{ll,tt,ss,2}.plvspctrm=mean(tmp3,4);
                freqhi_tMSasynch_comb{ll,tt,ss,2}.stdplv=std(tmp3,[],4);
                if numtests>=10
                  halftest=round(numtests/2);
                  freqhi_tMSasynch_comb{ll,tt,ss,3}=freqhi_tMSasynch_tmp{1};
                  freqhi_tMSasynch_comb{ll,tt,ss,4}=freqhi_tMSasynch_tmp{1};
                  freqhi_tMSasynch_comb{ll,tt,ss,3}.powspctrm=mean(tmp(:,:,:,1:halftest),4);
                  freqhi_tMSasynch_comb{ll,tt,ss,3}.stdpow=std(tmp(:,:,:,1:halftest),[],4);
                  freqhi_tMSasynch_comb{ll,tt,ss,3}.plvspctrm=mean(tmp3(:,:,:,1:halftest),4);
                  freqhi_tMSasynch_comb{ll,tt,ss,3}.stdplv=std(tmp3(:,:,:,1:halftest),[],4);
                  freqhi_tMSasynch_comb{ll,tt,ss,4}.powspctrm=mean(tmp(:,:,:,halftest+1:end),4);
                  freqhi_tMSasynch_comb{ll,tt,ss,4}.stdpow=std(tmp(:,:,:,halftest+1:end),[],4);
                  freqhi_tMSasynch_comb{ll,tt,ss,4}.plvspctrm=mean(tmp3(:,:,:,halftest+1:end),4);
                  freqhi_tMSasynch_comb{ll,tt,ss,4}.stdplv=std(tmp3(:,:,:,halftest+1:end),[],4);
                end
              end
            else
              freqhi_tMSasynch_comb{ll,tt,ss,2}=[];
            end
            tmphi_tMSasynch_comb2=tmp;
            tmphi_tMSasynch_comb2plv=tmp3;
            clear freqlo_tMSasynch_tmp freqhi_tMSasynch_tmp
            clear tmp tmp2 tmp3
            
            if ss==12 || ss==13
            else
              tlock_tacMSshift1{ll,tt,ss}=[]; % clearing to help save memory as running
              tlock_tacMSshift2{ll,tt,ss}=[]; % clearing to help save memory as running
            end
            
            % (TA0 + TA0_shifted70)
            freqlo_tMSsynch_tmp{1}=[];
            freqhi_tMSsynch_tmp{1}=[];
            for cc=1:numtests
              %               tmp=Shuffle(combindex(:));
              %               combuse=tmp(1:numt_trials(ll,tt,ss));
              combuse=Shuffle(1:numtsynch_trials(ll,tt,ss));
              
              tlock_fake=tlock_tacMSshift3{ll,tt,ss};
              for at=1:numtsynch_trials(ll,tt,ss)
                if cc==1  % do it as 'normal'
                  tind=at;aind=at;
                else
                  %                   [tind,aind]=find(combindex==combuse(at));
                  tind=at;
                  aind=combuse(at);
                end
                cfg=[];
                cfg.trials=tind;
                tmpt=ft_selectdata(cfg,tlock_tacMSshift3{ll,tt,ss});
                cfg.trials=aind;
                tmpa=ft_selectdata(cfg,tlock_tacMSshift4{ll,tt,ss});
                cfg=[];
                cfg.operation='add';
                cfg.parameter='trial';
                tmpsum=ft_math(cfg,tmpt,tmpa);
                if at==1
                  tlock_fake=tmpsum;
                end
                tlock_fake.trial(at,:,:)=tmpsum.trial(1,:,:);
              end
              cfg=[];
              cfg.method='mtmconvol';
              cfg.pad=4;
              cfg.foi=4:2:30;
              cfg.taper='hanning';
              cfg.toi=-1.2:timestepplv:1.3;
              cfg.t_ftimwin=4./cfg.foi;
              cfg.output='fourier';
              freqlo_tMSsynch_tmp{cc}=ft_freqanalysis(cfg,tlock_fake);
              freqlo_tMSsynch_tmp{cc}.powspctrm=squeeze(mean(abs(freqlo_tMSsynch_tmp{cc}.fourierspctrm).^2,1));
              freqlo_tMSsynch_tmp{cc}.plvspctrm=getplv(freqlo_tMSsynch_tmp{cc}.fourierspctrm);
              freqlo_tMSsynch_tmp{cc}=rmfield(freqlo_tMSsynch_tmp{cc},'fourierspctrm');
              cfg=[];
              cfg.method='mtmconvol';
              cfg.pad=4;
              cfg.foi=[30:5:45 55:5:80];
              cfg.taper='dpss';
              cfg.tapsmofrq=7*ones(1,length(cfg.foi));
              cfg.toi=-1.2:timestepplv:1.3;
              cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
              cfg.output='fourier';
              freqhi_tMSsynch_tmp{cc}=ft_freqanalysis(cfg,tlock_fake);
              freqhi_tMSsynch_tmp{cc}.powspctrm=squeeze(mean(abs(freqhi_tMSsynch_tmp{cc}.fourierspctrm).^2,1));
              freqhi_tMSsynch_tmp{cc}.plvspctrm=getplv(freqhi_tMSsynch_tmp{cc}.fourierspctrm);
              freqhi_tMSsynch_tmp{cc}=rmfield(freqhi_tMSsynch_tmp{cc},'fourierspctrm');
            end
            freqlo_tMSsynch_comb{ll,tt,ss,1}=freqlo_tMSsynch_tmp{1};
            clear tmp tmp2 tmp3
            if numtests>1
              for cc=1:numtests,tmp(:,:,:,cc)=freqlo_tMSsynch_tmp{cc}.powspctrm;end
              for cc=1:numtests,tmp3(:,:,:,cc)=freqlo_tMSsynch_tmp{cc}.plvspctrm;end
              if exist('tmp','var')
                freqlo_tMSsynch_comb{ll,tt,ss,2}=freqlo_tMSsynch_tmp{1};
                freqlo_tMSsynch_comb{ll,tt,ss,2}.powspctrm=mean(tmp,4);
                freqlo_tMSsynch_comb{ll,tt,ss,2}.stdpow=std(tmp,[],4);
                freqlo_tMSsynch_comb{ll,tt,ss,2}.plvspctrm=mean(tmp3,4);
                freqlo_tMSsynch_comb{ll,tt,ss,2}.stdplv=std(tmp3,[],4);
                if numtests>=10
                  halftest=round(numtests/2);
                  freqlo_tMSsynch_comb{ll,tt,ss,3}=freqlo_tMSsynch_tmp{1};
                  freqlo_tMSsynch_comb{ll,tt,ss,4}=freqlo_tMSsynch_tmp{1};
                  freqlo_tMSsynch_comb{ll,tt,ss,3}.powspctrm=mean(tmp(:,:,:,1:halftest),4);
                  freqlo_tMSsynch_comb{ll,tt,ss,3}.stdpow=std(tmp(:,:,:,1:halftest),[],4);
                  freqlo_tMSsynch_comb{ll,tt,ss,3}.plvspctrm=mean(tmp3(:,:,:,1:halftest),4);
                  freqlo_tMSsynch_comb{ll,tt,ss,3}.stdplv=std(tmp3(:,:,:,1:halftest),[],4);
                  freqlo_tMSsynch_comb{ll,tt,ss,4}.powspctrm=mean(tmp(:,:,:,halftest+1:end),4);
                  freqlo_tMSsynch_comb{ll,tt,ss,4}.stdpow=std(tmp(:,:,:,halftest+1:end),[],4);
                  freqlo_tMSsynch_comb{ll,tt,ss,4}.plvspctrm=mean(tmp3(:,:,:,halftest+1:end),4);
                  freqlo_tMSsynch_comb{ll,tt,ss,4}.stdplv=std(tmp3(:,:,:,halftest+1:end),[],4);
                end
              end
            else
              freqlo_tMSsynch_comb{ll,tt,ss,2}=[];
            end
            tmplo_tMSsynch_comb2=tmp;
            tmplo_tMSsynch_comb2plv=tmp3;
            freqhi_tMSsynch_comb{ll,tt,ss,1}=freqhi_tMSsynch_tmp{1};
            clear tmp tmp2 tmp3
            if numtests>1
              for cc=1:numtests,tmp(:,:,:,cc)=freqhi_tMSsynch_tmp{cc}.powspctrm;end
              for cc=1:numtests,tmp3(:,:,:,cc)=freqhi_tMSsynch_tmp{cc}.plvspctrm;end
              if exist('tmp','var')
                freqhi_tMSsynch_comb{ll,tt,ss,2}=freqhi_tMSsynch_tmp{1};
                freqhi_tMSsynch_comb{ll,tt,ss,2}.powspctrm=mean(tmp,4);
                freqhi_tMSsynch_comb{ll,tt,ss,2}.stdpow=std(tmp,[],4);
                freqhi_tMSsynch_comb{ll,tt,ss,2}.plvspctrm=mean(tmp3,4);
                freqhi_tMSsynch_comb{ll,tt,ss,2}.stdplv=std(tmp3,[],4);
                if numtests>=10
                  halftest=round(numtests/2);
                  freqhi_tMSsynch_comb{ll,tt,ss,3}=freqhi_tMSsynch_tmp{1};
                  freqhi_tMSsynch_comb{ll,tt,ss,4}=freqhi_tMSsynch_tmp{1};
                  freqhi_tMSsynch_comb{ll,tt,ss,3}.powspctrm=mean(tmp(:,:,:,1:halftest),4);
                  freqhi_tMSsynch_comb{ll,tt,ss,3}.stdpow=std(tmp(:,:,:,1:halftest),[],4);
                  freqhi_tMSsynch_comb{ll,tt,ss,3}.plvspctrm=mean(tmp3(:,:,:,1:halftest),4);
                  freqhi_tMSsynch_comb{ll,tt,ss,3}.stdplv=std(tmp3(:,:,:,1:halftest),[],4);
                  freqhi_tMSsynch_comb{ll,tt,ss,4}.powspctrm=mean(tmp(:,:,:,halftest+1:end),4);
                  freqhi_tMSsynch_comb{ll,tt,ss,4}.stdpow=std(tmp(:,:,:,halftest+1:end),[],4);
                  freqhi_tMSsynch_comb{ll,tt,ss,4}.plvspctrm=mean(tmp3(:,:,:,halftest+1:end),4);
                  freqhi_tMSsynch_comb{ll,tt,ss,4}.stdplv=std(tmp3(:,:,:,halftest+1:end),[],4);
                end
              end
            else
              freqhi_tMSsynch_comb{ll,tt,ss,2}=[];
            end
            tmphi_tMSsynch_comb2=tmp;
            tmphi_tMSsynch_comb2plv=tmp3;
            clear freqlo_tMSsynch_tmp freqhi_tMSsynch_tmp
            clear tmp tmp2 tmp3
            
            if ss==12 || ss==13
            else
              tlock_tacMSshift3{ll,tt,ss}=[]; % clearing to help save memory as running
              tlock_tacMSshift4{ll,tt,ss}=[]; % clearing to help save memory as running
            end
            
            % synch minus asynch
            for cc=1:size(tmplo_tMSsynch_comb2,1)
              for ff=1:size(tmplo_tMSsynch_comb2,2)
                for ttime=1:size(tmplo_tMSsynch_comb2,3)
                  if ~isnan(mean(squeeze(tmplo_tMSsynch_comb2(cc,ff,ttime,:)))) &&  ~isnan(mean(squeeze(tmplo_tMSasynch_comb2(cc,ff,ttime,:))))
                    [statsst{ll,tt,ss}.ttestlo_h(cc,ff,ttime),statsst{ll,tt,ss}.ttestlo_p(cc,ff,ttime),statsst{ll,tt,ss}.ttestlo_ci(cc,ff,ttime,1,:),statsst{ll,tt,ss}.ttestlo_stats{cc,ff,ttime,1}]=ttest2(squeeze(tmplo_tMSsynch_comb2(cc,ff,ttime,:)),squeeze(tmplo_tMSasynch_comb2(cc,ff,ttime,:)));
                    [statsst{ll,tt,ss}.kstestlo_h(cc,ff,ttime),statsst{ll,tt,ss}.kstestlo_p(cc,ff,ttime),statsst{ll,tt,ss}.kstestlo_stats{cc,ff,ttime,1}]=kstest2(squeeze(tmplo_tMSsynch_comb2(cc,ff,ttime,:)),squeeze(tmplo_tMSasynch_comb2(cc,ff,ttime,:)));
                    if numtests>=10
                      [statsst{ll,tt,ss}.ttestlo_h(cc,ff,ttime,2),statsst{ll,tt,ss}.ttestlo_p(cc,ff,ttime,2),statsst{ll,tt,ss}.ttestlo_ci(cc,ff,ttime,2,:),statsst{ll,tt,ss}.ttestlo_stats{cc,ff,ttime,2}]=  ttest2(squeeze(tmplo_tMSsynch_comb2(cc,ff,ttime,1:halftest)),squeeze(tmplo_tMSasynch_comb2(cc,ff,ttime,1:halftest)));
                      [statsst{ll,tt,ss}.kstestlo_h(cc,ff,ttime,2),statsst{ll,tt,ss}.kstestlo_p(cc,ff,ttime,2),statsst{ll,tt,ss}.kstestlo_stats{cc,ff,ttime,2}]=  kstest2(squeeze(tmplo_tMSsynch_comb2(cc,ff,ttime,1:halftest)),squeeze(tmplo_tMSasynch_comb2(cc,ff,ttime,1:halftest)));
                      [statsst{ll,tt,ss}.ttestlo_h(cc,ff,ttime,3),statsst{ll,tt,ss}.ttestlo_p(cc,ff,ttime,3),statsst{ll,tt,ss}.ttestlo_ci(cc,ff,ttime,3,:),statsst{ll,tt,ss}.ttestlo_stats{cc,ff,ttime,3}]=  ttest2(squeeze(tmplo_tMSsynch_comb2(cc,ff,ttime,halftest+1:end)),squeeze(tmplo_tMSasynch_comb2(cc,ff,ttime,halftest+1:end)));
                      [statsst{ll,tt,ss}.kstestlo_h(cc,ff,ttime,3),statsst{ll,tt,ss}.kstestlo_p(cc,ff,ttime,3),statsst{ll,tt,ss}.kstestlo_stats{cc,ff,ttime,3}]=  kstest2(squeeze(tmplo_tMSsynch_comb2(cc,ff,ttime,halftest+1:end)),squeeze(tmplo_tMSasynch_comb2(cc,ff,ttime,halftest+1:end)));
                    end
                    [statsstplv{ll,tt,ss}.ttestlo_h(cc,ff,ttime),statsstplv{ll,tt,ss}.ttestlo_p(cc,ff,ttime),statsstplv{ll,tt,ss}.ttestlo_ci(cc,ff,ttime,1,:),statsstplv{ll,tt,ss}.ttestlo_stats{cc,ff,ttime,1}]=ttest2(squeeze(tmplo_tMSsynch_comb2plv(cc,ff,ttime,:)),squeeze(tmplo_tMSasynch_comb2plv(cc,ff,ttime,:)));
                    [statsstplv{ll,tt,ss}.kstestlo_h(cc,ff,ttime),statsstplv{ll,tt,ss}.kstestlo_p(cc,ff,ttime),statsstplv{ll,tt,ss}.kstestlo_stats{cc,ff,ttime,1}]=kstest2(squeeze(tmplo_tMSsynch_comb2plv(cc,ff,ttime,:)),squeeze(tmplo_tMSasynch_comb2plv(cc,ff,ttime,:)));
                    if numtests>=10
                      [statsstplv{ll,tt,ss}.ttestlo_h(cc,ff,ttime,2),statsstplv{ll,tt,ss}.ttestlo_p(cc,ff,ttime,2),statsstplv{ll,tt,ss}.ttestlo_ci(cc,ff,ttime,2,:),statsstplv{ll,tt,ss}.ttestlo_stats{cc,ff,ttime,2}]=  ttest2(squeeze(tmplo_tMSsynch_comb2plv(cc,ff,ttime,1:halftest)),squeeze(tmplo_tMSasynch_comb2plv(cc,ff,ttime,1:halftest)));
                      [statsstplv{ll,tt,ss}.kstestlo_h(cc,ff,ttime,2),statsstplv{ll,tt,ss}.kstestlo_p(cc,ff,ttime,2),statsstplv{ll,tt,ss}.kstestlo_stats{cc,ff,ttime,2}]=  kstest2(squeeze(tmplo_tMSsynch_comb2plv(cc,ff,ttime,1:halftest)),squeeze(tmplo_tMSasynch_comb2plv(cc,ff,ttime,1:halftest)));
                      [statsstplv{ll,tt,ss}.ttestlo_h(cc,ff,ttime,3),statsstplv{ll,tt,ss}.ttestlo_p(cc,ff,ttime,3),statsstplv{ll,tt,ss}.ttestlo_ci(cc,ff,ttime,3,:),statsstplv{ll,tt,ss}.ttestlo_stats{cc,ff,ttime,3}]=  ttest2(squeeze(tmplo_tMSsynch_comb2plv(cc,ff,ttime,halftest+1:end)),squeeze(tmplo_tMSasynch_comb2plv(cc,ff,ttime,halftest+1:end)));
                      [statsstplv{ll,tt,ss}.kstestlo_h(cc,ff,ttime,3),statsstplv{ll,tt,ss}.kstestlo_p(cc,ff,ttime,3),statsstplv{ll,tt,ss}.kstestlo_stats{cc,ff,ttime,3}]=  kstest2(squeeze(tmplo_tMSsynch_comb2plv(cc,ff,ttime,halftest+1:end)),squeeze(tmplo_tMSasynch_comb2plv(cc,ff,ttime,halftest+1:end)));
                    end
                  end
                end
              end
              for ff=1:size(tmphi_tMSsynch_comb2,2)
                for ttime=1:size(tmphi_tMSsynch_comb2,3)
                  if ~isnan(mean(squeeze(tmphi_tMSsynch_comb2(cc,ff,ttime,:)))) &&  ~isnan(mean(squeeze(tmphi_tMSasynch_comb2(cc,ff,ttime,:))))
                    [statsst{ll,tt,ss}.ttesthi_h(cc,ff,ttime),statsst{ll,tt,ss}.ttesthi_p(cc,ff,ttime),statsst{ll,tt,ss}.ttesthi_ci(cc,ff,ttime,1,:),statsst{ll,tt,ss}.ttesthi_stats{cc,ff,ttime,1}]=ttest2(squeeze(tmphi_tMSsynch_comb2(cc,ff,ttime,:)),squeeze(tmphi_tMSasynch_comb2(cc,ff,ttime,:)));
                    [statsst{ll,tt,ss}.kstesthi_h(cc,ff,ttime),statsst{ll,tt,ss}.kstesthi_p(cc,ff,ttime),statsst{ll,tt,ss}.kstesthi_stats{cc,ff,ttime,1}]=kstest2(squeeze(tmphi_tMSsynch_comb2(cc,ff,ttime,:)),squeeze(tmphi_tMSasynch_comb2(cc,ff,ttime,:)));
                    if numtests>=10
                      [statsst{ll,tt,ss}.ttesthi_h(cc,ff,ttime,2),statsst{ll,tt,ss}.ttesthi_p(cc,ff,ttime,2),statsst{ll,tt,ss}.ttesthi_ci(cc,ff,ttime,2,:),statsst{ll,tt,ss}.ttesthi_stats{cc,ff,ttime,2}]=  ttest2(squeeze(tmphi_tMSsynch_comb2(cc,ff,ttime,1:halftest)),squeeze(tmphi_tMSasynch_comb2(cc,ff,ttime,1:halftest)));
                      [statsst{ll,tt,ss}.kstesthi_h(cc,ff,ttime,2),statsst{ll,tt,ss}.kstesthi_p(cc,ff,ttime,2),statsst{ll,tt,ss}.kstesthi_stats{cc,ff,ttime,2}]=  kstest2(squeeze(tmphi_tMSsynch_comb2(cc,ff,ttime,1:halftest)),squeeze(tmphi_tMSasynch_comb2(cc,ff,ttime,1:halftest)));
                      [statsst{ll,tt,ss}.ttesthi_h(cc,ff,ttime,3),statsst{ll,tt,ss}.ttesthi_p(cc,ff,ttime,3),statsst{ll,tt,ss}.ttesthi_ci(cc,ff,ttime,3,:),statsst{ll,tt,ss}.ttesthi_stats{cc,ff,ttime,3}]=  ttest2(squeeze(tmphi_tMSsynch_comb2(cc,ff,ttime,halftest+1:end)),squeeze(tmphi_tMSasynch_comb2(cc,ff,ttime,halftest+1:end)));
                      [statsst{ll,tt,ss}.kstesthi_h(cc,ff,ttime,3),statsst{ll,tt,ss}.kstesthi_p(cc,ff,ttime,3),statsst{ll,tt,ss}.kstesthi_stats{cc,ff,ttime,3}]=  kstest2(squeeze(tmphi_tMSsynch_comb2(cc,ff,ttime,halftest+1:end)),squeeze(tmphi_tMSasynch_comb2(cc,ff,ttime,halftest+1:end)));
                    end
                    [statsstplv{ll,tt,ss}.ttesthi_h(cc,ff,ttime),statsstplv{ll,tt,ss}.ttesthi_p(cc,ff,ttime),statsstplv{ll,tt,ss}.ttesthi_ci(cc,ff,ttime,1,:),statsstplv{ll,tt,ss}.ttesthi_stats{cc,ff,ttime,1}]=ttest2(squeeze(tmphi_tMSsynch_comb2plv(cc,ff,ttime,:)),squeeze(tmphi_tMSasynch_comb2plv(cc,ff,ttime,:)));
                    [statsstplv{ll,tt,ss}.kstesthi_h(cc,ff,ttime),statsstplv{ll,tt,ss}.kstesthi_p(cc,ff,ttime),statsstplv{ll,tt,ss}.kstesthi_stats{cc,ff,ttime,1}]=kstest2(squeeze(tmphi_tMSsynch_comb2plv(cc,ff,ttime,:)),squeeze(tmphi_tMSasynch_comb2plv(cc,ff,ttime,:)));
                    if numtests>=10
                      [statsstplv{ll,tt,ss}.ttesthi_h(cc,ff,ttime,2),statsstplv{ll,tt,ss}.ttesthi_p(cc,ff,ttime,2),statsstplv{ll,tt,ss}.ttesthi_ci(cc,ff,ttime,2,:),statsstplv{ll,tt,ss}.ttesthi_stats{cc,ff,ttime,2}]=  ttest2(squeeze(tmphi_tMSsynch_comb2plv(cc,ff,ttime,1:halftest)),squeeze(tmphi_tMSasynch_comb2plv(cc,ff,ttime,1:halftest)));
                      [statsstplv{ll,tt,ss}.kstesthi_h(cc,ff,ttime,2),statsstplv{ll,tt,ss}.kstesthi_p(cc,ff,ttime,2),statsstplv{ll,tt,ss}.kstesthi_stats{cc,ff,ttime,2}]=  kstest2(squeeze(tmphi_tMSsynch_comb2plv(cc,ff,ttime,1:halftest)),squeeze(tmphi_tMSasynch_comb2plv(cc,ff,ttime,1:halftest)));
                      [statsstplv{ll,tt,ss}.ttesthi_h(cc,ff,ttime,3),statsstplv{ll,tt,ss}.ttesthi_p(cc,ff,ttime,3),statsstplv{ll,tt,ss}.ttesthi_ci(cc,ff,ttime,3,:),statsstplv{ll,tt,ss}.ttesthi_stats{cc,ff,ttime,3}]=  ttest2(squeeze(tmphi_tMSsynch_comb2plv(cc,ff,ttime,halftest+1:end)),squeeze(tmphi_tMSasynch_comb2plv(cc,ff,ttime,halftest+1:end)));
                      [statsstplv{ll,tt,ss}.kstesthi_h(cc,ff,ttime,3),statsstplv{ll,tt,ss}.kstesthi_p(cc,ff,ttime,3),statsstplv{ll,tt,ss}.kstesthi_stats{cc,ff,ttime,3}]=  kstest2(squeeze(tmphi_tMSsynch_comb2plv(cc,ff,ttime,halftest+1:end)),squeeze(tmphi_tMSasynch_comb2plv(cc,ff,ttime,halftest+1:end)));
                    end
                  end
                end
              end
            end
            clear tmplo_* tmphi_*
            
          end % ll<5
          
          
          %           end % if ll<10
          disp(['ll_completed_' num2str(ll)]),toc
          if ll<9
            disp(['saving now temporary after ll ' num2str(ll) ' and ss ' num2str(ss)])
            if dofftadd
              save(['freq_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_tt' num2str(tt) '.mat'],'freq*comb','freq*fftadd','num*trials','-v7.3')
            else
              % NOTE: this below is with trialkc=0 implicit!!
              %                   if numtests>2
              %                     save(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_tt' num2str(tt) '.mat'],'freq*comb','stats*','num*trials','-v7.3')
              %                   else
              %                     save(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_tt' num2str(tt) '.mat'],'freq*comb','num*trials','-v7.3')
              %                   end
              if numtests>2
                save(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_tt' num2str(tt) '_trialkc' num2str(trialkc) '.mat'],'freq*comb','stats*','num*trials','-v7.3')
              else
%                 save(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_tt' num2str(tt) '_trialkc' num2str(trialkc) '.mat'],'freq*comb','num*trials','-v7.3')
                save(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'freq*comb','num*trials','-v7.3')
              end
            end
          end
        end % ll
        
      end % if exclude ss
    end % ss
    clear tlock_tac tlock_aud tlock_nul tr tlock*tlock tlock_tacMSshift*
    if ~phaset0
      if dofftadd
        save(['freq_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_tt' num2str(tt) '.mat'],'freq*comb','freq*fftadd','num*trials','-v7.3')
      else
        % NOTE: this below is with trialkc=0 implicit!!
        %             if numtests>2
        %               save(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_tt' num2str(tt) '.mat'],'freq*comb','stats*','num*trials','-v7.3')
        %             else
        %               save(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_tt' num2str(tt) '.mat'],'freq*comb','num*trials','-v7.3')
        %             end
        if numtests>2
          save(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_tt' num2str(tt) '_trialkc' num2str(trialkc) '.mat'],'freq*comb','stats*','num*trials','-v7.3')
        else
%           save(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_tt' num2str(tt) '_trialkc' num2str(trialkc) '.mat'],'freq*comb','num*trials','-v7.3')
          save(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'freq*comb','num*trials','-v7.3')
        end
      end
    end
    clear freq*comb freq*fftadd num*trials freq*time0
  end % iter
  %   end % tt
  
  continue
  
  %% Do AudTac
  clearvars -except ii sub *dir ii*use sleep minnumcomb raw_* hostname timestep* soades dofftadd statsa
  if 0 % temporarily comment out as this worked correctly.
    tacaud=0;
    %     if 0 % don't do this tacaud=0 for now
    for tt=[3]
      %       [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection2_wakeSleep_sepTacAud(ii,sleep,tt,tacaud);
      [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection2_wakeSleep_sepTacAud(ii,sleep,tt,tacaud,raw_tac, raw_aud, raw_nul);
      %               if ll==max(soalist)
      tlock_tac{10,tt,12}=[];
      tlock_tac{10,tt,13}=[];
      %               end
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
        if ~any([tr.a10trialkept{:,tt,ss}])
          ssuse=setdiff(ssuse,ss);
          continue
        else
          %   for tt=1:4
          
          try
            fsample=1/diff(tlock_tac{10,tt,ss}.time(1:2)); % sampling rate
          catch
            try
              fsample=1/diff(tlock_aud{10,tt,ss}.time(1:2)); % sampling rate
            catch
              fsample=1/diff(tlock_nul{10,tt,ss}.time(1:2)); % sampling rate
            end
          end
          
          %   for tt=1:4
          
          % Multisensory-shifted comparison
          for ll=[1 3 4] % need to do this before other loops as the fields get deleted as it loops through ll
            if ss==23
              if length(tr.amstrialkept1{ll,tt,12})>=2 && length(tr.amstrialkept1{ll,tt,13})<2
                cfg=[];
                cfg.trials=tr.amstrialkept1{ll,tt,12};
                tlock_audMSshift1{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{ll,tt,12});
                cfg=[];
                cfg.trials=tr.amstrialkept2{ll,tt,12};
                tlock_audMSshift2{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{10-ll,tt,12});
                cfg=[];
                cfg.trials=tr.amstrialkept3{ll,tt,12};
                tlock_audMSshift3{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{5,tt,12});
                cfg=[];
                cfg.trials=tr.amstrialkept4{ll,tt,12};
                tlock_audMSshift4{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{5,tt,12});
              elseif length(tr.amstrialkept1{ll,tt,12})<2 && length(tr.amstrialkept1{ll,tt,13})>=2
                cfg=[];
                cfg.trials=tr.amstrialkept1{ll,tt,13};
                tlock_audMSshift1{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{ll,tt,13});
                cfg=[];
                cfg.trials=tr.amstrialkept2{ll,tt,13};
                tlock_audMSshift2{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{10-ll,tt,13});
                cfg=[];
                cfg.trials=tr.amstrialkept3{ll,tt,13};
                tlock_audMSshift3{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{5,tt,13});
                cfg=[];
                cfg.trials=tr.amstrialkept4{ll,tt,13};
                tlock_audMSshift4{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{5,tt,13});
              elseif length(tr.amstrialkept1{ll,tt,12})>=2 && length(tr.amstrialkept1{ll,tt,13})>=2
                cfg=[];
                cfg.trials=tr.amstrialkept1{ll,tt,12};
                tmp12=ft_selectdata(cfg,tlock_aud{ll,tt,12});
                cfg=[];
                cfg.trials=tr.amstrialkept1{ll,tt,13};
                tmp13=ft_selectdata(cfg,tlock_aud{ll,tt,13});
                cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
                cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
                tlock_audMSshift1{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
                cfg=[];
                cfg.trials=tr.amstrialkept2{ll,tt,12};
                tmp12=ft_selectdata(cfg,tlock_aud{10-ll,tt,12});
                cfg=[];
                cfg.trials=tr.amstrialkept2{ll,tt,13};
                tmp13=ft_selectdata(cfg,tlock_aud{10-ll,tt,13});
                cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
                cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
                tlock_audMSshift2{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
                cfg=[];
                cfg.trials=tr.amstrialkept3{ll,tt,12};
                tmp12=ft_selectdata(cfg,tlock_aud{5,tt,12});
                cfg=[];
                cfg.trials=tr.amstrialkept3{ll,tt,13};
                tmp13=ft_selectdata(cfg,tlock_aud{5,tt,13});
                cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
                cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
                tlock_audMSshift3{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
                cfg=[];
                cfg.trials=tr.amstrialkept4{ll,tt,12};
                tmp12=ft_selectdata(cfg,tlock_aud{5,tt,12});
                cfg=[];
                cfg.trials=tr.amstrialkept4{ll,tt,13};
                tmp13=ft_selectdata(cfg,tlock_aud{5,tt,13});
                cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
                cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
                tlock_audMSshift4{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
              else
                %               tlock_aud{ll,tt,ss}=[];
                continue
              end
            else
              if length(tr.amstrialkept1{ll,tt,ss})>=2 && length(tr.amstrialkept2{ll,tt,ss})>=2 && length(tr.amstrialkept3{ll,tt,ss})>=2 && length(tr.amstrialkept4{ll,tt,ss})>=2
                cfg=[];
                cfg.trials=tr.amstrialkept1{ll,tt,ss};
                tlock_audMSshift1{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{ll,tt,ss});
                cfg=[];
                cfg.trials=tr.amstrialkept2{ll,tt,ss};
                tlock_audMSshift2{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{10-ll,tt,ss});
                cfg=[];
                cfg.trials=tr.amstrialkept3{ll,tt,ss};
                tlock_audMSshift3{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{5,tt,ss});
                cfg=[];
                cfg.trials=tr.amstrialkept4{ll,tt,ss};
                tlock_audMSshift4{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{5,tt,ss});
              else
                tlock_audMSshift1{ll,tt,ss}=[];
                tlock_audMSshift2{ll,tt,ss}=[];
                tlock_audMSshift3{ll,tt,ss}=[];
                tlock_audMSshift4{ll,tt,ss}=[];
              end
            end
            % do time-shifting here
            if ~isempty(tlock_audMSshift1{ll,tt,ss})
              warning off
              cfg=[];
              cfg.offset=round(fsample*(-soades(ll)));
              tlock_audMSshift2{ll,tt,ss}=ft_redefinetrial(cfg,tlock_audMSshift2{ll,tt,ss}); % Aud first shifted so that Aud at time 0
              tlock_audMSshift4{ll,tt,ss}=ft_redefinetrial(cfg,tlock_audMSshift4{ll,tt,ss}); % Simult shifted so that both at time of second stim of 'll' condition
              warning on
            end
          end %ll
          
          for ll=[soalist soalist+20 soalist+40]
            
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
                
                if length(tr.a10trialkept{ll,tt,12})<2 && length(tr.a10trialkept{ll,tt,13})>1
                  cfg=[];
                  cfg.trials=tr.a10trialkept{ll,tt,13};
                  tlock_aud{ll+20,tt,ss}=ft_selectdata(cfg,tlock_aud{10,tt,13}); % create +20 here as 10 will be different for every ll,tt
                elseif length(tr.a10trialkept{ll,tt,13})<2 && length(tr.a10trialkept{ll,tt,12})>1
                  cfg=[];
                  cfg.trials=tr.a10trialkept{ll,tt,12};
                  tlock_aud{ll+20,tt,ss}=ft_selectdata(cfg,tlock_aud{10,tt,12}); % create +20 here as 10 will be different for every ll,tt
                elseif length(tr.a10trialkept{ll,tt,13})>1 && length(tr.a10trialkept{ll,tt,12})>1
                  if ispc
                    S=memory; % memory doesn't exist on linux
                    gbavail=S.MemAvailableAllArrays/(1024^3);
                  else
                    [aa,bb]=system('head /proc/meminfo');
                    if ~isempty(strfind(hostname,'les')) % les-linux-fs3
                      gbavail=str2num(bb(strfind(bb,'MemFree')+8:strfind(bb,'MemFree')+23))/(1024^2);
                    elseif ~isempty(strfind(hostname,'LES')) % COLLES-151401
                      gb1=str2num(bb(strfind(bb,'Inactive:')+9:strfind(bb,'Inactive:')+23))/(1024^2);
                      gb2=str2num(bb(strfind(bb,'MemFree')+8:strfind(bb,'MemFree')+23))/(1024^2);
                      gbavail=max(gb1,gb2);
                    end
                  end
                  if gbavail < 1.5
                    save tmptlock.mat tlock_tac -v7.3
                    clear tlock_tac
                  end
                  cfg=[];
                  cfg.trials=tr.a10trialkept{ll,tt,12};
                  tmp12=ft_selectdata(cfg,tlock_aud{10,tt,12}); % create +20 here as 10 will be different for every ll,tt
                  cfg.trials=tr.a10trialkept{ll,tt,13};
                  tmp13=ft_selectdata(cfg,tlock_aud{10,tt,13}); % create +20 here as 10 will be different for every ll,tt
                  cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
                  cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
                  tlock_aud{ll+20,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
                else
                  numa_trials(ll,tt,ss)=0;
                  tlock_aud{ll+20,tt,ss}=[];
                  continue
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
                  if ispc
                    S=memory; % memory doesn't exist on linux
                    gbavail=S.MemAvailableAllArrays/(1024^3);
                  else
                    [aa,bb]=system('head /proc/meminfo');
                    if ~isempty(strfind(hostname,'les')) % les-linux-fs3
                      gbavail=str2num(bb(strfind(bb,'MemFree')+8:strfind(bb,'MemFree')+23))/(1024^2);
                    elseif ~isempty(strfind(hostname,'LES')) % COLLES-151401
                      gb1=str2num(bb(strfind(bb,'Inactive:')+9:strfind(bb,'Inactive:')+23))/(1024^2);
                      gb2=str2num(bb(strfind(bb,'MemFree')+8:strfind(bb,'MemFree')+23))/(1024^2);
                      gbavail=max(gb1,gb2);
                    end
                  end
                  if gbavail < 1.5
                    save tmptlock.mat tlock_tac -v7.3
                    clear tlock_tac
                  end
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
                  numa_trials(ll,tt,ss)=0;
                  continue
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
                  if ll==max(soalist) && ss<12 && ~phaset0 % we use it later in phaset0 for PCA
                    tlock_tac{10,tt,ss}=[];
                  end
                  
                  cfg=[];
                  cfg.trials=tr.a10trialkept{ll,tt,ss};
                  if length(cfg.trials)<2
                    tlock_aud{ll+20,tt,ss}=[]; % create +20 here as 10 will be different for every ll,tt
                  else
                    tlock_aud{ll+20,tt,ss}=ft_selectdata(cfg,tlock_aud{10,tt,ss}); % create +20 here as 10 will be different for every ll,tt
                  end
                  if ll==max(soalist) && ss<12 && ~phaset0 % we use it later in phaset0 for PCA
                    tlock_aud{10,tt,ss}=[];
                  end
                  
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
              if isempty(tlock_aud{ll+20,tt,ss})
                numa_trials(ll,tt,ss)=0;
              else
                numa_trials(ll,tt,ss)=size(tlock_aud{ll+20,tt,ss}.trial,1); % this should be the same for audll20, tacll40, audll, nulll60
              end
            end
            if ~exist('tlock_tac','var')
              load tmptlock.mat
              delete tmptlock.mat
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
                if ispc
                  S=memory; % memory doesn't exist on linux
                  gbavail=S.MemAvailableAllArrays/(1024^3);
                else
                  [aa,bb]=system('head /proc/meminfo');
                  if ~isempty(strfind(hostname,'les')) % les-linux-fs3
                    gbavail=str2num(bb(strfind(bb,'MemFree')+8:strfind(bb,'MemFree')+23))/(1024^2);
                  elseif ~isempty(strfind(hostname,'LES')) % COLLES-151401
                    gb1=str2num(bb(strfind(bb,'Inactive:')+9:strfind(bb,'Inactive:')+23))/(1024^2);
                    gb2=str2num(bb(strfind(bb,'MemFree')+8:strfind(bb,'MemFree')+23))/(1024^2);
                    gbavail=max(gb1,gb2);
                  end
                end
                if gbavail< 1.5
                  save tmptlock.mat tlock_nul -v7.3
                  clear tlock_nul
                end
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
                if ispc
                  S=memory; % memory doesn't exist on linux
                  gbavail=S.MemAvailableAllArrays/(1024^3);
                else
                  [aa,bb]=system('head /proc/meminfo');
                  if ~isempty(strfind(hostname,'les')) % les-linux-fs3
                    gbavail=str2num(bb(strfind(bb,'MemFree')+8:strfind(bb,'MemFree')+23))/(1024^2);
                  elseif ~isempty(strfind(hostname,'LES')) % COLLES-151401
                    gb1=str2num(bb(strfind(bb,'Inactive:')+9:strfind(bb,'Inactive:')+23))/(1024^2);
                    gb2=str2num(bb(strfind(bb,'MemFree')+8:strfind(bb,'MemFree')+23))/(1024^2);
                    gbavail=max(gb1,gb2);
                  end
                end
                if gbavail < 1.5
                  save tmptlock.mat tlock_nul -v7.3
                  clear tlock_nul
                end
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
            
            if ~exist('tlock_nul','var')
              load tmptlock.mat
              delete tmptlock.mat
            end
            
            % up to here same as for ERP
          end % ll
          
          %%
          
          % Do it two ways:
          
          if dofftadd
            % %%%%      % FFT then add    %%%%%%%%%%%%%
            % Way 1:  FFT+norm of each condition, then add conditions
            
            for ll=[soalist soalist+20]
              %           if 0
              if isfield(tlock_aud{ll,tt,ss},'trial') && size(tlock_aud{ll,tt,ss}.trial,1)
                %         cfg=[];
                %         cfg.lpfilter='yes';
                %         cfg.lpfreq=40;
                %         cfg.demean='yes';
                %         cfg.baselinewindow=[-1.7 -0.6];
                cfg=[];
                cfg.method='mtmconvol';
                cfg.pad=4;
                cfg.foi=4:2:30;
                cfg.taper='hanning';
                cfg.toi=-1.2:timestepfft:1.3;
                cfg.t_ftimwin=4./cfg.foi;
                cfg.output='powandcsd';
                %         cfg.keeptrials='yes';  % not necessary at this stage for keeping trials (yes later for source analysis and phase resetting)
                freqlo_aud{ll,tt,ss}=ft_freqanalysis(cfg,tlock_aud{ll,tt,ss});
                
                % Higher frequencies with dpss multitaper
                cfg=[];
                cfg.method='mtmconvol';
                cfg.pad=4;
                cfg.foi=[30:5:45 55:5:80];
                cfg.taper='dpss';
                cfg.tapsmofrq=7*ones(1,length(cfg.foi));
                cfg.toi=-1.2:timestepfft:1.3;
                cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
                cfg.output='powandcsd';
                freqhi_aud{ll,tt,ss}=ft_freqanalysis(cfg,tlock_aud{ll,tt,ss});
              else
                freqlo_aud{ll,tt,ss}=[];
                freqhi_aud{ll,tt,ss}=[];
              end
            end
            %       tlock_aud{ll,tt,ss}=[]; % clearing to help save memory as running
            for ll=soalist+40
              if numa_trials(ll-40,tt,ss) && isfield(tlock_tac{ll,tt,ss},'trial') && size(tlock_tac{ll,tt,ss}.trial,1)
                cfg=[];
                cfg.method='mtmconvol';
                cfg.pad=4;
                cfg.foi=4:2:30;
                cfg.taper='hanning';
                cfg.toi=-1.2:timestepfft:1.3;
                cfg.t_ftimwin=4./cfg.foi;
                cfg.output='powandcsd';
                freqlo_tac{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tac{ll,tt,ss});
                cfg=[];
                cfg.method='mtmconvol';
                cfg.pad=4;
                cfg.foi=[30:5:45 55:5:80];
                cfg.taper='dpss';
                cfg.tapsmofrq=7*ones(1,length(cfg.foi));
                cfg.toi=-1.2:timestepfft:1.3;
                cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
                cfg.output='powandcsd';
                freqhi_tac{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tac{ll,tt,ss});
              else
                freqlo_tac{ll,tt,ss}=[];
                freqhi_tac{ll,tt,ss}=[];
              end
              %       tlock_tac{ll,tt,ss}=[]; % clearing to help save memory as running
              %           end
            end
            for ll=[soalist+60]
              if numa_trials(ll-60,tt,ss) && isfield(tlock_nul{ll,tt,ss},'trial') && size(tlock_nul{ll,tt,ss}.trial,1)
                cfg=[];
                cfg.method='mtmconvol';
                cfg.pad=4;
                cfg.foi=4:2:30;
                cfg.taper='hanning';
                cfg.toi=-1.2:timestepfft:1.3;
                cfg.t_ftimwin=4./cfg.foi;
                cfg.output='powandcsd';
                freqlo_nul{ll,tt,ss}=ft_freqanalysis(cfg,tlock_nul{ll,tt,ss});
                cfg=[];
                cfg.method='mtmconvol';
                cfg.pad=4;
                cfg.foi=[30:5:45 55:5:80];
                cfg.taper='dpss';
                cfg.tapsmofrq=7*ones(1,length(cfg.foi));
                cfg.toi=-1.2:timestepfft:1.3;
                cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
                cfg.output='powandcsd';
                freqhi_nul{ll,tt,ss}=ft_freqanalysis(cfg,tlock_nul{ll,tt,ss});
              else
                freqlo_nul{ll,tt,ss}=[];
                freqhi_nul{ll,tt,ss}=[];
              end
              %       tlock_nul{ll,tt,ss}=[]; % clearing to help save memory as running
            end
            
            for ll=soalist
              
              % create sum of unisensory conditions
              if numa_trials(ll,tt,ss)
                cfg=[];
                cfg.operation='add';
                cfg.parameter='powspctrm';
                freqlo_audPtac_fftadd{ll,tt,ss}=ft_math(cfg,freqlo_aud{ll+20,tt,ss},freqlo_tac{ll+40,tt,ss}); % keeping aud centred but jitter tac (to match tlock_aud{X,tt,ss})
                freqhi_audPtac_fftadd{ll,tt,ss}=ft_math(cfg,freqhi_aud{ll+20,tt,ss},freqhi_tac{ll+40,tt,ss}); % keeping aud centred but jitter tac (to match tlock_aud{X,tt,ss})
                freqlo_aTacAlone_fftadd{ll,tt,ss}=freqlo_tac{ll+40,tt,ss};
                freqhi_aTacAlone_fftadd{ll,tt,ss}=freqhi_tac{ll+40,tt,ss};
                freqlo_aAudAlone_fftadd{ll,tt,ss}=freqlo_aud{ll+20,tt,ss};
                freqhi_aaudAlone_fftadd{ll,tt,ss}=freqhi_aud{ll+20,tt,ss};
                cfg.parameter='crsspctrm';
                tmp=ft_math(cfg,freqlo_aud{ll+20,tt,ss},freqlo_tac{ll+40,tt,ss});
                freqlo_audPtac_fftadd{ll,tt,ss}.crsspctrm=tmp.crsspctrm;
                freqlo_aAudAlone_fftadd{ll,tt,ss}.crsspctrm=freqlo_aud{ll+20,tt,ss}.crsspctrm;
                freqlo_aTacAlone_fftadd{ll,tt,ss}.crsspctrm=freqlo_tac{ll+40,tt,ss}.crsspctrm;
                tmp=ft_math(cfg,freqhi_aud{ll+20,tt,ss},freqhi_tac{ll+40,tt,ss});
                freqhi_audPtac_fftadd{ll,tt,ss}.crsspctrm=tmp.crsspctrm;
                freqhi_aAudAlone_fftadd{ll,tt,ss}.crsspctrm=freqhi_aud{ll+20,tt,ss}.crsspctrm;
                freqhi_aTacAlone_fftadd{ll,tt,ss}.crsspctrm=freqhi_tac{ll+40,tt,ss}.crsspctrm;
              else
                freqlo_audPtac_fftadd{ll,tt,ss}=[];
                freqhi_audPtac_fftadd{ll,tt,ss}=[];
                freqlo_aTacAlone_fftadd{ll,tt,ss}=[];
                freqhi_aTacAlone_fftadd{ll,tt,ss}=[];
                freqlo_aAudAlone_fftadd{ll,tt,ss}=[];
                freqhi_aaudAlone_fftadd{ll,tt,ss}=[];
              end
              
              
              % create sum of mulitsensory and null conditions
              %       freqlo_tacMSpN{ll,tt,ss}=ft_math(cfg,ft_selectdata(cfgavg,freqlo_tac{ll,tt,ss}),ft_selectdata(cfgavg,freqlo_nul{ll+40,tt,ss})); % keeping tac centred but jitter aud (to match tlock_tac{X,tt,ss})
              %       freqlo_audMSpN{ll,tt,ss}=ft_math(cfg,ft_selectdata(cfgavg,freqlo_aud{ll,tt,ss}),ft_selectdata(cfgavg,freqlo_nul{ll+40,tt,ss})); % keeping aud centred but jitter tac (to match tlock_aud{X,tt,ss})
              %       freqhi_tacMSpN{ll,tt,ss}=ft_math(cfg,ft_selectdata(cfgavg,freqhi_tac{ll,tt,ss}),ft_selectdata(cfgavg,freqhi_nul{ll+40,tt,ss})); % keeping tac centred but jitter aud (to match tlock_tac{X,tt,ss})
              %       freqhi_audMSpN{ll,tt,ss}=ft_math(cfg,ft_selectdata(cfgavg,freqhi_aud{ll,tt,ss}),ft_selectdata(cfgavg,freqhi_nul{ll+40,tt,ss})); % keeping aud centred but jitter tac (to match tlock_aud{X,tt,ss})
              if numa_trials(ll,tt,ss)
                cfg=[];
                cfg.operation='add';
                cfg.parameter='powspctrm';
                freqlo_audMSpN_fftadd{ll,tt,ss}=ft_math(cfg,freqlo_aud{ll,tt,ss},freqlo_nul{ll+60,tt,ss}); % keeping aud centred but jitter tac (to match tlock_aud{X,tt,ss})
                freqhi_audMSpN_fftadd{ll,tt,ss}=ft_math(cfg,freqhi_aud{ll,tt,ss},freqhi_nul{ll+60,tt,ss}); % keeping aud centred but jitter tac (to match tlock_aud{X,tt,ss})
                freqlo_aMSAlone_fftadd{ll,tt,ss}=freqlo_aud{ll,tt,ss};
                freqhi_aMSAlone_fftadd{ll,tt,ss}=freqhi_aud{ll,tt,ss};
                freqlo_aNulAlone_fftadd{ll,tt,ss}=freqlo_nul{ll+60,tt,ss};
                freqhi_aNulAlone_fftadd{ll,tt,ss}=freqhi_nul{ll+60,tt,ss};
                cfg.parameter='crsspctrm';
                tmp=ft_math(cfg,freqlo_aud{ll,tt,ss},freqlo_nul{ll+60,tt,ss}); % keeping aud centred but jitter tac (to match tlock_aud{X,tt,ss})
                freqlo_audMSpN_fftadd{ll,tt,ss}.crsspctrm=tmp.crsspctrm;
                freqlo_aMSAlone_fftadd{ll,tt,ss}.crsspctrm=freqlo_aud{ll,tt,ss}.crsspctrm;
                freqlo_aNulAlone_fftadd{ll,tt,ss}.crsspctrm=freqlo_nul{ll+60,tt,ss}.crsspctrm;
                tmp=ft_math(cfg,freqhi_aud{ll,tt,ss},freqhi_nul{ll+60,tt,ss}); % keeping aud centred but jitter tac (to match tlock_aud{X,tt,ss})
                freqhi_audMSpN_fftadd{ll,tt,ss}.crsspctrm=tmp.crsspctrm;
                freqhi_aMSAlone_fftadd{ll,tt,ss}.crsspctrm=freqhi_aud{ll,tt,ss}.crsspctrm;
                freqhi_aNulAlone_fftadd{ll,tt,ss}.crsspctrm=freqhi_nul{ll+60,tt,ss}.crsspctrm;
              else
                freqlo_audMSpN_fftadd{ll,tt,ss}=[];
                freqhi_audMSpN_fftadd{ll,tt,ss}=[];
                freqlo_aMSAlone_fftadd{ll,tt,ss}=[];
                freqhi_aMSAlone_fftadd{ll,tt,ss}=[];
                freqlo_aNulAlone_fftadd{ll,tt,ss}=[];
                freqhi_aNulAlone_fftadd{ll,tt,ss}=[];
              end
            end % ll
            
            for ll=[1 3 4]
              if ~isempty(tlock_audMSshift1{ll,tt,ss})&& size(tlock_audMSshift1{ll,tt,ss}.trial,1)
                cfg=[];
                cfg.method='mtmconvol';
                cfg.pad=4;
                cfg.foi=4:2:30;
                cfg.taper='hanning';
                cfg.toi=-1.2:timestepfft:1.3;
                cfg.t_ftimwin=4./cfg.foi;
                cfg.output='powandcsd';
                freqlo_audMSshift1{ll,tt,ss}=ft_freqanalysis(cfg,tlock_audMSshift1{ll,tt,ss});
                freqlo_audMSshift2{ll,tt,ss}=ft_freqanalysis(cfg,tlock_audMSshift2{ll,tt,ss});
                freqlo_audMSshift3{ll,tt,ss}=ft_freqanalysis(cfg,tlock_audMSshift3{ll,tt,ss});
                freqlo_audMSshift4{ll,tt,ss}=ft_freqanalysis(cfg,tlock_audMSshift4{ll,tt,ss});
                
                cfg=[];
                cfg.method='mtmconvol';
                cfg.pad=4;
                cfg.foi=[30:5:45 55:5:80];
                cfg.taper='dpss';
                cfg.tapsmofrq=7*ones(1,length(cfg.foi));
                cfg.toi=-1.2:timestepfft:1.3;
                cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
                cfg.output='powandcsd';
                freqhi_audMSshift1{ll,tt,ss}=ft_freqanalysis(cfg,tlock_audMSshift1{ll,tt,ss});
                freqhi_audMSshift2{ll,tt,ss}=ft_freqanalysis(cfg,tlock_audMSshift2{ll,tt,ss});
                freqhi_audMSshift3{ll,tt,ss}=ft_freqanalysis(cfg,tlock_audMSshift3{ll,tt,ss});
                freqhi_audMSshift4{ll,tt,ss}=ft_freqanalysis(cfg,tlock_audMSshift4{ll,tt,ss});
              else
                freqlo_audMSshift1{ll,tt,ss}=[];
                freqlo_audMSshift2{ll,tt,ss}=[];
                freqlo_audMSshift3{ll,tt,ss}=[];
                freqlo_audMSshift4{ll,tt,ss}=[];
                freqhi_audMSshift1{ll,tt,ss}=[];
                freqhi_audMSshift2{ll,tt,ss}=[];
                freqhi_audMSshift3{ll,tt,ss}=[];
                freqhi_audMSshift4{ll,tt,ss}=[];
              end
              
              if ~isempty(tlock_audMSshift1tlock{ll,tt,ss})
                % create sum of asynchronous multisensory conditions (e.g. AT70_shifted plus TA70)
                cfg=[];
                cfg.operation='add';
                cfg.parameter='powspctrm';
                freqlo_aMSasynch_fftadd{ll,tt,ss}=ft_math(cfg,freqlo_audMSshift1{ll,tt,ss},freqlo_audMSshift2{ll,tt,ss});
                freqhi_aMSasynch_fftadd{ll,tt,ss}=ft_math(cfg,freqhi_audMSshift1{ll,tt,ss},freqhi_audMSshift2{ll,tt,ss});
                cfg.parameter='crsspctrm';
                tmp=ft_math(cfg,freqlo_audMSshift1{ll,tt,ss},freqlo_audMSshift2{ll,tt,ss});
                freqlo_aMSasynch_fftadd{ll,tt,ss}.crsspctrm=tmp.crsspctrm;
                tmp=ft_math(cfg,freqhi_audMSshift1{ll,tt,ss},freqhi_audMSshift2{ll,tt,ss});
                freqhi_aMSasynch_fftadd{ll,tt,ss}.crsspctrm=tmp.crsspctrm;
                
                % create sum of synchronous multisensory conditions (e.g. TA0_shifted plus TA0)
                cfg=[];
                cfg.operation='add';
                cfg.parameter='powspctrm';
                freqlo_aMSsynch_fftadd{ll,tt,ss}=ft_math(cfg,freqlo_audMSshift3{ll,tt,ss},freqlo_audMSshift4{ll,tt,ss});
                freqhi_aMSsynch_fftadd{ll,tt,ss}=ft_math(cfg,freqhi_audMSshift3{ll,tt,ss},freqhi_audMSshift4{ll,tt,ss});
                cfg.parameter='crsspctrm';
                tmp=ft_math(cfg,freqlo_audMSshift3{ll,tt,ss},freqlo_audMSshift4{ll,tt,ss});
                freqlo_aMSsynch_fftadd{ll,tt,ss}.crsspctrm=tmp.crsspctrm;
                tmp=ft_math(cfg,freqhi_audMSshift3{ll,tt,ss},freqhi_audMSshift4{ll,tt,ss});
                freqhi_aMSsynch_fftadd{ll,tt,ss}.crsspctrm=tmp.crsspctrm;
              else
                freqlo_aMSasynch_fftadd{ll,tt,ss}=[];
                freqhi_aMSasynch_fftadd{ll,tt,ss}=[];
                freqlo_aMSsynch_fftadd{ll,tt,ss}=[];
                freqhi_aMSsynch_fftadd{ll,tt,ss}=[];
              end
            end % ll
          end % dofftadd
          
          
          
          %%
          % Way 2: Add first then FFT
          
          % Done two ways:
          %
          % 2a) freqlo_audPtac_comb{ll,tt,ss,1} is the result from adding
          % the timelock unaveraged data first, then abs(FFT), with trial
          % pairings as they are originally ordered.
          %
          % 2b) And freqlo_audPtac_comb{ll,tt,ss,2} is also created from adding
          % the timelock unaveraged data first, then abs(FFT) but with
          % trial pairings to be added randomised and the full distribution
          % computed, followed at last by averaging over all the
          % realisations of this distribution.
          %
          
          for ll=soalist
            % individual conditions first
            cfg=[];
            cfg.method='mtmconvol';
            cfg.pad=4;
            cfg.foi=4:2:30;
            cfg.taper='hanning';
            cfg.toi=-1.2:timestepplv:1.3;
            cfg.t_ftimwin=4./cfg.foi;
            %         cfg.keeptrials='yes';  % not necessary at this stage for keeping trials (yes later for source analysis and phase resetting)
            cfg.output='fourier';
            freqlo_aTacAlone_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tac{ll+40,tt,ss});
            freqlo_aAudAlone_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_aud{ll+20,tt,ss});
            freqlo_aMSAlone_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_aud{ll,tt,ss});
            freqlo_aNulAlone_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_nul{ll+60,tt,ss});
            
            freqlo_aTacAlone_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqlo_aTacAlone_comb{ll,tt,ss}.fourierspctrm).^2,1));
            freqlo_aTacAlone_comb{ll,tt,ss}.plvspctrm=getplv(freqlo_aTacAlone_comb{ll,tt,ss}.fourierspctrm);
            freqlo_aTacAlone_comb{ll,tt,ss}=rmfield(freqlo_aTacAlone_comb{ll,tt,ss},'fourierspctrm');
            freqlo_aAudAlone_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqlo_aAudAlone_comb{ll,tt,ss}.fourierspctrm).^2,1));
            freqlo_aAudAlone_comb{ll,tt,ss}.plvspctrm=getplv(freqlo_aAudAlone_comb{ll,tt,ss}.fourierspctrm);
            freqlo_aAudAlone_comb{ll,tt,ss}=rmfield(freqlo_aAudAlone_comb{ll,tt,ss},'fourierspctrm');
            freqlo_aMSAlone_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqlo_aMSAlone_comb{ll,tt,ss}.fourierspctrm).^2,1));
            freqlo_aMSAlone_comb{ll,tt,ss}.plvspctrm=getplv(freqlo_aMSAlone_comb{ll,tt,ss}.fourierspctrm);
            freqlo_aMSAlone_comb{ll,tt,ss}=rmfield(freqlo_aMSAlone_comb{ll,tt,ss},'fourierspctrm');
            freqlo_aNulAlone_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqlo_aNulAlone_comb{ll,tt,ss}.fourierspctrm).^2,1));
            freqlo_aNulAlone_comb{ll,tt,ss}.plvspctrm=getplv(freqlo_aNulAlone_comb{ll,tt,ss}.fourierspctrm);
            freqlo_aNulAlone_comb{ll,tt,ss}=rmfield(freqlo_aNulAlone_comb{ll,tt,ss},'fourierspctrm');
            
            
            cfg=[];
            cfg.method='mtmconvol';
            cfg.pad=4;
            cfg.foi=[30:5:45 55:5:80];
            cfg.taper='dpss';
            cfg.tapsmofrq=7*ones(1,length(cfg.foi));
            cfg.toi=-1.2:timestepplv:1.3;
            cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
            cfg.output='fourier';
            freqhi_aTacAlone_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tac{ll+40,tt,ss});
            freqhi_aAudAlone_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_aud{ll+20,tt,ss});
            freqhi_aMSAlone_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_aud{ll,tt,ss});
            freqhi_aNulAlone_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_nul{ll+60,tt,ss});
            
            freqhi_aTacAlone_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqhi_aTacAlone_comb{ll,tt,ss}.fourierspctrm).^2,1));
            freqhi_aTacAlone_comb{ll,tt,ss}.plvspctrm=getplv(freqhi_aTacAlone_comb{ll,tt,ss}.fourierspctrm);
            freqhi_aTacAlone_comb{ll,tt,ss}=rmfield(freqhi_aTacAlone_comb{ll,tt,ss},'fourierspctrm');
            freqhi_aAudAlone_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqhi_aAudAlone_comb{ll,tt,ss}.fourierspctrm).^2,1));
            freqhi_aAudAlone_comb{ll,tt,ss}.plvspctrm=getplv(freqhi_aAudAlone_comb{ll,tt,ss}.fourierspctrm);
            freqhi_aAudAlone_comb{ll,tt,ss}=rmfield(freqhi_aAudAlone_comb{ll,tt,ss},'fourierspctrm');
            freqhi_aMSAlone_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqhi_aMSAlone_comb{ll,tt,ss}.fourierspctrm).^2,1));
            freqhi_aMSAlone_comb{ll,tt,ss}.plvspctrm=getplv(freqhi_aMSAlone_comb{ll,tt,ss}.fourierspctrm);
            freqhi_aMSAlone_comb{ll,tt,ss}=rmfield(freqhi_aMSAlone_comb{ll,tt,ss},'fourierspctrm');
            freqhi_aNulAlone_comb{ll,tt,ss}.powspctrm=squeeze(mean(abs(freqhi_aNulAlone_comb{ll,tt,ss}.fourierspctrm).^2,1));
            freqhi_aNulAlone_comb{ll,tt,ss}.plvspctrm=getplv(freqhi_aNulAlone_comb{ll,tt,ss}.fourierspctrm);
            freqhi_aNulAlone_comb{ll,tt,ss}=rmfield(freqhi_aNulAlone_comb{ll,tt,ss},'fourierspctrm');
            
            %           if ll<10
            % audPtac
            numcomb=numa_trials(ll,tt,ss)^2;
            numtests=min(numcomb,minnumcomb);
            %               combindex=reshape(1:numcomb,numa_trials(ll,tt,ss),numa_trials(ll,tt,ss));
            freqlo_audPtac_tmp{1}=[];
            freqhi_audPtac_tmp{1}=[];
            for cc=1:numtests
              %                 tmp=Shuffle(combindex(:));
              %                 combuse=tmp(1:numa_trials(ll,tt,ss));
              combuse=Shuffle(1:numa_trials(ll,tt,ss));
              
              tlock_fake=tlock_aud{ll+20,tt,ss};
              for at=1:numa_trials(ll,tt,ss)
                if cc==1  % do it as 'normal'
                  tind=at;aind=at;
                else
                  %                     [aind,tind]=find(combindex==combuse(at));
                  aind=at;
                  tind=combuse(at);
                end
                cfg=[];
                cfg.trials=tind;
                tmpt=ft_selectdata(cfg,tlock_tac{ll+40,tt,ss});
                cfg.trials=aind;
                tmpa=ft_selectdata(cfg,tlock_aud{ll+20,tt,ss});
                cfg=[];
                cfg.operation='add';
                cfg.parameter='trial';
                tmpsum=ft_math(cfg,tmpt,tmpa);
                if at==1
                  tlock_fake=tmpsum;
                end
                tlock_fake.trial(at,:,:)=tmpsum.trial(1,:,:);
                %           tlock_fake.trial(at,:,:)=tlock_aud{ll+20,tt,ss}.trial(aind,:,:)+tlock_tac{ll+40,tt,ss}.trial(tind,:,:);
              end
              
              cfg=[];
              cfg.method='mtmconvol';
              cfg.pad=4;
              cfg.foi=4:2:30;
              cfg.taper='hanning';
              cfg.toi=-1.2:timestepplv:1.3;
              cfg.t_ftimwin=4./cfg.foi;
              %         cfg.keeptrials='yes';  % not necessary at this stage for keeping trials (yes later for source analysis and phase resetting)
              cfg.output='fourier';
              freqlo_audPtac_tmp{cc}=ft_freqanalysis(cfg,tlock_fake);
              freqlo_audPtac_tmp{cc}.powspctrm=squeeze(mean(abs(freqlo_audPtac_tmp{cc}.fourierspctrm).^2,1));
              freqlo_audPtac_tmp{cc}.plvspctrm=getplv(freqlo_audPtac_tmp{cc}.fourierspctrm);
              freqlo_audPtac_tmp{cc}=rmfield(freqlo_audPtac_tmp{cc},'fourierspctrm');
              cfg=[];
              cfg.method='mtmconvol';
              cfg.pad=4;
              cfg.foi=[30:5:45 55:5:80];
              cfg.taper='dpss';
              cfg.tapsmofrq=7*ones(1,length(cfg.foi));
              cfg.toi=-1.2:timestepplv:1.3;
              cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
              cfg.output='fourier';
              freqhi_audPtac_tmp{cc}=ft_freqanalysis(cfg,tlock_fake);
              freqhi_audPtac_tmp{cc}.powspctrm=squeeze(mean(abs(freqhi_audPtac_tmp{cc}.fourierspctrm).^2,1));
              freqhi_audPtac_tmp{cc}.plvspctrm=getplv(freqhi_audPtac_tmp{cc}.fourierspctrm);
              freqhi_audPtac_tmp{cc}=rmfield(freqhi_audPtac_tmp{cc},'fourierspctrm');
            end
            freqlo_audPtac_comb{ll,tt,ss,1}=freqlo_audPtac_tmp{1};
            clear tmp tmp2 tmp3
            if numtests>1
              for cc=1:numtests,tmp(:,:,:,cc)=freqlo_audPtac_tmp{cc}.powspctrm;end
              %            for cc=1:numtests,tmp2(:,:,:,cc)=freqlo_audPtac_tmp{cc}.crsspctrm;end
              for cc=1:numtests,tmp3(:,:,:,cc)=freqlo_audPtac_tmp{cc}.plvspctrm;end
              if exist('tmp','var')
                freqlo_audPtac_comb{ll,tt,ss,2}=freqlo_audPtac_tmp{1};
                freqlo_audPtac_comb{ll,tt,ss,2}.powspctrm=mean(tmp,4);
                freqlo_audPtac_comb{ll,tt,ss,2}.stdpow=std(tmp,[],4);
                %              freqlo_audPtac_comb{ll,tt,ss,2}.crsspctrm=mean(tmp2,4);
                %              freqlo_audPtac_comb{ll,tt,ss,2}.stdcrs=std(tmp2,[],4);
                freqlo_audPtac_comb{ll,tt,ss,2}.plvspctrm=mean(tmp3,4);
                freqlo_audPtac_comb{ll,tt,ss,2}.stdplv=std(tmp3,[],4);
                if numtests>=10
                  disp('must still code this part for saving means of lots combinations')
                  keyboard
                end
              end
            else
              freqlo_audPtac_comb{ll,tt,ss,2}=[];
            end
            tmplo_audPtac_comb2=tmp;
            tmplo_audPtac_comb2plv=tmp3;
            freqhi_audPtac_comb{ll,tt,ss,1}=freqhi_audPtac_tmp{1};
            clear tmp tmp2 tmp3
            if numtests>1
              for cc=1:numtests,tmp(:,:,:,cc)=freqhi_audPtac_tmp{cc}.powspctrm;end
              %            for cc=1:numtests,tmp2(:,:,:,cc)=freqhi_audPtac_tmp{cc}.crsspctrm;end
              for cc=1:numtests,tmp3(:,:,:,cc)=freqhi_audPtac_tmp{cc}.plvspctrm;end
              if exist('tmp','var')
                freqhi_audPtac_comb{ll,tt,ss,2}=freqhi_audPtac_tmp{1};
                freqhi_audPtac_comb{ll,tt,ss,2}.powspctrm=mean(tmp,4);
                freqhi_audPtac_comb{ll,tt,ss,2}.stdpow=std(tmp,[],4);
                %              freqhi_audPtac_comb{ll,tt,ss,2}.crsspctrm=mean(tmp2,4);
                %              freqhi_audPtac_comb{ll,tt,ss,2}.stdcrs=std(tmp2,[],4);
                freqhi_audPtac_comb{ll,tt,ss,2}.plvspctrm=mean(tmp3,4);
                freqhi_audPtac_comb{ll,tt,ss,2}.stdplv=std(tmp3,[],4);
              end
            else
              freqhi_audPtac_comb{ll,tt,ss,2}=[];
            end
            tmphi_audPtac_comb2=tmp;
            tmphi_audPtac_comb2plv=tmp3;
            clear freqlo_audPtac_tmp freqhi_audPtac_tmp
            clear tmp tmp2 tmp3
            if ss==12 || ss==13
            else
              tlock_tac{ll+40,tt,ss}=[]; % clearing to help save memory as running
              tlock_aud{ll+20,tt,ss}=[]; % clearing to help save memory as running
            end
            
            
            
            % audMSpN
            numcomb=numa_trials(ll,tt,ss)^2;
            numtests=min(numcomb,minnumcomb);
            %               combindex=reshape(1:numcomb,numa_trials(ll,tt,ss),numa_trials(ll,tt,ss));
            freqlo_audMSpN_tmp{1}=[];
            freqhi_audMSpN_tmp{1}=[];
            for cc=1:numtests
              %                 tmp=Shuffle(combindex(:));
              %                 combuse=tmp(1:numa_trials(ll,tt,ss));
              combuse=Shuffle(1:numa_trials(ll,tt,ss));
              
              tlock_fake=tlock_aud{ll,tt,ss};
              for at=1:numa_trials(ll,tt,ss)
                if cc==1  % do it as 'normal'
                  tind=at;aind=at;
                else
                  %                     [aind,tind]=find(combindex==combuse(at));
                  %                     aind=at;
                  tind=combuse(at);
                end
                cfg=[];
                cfg.trials=tind;
                tmpt=ft_selectdata(cfg,tlock_nul{ll+60,tt,ss});
                cfg.trials=aind;
                tmpa=ft_selectdata(cfg,tlock_aud{ll,tt,ss});
                cfg=[];
                cfg.operation='add';
                cfg.parameter='trial';
                tmpsum=ft_math(cfg,tmpt,tmpa);
                if at==1
                  tlock_fake=tmpsum;
                end
                tlock_fake.trial(at,:,:)=tmpsum.trial(1,:,:);
                %           tlock_fake.trial(at,:,:)=tlock_aud{ll,tt,ss}.trial(aind,:,:)+tlock_nul{ll+60,tt,ss}.trial(tind,:,:);
              end
              cfg=[];
              cfg.method='mtmconvol';
              cfg.pad=4;
              cfg.foi=4:2:30;
              cfg.taper='hanning';
              cfg.toi=-1.2:timestepplv:1.3;
              cfg.t_ftimwin=4./cfg.foi;
              %         cfg.keeptrials='yes';  % not necessary at this stage for keeping trials (yes later for source analysis and phase resetting)
              cfg.output='fourier';
              freqlo_audMSpN_tmp{cc}=ft_freqanalysis(cfg,tlock_fake);
              freqlo_audMSpN_tmp{cc}.powspctrm=squeeze(mean(abs(freqlo_audMSpN_tmp{cc}.fourierspctrm).^2,1));
              freqlo_audMSpN_tmp{cc}.plvspctrm=getplv(freqlo_audMSpN_tmp{cc}.fourierspctrm);
              freqlo_audMSpN_tmp{cc}=rmfield(freqlo_audMSpN_tmp{cc},'fourierspctrm');
              cfg=[];
              cfg.method='mtmconvol';
              cfg.pad=4;
              cfg.foi=[30:5:45 55:5:80];
              cfg.taper='dpss';
              cfg.tapsmofrq=7*ones(1,length(cfg.foi));
              cfg.toi=-1.2:timestepplv:1.3;
              cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
              cfg.output='fourier';
              freqhi_audMSpN_tmp{cc}=ft_freqanalysis(cfg,tlock_fake);
              freqhi_audMSpN_tmp{cc}.powspctrm=squeeze(mean(abs(freqhi_audMSpN_tmp{cc}.fourierspctrm).^2,1));
              freqhi_audMSpN_tmp{cc}.plvspctrm=getplv(freqhi_audMSpN_tmp{cc}.fourierspctrm);
              freqhi_audMSpN_tmp{cc}=rmfield(freqhi_audMSpN_tmp{cc},'fourierspctrm');
            end
            tmplo_audMSpN_comb2=tmp;
            freqlo_audMSpN_comb{ll,tt,ss,1}=freqlo_audMSpN_tmp{1};
            clear tmp tmp2 tmp3
            if numtests>1
              for cc=1:numtests,tmp(:,:,:,cc)=freqlo_audMSpN_tmp{cc}.powspctrm;end
              %            for cc=1:numtests,tmp2(:,:,:,cc)=freqlo_audMSpN_tmp{cc}.crsspctrm;end
              for cc=1:numtests,tmp3(:,:,:,cc)=freqlo_audMSpN_tmp{cc}.plvspctrm;end
              if exist('tmp','var')
                freqlo_audMSpN_comb{ll,tt,ss,2}=freqlo_audMSpN_tmp{1};
                freqlo_audMSpN_comb{ll,tt,ss,2}.powspctrm=mean(tmp,4);
                freqlo_audMSpN_comb{ll,tt,ss,2}.stdpow=std(tmp,[],4);
                %              freqlo_audMSpN_comb{ll,tt,ss,2}.crsspctrm=mean(tmp2,4);
                %              freqlo_audMSpN_comb{ll,tt,ss,2}.stdcrs=std(tmp2,[],4);
                freqlo_audMSpN_comb{ll,tt,ss,2}.plvspctrm=mean(tmp3,4);
                freqlo_audMSpN_comb{ll,tt,ss,2}.stdplv=std(tmp3,[],4);
              end
            else
              freqlo_audMSpN_comb{ll,tt,ss,2}=[];
            end
            tmphi_audMSpN_comb2=tmp;
            freqhi_audMSpN_comb{ll,tt,ss,1}=freqhi_audMSpN_tmp{1};
            clear tmp tmp2 tmp3
            if numtests>1
              for cc=1:numtests,tmp(:,:,:,cc)=freqhi_audMSpN_tmp{cc}.powspctrm;end
              %            for cc=1:numtests,tmp2(:,:,:,cc)=freqhi_audMSpN_tmp{cc}.crsspctrm;end
              for cc=1:numtests,tmp3(:,:,:,cc)=freqhi_audMSpN_tmp{cc}.plvspctrm;end
              if exist('tmp','var')
                freqhi_audMSpN_comb{ll,tt,ss,2}=freqhi_audMSpN_tmp{1};
                freqhi_audMSpN_comb{ll,tt,ss,2}.powspctrm=mean(tmp,4);
                freqhi_audMSpN_comb{ll,tt,ss,2}.stdpow=std(tmp,[],4);
                %              freqhi_audMSpN_comb{ll,tt,ss,2}.crsspctrm=mean(tmp2,4);
                %              freqhi_audMSpN_comb{ll,tt,ss,2}.stdcrs=std(tmp2,[],4);
                freqhi_audMSpN_comb{ll,tt,ss,2}.plvspctrm=mean(tmp3,4);
                freqhi_audMSpN_comb{ll,tt,ss,2}.stdplv=std(tmp3,[],4);
              end
            else
              freqhi_audMSpN_comb{ll,tt,ss,2}=[];
            end
            clear freqlo_audMSpN_tmp freqhi_audMSpN_tmp
            clear tmp tmp2 tmp3
            
            if ss==12 || ss==13
            else
              tlock_nul{ll+60,tt,ss}=[]; % clearing to help save memory as running
              tlock_aud{ll,tt,ss}=[]; % clearing to help save memory as running
            end
            
            
            for cc=1:size(tmplo_audPtac_comb2,1)
              for ff=1:size(tmplo_audPtac_comb2,2)
                for ttime=1:size(tmplo_audPtac_comb2,3)
                  [statsa{ll,tt,ss}.ttestlo_h(cc,ff,ttime),statsa{ll,tt,ss}.ttestlo_p(cc,ff,ttime),statsa{ll,tt,ss}.ttestlo_ci(cc,ff,ttime,1,:),statsa{ll,tt,ss}.ttestlo_stats{cc,ff,ttime,1}]=ttest2(squeeze(tmplo_audPtac_comb2(cc,ff,ttime,:)),squeeze(tmplo_audMSpN_comb2(cc,ff,ttime,:)));
                  [statsa{ll,tt,ss}.ttesthi_h(cc,ff,ttime),statsa{ll,tt,ss}.ttesthi_p(cc,ff,ttime),statsa{ll,tt,ss}.ttesthi_ci(cc,ff,ttime,1,:),statsa{ll,tt,ss}.ttesthi_stats{cc,ff,ttime,1}]=ttest2(squeeze(tmphi_audPtac_comb2(cc,ff,ttime,:)),squeeze(tmphi_audMSpN_comb2(cc,ff,ttime,:)));
                  [statsa{ll,tt,ss}.kstestlo_h(cc,ff,ttime),statsa{ll,tt,ss}.kstestlo_p(cc,ff,ttime),statsa{ll,tt,ss}.kstestlo_stats{cc,ff,ttime,1}]=kstest2(squeeze(tmplo_audPtac_comb2(cc,ff,ttime,:)),squeeze(tmplo_audMSpN_comb2(cc,ff,ttime,:)));
                  [statsa{ll,tt,ss}.kstesthi_h(cc,ff,ttime),statsa{ll,tt,ss}.kstesthi_p(cc,ff,ttime),statsa{ll,tt,ss}.kstesthi_stats{cc,ff,ttime,1}]=kstest2(squeeze(tmphi_audPtac_comb2(cc,ff,ttime,:)),squeeze(tmphi_audMSpN_comb2(cc,ff,ttime,:)));
                  if numtests<=10
                    [statsa{ll,tt,ss}.ttestlo_h(cc,ff,ttime,2),statsa{ll,tt,ss}.ttestlo_p(cc,ff,ttime,2),statsa{ll,tt,ss}.ttestlo_ci(cc,ff,ttime,2,:),statsa{ll,tt,ss}.ttestlo_stats{cc,ff,ttime,2}]=ttest2(squeeze(tmplo_audPtac_comb2(cc,ff,ttime,1:halftest)),squeeze(tmplo_audMSpN_comb2(cc,ff,ttime,1:halftest)));
                    [statsa{ll,tt,ss}.ttesthi_h(cc,ff,ttime,2),statsa{ll,tt,ss}.ttesthi_p(cc,ff,ttime,2),statsa{ll,tt,ss}.ttesthi_ci(cc,ff,ttime,2,:),statsa{ll,tt,ss}.ttesthi_stats{cc,ff,ttime,2}]=ttest2(squeeze(tmphi_audPtac_comb2(cc,ff,ttime,1:halftest)),squeeze(tmphi_audMSpN_comb2(cc,ff,ttime,1:halftest)));
                    [statsa{ll,tt,ss}.kstestlo_h(cc,ff,ttime,2),statsa{ll,tt,ss}.kstestlo_p(cc,ff,ttime,2),statsa{ll,tt,ss}.kstestlo_stats{cc,ff,ttime,2}]=kstest2(squeeze(tmplo_audPtac_comb2(cc,ff,ttime,1:halftest)),squeeze(tmplo_audMSpN_comb2(cc,ff,ttime,1:halftest)));
                    [statsa{ll,tt,ss}.kstesthi_h(cc,ff,ttime,2),statsa{ll,tt,ss}.kstesthi_p(cc,ff,ttime,2),statsa{ll,tt,ss}.kstesthi_stats{cc,ff,ttime,2}]=kstest2(squeeze(tmphi_audPtac_comb2(cc,ff,ttime,1:halftest)),squeeze(tmphi_audMSpN_comb2(cc,ff,ttime,1:halftest)));
                    [statsa{ll,tt,ss}.ttestlo_h(cc,ff,ttime,3),statsa{ll,tt,ss}.ttestlo_p(cc,ff,ttime,3),statsa{ll,tt,ss}.ttestlo_ci(cc,ff,ttime,3,:),statsa{ll,tt,ss}.ttestlo_stats{cc,ff,ttime,3}]=ttest2(squeeze(tmplo_audPtac_comb2(cc,ff,ttime,halftest+1:end)),squeeze(tmplo_audMSpN_comb2(cc,ff,ttime,halftest+1:end)));
                    [statsa{ll,tt,ss}.ttesthi_h(cc,ff,ttime,3),statsa{ll,tt,ss}.ttesthi_p(cc,ff,ttime,3),statsa{ll,tt,ss}.ttesthi_ci(cc,ff,ttime,3,:),statsa{ll,tt,ss}.ttesthi_stats{cc,ff,ttime,3}]=ttest2(squeeze(tmphi_audPtac_comb2(cc,ff,ttime,halftest+1:end)),squeeze(tmphi_audMSpN_comb2(cc,ff,ttime,halftest+1:end)));
                    [statsa{ll,tt,ss}.kstestlo_h(cc,ff,ttime,3),statsa{ll,tt,ss}.kstestlo_p(cc,ff,ttime,3),statsa{ll,tt,ss}.kstestlo_stats{cc,ff,ttime,3}]=kstest2(squeeze(tmplo_audPtac_comb2(cc,ff,ttime,halftest+1:end)),squeeze(tmplo_audMSpN_comb2(cc,ff,ttime,halftest+1:end)));
                    [statsa{ll,tt,ss}.kstesthi_h(cc,ff,ttime,3),statsa{ll,tt,ss}.kstesthi_p(cc,ff,ttime,3),statsa{ll,tt,ss}.kstesthi_stats{cc,ff,ttime,3}]=kstest2(squeeze(tmphi_audPtac_comb2(cc,ff,ttime,halftest+1:end)),squeeze(tmphi_audMSpN_comb2(cc,ff,ttime,halftest+1:end)));
                  end
                end
              end
            end
            clear tmplo_* tmphi_*
            
            %           end
            
            %           for ll=[soalist soalist+20 soalist+40]
            %           tlock_tac{ll,tt,ss}=[]; % clearing to help save memory as running
            %           tlock_aud{ll,tt,ss}=[]; % clearing to help save memory as running
            %           tlock_nul{ll,tt,ss}=[]; % clearing to help save memory as running
            
            
            
            if ll<5
              numasynch_trials(ll,tt,ss)=size(tlock_audMSshift1{ll,tt,ss}.trialinfo,1);
              % (AT70 + TA70)
              numcomb=numasynch_trials(ll,tt,ss)^2;
              numtests=min(numcomb,minnumcomb);
              %             combindex=reshape(1:numcomb,numt_trials(ll,tt,ss),numt_trials(ll,tt,ss));
              freqlo_aMSasynch_tmp{1}=[];
              freqhi_aMSasynch_tmp{1}=[];
              for cc=1:numtests
                %               tmp=Shuffle(combindex(:));
                %               combuse=tmp(1:numt_trials(ll,tt,ss));
                combuse=Shuffle(1:numasynch_trials(ll,tt,ss));
                
                tlock_fake=tlock_audMSshift1{ll,tt,ss};
                for at=1:numasynch_trials(ll,tt,ss)
                  if cc==1  % do it as 'normal'
                    tind=at;aind=at;
                  else
                    %                   [tind,aind]=find(combindex==combuse(at));
                    tind=at;
                    aind=combuse(at);
                  end
                  cfg=[];
                  cfg.trials=tind;
                  tmpt=ft_selectdata(cfg,tlock_audMSshift1{ll,tt,ss});
                  cfg.trials=aind;
                  tmpa=ft_selectdata(cfg,tlock_audMSshift2{ll,tt,ss});
                  cfg=[];
                  cfg.operation='add';
                  cfg.parameter='trial';
                  tmpsum=ft_math(cfg,tmpt,tmpa);
                  if at==1
                    tlock_fake=tmpsum;
                  end
                  tlock_fake.trial(at,:,:)=tmpsum.trial(1,:,:);
                end
                cfg=[];
                cfg.method='mtmconvol';
                cfg.pad=4;
                cfg.foi=4:2:30;
                cfg.taper='hanning';
                cfg.toi=-1.2:timestepplv:1.3;
                cfg.t_ftimwin=4./cfg.foi;
                cfg.output='fourier';
                freqlo_aMSasynch_tmp{cc}=ft_freqanalysis(cfg,tlock_fake);
                freqlo_aMSasynch_tmp{cc}.powspctrm=squeeze(mean(abs(freqlo_aMSasynch_tmp{cc}.fourierspctrm).^2,1));
                freqlo_aMSasynch_tmp{cc}.plvspctrm=getplv(freqlo_aMSasynch_tmp{cc}.fourierspctrm);
                freqlo_aMSasynch_tmp{cc}=rmfield(freqlo_aMSasynch_tmp{cc},'fourierspctrm');
                cfg=[];
                cfg.method='mtmconvol';
                cfg.pad=4;
                cfg.foi=[30:5:45 55:5:80];
                cfg.taper='dpss';
                cfg.tapsmofrq=7*ones(1,length(cfg.foi));
                cfg.toi=-1.2:timestepplv:1.3;
                cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
                cfg.output='fourier';
                freqhi_aMSasynch_tmp{cc}=ft_freqanalysis(cfg,tlock_fake);
                freqhi_aMSasynch_tmp{cc}.powspctrm=squeeze(mean(abs(freqhi_aMSasynch_tmp{cc}.fourierspctrm).^2,1));
                freqhi_aMSasynch_tmp{cc}.plvspctrm=getplv(freqhi_aMSasynch_tmp{cc}.fourierspctrm);
                freqhi_aMSasynch_tmp{cc}=rmfield(freqhi_aMSasynch_tmp{cc},'fourierspctrm');
              end
              freqlo_aMSasynch_comb{ll,tt,ss,1}=freqlo_aMSasynch_tmp{1};
              clear tmp tmp2 tmp3
              if numtests>1
                for cc=1:numtests,tmp(:,:,:,cc)=freqlo_aMSasynch_tmp{cc}.powspctrm;end
                for cc=1:numtests,tmp3(:,:,:,cc)=freqlo_aMSasynch_tmp{cc}.plvspctrm;end
                if exist('tmp','var')
                  freqlo_aMSasynch_comb{ll,tt,ss,2}=freqlo_aMSasynch_tmp{1};
                  freqlo_aMSasynch_comb{ll,tt,ss,2}.powspctrm=mean(tmp,4);
                  freqlo_aMSasynch_comb{ll,tt,ss,2}.stdpow=std(tmp,[],4);
                  freqlo_aMSasynch_comb{ll,tt,ss,2}.crsspctrm=mean(tmp3,4);
                  freqlo_aMSasynch_comb{ll,tt,ss,2}.stdcrs=std(tmp3,[],4);
                end
              else
                freqlo_aMSasynch_comb{ll,tt,ss,2}=[];
              end
              freqhi_aMSasynch_comb{ll,tt,ss,1}=freqhi_aMSasynch_tmp{1};
              clear tmp tmp2 tmp3
              if numtests>1
                for cc=1:numtests,tmp(:,:,:,cc)=freqhi_aMSasynch_tmp{cc}.powspctrm;end
                for cc=1:numtests,tmp3(:,:,:,cc)=freqhi_aMSasynch_tmp{cc}.plvspctrm;end
                if exist('tmp','var')
                  freqhi_aMSasynch_comb{ll,tt,ss,2}=freqhi_aMSasynch_tmp{1};
                  freqhi_aMSasynch_comb{ll,tt,ss,2}.powspctrm=mean(tmp,4);
                  freqhi_aMSasynch_comb{ll,tt,ss,2}.stdpow=std(tmp,[],4);
                  freqhi_aMSasynch_comb{ll,tt,ss,2}.crsspctrm=mean(tmp3,4);
                  freqhi_aMSasynch_comb{ll,tt,ss,2}.stdcrs=std(tmp3,[],4);
                end
              else
                freqhi_aMSasynch_comb{ll,tt,ss,2}=[];
              end
              clear freqlo_aMSasynch_tmp freqhi_aMSasynch_tmp
              clear tmp tmp2 tmp3
              
              if ss==12 || ss==13
              else
                tlock_audMSshift1{ll,tt,ss}=[]; % clearing to help save memory as running
                tlock_audMSshift2{ll,tt,ss}=[]; % clearing to help save memory as running
              end
              
              % (TA0 + TA0_shifted70)
              freqlo_aMSsynch_tmp{1}=[];
              freqhi_aMSsynch_tmp{1}=[];
              for cc=1:numtests
                combuse=Shuffle(1:numasynch_trials(ll,tt,ss));
                
                tlock_fake=tlock_audMSshift3{ll,tt,ss};
                for at=1:numasynch_trials(ll,tt,ss)
                  if cc==1  % do it as 'normal'
                    tind=at;aind=at;
                  else
                    tind=at;
                    aind=combuse(at);
                  end
                  cfg=[];
                  cfg.trials=tind;
                  tmpt=ft_selectdata(cfg,tlock_audMSshift3{ll,tt,ss});
                  cfg.trials=aind;
                  tmpa=ft_selectdata(cfg,tlock_audMSshift4{ll,tt,ss});
                  cfg=[];
                  cfg.operation='add';
                  cfg.parameter='trial';
                  tmpsum=ft_math(cfg,tmpt,tmpa);
                  if at==1
                    tlock_fake=tmpsum;
                  end
                  tlock_fake.trial(at,:,:)=tmpsum.trial(1,:,:);
                end
                cfg=[];
                cfg.method='mtmconvol';
                cfg.pad=4;
                cfg.foi=4:2:30;
                cfg.taper='hanning';
                cfg.toi=-1.2:timestepplv:1.3;
                cfg.t_ftimwin=4./cfg.foi;
                cfg.output='fourier';
                freqlo_aMSsynch_tmp{cc}=ft_freqanalysis(cfg,tlock_fake);
                freqlo_aMSsynch_tmp{cc}.powspctrm=squeeze(mean(abs(freqlo_aMSsynch_tmp{cc}.fourierspctrm).^2,1));
                freqlo_aMSsynch_tmp{cc}.plvspctrm=getplv(freqlo_aMSsynch_tmp{cc}.fourierspctrm);
                freqlo_aMSsynch_tmp{cc}=rmfield(freqlo_aMSsynch_tmp{cc},'fourierspctrm');
                cfg=[];
                cfg.method='mtmconvol';
                cfg.pad=4;
                cfg.foi=[30:5:45 55:5:80];
                cfg.taper='dpss';
                cfg.tapsmofrq=7*ones(1,length(cfg.foi));
                cfg.toi=-1.2:timestepplv:1.3;
                cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
                cfg.output='fourier';
                freqhi_aMSsynch_tmp{cc}=ft_freqanalysis(cfg,tlock_fake);
                freqhi_aMSsynch_tmp{cc}.powspctrm=squeeze(mean(abs(freqhi_aMSsynch_tmp{cc}.fourierspctrm).^2,1));
                freqhi_aMSsynch_tmp{cc}.plvspctrm=getplv(freqhi_aMSsynch_tmp{cc}.fourierspctrm);
                freqhi_aMSsynch_tmp{cc}=rmfield(freqhi_aMSsynch_tmp{cc},'fourierspctrm');
              end
              freqlo_aMSsynch_comb{ll,tt,ss,1}=freqlo_aMSsynch_tmp{1};
              clear tmp tmp2 tmp3
              if numtests>1
                for cc=1:numtests,tmp(:,:,:,cc)=freqlo_aMSsynch_tmp{cc}.powspctrm;end
                for cc=1:numtests,tmp3(:,:,:,cc)=freqlo_aMSsynch_tmp{cc}.plvspctrm;end
                if exist('tmp','var')
                  freqlo_aMSsynch_comb{ll,tt,ss,2}=freqlo_aMSsynch_tmp{1};
                  freqlo_aMSsynch_comb{ll,tt,ss,2}.powspctrm=mean(tmp,4);
                  freqlo_aMSsynch_comb{ll,tt,ss,2}.stdpow=std(tmp,[],4);
                  freqlo_aMSsynch_comb{ll,tt,ss,2}.crsspctrm=mean(tmp3,4);
                  freqlo_aMSsynch_comb{ll,tt,ss,2}.stdcrs=std(tmp3,[],4);
                end
              else
                freqlo_aMSsynch_comb{ll,tt,ss,2}=[];
              end
              freqhi_aMSsynch_comb{ll,tt,ss,1}=freqhi_aMSsynch_tmp{1};
              clear tmp tmp2 tmp3
              if numtests>1
                for cc=1:numtests,tmp(:,:,:,cc)=freqhi_aMSsynch_tmp{cc}.powspctrm;end
                for cc=1:numtests,tmp3(:,:,:,cc)=freqhi_aMSsynch_tmp{cc}.plvspctrm;end
                if exist('tmp','var')
                  freqhi_aMSsynch_comb{ll,tt,ss,2}=freqhi_aMSsynch_tmp{1};
                  freqhi_aMSsynch_comb{ll,tt,ss,2}.powspctrm=mean(tmp,4);
                  freqhi_aMSsynch_comb{ll,tt,ss,2}.stdpow=std(tmp,[],4);
                  freqhi_aMSsynch_comb{ll,tt,ss,2}.crsspctrm=mean(tmp3,4);
                  freqhi_aMSsynch_comb{ll,tt,ss,2}.stdcrs=std(tmp3,[],4);
                end
              else
                freqhi_aMSsynch_comb{ll,tt,ss,2}=[];
              end
              clear freqlo_aMSsynch_tmp freqhi_aMSsynch_tmp
              clear tmp tmp2 tmp3
              
              if ss==12 || ss==13
              else
                tlock_audMSshift3{ll,tt,ss}=[]; % clearing to help save memory as running
                tlock_audMSshift4{ll,tt,ss}=[]; % clearing to help save memory as running
              end
              
              for cc=1:size(tmplo_audPtac_comb2,1)
                for ff=1:size(tmplo_audPtac_comb2,2)
                  for ttime=1:size(tmplo_audPtac_comb2,3)
                    [statsa{ll,tt,ss}.ttestlo_h(cc,ff,ttime),statsa{ll,tt,ss}.ttestlo_p(cc,ff,ttime),statsa{ll,tt,ss}.ttestlo_ci(cc,ff,ttime,1,:),statsa{ll,tt,ss}.ttestlo_stats{cc,ff,ttime,1}]=ttest2(squeeze(tmplo_audPtac_comb2(cc,ff,ttime,:)),squeeze(tmplo_audMSpN_comb2(cc,ff,ttime,:)));
                    [statsa{ll,tt,ss}.ttesthi_h(cc,ff,ttime),statsa{ll,tt,ss}.ttesthi_p(cc,ff,ttime),statsa{ll,tt,ss}.ttesthi_ci(cc,ff,ttime,1,:),statsa{ll,tt,ss}.ttesthi_stats{cc,ff,ttime,1}]=ttest2(squeeze(tmphi_audPtac_comb2(cc,ff,ttime,:)),squeeze(tmphi_audMSpN_comb2(cc,ff,ttime,:)));
                    [statsa{ll,tt,ss}.kstestlo_h(cc,ff,ttime),statsa{ll,tt,ss}.kstestlo_p(cc,ff,ttime),statsa{ll,tt,ss}.kstestlo_stats{cc,ff,ttime,1}]=kstest2(squeeze(tmplo_audPtac_comb2(cc,ff,ttime,:)),squeeze(tmplo_audMSpN_comb2(cc,ff,ttime,:)));
                    [statsa{ll,tt,ss}.kstesthi_h(cc,ff,ttime),statsa{ll,tt,ss}.kstesthi_p(cc,ff,ttime),statsa{ll,tt,ss}.kstesthi_stats{cc,ff,ttime,1}]=kstest2(squeeze(tmphi_audPtac_comb2(cc,ff,ttime,:)),squeeze(tmphi_audMSpN_comb2(cc,ff,ttime,:)));
                    if numtests<=10
                      [statsa{ll,tt,ss}.ttestlo_h(cc,ff,ttime,2),statsa{ll,tt,ss}.ttestlo_p(cc,ff,ttime,2),statsa{ll,tt,ss}.ttestlo_ci(cc,ff,ttime,2,:),statsa{ll,tt,ss}.ttestlo_stats{cc,ff,ttime,2}]=ttest2(squeeze(tmplo_audPtac_comb2(cc,ff,ttime,1:halftest)),squeeze(tmplo_audMSpN_comb2(cc,ff,ttime,1:halftest)));
                      [statsa{ll,tt,ss}.ttesthi_h(cc,ff,ttime,2),statsa{ll,tt,ss}.ttesthi_p(cc,ff,ttime,2),statsa{ll,tt,ss}.ttesthi_ci(cc,ff,ttime,2,:),statsa{ll,tt,ss}.ttesthi_stats{cc,ff,ttime,2}]=ttest2(squeeze(tmphi_audPtac_comb2(cc,ff,ttime,1:halftest)),squeeze(tmphi_audMSpN_comb2(cc,ff,ttime,1:halftest)));
                      [statsa{ll,tt,ss}.kstestlo_h(cc,ff,ttime,2),statsa{ll,tt,ss}.kstestlo_p(cc,ff,ttime,2),statsa{ll,tt,ss}.kstestlo_stats{cc,ff,ttime,2}]=kstest2(squeeze(tmplo_audPtac_comb2(cc,ff,ttime,1:halftest)),squeeze(tmplo_audMSpN_comb2(cc,ff,ttime,1:halftest)));
                      [statsa{ll,tt,ss}.kstesthi_h(cc,ff,ttime,2),statsa{ll,tt,ss}.kstesthi_p(cc,ff,ttime,2),statsa{ll,tt,ss}.kstesthi_stats{cc,ff,ttime,2}]=kstest2(squeeze(tmphi_audPtac_comb2(cc,ff,ttime,1:halftest)),squeeze(tmphi_audMSpN_comb2(cc,ff,ttime,1:halftest)));
                      [statsa{ll,tt,ss}.ttestlo_h(cc,ff,ttime,3),statsa{ll,tt,ss}.ttestlo_p(cc,ff,ttime,3),statsa{ll,tt,ss}.ttestlo_ci(cc,ff,ttime,3,:),statsa{ll,tt,ss}.ttestlo_stats{cc,ff,ttime,3}]=ttest2(squeeze(tmplo_audPtac_comb2(cc,ff,ttime,halftest+1:end)),squeeze(tmplo_audMSpN_comb2(cc,ff,ttime,halftest+1:end)));
                      [statsa{ll,tt,ss}.ttesthi_h(cc,ff,ttime,3),statsa{ll,tt,ss}.ttesthi_p(cc,ff,ttime,3),statsa{ll,tt,ss}.ttesthi_ci(cc,ff,ttime,3,:),statsa{ll,tt,ss}.ttesthi_stats{cc,ff,ttime,3}]=ttest2(squeeze(tmphi_audPtac_comb2(cc,ff,ttime,halftest+1:end)),squeeze(tmphi_audMSpN_comb2(cc,ff,ttime,halftest+1:end)));
                      [statsa{ll,tt,ss}.kstestlo_h(cc,ff,ttime,3),statsa{ll,tt,ss}.kstestlo_p(cc,ff,ttime,3),statsa{ll,tt,ss}.kstestlo_stats{cc,ff,ttime,3}]=kstest2(squeeze(tmplo_audPtac_comb2(cc,ff,ttime,halftest+1:end)),squeeze(tmplo_audMSpN_comb2(cc,ff,ttime,halftest+1:end)));
                      [statsa{ll,tt,ss}.kstesthi_h(cc,ff,ttime,3),statsa{ll,tt,ss}.kstesthi_p(cc,ff,ttime,3),statsa{ll,tt,ss}.kstesthi_stats{cc,ff,ttime,3}]=kstest2(squeeze(tmphi_audPtac_comb2(cc,ff,ttime,halftest+1:end)),squeeze(tmphi_audMSpN_comb2(cc,ff,ttime,halftest+1:end)));
                    end
                  end
                end
              end
              clear tmplo_* tmphi_*
              
            end % if l<5
            
            
            
          end % ll
        end
      end % ss
      
      % Question: do baseline correction prior to enter in to stats?
      
      clear tlock_tac tlock_aud tlock_nul tr tlock*tlock tlock_tacMSshift*
      
      if dofftadd
        save(['freq_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_tt' num2str(tt) '.mat'],'freq*comb','freq*fftadd','num*trials','-v7.3')
      else
        % NOTE: this below is with trialkc=0 implicit!!
        %           save(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_tt' num2str(tt) '.mat'],'freq*comb','statsa','num*trials','-v7.3')
        save(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_tt' num2str(tt) '_trialkc' num2str(trialkc) '.mat'],'freq*comb','statsa','num*trials','-v7.3')
      end
      clear freq*comb freq*fftadd num*trials
    end % tt
  end
  
  %     end
  
end % ii
% end % sleep

return






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sensor level group stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checking if Null, Tac-alone, and Aud-alone from each of 7 SOA conditions not signficantly diff from each other.
load elec1010_neighb.mat

% allcond_sameN=1;
soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];
tacbasemax=min(-.15,soades-.15);
audbasemax=min(-.15,fliplr(soades)-.15);
msstart=max(0,soades);

ylimlo=[4 7; 8 12; 14 30];
ylimhi=[30 45; 55 80];

soadesc={'Aud first by 500ms' '' 'Aud first by 70ms' 'Aud first by 20ms' 'Simultaneous' 'Tac first by 20ms' 'Tac first by 70ms' '' 'Tac first by 500ms'};
plotflag=1;
printflag=1;
statspowflag=1;
statsplvflag=1;
audtacflag=0;
fftaddflag=0;
synchasynch=0;

soalist=[1 3 4 5 6 7 9];
% soalist=[3 4 5 6 7 9];

chanuse_sleep0={'all' '-F4'};
chanuse_sleep1={'all' '-AF7' '-AF3' '-Fp1'};
for sleep=[1]
  for tt=3
    
    clearvars -except ll tt sub *dir ii*use sleep *flag figind soadesc soalist chanuse* ylim* neigh* *basemax
    if sleep
      chanuse=chanuse_sleep1;
      subuseall=iiBuse;
      iteruse=11;
      trialkc=0;
      %       ssuse=12;
      sleepcond='Sleep N2';
    else
      chanuse=chanuse_sleep0;
      subuseall=setdiff(iiSuse,[]);
      iteruse=27;
      trialkc=-1;
      %       ssuse=10; % awake
      sleepcond='Awake W';
    end
    %     ss=ssuse;
    
    submin=subuseall(1)-1;
    subuseind=0;
    usetr=1;
    iteruse
    
    for ii=subuseall
      cd([edir sub{ii} ])
      %       load(['freq_diffs_averef_' sub{ii} '.mat']);
      %       try
      %         load(['freq_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_tt' num2str(tt) '.mat'])
      %         load(['freq_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(0) '_tt' num2str(tt) '.mat'])
      %       catch
      %         if tt==2
      %           load(['freq_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '.mat'])
      %           load(['freq_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(0) '.mat'])
      %         end
      %       end
      
      % NOTE: this below is with trialkc=0 implicit!!
      %         load(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_tt' num2str(tt) '.mat'])
      try
        load(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_tt' num2str(tt) '_trialkc' num2str(trialkc) '.mat'])
      catch
        if trialkc==0
          load(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_tt' num2str(tt) '.mat'])
        end
      end
      
      
      %       try
      %         tkt=load(['trialkeptTFR_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iteruse) '_usetr' num2str(usetr) '_trialkc' num2str(trialkc) '.mat']);
      %       catch
      %         tkt=load(['trialkeptTFR_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '.mat']);
      %       end
      %         for ss=ssuse
      %       subuse=subuseall; % reset to all for each sleep stage
      subuseind=subuseind+1;
      
      for ss=10:12
        for ll=soalist
          
          numtrt(ll,tt,ss,subuseind)=numt_trials(ll,tt,ss);
          if numtrt(ll,tt,ss,subuseind)<20
            subuse(ll,ss,subuseind)=nan;
          else
            subuse(ll,ss,subuseind)=1;
          end
          
          if ~isnan(subuse(ll,ss,subuseind))
            freqhiall_tNulAlone_comb{ll,ss}{subuseind}=freqhi_tNulAlone_comb{ll,tt,ss};
            freqloall_tNulAlone_comb{ll,ss}{subuseind}=freqlo_tNulAlone_comb{ll,tt,ss};
            freqhiall_tTacAlone_comb{ll,ss}{subuseind}=freqhi_tTacAlone_comb{ll,tt,ss};
            freqloall_tTacAlone_comb{ll,ss}{subuseind}=freqlo_tTacAlone_comb{ll,tt,ss};
            freqhiall_tAudAlone_comb{ll,ss}{subuseind}=freqhi_tAudAlone_comb{ll,tt,ss};
            freqloall_tAudAlone_comb{ll,ss}{subuseind}=freqlo_tAudAlone_comb{ll,tt,ss};
            freqhiall_tMSAlone_comb{ll,ss}{subuseind} =freqhi_tMSAlone_comb{ll,tt,ss};
            freqloall_tMSAlone_comb{ll,ss}{subuseind} =freqlo_tMSAlone_comb{ll,tt,ss};
            freqloall_tNulAlone_comb{ll,ss}{subuseind}.dimord='chan_freq_time';
            freqloall_tTacAlone_comb{ll,ss}{subuseind}.dimord='chan_freq_time';
            freqloall_tAudAlone_comb{ll,ss}{subuseind}.dimord='chan_freq_time';
            freqloall_tMSAlone_comb{ll,ss}{subuseind}.dimord='chan_freq_time';
            freqhiall_tNulAlone_comb{ll,ss}{subuseind}.dimord='chan_freq_time';
            freqhiall_tTacAlone_comb{ll,ss}{subuseind}.dimord='chan_freq_time';
            freqhiall_tAudAlone_comb{ll,ss}{subuseind}.dimord='chan_freq_time';
            freqhiall_tMSAlone_comb{ll,ss}{subuseind}.dimord='chan_freq_time';
          end
        end %ll
      end
    end %ii
    subuseindfinal=subuseind;
    
    
    for ss=10:12
      
      cfg=[];
      cfg.keepindividual='yes';
      cfg.parameter={'powspctrm' 'plvspctrm'};
      for ll=soalist
        usesub=~isnan(subuse(ll,ss,:));
        grindlo_tNulAlone{ll,ss}=ft_freqgrandaverage(cfg,freqloall_tNulAlone_comb{ll,ss}{usesub});
        grindhi_tNulAlone{ll,ss}=ft_freqgrandaverage(cfg,freqhiall_tNulAlone_comb{ll,ss}{usesub});
        grindlo_tTacAlone{ll,ss}=ft_freqgrandaverage(cfg,freqloall_tTacAlone_comb{ll,ss}{usesub});
        grindhi_tTacAlone{ll,ss}=ft_freqgrandaverage(cfg,freqhiall_tTacAlone_comb{ll,ss}{usesub});
        grindlo_tAudAlone{ll,ss}=ft_freqgrandaverage(cfg,freqloall_tAudAlone_comb{ll,ss}{usesub});
        grindhi_tAudAlone{ll,ss}=ft_freqgrandaverage(cfg,freqhiall_tAudAlone_comb{ll,ss}{usesub});
        grindlo_tMSAlone{ll,ss} =ft_freqgrandaverage(cfg,freqloall_tMSAlone_comb{ll,ss}{usesub});
        grindhi_tMSAlone{ll,ss} =ft_freqgrandaverage(cfg,freqhiall_tMSAlone_comb{ll,ss}{usesub});
      end % ll
      
      cfg=[];
      cfg.neighbours=neighbours;
      cfg.method='montecarlo';
      cfg.numrandomization=2000;
      cfg.correctm='cluster';
      cfg.clusteralpha = 0.05;
      cfg.clusterstatistic = 'maxsum';
      cfg.minnbchan = 2;
      cfg.ivar=1;
      cfg.uvar=2;
      
      if statspowflag
        cfg.statistic='depsamplesT';
        cfg.parameter='powspctrm';
        for ll=soalist
          nsub=size(grindlo_tNulAlone{ll,ss}.powspctrm,1);
          cfg.design=zeros(2,2*nsub);
          cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
          cfg.design(2,:)=[1:nsub 1:nsub];
          cfg.latency=[-.15-audbasemax(5); .55-audbasemax(5)];
          stattl_mc_TacVsNul{ll,ss}=ft_freqstatistics(cfg, grindlo_tTacAlone{ll,ss}, grindlo_tNulAlone{ll,ss});
          if ll>=5
            cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
          else
            cfg.latency=[.15+tacbasemax(ll); .85+tacbasemax(ll)];
          end
          stattl_mc_AudVsNul{ll,ss}=ft_freqstatistics(cfg, grindlo_tAudAlone{ll,ss}, grindlo_tNulAlone{ll,ss});
          cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
          stattl_mc_MSVsNul{ll,ss} =ft_freqstatistics(cfg, grindlo_tMSAlone{ll,ss}, grindlo_tNulAlone{ll,ss});
        end
        
      end
      
      if statsplvflag
        cfg.statistic='diff_itc';
        cfg.parameter='plvspctrm';
        cfg.complex='diffabs'; % default; not sensitive to phase differences between conditions
        cfg.clusterthreshold = 'nonparametric_individual'; % or 'nonparametric_individual' see clusterstat.m
        for ll=soalist
          nsub=size(grindlo_tNulAlone{ll,ss}.powspctrm,1);
          cfg.design=zeros(2,2*nsub);
          cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
          cfg.design(2,:)=[1:nsub 1:nsub];
          cfg.latency=[-.15-audbasemax(5); .55-audbasemax(5)];
          stattl_mcplv_TacVsNul{ll,ss}=ft_freqstatistics(cfg, grindlo_tTacAlone{ll,ss}, grindlo_tNulAlone{ll,ss});
          if ll>=5
            cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
          else
            cfg.latency=[.15+tacbasemax(ll); .85+tacbasemax(ll)];
          end
          stattl_mcplv_AudVsNul{ll,ss}=ft_freqstatistics(cfg, grindlo_tAudAlone{ll,ss}, grindlo_tNulAlone{ll,ss});
          cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
          stattl_mcplv_MSVsNul{ll,ss} =ft_freqstatistics(cfg, grindlo_tMSAlone{ll,ss}, grindlo_tNulAlone{ll,ss});
        end
      end
      
      
      % NOTE: without trialkc label, it refers to trialkc=0;
      %     if ~exist([edir 'stats_TFR_UniMSNul_' num2str(tt) num2str(ss) num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'file')
      %       %     save([edir 'statsgrave_TFR_UniMSNul_' num2str(tt) num2str(ss) num2str(sleep) '.mat'],'stat*','grave*');
      %       save([edir 'stats_TFR_UniMSNul_' num2str(tt) num2str(ss) num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'stat*');
      %     else
      %       save([edir 'stats_TFR_UniMSNul_' num2str(tt) num2str(ss) num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'stat*','-append');
      %     end
      if ~exist([edir 'stats_TFR_UniMSNul_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'file')
        %     save([edir 'statsgrave_TFR_UniMSNul_' num2str(tt) num2str(ss) num2str(sleep) '.mat'],'stat*','grave*');
        save([edir 'stats_TFR_UniMSNul_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'stat*','subuse');
      else
        save([edir 'stats_TFR_UniMSNul_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'stat*','subuse','-append');
      end
      
      load([edir 'sortTacWP200.mat'],'iiuse_*halfWwideP2N1');
      
      for ll=soalist
        subind_top=dsearchn(iiBuse',iiuse_tophalfWwideP2N1');
        subind_bot=dsearchn(iiBuse',iiuse_bothalfWwideP2N1');
        usesub=find(squeeze(~isnan(subuse(ll,ss,:))));
        
        grindlo_tTacAlone_top{ll,ss}=grindlo_tTacAlone{ll,ss};
        grindlo_tTacAlone_top{ll,ss}.powspctrm=grindlo_tTacAlone{ll,ss}.powspctrm(dsearchn(usesub,intersect(usesub,subind_top)),:,:,:);
        grindlo_tTacAlone_top{ll,ss}.plvspctrm=grindlo_tTacAlone{ll,ss}.plvspctrm(dsearchn(usesub,intersect(usesub,subind_top)),:,:,:);
        grindlo_tTacAlone_bot{ll,ss}=grindlo_tTacAlone{ll,ss};
        grindlo_tTacAlone_bot{ll,ss}.powspctrm=grindlo_tTacAlone{ll,ss}.powspctrm(dsearchn(usesub,intersect(usesub,subind_bot)),:,:,:);
        grindlo_tTacAlone_bot{ll,ss}.plvspctrm=grindlo_tTacAlone{ll,ss}.plvspctrm(dsearchn(usesub,intersect(usesub,subind_bot)),:,:,:);
        grindlo_tNulAlone_top{ll,ss}=grindlo_tNulAlone{ll,ss};
        grindlo_tNulAlone_top{ll,ss}.powspctrm=grindlo_tNulAlone{ll,ss}.powspctrm(dsearchn(usesub,intersect(usesub,subind_top)),:,:,:);
        grindlo_tNulAlone_top{ll,ss}.plvspctrm=grindlo_tNulAlone{ll,ss}.plvspctrm(dsearchn(usesub,intersect(usesub,subind_top)),:,:,:);
        grindlo_tNulAlone_bot{ll,ss}=grindlo_tNulAlone{ll,ss};
        grindlo_tNulAlone_bot{ll,ss}.powspctrm=grindlo_tNulAlone{ll,ss}.powspctrm(dsearchn(usesub,intersect(usesub,subind_bot)),:,:,:);
        grindlo_tNulAlone_bot{ll,ss}.plvspctrm=grindlo_tNulAlone{ll,ss}.plvspctrm(dsearchn(usesub,intersect(usesub,subind_bot)),:,:,:);
        grindlo_tAudAlone_top{ll,ss}=grindlo_tAudAlone{ll,ss};
        grindlo_tAudAlone_top{ll,ss}.powspctrm=grindlo_tAudAlone{ll,ss}.powspctrm(dsearchn(usesub,intersect(usesub,subind_top)),:,:,:);
        grindlo_tAudAlone_top{ll,ss}.plvspctrm=grindlo_tAudAlone{ll,ss}.plvspctrm(dsearchn(usesub,intersect(usesub,subind_top)),:,:,:);
        grindlo_tAudAlone_bot{ll,ss}=grindlo_tAudAlone{ll,ss};
        grindlo_tAudAlone_bot{ll,ss}.powspctrm=grindlo_tAudAlone{ll,ss}.powspctrm(dsearchn(usesub,intersect(usesub,subind_bot)),:,:,:);
        grindlo_tAudAlone_bot{ll,ss}.plvspctrm=grindlo_tAudAlone{ll,ss}.plvspctrm(dsearchn(usesub,intersect(usesub,subind_bot)),:,:,:);
        grindlo_tMSAlone_top{ll,ss}=grindlo_tMSAlone{ll,ss};
        grindlo_tMSAlone_top{ll,ss}.powspctrm=grindlo_tMSAlone{ll,ss}.powspctrm(dsearchn(usesub,intersect(usesub,subind_top)),:,:,:);
        grindlo_tMSAlone_top{ll,ss}.plvspctrm=grindlo_tMSAlone{ll,ss}.plvspctrm(dsearchn(usesub,intersect(usesub,subind_top)),:,:,:);
        grindlo_tMSAlone_bot{ll,ss}=grindlo_tMSAlone{ll,ss};
        grindlo_tMSAlone_bot{ll,ss}.powspctrm=grindlo_tMSAlone{ll,ss}.powspctrm(dsearchn(usesub,intersect(usesub,subind_bot)),:,:,:);
        grindlo_tMSAlone_bot{ll,ss}.plvspctrm=grindlo_tMSAlone{ll,ss}.plvspctrm(dsearchn(usesub,intersect(usesub,subind_bot)),:,:,:);
      end
      
      cfg=[];
      cfg.neighbours=neighbours;
      cfg.method='montecarlo';
      cfg.numrandomization=2000;
      % cfg.correctm='holm';
      cfg.correctm='cluster';
      cfg.clusteralpha = 0.05;
      cfg.clusterstatistic = 'maxsum';
      cfg.minnbchan = 2;
      cfg.ivar=1;
      cfg.uvar=2;
      
      
      if statspowflag
        cfg.parameter='powspctrm';
        for ll=soalist
          % top
          cfg.statistic='depsamplesT';
          cfg.uvar=2;
          nsub=size(grindlo_tNulAlone_top{ll,ss}.powspctrm,1);
          if nsub>1
            cfg.design=zeros(2,2*nsub);
            cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
            cfg.design(2,:)=[1:nsub 1:nsub];
            
            cfg.latency=[-.15-audbasemax(5); .55-audbasemax(5)];
            stattl_mc_TacVsNul_top{ll,ss}=ft_freqstatistics(cfg, grindlo_tTacAlone_top{ll,ss}, grindlo_tNulAlone_top{ll,ss});
            if ll>=5
              cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            else
              cfg.latency=[.15+tacbasemax(ll); .85+tacbasemax(ll)];
            end
            stattl_mc_AudVsNul_top{ll,ss}=ft_freqstatistics(cfg, grindlo_tAudAlone_top{ll,ss}, grindlo_tNulAlone_top{ll,ss});
            cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            stattl_mc_MSVsNul_top{ll,ss} =ft_freqstatistics(cfg, grindlo_tMSAlone_top{ll,ss}, grindlo_tNulAlone_top{ll,ss});
          else
            stattl_mc_TacVsNul_top{ll,ss}=[];
            stattl_mc_AudVsNul_top{ll,ss}=[];
            stattl_mc_MSVsNul_top{ll,ss} =[];
          end
          
          
          % bottom
          cfg.statistic='depsamplesT';
          cfg.uvar=2;
          nsub=size(grindlo_tNulAlone_bot{ll,ss}.powspctrm,1);
          if nsub>1
            cfg.design=zeros(2,2*nsub);
            cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
            cfg.design(2,:)=[1:nsub 1:nsub];
            
            cfg.latency=[-.15-audbasemax(5); .55-audbasemax(5)];
            stattl_mc_TacVsNul_bot{ll,ss}=ft_freqstatistics(cfg, grindlo_tTacAlone_bot{ll,ss}, grindlo_tNulAlone_bot{ll,ss});
            if ll>=5
              cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            else
              cfg.latency=[.15+tacbasemax(ll); .85+tacbasemax(ll)];
            end
            stattl_mc_AudVsNul_bot{ll,ss}=ft_freqstatistics(cfg, grindlo_tAudAlone_bot{ll,ss}, grindlo_tNulAlone_bot{ll,ss});
            cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            stattl_mc_MSVsNul_bot{ll,ss} =ft_freqstatistics(cfg, grindlo_tMSAlone_bot{ll,ss}, grindlo_tNulAlone_bot{ll,ss});
          else
            stattl_mc_TacVsNul_bot{ll,ss}=[];
            stattl_mc_AudVsNul_bot{ll,ss}=[];
            stattl_mc_MSVsNul_bot{ll,ss} =[];
          end
          
          % top vs bottom
          cfg.statistic='indepsamplesT';
          cfg.uvar=[];
          nsubt=size(grindlo_tNulAlone_top{ll,ss}.powspctrm,1);
          nsubb=size(grindlo_tNulAlone_bot{ll,ss}.powspctrm,1);
          if nsubb>1 && nsubt>1
            cfg.design=zeros(2,nsubb + nsubt);
            cfg.design(1,:)=[ones(1,nsubt) 2*ones(1,nsubb)];
            cfg.design(2,:)=[1:nsubt 1:nsubb];
            
            cfg.latency=[-.15-audbasemax(5); .55-audbasemax(5)];
            stattl_mc_TacTVsTacB{ll,ss}=ft_freqstatistics(cfg, grindlo_tTacAlone_top{ll,ss}, grindlo_tTacAlone_bot{ll,ss});
            if ll>=5
              cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            else
              cfg.latency=[.15+tacbasemax(ll); .85+tacbasemax(ll)];
            end
            stattl_mc_AudTVsAudB{ll,ss}=ft_freqstatistics(cfg, grindlo_tAudAlone_top{ll,ss}, grindlo_tAudAlone_bot{ll,ss});
            cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            stattl_mc_MSTVsMSB{ll,ss}=ft_freqstatistics(cfg, grindlo_tMSAlone_top{ll,ss}, grindlo_tMSAlone_bot{ll,ss});
          else
            stattl_mc_TacTVsTacB{ll,ss}=[];
            stattl_mc_AudTVsAudB{ll,ss}=[];
            stattl_mc_MSTVsMSB{ll,ss}=[];
          end
        end %ll
      end
      
      if statsplvflag
        cfg.statistic='diff_itc';
        cfg.uvar=[];
        cfg.parameter='plvspctrm';
        cfg.complex='diffabs'; % default; not sensitive to phase differences between conditions
        cfg.clusterthreshold = 'nonparametric_individual'; % or 'nonparametric_individual' see clusterstat.m
        for ll=soalist
          
          % top
          nsub=size(grindlo_tNulAlone_top{ll,ss}.powspctrm,1);
          if nsub>1
            cfg.design=zeros(2,2*nsub);
            cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
            cfg.design(2,:)=[1:nsub 1:nsub];
            cfg.latency=[-.15-audbasemax(5); .55-audbasemax(5)];
            stattl_mcplv_TacVsNul_top{ll,ss}=ft_freqstatistics(cfg, grindlo_tTacAlone_top{ll,ss}, grindlo_tNulAlone_top{ll,ss});
            if ll>=5
              cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            else
              cfg.latency=[.15+tacbasemax(ll); .85+tacbasemax(ll)];
            end
            stattl_mcplv_AudVsNul_top{ll,ss}=ft_freqstatistics(cfg, grindlo_tAudAlone_top{ll,ss}, grindlo_tNulAlone_top{ll,ss});
            cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            stattl_mcplv_MSVsNul_top{ll,ss} =ft_freqstatistics(cfg, grindlo_tMSAlone_top{ll,ss}, grindlo_tNulAlone_top{ll,ss});
          else
            stattl_mcplv_TacVsNul_top{ll,ss}=[];
            stattl_mcplv_AudVsNul_top{ll,ss}=[];
            stattl_mcplv_MSVsNul_top{ll,ss} =[];
          end
          
          % bottom
          nsub=size(grindlo_tNulAlone_bot{ll,ss}.powspctrm,1);
          if nsub>1
            cfg.design=zeros(2,2*nsub);
            cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
            cfg.design(2,:)=[1:nsub 1:nsub];
            cfg.latency=[-.15-audbasemax(5); .55-audbasemax(5)];
            stattl_mcplv_TacVsNul_bot{ll,ss}=ft_freqstatistics(cfg, grindlo_tTacAlone_bot{ll,ss}, grindlo_tNulAlone_bot{ll,ss});
            if ll>=5
              cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            else
              cfg.latency=[.15+tacbasemax(ll); .85+tacbasemax(ll)];
            end
            stattl_mcplv_AudVsNul_bot{ll,ss}=ft_freqstatistics(cfg, grindlo_tAudAlone_bot{ll,ss}, grindlo_tNulAlone_bot{ll,ss});
            cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            stattl_mcplv_MSVsNul_bot{ll,ss} =ft_freqstatistics(cfg, grindlo_tMSAlone_bot{ll,ss}, grindlo_tNulAlone_bot{ll,ss});
          else
            stattl_mcplv_TacVsNul_bot{ll,ss}=[];
            stattl_mcplv_AudVsNul_bot{ll,ss}=[];
            stattl_mcplv_MSVsNul_bot{ll,ss} =[];
          end
          
          % top vs bottom
          nsubt=size(grindlo_tNulAlone_top{ll,ss}.powspctrm,1);
          nsubb=size(grindlo_tNulAlone_bot{ll,ss}.powspctrm,1);
          if nsubt>1 && nsubb>1
            cfg.design=zeros(2,nsubb + nsubt);
            cfg.design(1,:)=[ones(1,nsubt) 2*ones(1,nsubb)];
            cfg.design(2,:)=[1:nsubt 1:nsubb];
            cfg.latency=[-.15-audbasemax(5); .55-audbasemax(5)];
            stattl_mcplv_TacTVsTacB{ll,ss}=ft_freqstatistics(cfg, grindlo_tTacAlone_top{ll,ss}, grindlo_tTacAlone_bot{ll,ss});
            if ll>=5
              cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            else
              cfg.latency=[.15+tacbasemax(ll); .85+tacbasemax(ll)];
            end
            stattl_mcplv_AudTVsAudB{ll,ss}=ft_freqstatistics(cfg, grindlo_tAudAlone_top{ll,ss}, grindlo_tAudAlone_bot{ll,ss});
            cfg.latency=[-.15-audbasemax(ll); .55-audbasemax(ll)];
            stattl_mcplv_MSTVsMSB{ll,ss}=ft_freqstatistics(cfg, grindlo_tMSAlone_top{ll,ss}, grindlo_tMSAlone_bot{ll,ss});
          else
          end
        end %ll
        
      end
      
      % NOTE: without trialkc label, it refers to trialkc=0;
      %     save([edir 'stats_TFR_UniMSNul_' num2str(tt) num2str(ss) num2str(sleep) '.mat'],'stat*');
      %     save([edir 'stats_TFR_UniMSNul_' num2str(tt) num2str(ss) num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'stat*','-append');
      save([edir 'stats_TFR_UniMSNul_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'stat*','-append');
    end % ss
    
  end % tt
end % sleep


%%  Main stats

load elec1010_neighb.mat

% allcond_sameN=1;
soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];
tacbasemax=min(-.15,soades-.15);
audbasemax=min(-.15,fliplr(soades)-.15);

ylimlo=[4 7; 8 12; 14 30];
ylimhi=[30 45; 55 80];

soadesc={'Aud first by 500ms' '' 'Aud first by 70ms' 'Aud first by 20ms' 'Simultaneous' 'Tac first by 20ms' 'Tac first by 70ms' '' 'Tac first by 500ms'};
plotflag=0;
printflag=0;
audtacflag=0;
comb2flag=1;
fftaddflag=0;
synchasynch=0;
mcseed=13;  % montecarlo cfg.randomseed
usetr=3;

soalist=[1 3 4 5 6 7 9];
% soalist=[3 4 5 6 7 9];

chanuse_sleep0={'all' '-F4'};
chanuse_sleep1={'all' '-AF7' '-AF3' '-Fp1'};

tt=3;
sleep=0;
% for sleep=[1]
if sleep
  chanuse=chanuse_sleep1;
else
  chanuse=chanuse_sleep0;
end
%   for tt=[3]
figind=1;


for ll=soalist
% for ll=[4 5 6 7 9];
  clearvars -except ll tt sub *dir ii*use sleep *flag figind soadesc soalist chanuse* ylim* neigh* *basemax synch* mcseed usetr
  
  statsflag=1;

if sleep
    subuseall=iiBuse;
    iter=11;
    trialkc=0;
  else
    subuseall=setdiff(iiSuse,[])
%     iter=27;
    iter=31;
    trialkc=-1;
  end
  
  submin=subuseall(1)-1;
  subuseind=0;
  
  % Baseline correct each participant prior to entering to stats???? NO
  for ii=subuseall
%   for ii=8:9
    cd([edir sub{ii} ])
    %       load(['freq_diffs_averef_' sub{ii} '.mat']);
    try
      if fftaddflag
        load(['freq_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_tt' num2str(tt) '.mat'])
      else
%         load(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_tt' num2str(tt) '.mat'])
        load(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'])
      end
      if audtacflag
        load(['freq_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(0) '_tt' num2str(tt) '.mat'])
      end
    catch
      if tt==2
        load(['freq_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '.mat'])
        if audtacflag
          load(['freq_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(0) '.mat'])
        end
      end
    end
    if audtacflag
      tka=load(['trialkeptTFR_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(0) '_iter' num2str(iter) '_usetr' num2str(usetr) '.mat']);
    end
    try
      if usetr==2
        % load usetr=0 here; then later down load usetr=2
        tkt=load(['trialkeptTFR_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter) '_usetr' num2str(3) '_trialkc' num2str(trialkc) '.mat']);
      else
        tkt=load(['trialkeptTFR_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter) '_usetr' num2str(usetr) '.mat']);
      end
    catch
      tkt=load(['trialkeptTFR_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter) '_usetr' num2str(usetr) '_trialkc' num2str(trialkc) '.mat']);
    end
    
    if 1
      % THis is preliminary...how best to included all stages later on?
      if sleep==0
        ssuse=10; % awake
        sleepcond='Awake W';
      elseif sleep==1
        ssuse=12; % this is concatenation of N2 and N3
        sleepcond='Sleep N2';
        %             ssuse=23; % this is concatenation of N2 and N3
        %             sleepcond='Sleep (N2+N3)';
      end
    else
      if sleep
        ssuse=tkt.tr.stageuse;
      else
        ssuse=tka.tr.stageuse;
      end
    end
    
    ss=ssuse;
    %         for ss=ssuse
    subuse=subuseall; % reset to all for each sleep stage
    numtrt(ll,tt,ss,ii-submin)=tkt.numcondtfinal(ll,tt,ss);
    if audtacflag
      numtra(ll,tt,ss,ii-submin)=tka.numcondafinal(ll,tt,ss);
    end
    
    
    if min(numtrt(ll,tt,ss,ii-submin),[],1)<20 % what is best number to use here?
      subuse=setdiff(subuse,ii);
    else
      subuseind=subuseind+1;
      freqloall_tacPaud_comb1{subuseind,1}=freqlo_tacPaud_comb{ll,tt,ss,1};
      freqloall_tacMSpN_comb1{subuseind,1}=freqlo_tacMSpN_comb{ll,tt,ss,1};
      freqhiall_tacPaud_comb1{subuseind,1}=freqhi_tacPaud_comb{ll,tt,ss,1};
      freqhiall_tacMSpN_comb1{subuseind,1}=freqhi_tacMSpN_comb{ll,tt,ss,1};
      freqhiall_tNulAlone_comb{subuseind,1}=freqhi_tNulAlone_comb{ll,tt,ss};
      freqloall_tNulAlone_comb{subuseind,1}=freqlo_tNulAlone_comb{ll,tt,ss};
      freqhiall_tTacAlone_comb{subuseind,1}=freqhi_tTacAlone_comb{ll,tt,ss};
      freqloall_tTacAlone_comb{subuseind,1}=freqlo_tTacAlone_comb{ll,tt,ss};
      freqhiall_tAudAlone_comb{subuseind,1}=freqhi_tAudAlone_comb{ll,tt,ss};
      freqloall_tAudAlone_comb{subuseind,1}=freqlo_tAudAlone_comb{ll,tt,ss};
      freqhiall_tMSAlone_comb{subuseind,1}=freqhi_tMSAlone_comb{ll,tt,ss};
      freqloall_tMSAlone_comb{subuseind,1}=freqlo_tMSAlone_comb{ll,tt,ss};
      freqloall_tacPaud_comb1{subuseind,1}.dimord='chan_freq_time';
      freqloall_tacMSpN_comb1{subuseind,1}.dimord='chan_freq_time';
      freqhiall_tacPaud_comb1{subuseind,1}.dimord='chan_freq_time';
      freqhiall_tacMSpN_comb1{subuseind,1}.dimord='chan_freq_time';
      if synchasynch && ll<5
        freqloall_tMSasynch_comb1{subuseind}=freqlo_tMSasynch_comb{ll,tt,ss};
        freqloall_tMSsynch_comb1{subuseind}=freqlo_tMSsynch_comb{ll,tt,ss};
        freqhiall_tMSasynch_comb1{subuseind}=freqhi_tMSasynch_comb{ll,tt,ss};
        freqhiall_tMSsynch_comb1{subuseind}=freqhi_tMSsynch_comb{ll,tt,ss};
        freqloall_tMSasynch_comb1{subuseind}.dimord='chan_freq_time';
        freqloall_tMSsynch_comb1{subuseind}.dimord='chan_freq_time';
        freqhiall_tMSasynch_comb1{subuseind}.dimord='chan_freq_time';
        freqhiall_tMSsynch_comb1{subuseind}.dimord='chan_freq_time';
      end
      if comb2flag
        freqloall_tacPaud_comb2{subuseind,1}=freqlo_tacPaud_comb{ll,tt,ss,2};
        freqloall_tacMSpN_comb2{subuseind,1}=freqlo_tacMSpN_comb{ll,tt,ss,2};
        freqhiall_tacPaud_comb2{subuseind,1}=freqhi_tacPaud_comb{ll,tt,ss,2};
        freqhiall_tacMSpN_comb2{subuseind,1}=freqhi_tacMSpN_comb{ll,tt,ss,2};
        freqloall_tacPaud_comb2{subuseind,1}.dimord='chan_freq_time';
        freqloall_tacMSpN_comb2{subuseind,1}.dimord='chan_freq_time';
        freqhiall_tacPaud_comb2{subuseind,1}.dimord='chan_freq_time';
        freqhiall_tacMSpN_comb2{subuseind,1}.dimord='chan_freq_time';
        if synchasynch ll<5
          freqloall_tMSasynch_comb2{subuseind}=freqlo_tMSasynch_comb{ll,tt,ss};
          freqloall_tMSsynch_comb2{subuseind}=freqlo_tMSsynch_comb{ll,tt,ss};
          freqhiall_tMSasynch_comb2{subuseind}=freqhi_tMSasynch_comb{ll,tt,ss};
          freqhiall_tMSsynch_comb2{subuseind}=freqhi_tMSsynch_comb{ll,tt,ss};
          freqloall_tMSasynch_comb2{subuseind}.dimord='chan_freq_time';
          freqloall_tMSsynch_comb2{subuseind}.dimord='chan_freq_time';
          freqhiall_tMSasynch_comb2{subuseind}.dimord='chan_freq_time';
          freqhiall_tMSsynch_comb2{subuseind}.dimord='chan_freq_time';
        end
      end
      freqhiall_tNulAlone_comb{subuseind,1}.dimord='chan_freq_time';
      freqloall_tNulAlone_comb{subuseind,1}.dimord='chan_freq_time';
      freqhiall_tTacAlone_comb{subuseind,1}.dimord='chan_freq_time';
      freqloall_tTacAlone_comb{subuseind,1}.dimord='chan_freq_time';
      freqhiall_tAudAlone_comb{subuseind,1}.dimord='chan_freq_time';
      freqloall_tAudAlone_comb{subuseind,1}.dimord='chan_freq_time';
      freqhiall_tMSAlone_comb{subuseind,1}.dimord='chan_freq_time';
      freqloall_tMSAlone_comb{subuseind,1}.dimord='chan_freq_time';
      if fftaddflag
        freqloall_tacPaud_fftadd{subuseind}=freqlo_tacPaud_fftadd{ll,tt,ss};
        freqloall_tacMSpN_fftadd{subuseind}=freqlo_tacMSpN_fftadd{ll,tt,ss};
        freqhiall_tacPaud_fftadd{subuseind}=freqhi_tacPaud_fftadd{ll,tt,ss};
        freqhiall_tacMSpN_fftadd{subuseind}=freqhi_tacMSpN_fftadd{ll,tt,ss};
        freqhiall_tNulAlone_fftadd{subuseind}=freqhi_tNulAlone_fftadd{ll,tt,ss};
        freqloall_tNulAlone_fftadd{subuseind}=freqlo_tNulAlone_fftadd{ll,tt,ss};
        freqhiall_tTacAlone_fftadd{subuseind}=freqhi_tTacAlone_fftadd{ll,tt,ss};
        freqloall_tTacAlone_fftadd{subuseind}=freqlo_tTacAlone_fftadd{ll,tt,ss};
        freqhiall_tAudAlone_fftadd{subuseind}=freqhi_tAudAlone_fftadd{ll,tt,ss};
        freqloall_tAudAlone_fftadd{subuseind}=freqlo_tAudAlone_fftadd{ll,tt,ss};
        freqhiall_tMSAlone_fftadd{subuseind}=freqhi_tMSAlone_fftadd{ll,tt,ss};
        freqloall_tMSAlone_fftadd{subuseind}=freqlo_tMSAlone_fftadd{ll,tt,ss};
        freqloall_tacPaud_fftadd{subuseind}.dimord='chan_freq_time';
        freqloall_tacMSpN_fftadd{subuseind}.dimord='chan_freq_time';
        freqhiall_tacPaud_fftadd{subuseind}.dimord='chan_freq_time';
        freqhiall_tacMSpN_fftadd{subuseind}.dimord='chan_freq_time';
        freqhiall_tNulAlone_fftadd{subuseind}.dimord='chan_freq_time';
        freqloall_tNulAlone_fftadd{subuseind}.dimord='chan_freq_time';
        freqhiall_tTacAlone_fftadd{subuseind}.dimord='chan_freq_time';
        freqloall_tTacAlone_fftadd{subuseind}.dimord='chan_freq_time';
        freqhiall_tAudAlone_fftadd{subuseind}.dimord='chan_freq_time';
        freqloall_tAudAlone_fftadd{subuseind}.dimord='chan_freq_time';
        freqhiall_tMSAlone_fftadd{subuseind}.dimord='chan_freq_time';
        freqloall_tMSAlone_fftadd{subuseind}.dimord='chan_freq_time';
        freqloall_tacPaud_fftadd{subuseind}.crsspctrm=[];
        freqloall_tacMSpN_fftadd{subuseind}.crsspctrm=[];
        freqhiall_tacPaud_fftadd{subuseind}.crsspctrm=[];
        freqhiall_tacMSpN_fftadd{subuseind}.crsspctrm=[];
        freqhiall_tNulAlone_fftadd{subuseind}.crsspctrm=[];
        freqloall_tNulAlone_fftadd{subuseind}.crsspctrm=[];
        freqhiall_tTacAlone_fftadd{subuseind}.crsspctrm=[];
        freqloall_tTacAlone_fftadd{subuseind}.crsspctrm=[];
        freqhiall_tAudAlone_fftadd{subuseind}.crsspctrm=[];
        freqloall_tAudAlone_fftadd{subuseind}.crsspctrm=[];
        freqhiall_tMSAlone_fftadd{subuseind}.crsspctrm=[];
        freqloall_tMSAlone_fftadd{subuseind}.crsspctrm=[];
      end
      
      if audtacflag
        freqloall_audPtac_comb1{subuseind}=freqlo_audPtac_comb{ll,tt,ss,1};
        freqloall_audMSpN_comb1{subuseind}=freqlo_audMSpN_comb{ll,tt,ss,1};
        freqhiall_audPtac_comb1{subuseind}=freqhi_audPtac_comb{ll,tt,ss,1};
        freqhiall_audMSpN_comb1{subuseind}=freqhi_audMSpN_comb{ll,tt,ss,1};
        freqhiall_aNulAlone_comb{subuseind}=freqhi_aNulAlone_comb{ll,tt,ss};
        freqloall_aNulAlone_comb{subuseind}=freqlo_aNulAlone_comb{ll,tt,ss};
        freqhiall_aTacAlone_comb{subuseind}=freqhi_aTacAlone_comb{ll,tt,ss};
        freqloall_aTacAlone_comb{subuseind}=freqlo_aTacAlone_comb{ll,tt,ss};
        freqhiall_aAudAlone_comb{subuseind}=freqhi_aAudAlone_comb{ll,tt,ss};
        freqloall_aAudAlone_comb{subuseind}=freqlo_aAudAlone_comb{ll,tt,ss};
        freqhiall_aMSAlone_comb{subuseind}=freqhi_aMSAlone_comb{ll,tt,ss};
        freqloall_aMSAlone_comb{subuseind}=freqlo_aMSAlone_comb{ll,tt,ss};
        freqloall_audPtac_comb1{subuseind}.dimord='chan_freq_time';
        freqloall_audMSpN_comb1{subuseind}.dimord='chan_freq_time';
        freqhiall_audPtac_comb1{subuseind}.dimord='chan_freq_time';
        freqhiall_audMSpN_comb1{subuseind}.dimord='chan_freq_time';
        if comb2flag
          freqloall_audPtac_comb2{subuseind}=freqlo_audPtac_comb{ll,tt,ss,2};
          freqloall_audMSpN_comb2{subuseind}=freqlo_audMSpN_comb{ll,tt,ss,2};
          freqhiall_audPtac_comb2{subuseind}=freqhi_audPtac_comb{ll,tt,ss,2};
          freqhiall_audMSpN_comb2{subuseind}=freqhi_audMSpN_comb{ll,tt,ss,2};
          freqloall_audPtac_comb2{subuseind}.dimord='chan_freq_time';
          freqloall_audMSpN_comb2{subuseind}.dimord='chan_freq_time';
          freqhiall_audPtac_comb2{subuseind}.dimord='chan_freq_time';
          freqhiall_audMSpN_comb2{subuseind}.dimord='chan_freq_time';
        end
        freqhiall_aNulAlone_comb{subuseind}.dimord='chan_freq_time';
        freqloall_aNulAlone_comb{subuseind}.dimord='chan_freq_time';
        freqhiall_aTacAlone_comb{subuseind}.dimord='chan_freq_time';
        freqloall_aTacAlone_comb{subuseind}.dimord='chan_freq_time';
        freqhiall_aAudAlone_comb{subuseind}.dimord='chan_freq_time';
        freqloall_aAudAlone_comb{subuseind}.dimord='chan_freq_time';
        freqhiall_aMSAlone_comb{subuseind}.dimord='chan_freq_time';
        freqloall_aMSAlone_comb{subuseind}.dimord='chan_freq_time';
        try
          freqloall_audPtac_fftadd{subuseind}=freqlo_audPtac_fftadd{ll,tt,ss};
          freqloall_audMSpN_fftadd{subuseind}=freqlo_audMSpN_fftadd{ll,tt,ss};
          freqhiall_audPtac_fftadd{subuseind}=freqhi_audPtac_fftadd{ll,tt,ss};
          freqhiall_audMSpN_fftadd{subuseind}=freqhi_audMSpN_fftadd{ll,tt,ss};
          freqhiall_aNulAlone_fftadd{subuseind}=freqhi_aNulAlone_fftadd{ll,tt,ss};
          freqloall_aNulAlone_fftadd{subuseind}=freqlo_aNulAlone_fftadd{ll,tt,ss};
          freqhiall_aTacAlone_fftadd{subuseind}=freqhi_aTacAlone_fftadd{ll,tt,ss};
          freqloall_aTacAlone_fftadd{subuseind}=freqlo_aTacAlone_fftadd{ll,tt,ss};
          freqhiall_aAudAlone_fftadd{subuseind}=freqhi_aAudAlone_fftadd{ll,tt,ss};
          freqloall_aAudAlone_fftadd{subuseind}=freqlo_aAudAlone_fftadd{ll,tt,ss};
          freqhiall_aMSAlone_fftadd{subuseind}=freqhi_aMSAlone_fftadd{ll,tt,ss};
          freqloall_aMSAlone_fftadd{subuseind}=freqlo_aMSAlone_fftadd{ll,tt,ss};
          freqloall_audPtac_fftadd{subuseind}.dimord='chan_freq_time';
          freqloall_audMSpN_fftadd{subuseind}.dimord='chan_freq_time';
          freqhiall_audPtac_fftadd{subuseind}.dimord='chan_freq_time';
          freqhiall_audMSpN_fftadd{subuseind}.dimord='chan_freq_time';
          freqhiall_aNulAlone_fftadd{subuseind}.dimord='chan_freq_time';
          freqloall_aNulAlone_fftadd{subuseind}.dimord='chan_freq_time';
          freqhiall_aTacAlone_fftadd{subuseind}.dimord='chan_freq_time';
          freqloall_aTacAlone_fftadd{subuseind}.dimord='chan_freq_time';
          freqhiall_aAudAlone_fftadd{subuseind}.dimord='chan_freq_time';
          freqloall_aAudAlone_fftadd{subuseind}.dimord='chan_freq_time';
          freqhiall_aMSAlone_fftadd{subuseind}.dimord='chan_freq_time';
          freqloall_aMSAlone_fftadd{subuseind}.dimord='chan_freq_time';
          freqloall_audPtac_fftadd{subuseind}.crsspctrm=[];
          freqloall_audMSpN_fftadd{subuseind}.crsspctrm=[];
          freqhiall_audPtac_fftadd{subuseind}.crsspctrm=[];
          freqhiall_audMSpN_fftadd{subuseind}.crsspctrm=[];
          freqhiall_aNulAlone_fftadd{subuseind}.crsspctrm=[];
          freqloall_aNulAlone_fftadd{subuseind}.crsspctrm=[];
          freqhiall_aTacAlone_fftadd{subuseind}.crsspctrm=[];
          freqloall_aTacAlone_fftadd{subuseind}.crsspctrm=[];
          freqhiall_aAudAlone_fftadd{subuseind}.crsspctrm=[];
          freqloall_aAudAlone_fftadd{subuseind}.crsspctrm=[];
          freqhiall_aMSAlone_fftadd{subuseind}.crsspctrm=[];
          freqloall_aMSAlone_fftadd{subuseind}.crsspctrm=[];
        end
      end
    end
    %         end % ss
    clear freqlo_* freqhi_*
    
    if usetr==2
      load(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter+1) '_trialkc' num2str(trialkc) '.mat'])
%       tkt=load(['trialkeptTFR_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter+1) '_usetr' num2str(2) '_trialkc' num2str(trialkc) '.mat']);

      freqloall_tacPaud_comb1{subuseind,2}=freqlo_tacPaud_comb{ll,tt,ss,1};
      freqloall_tacMSpN_comb1{subuseind,2}=freqlo_tacMSpN_comb{ll,tt,ss,1};
      freqloall_tNulAlone_comb{subuseind,2}=freqlo_tNulAlone_comb{ll,tt,ss};
      freqloall_tTacAlone_comb{subuseind,2}=freqlo_tTacAlone_comb{ll,tt,ss};
      freqloall_tAudAlone_comb{subuseind,2}=freqlo_tAudAlone_comb{ll,tt,ss};
      freqloall_tMSAlone_comb{subuseind,2}=freqlo_tMSAlone_comb{ll,tt,ss};
      freqloall_tacPaud_comb1{subuseind,2}.dimord='chan_freq_time';
      freqloall_tacMSpN_comb1{subuseind,2}.dimord='chan_freq_time';
      freqloall_tNulAlone_comb{subuseind,2}.dimord='chan_freq_time';
      freqloall_tTacAlone_comb{subuseind,2}.dimord='chan_freq_time';
      freqloall_tAudAlone_comb{subuseind,2}.dimord='chan_freq_time';
      freqloall_tMSAlone_comb{subuseind,2}.dimord='chan_freq_time';
      if comb2flag
        freqloall_tacPaud_comb2{subuseind,2}=freqlo_tacPaud_comb{ll,tt,ss,2};
        freqloall_tacMSpN_comb2{subuseind,2}=freqlo_tacMSpN_comb{ll,tt,ss,2};
        freqloall_tacPaud_comb2{subuseind,2}.dimord='chan_freq_time';
        freqloall_tacMSpN_comb2{subuseind,2}.dimord='chan_freq_time';
      end
      freqhiall_tacPaud_comb1{subuseind,2}=freqhi_tacPaud_comb{ll,tt,ss,1};
      freqhiall_tacMSpN_comb1{subuseind,2}=freqhi_tacMSpN_comb{ll,tt,ss,1};
      freqhiall_tNulAlone_comb{subuseind,2}=freqhi_tNulAlone_comb{ll,tt,ss};
      freqhiall_tTacAlone_comb{subuseind,2}=freqhi_tTacAlone_comb{ll,tt,ss};
      freqhiall_tAudAlone_comb{subuseind,2}=freqhi_tAudAlone_comb{ll,tt,ss};
      freqhiall_tMSAlone_comb{subuseind,2}=freqhi_tMSAlone_comb{ll,tt,ss};
      freqhiall_tacPaud_comb1{subuseind,2}.dimord='chan_freq_time';
      freqhiall_tacMSpN_comb1{subuseind,2}.dimord='chan_freq_time';
      freqhiall_tNulAlone_comb{subuseind,2}.dimord='chan_freq_time';
      freqhiall_tTacAlone_comb{subuseind,2}.dimord='chan_freq_time';
      freqhiall_tAudAlone_comb{subuseind,2}.dimord='chan_freq_time';
      freqhiall_tMSAlone_comb{subuseind,2}.dimord='chan_freq_time';
      if comb2flag
        freqhiall_tacPaud_comb2{subuseind,2}=freqhi_tacPaud_comb{ll,tt,ss,2};
        freqhiall_tacMSpN_comb2{subuseind,2}=freqhi_tacMSpN_comb{ll,tt,ss,2};
        freqhiall_tacPaud_comb2{subuseind,2}.dimord='chan_freq_time';
        freqhiall_tacMSpN_comb2{subuseind,2}.dimord='chan_freq_time';
      end
    end
  end % ii
  
  
  
  for ii=1:subuseind
    cfg=[];
    cfg.operation='subtract';
    cfg.parameter='powspctrm';
    if fftaddflag
      freqloall_TPA_MSPN_fftadd{ii}=ft_math(cfg,freqloall_tacPaud_fftadd{ii},freqloall_tacMSpN_fftadd{ii});
      freqhiall_TPA_MSPN_fftadd{ii}=ft_math(cfg,freqhiall_tacPaud_fftadd{ii},freqhiall_tacMSpN_fftadd{ii});
      if 0
        freqloall_TPA_MSPN_fftaddbn{ii}=ft_math(cfg,freqloall_tacPaud_fftaddbn{ii},freqloall_tacMSpN_fftaddbn{ii});
        freqhiall_TPA_MSPN_fftaddbn{ii}=ft_math(cfg,freqhiall_tacPaud_fftaddbn{ii},freqhiall_tacMSpN_fftaddbn{ii});
        freqloall_TPA_MSPN_comb1bn{ii}=ft_math(cfg,freqloall_tacPaud_comb1bn{ii},freqloall_tacMSpN_comb1bn{ii});
        freqhiall_TPA_MSPN_comb1bn{ii}=ft_math(cfg,freqhiall_tacPaud_comb1bn{ii},freqhiall_tacMSpN_comb1bn{ii});
      end
      
      if audtacflag
        freqloall_APT_MSPN_fftadd{ii}=ft_math(cfg,freqloall_audPtac_fftadd{ii},freqloall_audMSpN_fftadd{ii});
        freqhiall_APT_MSPN_fftadd{ii}=ft_math(cfg,freqhiall_audPtac_fftadd{ii},freqhiall_audMSpN_fftadd{ii});
        if 0
          freqloall_APT_MSPN_fftaddbn{ii}=ft_math(cfg,freqloall_audPtac_fftaddbn{ii},freqloall_audMSpN_fftaddbn{ii});
          freqhiall_APT_MSPN_fftaddbn{ii}=ft_math(cfg,freqhiall_audPtac_fftaddbn{ii},freqhiall_audMSpN_fftaddbn{ii});
          freqloall_APT_MSPN_comb1bn{ii}=ft_math(cfg,freqloall_audPtac_comb1bn{ii},freqloall_audMSpN_comb1bn{ii});
          freqhiall_APT_MSPN_comb1bn{ii}=ft_math(cfg,freqhiall_audPtac_comb1bn{ii},freqhiall_audMSpN_comb1bn{ii});
        end
      end
    end
    cfg.parameter={'powspctrm' 'plvspctrm'};
    if usetr~=2
      iterinduse=1;
    else
      iterinduse=2;
    end
    for iterind=1:iterinduse
      freqloall_TPA_MSPN_comb1{ii,iterind}=ft_math(cfg,freqloall_tacPaud_comb1{ii,iterind},freqloall_tacMSpN_comb1{ii,iterind});
      freqhiall_TPA_MSPN_comb1{ii,iterind}=ft_math(cfg,freqhiall_tacPaud_comb1{ii,iterind},freqhiall_tacMSpN_comb1{ii,iterind});
    end
    if synchasynch && ll<5
      freqloall_TMSs_TMSa_comb1{ii}=ft_math(cfg,freqloall_tMSsynch_comb1{ii},freqloall_tMSasynch_comb1{ii});
      freqhiall_TMSs_TMSa_comb1{ii}=ft_math(cfg,freqhiall_tMSsynch_comb1{ii},freqhiall_tMSasynch_comb1{ii});
    end
    if comb2flag
      for iterind=1:iterinduse
        freqloall_TPA_MSPN_comb2{ii,iterind}=ft_math(cfg,freqloall_tacPaud_comb2{ii,iterind},freqloall_tacMSpN_comb2{ii,iterind});
        freqhiall_TPA_MSPN_comb2{ii,iterind}=ft_math(cfg,freqhiall_tacPaud_comb2{ii,iterind},freqhiall_tacMSpN_comb2{ii,iterind});
      end
      if synchasynch && ll<5
        freqloall_TMSs_TMSa_comb2{ii}=ft_math(cfg,freqloall_tMSsynch_comb2{ii},freqloall_tMSasynch_comb2{ii});
        freqhiall_TMSs_TMSa_comb2{ii}=ft_math(cfg,freqhiall_tMSsynch_comb2{ii},freqhiall_tMSasynch_comb2{ii});
      end
    end
    
    if audtacflag
      freqloall_APT_MSPN_comb1{ii}=ft_math(cfg,freqloall_audPtac_comb1{ii},freqloall_audMSpN_comb1{ii});
      freqhiall_APT_MSPN_comb1{ii}=ft_math(cfg,freqhiall_audPtac_comb1{ii},freqhiall_audMSpN_comb1{ii});
      if comb2flag
        freqloall_APT_MSPN_comb2{ii}=ft_math(cfg,freqloall_audPtac_comb2{ii},freqloall_audMSpN_comb2{ii});
        freqhiall_APT_MSPN_comb2{ii}=ft_math(cfg,freqhiall_audPtac_comb2{ii},freqhiall_audMSpN_comb2{ii});
      end
    end
    if usetr==2
      cfg.operation='(x1+x2)/2';
      freqloall_TPA_MSPN_comb1_usetr2{ii}=ft_math(cfg,freqloall_TPA_MSPN_comb1{ii,:});
      freqhiall_TPA_MSPN_comb1_usetr2{ii}=ft_math(cfg,freqhiall_TPA_MSPN_comb1{ii,:});
      freqloall_TPA_MSPN_comb2_usetr2{ii}=ft_math(cfg,freqloall_TPA_MSPN_comb2{ii,:});
      freqhiall_TPA_MSPN_comb2_usetr2{ii}=ft_math(cfg,freqhiall_TPA_MSPN_comb2{ii,:});
      freqloall_tacPaud_comb1tmp{ii}=ft_math(cfg,freqloall_tacPaud_comb1{ii,:});
      freqloall_tacPaud_comb2tmp{ii}=ft_math(cfg,freqloall_tacPaud_comb2{ii,:});
      freqhiall_tacPaud_comb1tmp{ii}=ft_math(cfg,freqhiall_tacPaud_comb1{ii,:});
      freqhiall_tacPaud_comb2tmp{ii}=ft_math(cfg,freqhiall_tacPaud_comb2{ii,:});
      freqloall_tacMSpN_comb1tmp{ii}=ft_math(cfg,freqloall_tacMSpN_comb1{ii,:});
      freqloall_tacMSpN_comb2tmp{ii}=ft_math(cfg,freqloall_tacMSpN_comb2{ii,:});
      freqhiall_tacMSpN_comb1tmp{ii}=ft_math(cfg,freqhiall_tacMSpN_comb1{ii,:});
      freqhiall_tacMSpN_comb2tmp{ii}=ft_math(cfg,freqhiall_tacMSpN_comb2{ii,:});
      freqloall_tTacAlone_combtmp{ii}=ft_math(cfg,freqloall_tTacAlone_comb{ii,:});
      freqloall_tAudAlone_combtmp{ii}=ft_math(cfg,freqloall_tAudAlone_comb{ii,:});
      freqloall_tNulAlone_combtmp{ii}=ft_math(cfg,freqloall_tNulAlone_comb{ii,:});
      freqloall_tMSAlone_combtmp{ii}=ft_math(cfg,freqloall_tMSAlone_comb{ii,:});
      freqhiall_tTacAlone_combtmp{ii}=ft_math(cfg,freqhiall_tTacAlone_comb{ii,:});
      freqhiall_tAudAlone_combtmp{ii}=ft_math(cfg,freqhiall_tAudAlone_comb{ii,:});
      freqhiall_tNulAlone_combtmp{ii}=ft_math(cfg,freqhiall_tNulAlone_comb{ii,:});
      freqhiall_tMSAlone_combtmp{ii}=ft_math(cfg,freqhiall_tMSAlone_comb{ii,:});
    end
  end % end ii
  if usetr==2
    freqloall_tacPaud_comb1=freqloall_tacPaud_comb1tmp;
    freqloall_tacPaud_comb2=freqloall_tacPaud_comb2tmp;
    freqhiall_tacPaud_comb1=freqhiall_tacPaud_comb1tmp;
    freqhiall_tacPaud_comb2=freqhiall_tacPaud_comb2tmp;
    freqloall_tacMSpN_comb1=freqloall_tacMSpN_comb1tmp;
    freqloall_tacMSpN_comb2=freqloall_tacMSpN_comb2tmp;
    freqhiall_tacMSpN_comb1=freqhiall_tacMSpN_comb1tmp;
    freqhiall_tacMSpN_comb2=freqhiall_tacMSpN_comb2tmp;
    freqloall_tTacAlone_comb=freqloall_tTacAlone_combtmp;
    freqloall_tAudAlone_comb=freqloall_tAudAlone_combtmp;
    freqloall_tNulAlone_comb=freqloall_tNulAlone_combtmp;
    freqloall_tMSAlone_comb=freqloall_tMSAlone_combtmp;
    freqhiall_tTacAlone_comb=freqhiall_tTacAlone_combtmp;
    freqhiall_tAudAlone_comb=freqhiall_tAudAlone_combtmp;
    freqhiall_tNulAlone_comb=freqhiall_tNulAlone_combtmp;
    freqhiall_tMSAlone_comb=freqhiall_tMSAlone_combtmp;
    clear freq*tmp
  end
  
  cfg=[];
  if fftaddflag
    gravelo_TPA_MSPN_fftadd=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_fftadd{:});
    gravehi_TPA_MSPN_fftadd=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_fftadd{:});
    if 0
      gravelo_TPA_MSPN_fftaddbn=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_fftaddbn{:});
      gravehi_TPA_MSPN_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_fftaddbn{:});
      gravelo_TPA_MSPN_comb1bn=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb1bn{:});
      gravehi_TPA_MSPN_comb1bn=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb1bn{:});
    end
    
    if audtacflag
      gravelo_APT_MSPN_fftadd=ft_freqgrandaverage(cfg,freqloall_APT_MSPN_fftadd{:});
      gravehi_APT_MSPN_fftadd=ft_freqgrandaverage(cfg,freqhiall_APT_MSPN_fftadd{:});
      if 0
        gravelo_APT_MSPN_fftaddbn=ft_freqgrandaverage(cfg,freqloall_APT_MSPN_fftaddbn{:});
        gravehi_APT_MSPN_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_APT_MSPN_fftaddbn{:});
        gravelo_APT_MSPN_comb1bn=ft_freqgrandaverage(cfg,freqloall_APT_MSPN_comb1bn{:});
        gravehi_APT_MSPN_comb1bn=ft_freqgrandaverage(cfg,freqhiall_APT_MSPN_comb1bn{:});
      end
    end
  end
  cfg.parameter={'powspctrm' 'plvspctrm'};
  if usetr==2
    gravelo_TPA_MSPN_comb1=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb1_usetr2{:});
    gravehi_TPA_MSPN_comb1=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb1_usetr2{:});
  else
    gravelo_TPA_MSPN_comb1=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb1{:});
    gravehi_TPA_MSPN_comb1=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb1{:});
  end
  if synchasynch && ll<5
    gravelo_TMSs_TMSa_comb1=ft_freqgrandaverage(cfg,freqloall_TMSs_TMSa_comb1{:});
    gravehi_TMSs_TMSa_comb1=ft_freqgrandaverage(cfg,freqhiall_TMSs_TMSa_comb1{:});
  end
  if comb2flag
    if usetr==2
      gravelo_TPA_MSPN_comb2=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb2_usetr2{:});
      gravehi_TPA_MSPN_comb2=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb2_usetr2{:});
    else
      gravelo_TPA_MSPN_comb2=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb2{:});
      gravehi_TPA_MSPN_comb2=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb2{:});
    end
    if synchasynch && ll<5
      gravelo_TMSs_TMSa_comb2=ft_freqgrandaverage(cfg,freqloall_TMSs_TMSa_comb2{:});
      gravehi_TMSs_TMSa_comb2=ft_freqgrandaverage(cfg,freqhiall_TMSs_TMSa_comb2{:});
    end
  end
  if audtacflag
    gravelo_APT_MSPN_comb1=ft_freqgrandaverage(cfg,freqloall_APT_MSPN_comb1{:});
    gravehi_APT_MSPN_comb1=ft_freqgrandaverage(cfg,freqhiall_APT_MSPN_comb1{:});
    if comb2flag
      gravelo_APT_MSPN_comb2=ft_freqgrandaverage(cfg,freqloall_APT_MSPN_comb2{:});
      gravehi_APT_MSPN_comb2=ft_freqgrandaverage(cfg,freqhiall_APT_MSPN_comb2{:});
    end
  end
  
  if 0
    figure;
    cfg=[];
    cfg.layout='elec1010.lay';
    cfg.baselinetype='relchange';
    cfg.baseline=[-1.2 tacbasemax(ll)];
    cfg.zlim='maxabs';
    cfg.parameter='powspctrm';
    ft_multiplotTFR(cfg,gravehi_TPA_MSPN_comb1);
    ft_multiplotTFR(cfg,gravelo_TPA_MSPN_comb1);
    cfg.parameter='plvavgabs';
    ft_multiplotTFR(cfg,gravehi_TPA_MSPN_comb1);
    ft_multiplotTFR(cfg,gravelo_TPA_MSPN_comb1);
  end
  
  
  
  % assessing unisensory
  if 0
    if plotflag
      
      for ii=1:length(subuse)
        figure;
        cfg=[];
        cfg.layout='elec1010.lay';
        cfg.baselinetype='relchange';
        cfg.baseline=[-1.2 tacbasemax(ll)];
        cfg.zlim=[-1.5 1.5];
        cfg.ylim=[4 24];
        subplot(2,2,1);ft_multiplotTFR(cfg,freqloall_tNulAlone_comb{ii});
        subplot(2,2,2);ft_multiplotTFR(cfg,freqloall_tTacAlone_comb{ii});
        subplot(2,2,3);ft_multiplotTFR(cfg,freqloall_tAudAlone_comb{ii});
        subplot(2,2,4);ft_multiplotTFR(cfg,freqloall_tMSAlone_comb{ii});
      end
    end
  end
  
  
  
  % assessing combination of conditions
  if 0
    if plotflag
      figure(20);
      for ii=1:length(subuse)
        cfg=[];
        cfg.parameter='powspctrm';
        cfg.operation='subtract';
        diff{ii}=ft_math(cfg,freqloall_tacPaud_comb1{ii},freqloall_tacMSpN_comb1{ii});
        subplot(3,length(subuse),ii);imagesc(freqloall_tacPaud_comb1{ii}.time,freqloall_tacPaud_comb1{ii}.freq,squeeze(mean(freqloall_tacPaud_comb1{ii}.powspctrm,1)));caxis([-6 6])
        subplot(3,length(subuse),length(subuse)+ii);imagesc(freqloall_tacMSpN_comb1{ii}.time,freqloall_tacPaud_comb1{ii}.freq,squeeze(mean(freqloall_tacMSpN_comb1{ii}.powspctrm,1)));caxis([-6 6])
        subplot(3,length(subuse),2*length(subuse)+ii);imagesc(diff{ii}.time,freqloall_tacPaud_comb1{ii}.freq,squeeze(mean(diff{ii}.powspctrm,1)));caxis([-6 6])
      end
      if printflag
        print(20,[fdir 'sumuni_mspn_fl_diff_ta_cond' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
      end
      figure(22);
      for ii=1:length(subuse)
        cfg=[];
        cfg.parameter='powspctrm';
        cfg.operation='subtract';
        diff{ii}=ft_math(cfg,freqhiall_tacPaud_comb1{ii},freqhiall_tacMSpN_comb1{ii});
        subplot(3,length(subuse),ii);imagesc(freqhiall_tacPaud_comb1{ii}.time,freqhiall_tacPaud_comb1{ii}.freq,squeeze(mean(freqhiall_tacPaud_comb1{ii}.powspctrm,1)));caxis([-1 1])
        subplot(3,length(subuse),length(subuse)+ii);imagesc(freqhiall_tacMSpN_comb1{ii}.time,freqhiall_tacPaud_comb1{ii}.freq,squeeze(mean(freqhiall_tacMSpN_comb1{ii}.powspctrm,1)));caxis([-1 1])
        subplot(3,length(subuse),2*length(subuse)+ii);imagesc(diff{ii}.time,freqhiall_tacPaud_comb1{ii}.freq,squeeze(mean(diff{ii}.powspctrm,1)));caxis([-.2 .2])
      end
      if printflag
        print(22,[fdir 'sumuni_mspn_fh_diff_ta_cond' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
      end
      if audtacflag
        figure(21);
        for ii=1:length(subuse)
          cfg=[];
          cfg.parameter='powspctrm';
          cfg.operation='subtract';
          diff{ii}=ft_math(cfg,freqloall_audPtac_comb1{ii},freqloall_audMSpN_comb1{ii});
          subplot(3,length(subuse),ii);imagesc(freqloall_audPtac_comb1{ii}.time,freqloall_tacPaud_comb1{ii}.freq,squeeze(mean(freqloall_audPtac_comb1{ii}.powspctrm,1)));caxis([-6 6])
          subplot(3,length(subuse),length(subuse)+ii);imagesc(freqloall_audMSpN_comb1{ii}.time,freqloall_tacPaud_comb1{ii}.freq,squeeze(mean(freqloall_audMSpN_comb1{ii}.powspctrm,1)));caxis([-6 6])
          subplot(3,length(subuse),2*length(subuse)+ii);imagesc(diff{ii}.time,freqloall_tacPaud_comb1{ii}.freq,squeeze(mean(diff{ii}.powspctrm,1)));caxis([-6 6])
        end
        if printflag
          print(21,[fdir 'sumuni_mspn_fl_diff_at_cond' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
        end
        figure(23);
        for ii=1:length(subuse)
          cfg=[];
          cfg.parameter='powspctrm';
          cfg.operation='subtract';
          diff{ii}=ft_math(cfg,freqhiall_audPtac_comb1{ii},freqhiall_audMSpN_comb1{ii});
          subplot(3,length(subuse),ii);imagesc(freqhiall_audPtac_comb1{ii}.time,freqhiall_tacPaud_comb1{ii}.freq,squeeze(mean(freqhiall_audPtac_comb1{ii}.powspctrm,1)));caxis([-1 1])
          subplot(3,length(subuse),length(subuse)+ii);imagesc(freqhiall_audMSpN_comb1{ii}.time,freqhiall_tacPaud_comb1{ii}.freq,squeeze(mean(freqhiall_audMSpN_comb1{ii}.powspctrm,1)));caxis([-1 1])
          subplot(3,length(subuse),2*length(subuse)+ii);imagesc(diff{ii}.time,freqhiall_tacPaud_comb1{ii}.freq,squeeze(mean(diff{ii}.powspctrm,1)));caxis([-.2 .2])
        end
        if printflag
          print(23,[fdir 'sumuni_mspn_fh_diff_at_cond' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
        end
      end
    end
  end
  
  
  if fftaddflag
    for ii=1:length(subuse)
      
      cfg=[];
      cfg.baselinetype='relchange';
      cfg.baseline=[-1.3 tacbasemax(ll)];
      freqloall_tacPaud_fftaddbn{ii}=ft_freqbaseline(cfg,freqloall_tacPaud_fftadd{ii});
      freqloall_tacMSpN_fftaddbn{ii}=ft_freqbaseline(cfg,freqloall_tacMSpN_fftadd{ii});
      freqhiall_tacPaud_fftaddbn{ii}=ft_freqbaseline(cfg,freqhiall_tacPaud_fftadd{ii});
      freqhiall_tacMSpN_fftaddbn{ii}=ft_freqbaseline(cfg,freqhiall_tacMSpN_fftadd{ii});
      freqloall_tacPaud_comb1bn{ii}=ft_freqbaseline(cfg,freqloall_tacPaud_comb1{ii});
      freqloall_tacMSpN_comb1bn{ii}=ft_freqbaseline(cfg,freqloall_tacMSpN_comb1{ii});
      freqhiall_tacPaud_comb1bn{ii}=ft_freqbaseline(cfg,freqhiall_tacPaud_comb1{ii});
      freqhiall_tacMSpN_comb1bn{ii}=ft_freqbaseline(cfg,freqhiall_tacMSpN_comb1{ii});
      
      if audtacflag
        cfg.baseline=[-1.3 audbasemax(ll)];
        freqloall_audPtac_fftaddbn{ii}=ft_freqbaseline(cfg,freqloall_audPtac_fftadd{ii});
        freqloall_audMSpN_fftaddbn{ii}=ft_freqbaseline(cfg,freqloall_audMSpN_fftadd{ii});
        freqhiall_audPtac_fftaddbn{ii}=ft_freqbaseline(cfg,freqhiall_audPtac_fftadd{ii});
        freqhiall_audMSpN_fftaddbn{ii}=ft_freqbaseline(cfg,freqhiall_audMSpN_fftadd{ii});
        freqloall_audPtac_comb1bn{ii}=ft_freqbaseline(cfg,freqloall_audPtac_comb1{ii});
        freqloall_audMSpN_comb1bn{ii}=ft_freqbaseline(cfg,freqloall_audMSpN_comb1{ii});
        freqhiall_audPtac_comb1bn{ii}=ft_freqbaseline(cfg,freqhiall_audPtac_comb1{ii});
        freqhiall_audMSpN_comb1bn{ii}=ft_freqbaseline(cfg,freqhiall_audMSpN_comb1{ii});
      end
      
      %       baselinetime=[dsearchn(freqloall_tacPaud_fftadd{ii}.time',-1.2) dsearchn(freqloall_tacPaud_fftadd{ii}.time',-0.7)];
      %       freqloall_tacPaud_fabn{ii}=freqloall_tacPaud_fftadd{ii};
      %       freqloall_tacPaud_fabn{ii}.powspctrm=freqloall_tacPaud_fabn{ii}.powspctrm./repmat(nanmean(freqloall_tacPaud_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_fabn{ii}.time)]);
      %       freqloall_tacMSpN_fabn{ii}=freqloall_tacMSpN_fftadd{ii};
      %       freqloall_tacMSpN_fabn{ii}.powspctrm=freqloall_tacMSpN_fabn{ii}.powspctrm./repmat(nanmean(freqloall_tacMSpN_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_fabn{ii}.time)]);
      %       freqloall_audPtac_fabn{ii}=freqloall_audPtac_fftadd{ii};
      %       freqloall_audPtac_fabn{ii}.powspctrm=freqloall_audPtac_fabn{ii}.powspctrm./repmat(nanmean(freqloall_audPtac_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_fabn{ii}.time)]);
      %       freqloall_audMSpN_fabn{ii}=freqloall_audMSpN_fftadd{ii};
      %       freqloall_audMSpN_fabn{ii}.powspctrm=freqloall_audMSpN_fabn{ii}.powspctrm./repmat(nanmean(freqloall_audMSpN_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_fabn{ii}.time)]);
      %       freqhiall_tacPaud_fabn{ii}=freqhiall_tacPaud_fftadd{ii};
      %       freqhiall_tacPaud_fabn{ii}.powspctrm=freqhiall_tacPaud_fabn{ii}.powspctrm./repmat(nanmean(freqhiall_tacPaud_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_fabn{ii}.time)]);
      %       freqhiall_tacMSpN_fabn{ii}=freqhiall_tacMSpN_fftadd{ii};
      %       freqhiall_tacMSpN_fabn{ii}.powspctrm=freqhiall_tacMSpN_fabn{ii}.powspctrm./repmat(nanmean(freqhiall_tacMSpN_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_fabn{ii}.time)]);
      %       freqhiall_audPtac_fabn{ii}=freqhiall_audPtac_fftadd{ii};
      %       freqhiall_audPtac_fabn{ii}.powspctrm=freqhiall_audPtac_fabn{ii}.powspctrm./repmat(nanmean(freqhiall_audPtac_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_fabn{ii}.time)]);
      %       freqhiall_audMSpN_fabn{ii}=freqhiall_audMSpN_fftadd{ii};
      %       freqhiall_audMSpN_fabn{ii}.powspctrm=freqhiall_audMSpN_fabn{ii}.powspctrm./repmat(nanmean(freqhiall_audMSpN_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_fabn{ii}.time)]);
      %       baselinetime=[dsearchn(freqloall_tacPaud_comb1{ii}.time',-1.2) dsearchn(freqloall_tacPaud_comb1{ii}.time',-0.7)];
      %       freqloall_tacPaud_cbbn{ii}=freqloall_tacPaud_comb1{ii};
      %       freqloall_tacPaud_cbbn{ii}.powspctrm=freqloall_tacPaud_cbbn{ii}.powspctrm./repmat(nanmean(freqloall_tacPaud_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_cbbn{ii}.time)]);
      %       freqloall_tacMSpN_cbbn{ii}=freqloall_tacMSpN_comb1{ii};
      %       freqloall_tacMSpN_cbbn{ii}.powspctrm=freqloall_tacMSpN_cbbn{ii}.powspctrm./repmat(nanmean(freqloall_tacMSpN_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_cbbn{ii}.time)]);
      %       freqloall_audPtac_cbbn{ii}=freqloall_audPtac_comb1{ii};
      %       freqloall_audPtac_cbbn{ii}.powspctrm=freqloall_audPtac_cbbn{ii}.powspctrm./repmat(nanmean(freqloall_audPtac_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_cbbn{ii}.time)]);
      %       freqloall_audMSpN_cbbn{ii}=freqloall_audMSpN_comb1{ii};
      %       freqloall_audMSpN_cbbn{ii}.powspctrm=freqloall_audMSpN_cbbn{ii}.powspctrm./repmat(nanmean(freqloall_audMSpN_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_cbbn{ii}.time)]);
      %       freqhiall_tacPaud_cbbn{ii}=freqhiall_tacPaud_comb1{ii};
      %       freqhiall_tacPaud_cbbn{ii}.powspctrm=freqhiall_tacPaud_cbbn{ii}.powspctrm./repmat(nanmean(freqhiall_tacPaud_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_cbbn{ii}.time)]);
      %       freqhiall_tacMSpN_cbbn{ii}=freqhiall_tacMSpN_comb1{ii};
      %       freqhiall_tacMSpN_cbbn{ii}.powspctrm=freqhiall_tacMSpN_cbbn{ii}.powspctrm./repmat(nanmean(freqhiall_tacMSpN_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_cbbn{ii}.time)]);
      %       freqhiall_audPtac_cbbn{ii}=freqhiall_audPtac_comb1{ii};
      %       freqhiall_audPtac_cbbn{ii}.powspctrm=freqhiall_audPtac_cbbn{ii}.powspctrm./repmat(nanmean(freqhiall_audPtac_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_cbbn{ii}.time)]);
      %       freqhiall_audMSpN_cbbn{ii}=freqhiall_audMSpN_comb1{ii};
      %       freqhiall_audMSpN_cbbn{ii}.powspctrm=freqhiall_audMSpN_cbbn{ii}.powspctrm./repmat(nanmean(freqhiall_audMSpN_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_cbbn{ii}.time)]);
    end
  end
  
  cfg=[];
  cfg.keepindividual='yes';
  if fftaddflag
    grindlo_tacPaud_fftadd=ft_freqgrandaverage(cfg,freqloall_tacPaud_fftadd{:});
    grindlo_tacMSpN_fftadd=ft_freqgrandaverage(cfg,freqloall_tacMSpN_fftadd{:});
    grindhi_tacPaud_fftadd=ft_freqgrandaverage(cfg,freqhiall_tacPaud_fftadd{:});
    grindhi_tacMSpN_fftadd=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_fftadd{:});
    grindlo_tNulAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_tNulAlone_fftadd{:});
    grindhi_tNulAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_tNulAlone_fftadd{:});
    grindlo_tTacAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_tTacAlone_fftadd{:});
    grindhi_tTacAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_tTacAlone_fftadd{:});
    grindlo_tAudAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_tAudAlone_fftadd{:});
    grindhi_tAudAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_tAudAlone_fftadd{:});
    grindlo_tMSAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_tMSAlone_fftadd{:});
    grindhi_tMSAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_tMSAlone_fftadd{:});
    grindlo_tacPaud_fftaddbn=ft_freqgrandaverage(cfg,freqloall_tacPaud_fftaddbn{:});
    grindlo_tacMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqloall_tacMSpN_fftaddbn{:});
    grindhi_tacPaud_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_tacPaud_fftaddbn{:});
    grindhi_tacMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_fftaddbn{:});
    grindlo_tacPaud_comb1bn=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb1bn{:});
    grindlo_tacMSpN_comb1bn=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb1bn{:});
    grindhi_tacPaud_comb1bn=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb1bn{:});
    grindhi_tacMSpN_comb1bn=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb1bn{:});
    
    if audtacflag
      grindlo_audPtac_fftadd=ft_freqgrandaverage(cfg,freqloall_audPtac_fftadd{:});
      grindlo_audMSpN_fftadd=ft_freqgrandaverage(cfg,freqloall_audMSpN_fftadd{:});
      grindhi_audPtac_fftadd=ft_freqgrandaverage(cfg,freqhiall_audPtac_fftadd{:});
      grindhi_audMSpN_fftadd=ft_freqgrandaverage(cfg,freqhiall_audMSpN_fftadd{:});
      grindlo_audPtac_fftaddbn=ft_freqgrandaverage(cfg,freqloall_audPtac_fftaddbn{:});
      grindlo_audMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqloall_audMSpN_fftaddbn{:});
      grindhi_audPtac_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_audPtac_fftaddbn{:});
      grindhi_audMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_audMSpN_fftaddbn{:});
      grindlo_audPtac_comb1bn=ft_freqgrandaverage(cfg,freqloall_audPtac_comb1bn{:});
      grindlo_audMSpN_comb1bn=ft_freqgrandaverage(cfg,freqloall_audMSpN_comb1bn{:});
      grindhi_audPtac_comb1bn=ft_freqgrandaverage(cfg,freqhiall_audPtac_comb1bn{:});
      grindhi_audMSpN_comb1bn=ft_freqgrandaverage(cfg,freqhiall_audMSpN_comb1bn{:});
      grindlo_aNulAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_aNulAlone_fftadd{:});
      grindhi_aNulAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_aNulAlone_fftadd{:});
      grindlo_aTacAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_aTacAlone_fftadd{:});
      grindhi_aTacAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_aTacAlone_fftadd{:});
      grindlo_aAudAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_aAudAlone_fftadd{:});
      grindhi_aAudAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_aAudAlone_fftadd{:});
      grindlo_aMSAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_aMSAlone_fftadd{:});
      grindhi_aMSAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_aMSAlone_fftadd{:});
    end
  end
  
  cfg.parameter={'powspctrm' 'plvspctrm'};
  grindlo_tacPaud_comb1=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb1{:});
  grindlo_tacMSpN_comb1=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb1{:});
  grindhi_tacPaud_comb1=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb1{:});
  grindhi_tacMSpN_comb1=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb1{:});
  if synchasynch && ll<5
    grindlo_tMSsynch_comb1=ft_freqgrandaverage(cfg,freqloall_tMSsynch_comb1{:});
    grindlo_tMSasynch_comb1=ft_freqgrandaverage(cfg,freqloall_tMSasynch_comb1{:});
    grindhi_tMSsynch_comb1=ft_freqgrandaverage(cfg,freqhiall_tMSsynch_comb1{:});
    grindhi_tMSasynch_comb1=ft_freqgrandaverage(cfg,freqhiall_tMSasynch_comb1{:});
  end
  if comb2flag
    grindlo_tacPaud_comb2=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb2{:});
    grindlo_tacMSpN_comb2=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb2{:});
    grindhi_tacPaud_comb2=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb2{:});
    grindhi_tacMSpN_comb2=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb2{:});
    if synchasynch && ll<5
      grindlo_tMSsynch_comb2=ft_freqgrandaverage(cfg,freqloall_tMSsynch_comb2{:});
      grindlo_tMSasynch_comb2=ft_freqgrandaverage(cfg,freqloall_tMSasynch_comb2{:});
      grindhi_tMSsynch_comb2=ft_freqgrandaverage(cfg,freqhiall_tMSsynch_comb2{:});
      grindhi_tMSasynch_comb2=ft_freqgrandaverage(cfg,freqhiall_tMSasynch_comb2{:});
    end
  end
  grindlo_tNulAlone_comb=ft_freqgrandaverage(cfg,freqloall_tNulAlone_comb{:});
  grindhi_tNulAlone_comb=ft_freqgrandaverage(cfg,freqhiall_tNulAlone_comb{:});
  grindlo_tTacAlone_comb=ft_freqgrandaverage(cfg,freqloall_tTacAlone_comb{:});
  grindhi_tTacAlone_comb=ft_freqgrandaverage(cfg,freqhiall_tTacAlone_comb{:});
  grindlo_tAudAlone_comb=ft_freqgrandaverage(cfg,freqloall_tAudAlone_comb{:});
  grindhi_tAudAlone_comb=ft_freqgrandaverage(cfg,freqhiall_tAudAlone_comb{:});
  grindlo_tMSAlone_comb=ft_freqgrandaverage(cfg,freqloall_tMSAlone_comb{:});
  grindhi_tMSAlone_comb=ft_freqgrandaverage(cfg,freqhiall_tMSAlone_comb{:});
  if usetr==2
    grindlo_TPA_MSPN_comb1=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb1_usetr2{:});
    grindhi_TPA_MSPN_comb1=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb1_usetr2{:});
    if comb2flag
      grindlo_TPA_MSPN_comb2=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb2_usetr2{:});
      grindhi_TPA_MSPN_comb2=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb2_usetr2{:});
    end
  end
  
  if audtacflag
    grindlo_audPtac_comb1=ft_freqgrandaverage(cfg,freqloall_audPtac_comb1{:});
    grindlo_audMSpN_comb1=ft_freqgrandaverage(cfg,freqloall_audMSpN_comb1{:});
    grindhi_audPtac_comb1=ft_freqgrandaverage(cfg,freqhiall_audPtac_comb1{:});
    grindhi_audMSpN_comb1=ft_freqgrandaverage(cfg,freqhiall_audMSpN_comb1{:});
    if comb2flag
      grindlo_audPtac_comb2=ft_freqgrandaverage(cfg,freqloall_audPtac_comb2{:});
      grindlo_audMSpN_comb2=ft_freqgrandaverage(cfg,freqloall_audMSpN_comb2{:});
      grindhi_audPtac_comb2=ft_freqgrandaverage(cfg,freqhiall_audPtac_comb2{:});
      grindhi_audMSpN_comb2=ft_freqgrandaverage(cfg,freqhiall_audMSpN_comb2{:});
    end
    grindlo_aNulAlone_comb=ft_freqgrandaverage(cfg,freqloall_aNulAlone_comb{:});
    grindhi_aNulAlone_comb=ft_freqgrandaverage(cfg,freqhiall_aNulAlone_comb{:});
    grindlo_aTacAlone_comb=ft_freqgrandaverage(cfg,freqloall_aTacAlone_comb{:});
    grindhi_aTacAlone_comb=ft_freqgrandaverage(cfg,freqhiall_aTacAlone_comb{:});
    grindlo_aAudAlone_comb=ft_freqgrandaverage(cfg,freqloall_aAudAlone_comb{:});
    grindhi_aAudAlone_comb=ft_freqgrandaverage(cfg,freqhiall_aAudAlone_comb{:});
    grindlo_aMSAlone_comb=ft_freqgrandaverage(cfg,freqloall_aMSAlone_comb{:});
    grindhi_aMSAlone_comb=ft_freqgrandaverage(cfg,freqhiall_aMSAlone_comb{:});
  end
  
  
  cfg=[];
  if fftaddflag
    gravelo_tacPaud_fftadd=ft_freqgrandaverage(cfg,freqloall_tacPaud_fftadd{:});
    gravelo_tacMSpN_fftadd=ft_freqgrandaverage(cfg,freqloall_tacMSpN_fftadd{:});
    gravehi_tacPaud_fftadd=ft_freqgrandaverage(cfg,freqhiall_tacPaud_fftadd{:});
    gravehi_tacMSpN_fftadd=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_fftadd{:});
    gravelo_tNulAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_tNulAlone_fftadd{:});
    gravehi_tNulAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_tNulAlone_fftadd{:});
    gravelo_tTacAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_tTacAlone_fftadd{:});
    gravehi_tTacAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_tTacAlone_fftadd{:});
    gravelo_tAudAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_tAudAlone_fftadd{:});
    gravehi_tAudAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_tAudAlone_fftadd{:});
    gravelo_tMSAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_tMSAlone_fftadd{:});
    gravehi_tMSAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_tMSAlone_fftadd{:});
    gravelo_tacPaud_fftaddbn=ft_freqgrandaverage(cfg,freqloall_tacPaud_fftaddbn{:});
    gravelo_tacMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqloall_tacMSpN_fftaddbn{:});
    gravehi_tacPaud_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_tacPaud_fftaddbn{:});
    gravehi_tacMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_fftaddbn{:});
    gravelo_tacPaud_comb1bn=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb1bn{:});
    gravelo_tacMSpN_comb1bn=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb1bn{:});
    gravehi_tacPaud_comb1bn=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb1bn{:});
    gravehi_tacMSpN_comb1bn=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb1bn{:});
    
    if audtacflag
      gravelo_audPtac_fftadd=ft_freqgrandaverage(cfg,freqloall_audPtac_fftadd{:});
      gravelo_audMSpN_fftadd=ft_freqgrandaverage(cfg,freqloall_audMSpN_fftadd{:});
      gravehi_audPtac_fftadd=ft_freqgrandaverage(cfg,freqhiall_audPtac_fftadd{:});
      gravehi_audMSpN_fftadd=ft_freqgrandaverage(cfg,freqhiall_audMSpN_fftadd{:});
      gravelo_audPtac_fftaddbn=ft_freqgrandaverage(cfg,freqloall_audPtac_fftaddbn{:});
      gravelo_audMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqloall_audMSpN_fftaddbn{:});
      gravehi_audPtac_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_audPtac_fftaddbn{:});
      gravehi_audMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_audMSpN_fftaddbn{:});
      gravelo_audPtac_comb1bn=ft_freqgrandaverage(cfg,freqloall_audPtac_comb1bn{:});
      gravelo_audMSpN_comb1bn=ft_freqgrandaverage(cfg,freqloall_audMSpN_comb1bn{:});
      gravehi_audPtac_comb1bn=ft_freqgrandaverage(cfg,freqhiall_audPtac_comb1bn{:});
      gravehi_audMSpN_comb1bn=ft_freqgrandaverage(cfg,freqhiall_audMSpN_comb1bn{:});
      gravelo_aNulAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_aNulAlone_fftadd{:});
      gravehi_aNulAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_aNulAlone_fftadd{:});
      gravelo_aTacAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_aTacAlone_fftadd{:});
      gravehi_aTacAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_aTacAlone_fftadd{:});
      gravelo_aAudAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_aAudAlone_fftadd{:});
      gravehi_aAudAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_aAudAlone_fftadd{:});
      gravelo_aMSAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_aMSAlone_fftadd{:});
      gravehi_aMSAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_aMSAlone_fftadd{:});
    end
  end
  
  
  cfg.parameter={'powspctrm' 'plvspctrm'};
  gravelo_tacPaud_comb1=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb1{:});
  gravelo_tacMSpN_comb1=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb1{:});
  gravehi_tacPaud_comb1=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb1{:});
  gravehi_tacMSpN_comb1=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb1{:});
  if synchasynch && ll<5
    gravelo_tMSsynch_comb1=ft_freqgrandaverage(cfg,freqloall_tMSsynch_comb1{:});
    gravelo_tMSasynch_comb1=ft_freqgrandaverage(cfg,freqloall_tMSasynch_comb1{:});
    gravehi_tMSsynch_comb1=ft_freqgrandaverage(cfg,freqhiall_tMSsynch_comb1{:});
    gravehi_tMSasynch_comb1=ft_freqgrandaverage(cfg,freqhiall_tMSasynch_comb1{:});
  end
  if comb2flag
    gravelo_tacPaud_comb2=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb2{:});
    gravelo_tacMSpN_comb2=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb2{:});
    gravehi_tacPaud_comb2=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb2{:});
    gravehi_tacMSpN_comb2=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb2{:});
    if synchasynch && ll<5
      gravelo_tMSsynch_comb2=ft_freqgrandaverage(cfg,freqloall_tMSsynch_comb2{:});
      gravelo_tMSasynch_comb2=ft_freqgrandaverage(cfg,freqloall_tMSasynch_comb2{:});
      gravehi_tMSsynch_comb2=ft_freqgrandaverage(cfg,freqhiall_tMSsynch_comb2{:});
      gravehi_tMSasynch_comb2=ft_freqgrandaverage(cfg,freqhiall_tMSasynch_comb2{:});
    end
  end
  gravelo_tNulAlone_comb=ft_freqgrandaverage(cfg,freqloall_tNulAlone_comb{:});
  gravehi_tNulAlone_comb=ft_freqgrandaverage(cfg,freqhiall_tNulAlone_comb{:});
  gravelo_tTacAlone_comb=ft_freqgrandaverage(cfg,freqloall_tTacAlone_comb{:});
  gravehi_tTacAlone_comb=ft_freqgrandaverage(cfg,freqhiall_tTacAlone_comb{:});
  gravelo_tAudAlone_comb=ft_freqgrandaverage(cfg,freqloall_tAudAlone_comb{:});
  gravehi_tAudAlone_comb=ft_freqgrandaverage(cfg,freqhiall_tAudAlone_comb{:});
  gravelo_tMSAlone_comb=ft_freqgrandaverage(cfg,freqloall_tMSAlone_comb{:});
  gravehi_tMSAlone_comb=ft_freqgrandaverage(cfg,freqhiall_tMSAlone_comb{:});
  
  if audtacflag
    gravelo_audPtac_comb1=ft_freqgrandaverage(cfg,freqloall_audPtac_comb1{:});
    gravelo_audMSpN_comb1=ft_freqgrandaverage(cfg,freqloall_audMSpN_comb1{:});
    gravehi_audPtac_comb1=ft_freqgrandaverage(cfg,freqhiall_audPtac_comb1{:});
    gravehi_audMSpN_comb1=ft_freqgrandaverage(cfg,freqhiall_audMSpN_comb1{:});
    if comb2flag
      gravelo_audPtac_comb2=ft_freqgrandaverage(cfg,freqloall_audPtac_comb2{:});
      gravelo_audMSpN_comb2=ft_freqgrandaverage(cfg,freqloall_audMSpN_comb2{:});
      gravehi_audPtac_comb2=ft_freqgrandaverage(cfg,freqhiall_audPtac_comb2{:});
      gravehi_audMSpN_comb2=ft_freqgrandaverage(cfg,freqhiall_audMSpN_comb2{:});
    end
    gravelo_aNulAlone_comb=ft_freqgrandaverage(cfg,freqloall_aNulAlone_comb{:});
    gravehi_aNulAlone_comb=ft_freqgrandaverage(cfg,freqhiall_aNulAlone_comb{:});
    gravelo_aTacAlone_comb=ft_freqgrandaverage(cfg,freqloall_aTacAlone_comb{:});
    gravehi_aTacAlone_comb=ft_freqgrandaverage(cfg,freqhiall_aTacAlone_comb{:});
    gravelo_aAudAlone_comb=ft_freqgrandaverage(cfg,freqloall_aAudAlone_comb{:});
    gravehi_aAudAlone_comb=ft_freqgrandaverage(cfg,freqhiall_aAudAlone_comb{:});
    gravelo_aMSAlone_comb=ft_freqgrandaverage(cfg,freqloall_aMSAlone_comb{:});
    gravehi_aMSAlone_comb=ft_freqgrandaverage(cfg,freqhiall_aMSAlone_comb{:});
  end
  
  if plotflag
    powalonemaxz=[0.5 0.4 0.3 0.2 0.2]; % one per frequency band
    plvabsmaxz=[.6 .4 .2 .1 .1];
    % individual conditions first
    
    % fftadd first
    if 0 % this is redundant with comb for non-added-together conditions
      topoplotTFR_highlight(61,gravelo_tNulAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(62,gravelo_tTacAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(63,gravelo_tAudAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(64,gravelo_tMSAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(71,gravelo_tNulAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(72,gravelo_tTacAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(73,gravelo_tAudAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(74,gravelo_tMSAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(81,gravelo_tNulAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(82,gravelo_tTacAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(83,gravelo_tAudAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(84,gravelo_tMSAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(66,gravehi_tNulAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(67,gravehi_tTacAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(68,gravehi_tAudAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(69,gravehi_tMSAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(76,gravehi_tNulAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(77,gravehi_tTacAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(78,gravehi_tAudAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(79,gravehi_tMSAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05);
      if printflag
        ylim=ylimlo(1,:);
        print(61,[fdir 'gravelo_tNulAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(62,[fdir 'gravelo_tTacAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(63,[fdir 'gravelo_tAudAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(64,[fdir 'gravelo_tMSAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        ylim=ylimlo(2,:);
        print(71,[fdir 'gravelo_tNulAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(72,[fdir 'gravelo_tTacAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(73,[fdir 'gravelo_tAudAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(74,[fdir 'gravelo_tMSAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        ylim=ylimlo(3,:);
        print(81,[fdir 'gravelo_tNulAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(82,[fdir 'gravelo_tTacAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(83,[fdir 'gravelo_tAudAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(84,[fdir 'gravelo_tMSAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        ylim=ylimhi(1,:);
        print(66,[fdir 'gravehi_tNulAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(67,[fdir 'gravehi_tTacAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(68,[fdir 'gravehi_tAudAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(69,[fdir 'gravehi_tMSAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        ylim=ylimhi(2,:);
        print(76,[fdir 'gravehi_tNulAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(77,[fdir 'gravehi_tTacAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(78,[fdir 'gravehi_tAudAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(79,[fdir 'gravehi_tMSAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      end
      close all
    end
    
    % comb second
    topoplotTFR_highlight(61,gravelo_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(1) powalonemaxz(1)]);
    topoplotTFR_highlight(62,gravelo_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(1) powalonemaxz(1)]);
    topoplotTFR_highlight(63,gravelo_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(1) powalonemaxz(1)]);
    topoplotTFR_highlight(64,gravelo_tMSAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(1) powalonemaxz(1)]);
    topoplotTFR_highlight(71,gravelo_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(2) powalonemaxz(2)]);
    topoplotTFR_highlight(72,gravelo_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(2) powalonemaxz(2)]);
    topoplotTFR_highlight(73,gravelo_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(2) powalonemaxz(2)]);
    topoplotTFR_highlight(74,gravelo_tMSAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(2) powalonemaxz(2)]);
    topoplotTFR_highlight(81,gravelo_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(3) powalonemaxz(3)]);
    topoplotTFR_highlight(82,gravelo_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(3) powalonemaxz(3)]);
    topoplotTFR_highlight(83,gravelo_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(3) powalonemaxz(3)]);
    topoplotTFR_highlight(84,gravelo_tMSAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(3) powalonemaxz(3)]);
    topoplotTFR_highlight(66,gravehi_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(4) powalonemaxz(4)]);
    topoplotTFR_highlight(67,gravehi_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(4) powalonemaxz(4)]);
    topoplotTFR_highlight(68,gravehi_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(4) powalonemaxz(4)]);
    topoplotTFR_highlight(69,gravehi_tMSAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(4) powalonemaxz(4)]);
    topoplotTFR_highlight(76,gravehi_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(5) powalonemaxz(5)]);
    topoplotTFR_highlight(77,gravehi_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(5) powalonemaxz(5)]);
    topoplotTFR_highlight(78,gravehi_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(5) powalonemaxz(5)]);
    topoplotTFR_highlight(79,gravehi_tMSAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(5) powalonemaxz(5)]);
    if printflag
      ylim=ylimlo(1,:);
      print(61,[fdir 'gravelo_tNulAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(62,[fdir 'gravelo_tTacAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(63,[fdir 'gravelo_tAudAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(64,[fdir 'gravelo_tMSAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(2,:);
      print(71,[fdir 'gravelo_tNulAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(72,[fdir 'gravelo_tTacAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(73,[fdir 'gravelo_tAudAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(74,[fdir 'gravelo_tMSAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(3,:);
      print(81,[fdir 'gravelo_tNulAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(82,[fdir 'gravelo_tTacAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(83,[fdir 'gravelo_tAudAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(84,[fdir 'gravelo_tMSAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(1,:);
      print(66,[fdir 'gravehi_tNulAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(67,[fdir 'gravehi_tTacAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(68,[fdir 'gravehi_tAudAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(69,[fdir 'gravehi_tMSAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(2,:);
      print(76,[fdir 'gravehi_tNulAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(77,[fdir 'gravehi_tTacAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(78,[fdir 'gravehi_tAudAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(79,[fdir 'gravehi_tMSAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
    end
    close all
    
    % comb with PLV
    topoplotTFR_highlight(141,gravelo_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(142,gravelo_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(143,gravelo_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(144,gravelo_tMSAlone_comb ,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(151,gravelo_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(152,gravelo_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(153,gravelo_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(154,gravelo_tMSAlone_comb ,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(161,gravelo_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(162,gravelo_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(163,gravelo_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(164,gravelo_tMSAlone_comb ,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(146,gravehi_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(147,gravehi_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(148,gravehi_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(149,gravehi_tMSAlone_comb ,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(156,gravehi_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(5)]);
    topoplotTFR_highlight(157,gravehi_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(5)]);
    topoplotTFR_highlight(158,gravehi_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(5)]);
    topoplotTFR_highlight(159,gravehi_tMSAlone_comb ,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(5)]);
    
    topoplotTFR_highlight(1141,gravelo_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(1142,gravelo_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(1143,gravelo_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(1144,gravelo_tMSAlone_comb ,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(1151,gravelo_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(1152,gravelo_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(1153,gravelo_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(1154,gravelo_tMSAlone_comb ,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(1161,gravelo_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(1162,gravelo_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(1163,gravelo_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(1164,gravelo_tMSAlone_comb ,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(1146,gravehi_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(1147,gravehi_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(1148,gravehi_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(1149,gravehi_tMSAlone_comb ,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(1156,gravehi_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(5)]);
    topoplotTFR_highlight(1157,gravehi_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(5)]);
    topoplotTFR_highlight(1158,gravehi_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(5)]);
    topoplotTFR_highlight(1159,gravehi_tMSAlone_comb ,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(5)]);
    if printflag
      ylim=ylimlo(1,:);
      print(141,[fdir 'gravelo_tNulAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(142,[fdir 'gravelo_tTacAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(143,[fdir 'gravelo_tAudAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(144,[fdir 'gravelo_tMSAlone_comb_plvabsavg_topoOverTime_ta_cond'  num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(2,:);
      print(151,[fdir 'gravelo_tNulAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(152,[fdir 'gravelo_tTacAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(153,[fdir 'gravelo_tAudAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(154,[fdir 'gravelo_tMSAlone_comb_plvabsavg_topoOverTime_ta_cond'  num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(3,:);
      print(161,[fdir 'gravelo_tNulAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(162,[fdir 'gravelo_tTacAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(163,[fdir 'gravelo_tAudAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(164,[fdir 'gravelo_tMSAlone_comb_plvabsavg_topoOverTime_ta_cond'  num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(1,:);
      print(146,[fdir 'gravehi_tNulAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(147,[fdir 'gravehi_tTacAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(148,[fdir 'gravehi_tAudAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(149,[fdir 'gravehi_tMSAlone_comb_plvabsavg_topoOverTime_ta_cond'  num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(2,:);
      print(156,[fdir 'gravehi_tNulAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(157,[fdir 'gravehi_tTacAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(158,[fdir 'gravehi_tAudAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(159,[fdir 'gravehi_tMSAlone_comb_plvabsavg_topoOverTime_ta_cond'  num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      
      ylim=ylimlo(1,:);
      print(1141,[fdir 'gravelo_tNulAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1142,[fdir 'gravelo_tTacAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1143,[fdir 'gravelo_tAudAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1144,[fdir 'gravelo_tMSAlone_comb_plvavgabs_topoOverTime_ta_cond'  num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(2,:);
      print(1151,[fdir 'gravelo_tNulAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1152,[fdir 'gravelo_tTacAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1153,[fdir 'gravelo_tAudAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1154,[fdir 'gravelo_tMSAlone_comb_plvavgabs_topoOverTime_ta_cond'  num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(3,:);
      print(1161,[fdir 'gravelo_tNulAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1162,[fdir 'gravelo_tTacAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1163,[fdir 'gravelo_tAudAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1164,[fdir 'gravelo_tMSAlone_comb_plvavgabs_topoOverTime_ta_cond'  num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(1,:);
      print(1146,[fdir 'gravehi_tNulAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1147,[fdir 'gravehi_tTacAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1148,[fdir 'gravehi_tAudAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1149,[fdir 'gravehi_tMSAlone_comb_plvavgabs_topoOverTime_ta_cond'  num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(2,:);
      print(1156,[fdir 'gravehi_tNulAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1157,[fdir 'gravehi_tTacAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1158,[fdir 'gravehi_tAudAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1159,[fdir 'gravehi_tMSAlone_comb_plvavgabs_topoOverTime_ta_cond'  num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
    end
    close all
    %         keyboard % single condition alone plvavgang and plvangavg ??
    
    
    % contrast of conditions second
    % power first
    if fftaddflag
      topoplotTFR_highlight(110,gravelo_tacPaud_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(1) powalonemaxz(1)]);
      topoplotTFR_highlight(111,gravelo_tacMSpN_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(1) powalonemaxz(1)]);
      topoplotTFR_highlight(112,gravelo_tacPaud_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(2) powalonemaxz(2)]);
      topoplotTFR_highlight(113,gravelo_tacMSpN_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(2) powalonemaxz(2)]);
      topoplotTFR_highlight(114,gravelo_tacPaud_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(3) powalonemaxz(3)]);
      topoplotTFR_highlight(115,gravelo_tacMSpN_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(3) powalonemaxz(3)]);
      topoplotTFR_highlight(116,gravehi_tacPaud_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(4) powalonemaxz(4)]);
      topoplotTFR_highlight(117,gravehi_tacMSpN_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(4) powalonemaxz(4)]);
      topoplotTFR_highlight(118,gravehi_tacPaud_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(5) powalonemaxz(5)]);
      topoplotTFR_highlight(119,gravehi_tacMSpN_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(5) powalonemaxz(5)]);
    end
    
    topoplotTFR_highlight(120,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(1) powalonemaxz(1)]);
    topoplotTFR_highlight(121,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(1) powalonemaxz(1)]);
    topoplotTFR_highlight(122,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(2) powalonemaxz(2)]);
    topoplotTFR_highlight(123,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(2) powalonemaxz(2)]);
    topoplotTFR_highlight(124,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(3) powalonemaxz(3)]);
    topoplotTFR_highlight(125,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(3) powalonemaxz(3)]);
    topoplotTFR_highlight(126,gravehi_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(4) powalonemaxz(4)]);
    topoplotTFR_highlight(127,gravehi_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(4) powalonemaxz(4)]);
    topoplotTFR_highlight(128,gravehi_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(5) powalonemaxz(5)]);
    topoplotTFR_highlight(129,gravehi_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(5) powalonemaxz(5)]);
    
    if comb2flag
      topoplotTFR_highlight(130,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(1) powalonemaxz(1)]);
      topoplotTFR_highlight(131,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(1) powalonemaxz(1)]);
      topoplotTFR_highlight(132,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(2) powalonemaxz(2)]);
      topoplotTFR_highlight(133,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(2) powalonemaxz(2)]);
      topoplotTFR_highlight(134,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(3) powalonemaxz(3)]);
      topoplotTFR_highlight(135,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(3) powalonemaxz(3)]);
      topoplotTFR_highlight(136,gravehi_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(4) powalonemaxz(4)]);
      topoplotTFR_highlight(137,gravehi_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(4) powalonemaxz(4)]);
      topoplotTFR_highlight(138,gravehi_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(5) powalonemaxz(5)]);
      topoplotTFR_highlight(139,gravehi_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(5) powalonemaxz(5)]);
    end
    
    if printflag
      ylim=ylimlo(1,:);
      if fftaddflag
        print(110,[fdir 'gravelo_tacPaud_fftadd_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(111,[fdir 'gravelo_tacMSpN_fftadd_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      end
      print(120,[fdir 'gravelo_tacPaud_comb1_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(121,[fdir 'gravelo_tacMSpN_comb1_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(130,[fdir 'gravelo_tacPaud_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(131,[fdir 'gravelo_tacMSpN_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(2,:);
      if fftaddflag
        print(112,[fdir 'gravelo_tacPaud_fftadd_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(113,[fdir 'gravelo_tacMSpN_fftadd_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      end
      print(122,[fdir 'gravelo_tacPaud_comb1_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(123,[fdir 'gravelo_tacMSpN_comb1_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(132,[fdir 'gravelo_tacPaud_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(133,[fdir 'gravelo_tacMSpN_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(3,:);
      if fftaddflag
        print(114,[fdir 'gravelo_tacPaud_fftadd_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(115,[fdir 'gravelo_tacMSpN_fftadd_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      end
      print(124,[fdir 'gravelo_tacPaud_comb1_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(125,[fdir 'gravelo_tacMSpN_comb1_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(134,[fdir 'gravelo_tacPaud_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(135,[fdir 'gravelo_tacMSpN_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(1,:);
      if fftaddflag
        print(116,[fdir 'gravehi_tacPaud_fftadd_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(117,[fdir 'gravehi_tacMSpN_fftadd_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      end
      print(126,[fdir 'gravehi_tacPaud_comb1_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(127,[fdir 'gravehi_tacMSpN_comb1_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(136,[fdir 'gravehi_tacPaud_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(137,[fdir 'gravehi_tacMSpN_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(2,:);
      if fftaddflag
        print(118,[fdir 'gravehi_tacPaud_fftadd_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(119,[fdir 'gravehi_tacMSpN_fftadd_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      end
      print(128,[fdir 'gravehi_tacPaud_comb1_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(129,[fdir 'gravehi_tacMSpN_comb1_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(138,[fdir 'gravehi_tacPaud_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(139,[fdir 'gravehi_tacMSpN_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
    end
    close all
    
    % plv second
    topoplotTFR_highlight(220,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(221,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(222,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(223,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(224,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(225,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(226,gravehi_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(227,gravehi_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(228,gravehi_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(5)]);
    topoplotTFR_highlight(229,gravehi_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(5)]);
    topoplotTFR_highlight(1220,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(1221,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(1222,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(1223,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(1224,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(1225,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(1226,gravehi_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(1227,gravehi_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(1228,gravehi_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(5)]);
    topoplotTFR_highlight(1229,gravehi_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(5)]);
    
    topoplotTFR_highlight(250,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvangavg',[-pi pi]);
    topoplotTFR_highlight(250,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvangavg',[-pi pi]);
    topoplotTFR_highlight(251,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvangavg',[-pi pi]);
    topoplotTFR_highlight(252,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvangavg',[-pi pi]);
    topoplotTFR_highlight(253,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvangavg',[-pi pi]);
    topoplotTFR_highlight(254,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvangavg',[-pi pi]);
    topoplotTFR_highlight(255,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvangavg',[-pi pi]);
    topoplotTFR_highlight(256,gravehi_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvangavg',[-pi pi]);
    topoplotTFR_highlight(257,gravehi_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvangavg',[-pi pi]);
    topoplotTFR_highlight(258,gravehi_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvangavg',[-pi pi]);
    topoplotTFR_highlight(259,gravehi_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvangavg',[-pi pi]);
    topoplotTFR_highlight(1250,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgang',[-pi pi]);
    topoplotTFR_highlight(1251,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgang',[-pi pi]);
    topoplotTFR_highlight(1252,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgang',[-pi pi]);
    topoplotTFR_highlight(1253,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgang',[-pi pi]);
    topoplotTFR_highlight(1254,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgang',[-pi pi]);
    topoplotTFR_highlight(1255,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgang',[-pi pi]);
    topoplotTFR_highlight(1256,gravehi_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgang',[-pi pi]);
    topoplotTFR_highlight(1257,gravehi_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgang',[-pi pi]);
    topoplotTFR_highlight(1258,gravehi_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgang',[-pi pi]);
    topoplotTFR_highlight(1259,gravehi_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgang',[-pi pi]);
    
    if comb2flag
      topoplotTFR_highlight(230,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(1)]);
      topoplotTFR_highlight(231,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(1)]);
      topoplotTFR_highlight(232,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(2)]);
      topoplotTFR_highlight(233,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(2)]);
      topoplotTFR_highlight(234,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(3)]);
      topoplotTFR_highlight(235,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(3)]);
      topoplotTFR_highlight(236,gravehi_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(4)]);
      topoplotTFR_highlight(237,gravehi_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(4)]);
      topoplotTFR_highlight(238,gravehi_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(5)]);
      topoplotTFR_highlight(239,gravehi_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(5)]);
      topoplotTFR_highlight(1230,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(1)]);
      topoplotTFR_highlight(1231,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(1)]);
      topoplotTFR_highlight(1232,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(2)]);
      topoplotTFR_highlight(1233,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(2)]);
      topoplotTFR_highlight(1234,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(3)]);
      topoplotTFR_highlight(1235,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(3)]);
      topoplotTFR_highlight(1236,gravehi_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(4)]);
      topoplotTFR_highlight(1237,gravehi_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(4)]);
      topoplotTFR_highlight(1238,gravehi_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(5)]);
      topoplotTFR_highlight(1239,gravehi_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(5)]);
      
      topoplotTFR_highlight(260,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvangavg',[-pi pi]);
      topoplotTFR_highlight(261,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvangavg',[-pi pi]);
      topoplotTFR_highlight(262,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvangavg',[-pi pi]);
      topoplotTFR_highlight(263,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvangavg',[-pi pi]);
      topoplotTFR_highlight(264,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvangavg',[-pi pi]);
      topoplotTFR_highlight(265,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvangavg',[-pi pi]);
      topoplotTFR_highlight(266,gravehi_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvangavg',[-pi pi]);
      topoplotTFR_highlight(267,gravehi_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvangavg',[-pi pi]);
      topoplotTFR_highlight(268,gravehi_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvangavg',[-pi pi]);
      topoplotTFR_highlight(269,gravehi_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvangavg',[-pi pi]);
      topoplotTFR_highlight(1260,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgang',[-pi pi]);
      topoplotTFR_highlight(1261,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgang',[-pi pi]);
      topoplotTFR_highlight(1262,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgang',[-pi pi]);
      topoplotTFR_highlight(1263,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgang',[-pi pi]);
      topoplotTFR_highlight(1264,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgang',[-pi pi]);
      topoplotTFR_highlight(1265,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgang',[-pi pi]);
      topoplotTFR_highlight(1266,gravehi_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgang',[-pi pi]);
      topoplotTFR_highlight(1267,gravehi_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgang',[-pi pi]);
      topoplotTFR_highlight(1268,gravehi_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgang',[-pi pi]);
      topoplotTFR_highlight(1269,gravehi_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgang',[-pi pi]);
    end
    
    if printflag
      ylim=ylimlo(1,:);
      print(220,[fdir 'gravelo_tacPaud_comb1_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(221,[fdir 'gravelo_tacMSpN_comb1_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(250,[fdir 'gravelo_tacPaud_comb1_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(251,[fdir 'gravelo_tacMSpN_comb1_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(2,:);
      print(222,[fdir 'gravelo_tacPaud_comb1_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(223,[fdir 'gravelo_tacMSpN_comb1_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(252,[fdir 'gravelo_tacPaud_comb1_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(253,[fdir 'gravelo_tacMSpN_comb1_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(3,:);
      print(224,[fdir 'gravelo_tacPaud_comb1_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(225,[fdir 'gravelo_tacMSpN_comb1_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(254,[fdir 'gravelo_tacPaud_comb1_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(255,[fdir 'gravelo_tacMSpN_comb1_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(1,:);
      print(226,[fdir 'gravehi_tacPaud_comb1_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(227,[fdir 'gravehi_tacMSpN_comb1_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(256,[fdir 'gravehi_tacPaud_comb1_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(257,[fdir 'gravehi_tacMSpN_comb1_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(2,:);
      print(228,[fdir 'gravehi_tacPaud_comb1_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(229,[fdir 'gravehi_tacMSpN_comb1_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(258,[fdir 'gravehi_tacPaud_comb1_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(259,[fdir 'gravehi_tacMSpN_comb1_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      
      ylim=ylimlo(1,:);
      print(1220,[fdir 'gravelo_tacPaud_comb1_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1221,[fdir 'gravelo_tacMSpN_comb1_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1250,[fdir 'gravelo_tacPaud_comb1_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1251,[fdir 'gravelo_tacMSpN_comb1_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(2,:);
      print(1222,[fdir 'gravelo_tacPaud_comb1_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1223,[fdir 'gravelo_tacMSpN_comb1_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1252,[fdir 'gravelo_tacPaud_comb1_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1253,[fdir 'gravelo_tacMSpN_comb1_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(3,:);
      print(1224,[fdir 'gravelo_tacPaud_comb1_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1225,[fdir 'gravelo_tacMSpN_comb1_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1254,[fdir 'gravelo_tacPaud_comb1_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1255,[fdir 'gravelo_tacMSpN_comb1_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(1,:);
      print(1226,[fdir 'gravehi_tacPaud_comb1_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1227,[fdir 'gravehi_tacMSpN_comb1_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1256,[fdir 'gravehi_tacPaud_comb1_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1257,[fdir 'gravehi_tacMSpN_comb1_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(2,:);
      print(1228,[fdir 'gravehi_tacPaud_comb1_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1229,[fdir 'gravehi_tacMSpN_comb1_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1258,[fdir 'gravehi_tacPaud_comb1_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1259,[fdir 'gravehi_tacMSpN_comb1_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      
      if comb2flag
        ylim=ylimlo(1,:);
        print(230,[fdir 'gravelo_tacPaud_comb2_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(231,[fdir 'gravelo_tacMSpN_comb2_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(260,[fdir 'gravelo_tacPaud_comb2_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(261,[fdir 'gravelo_tacMSpN_comb2_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        ylim=ylimlo(2,:);
        print(232,[fdir 'gravelo_tacPaud_comb2_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(233,[fdir 'gravelo_tacMSpN_comb2_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(262,[fdir 'gravelo_tacPaud_comb2_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(263,[fdir 'gravelo_tacMSpN_comb2_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        ylim=ylimlo(3,:);
        print(234,[fdir 'gravelo_tacPaud_comb2_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(235,[fdir 'gravelo_tacMSpN_comb2_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(264,[fdir 'gravelo_tacPaud_comb2_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(265,[fdir 'gravelo_tacMSpN_comb2_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        ylim=ylimhi(1,:);
        print(236,[fdir 'gravehi_tacPaud_comb2_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(237,[fdir 'gravehi_tacMSpN_comb2_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(266,[fdir 'gravehi_tacPaud_comb2_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(267,[fdir 'gravehi_tacMSpN_comb2_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        ylim=ylimhi(2,:);
        print(238,[fdir 'gravehi_tacPaud_comb2_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(239,[fdir 'gravehi_tacMSpN_comb2_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(268,[fdir 'gravehi_tacPaud_comb2_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(269,[fdir 'gravehi_tacMSpN_comb2_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        
        ylim=ylimlo(1,:);
        print(1230,[fdir 'gravelo_tacPaud_comb2_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(1231,[fdir 'gravelo_tacMSpN_comb2_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(1260,[fdir 'gravelo_tacPaud_comb2_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(1261,[fdir 'gravelo_tacMSpN_comb2_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        ylim=ylimlo(2,:);
        print(1232,[fdir 'gravelo_tacPaud_comb2_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(1233,[fdir 'gravelo_tacMSpN_comb2_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(1262,[fdir 'gravelo_tacPaud_comb2_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(1263,[fdir 'gravelo_tacMSpN_comb2_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        ylim=ylimlo(3,:);
        print(1234,[fdir 'gravelo_tacPaud_comb2_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(1235,[fdir 'gravelo_tacMSpN_comb2_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(1264,[fdir 'gravelo_tacPaud_comb2_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(1265,[fdir 'gravelo_tacMSpN_comb2_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        ylim=ylimhi(1,:);
        print(1236,[fdir 'gravehi_tacPaud_comb2_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(1237,[fdir 'gravehi_tacMSpN_comb2_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(1266,[fdir 'gravehi_tacPaud_comb2_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(1267,[fdir 'gravehi_tacMSpN_comb2_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        ylim=ylimhi(2,:);
        print(1238,[fdir 'gravehi_tacPaud_comb2_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(1239,[fdir 'gravehi_tacMSpN_comb2_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(1268,[fdir 'gravehi_tacPaud_comb2_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(1269,[fdir 'gravehi_tacMSpN_comb2_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      end
    end
    close all
    
    
    
    %       if 0
    %       topoplotTFR_highlight(31,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       topoplotTFR_highlight(32,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       topoplotTFR_highlight(33,gravehi_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       topoplotTFR_highlight(34,gravehi_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       topoplotTFR_highlight(41,gravelo_tacPaud_fftaddbn,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       topoplotTFR_highlight(42,gravelo_tacMSpN_fftaddbn,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       topoplotTFR_highlight(43,gravehi_tacPaud_fftaddbn,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       topoplotTFR_highlight(44,gravehi_tacMSpN_fftaddbn,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       topoplotTFR_highlight(51,gravelo_tacPaud_comb1bn,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       topoplotTFR_highlight(52,gravelo_tacMSpN_comb1bn,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       topoplotTFR_highlight(53,gravehi_tacPaud_comb1bn,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       topoplotTFR_highlight(54,gravehi_tacMSpN_comb1bn,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       end
    
    if audtacflag
      topoplotTFR_highlight(15,gravelo_audPtac_fftadd,[audbasemax(ll); .55-tacbasemax(ll)],ylimlo(2,:),[-1.2 audbasemax(ll)],[],0.05);
      topoplotTFR_highlight(16,gravelo_audMSpN_fftadd,[audbasemax(ll); .55-tacbasemax(ll)],ylimlo(2,:),[-1.2 audbasemax(ll)],[],0.05);
      topoplotTFR_highlight(17,gravehi_audPtac_fftadd,[audbasemax(ll); .55-tacbasemax(ll)],ylimhi(2,:),[-1.2 audbasemax(ll)],[],0.05);
      topoplotTFR_highlight(18,gravehi_audMSpN_fftadd,[audbasemax(ll); .55-tacbasemax(ll)],ylimhi(2,:),[-1.2 audbasemax(ll)],[],0.05);
      
      topoplotTFR_highlight(25,gravelo_audPtac_comb1,[audbasemax(ll); .55-tacbasemax(ll)],ylimlo(2,:),[-1.2 audbasemax(ll)],[],0.05);
      topoplotTFR_highlight(26,gravelo_audMSpN_comb1,[audbasemax(ll); .55-tacbasemax(ll)],ylimlo(2,:),[-1.2 audbasemax(ll)],[],0.05);
      topoplotTFR_highlight(27,gravehi_audPtac_comb1,[audbasemax(ll); .55-tacbasemax(ll)],ylimhi(2,:),[-1.2 audbasemax(ll)],[],0.05);
      topoplotTFR_highlight(28,gravehi_audMSpN_comb1,[audbasemax(ll); .55-tacbasemax(ll)],ylimhi(2,:),[-1.2 audbasemax(ll)],[],0.05);
      
      if comb2flag
        topoplotTFR_highlight(35,gravelo_audPtac_comb2,[audbasemax(ll); .55-tacbasemax(ll)],ylimlo(2,:),[-1.2 audbasemax(ll)],[],0.05);
        topoplotTFR_highlight(36,gravelo_audMSpN_comb2,[audbasemax(ll); .55-tacbasemax(ll)],ylimlo(2,:),[-1.2 audbasemax(ll)],[],0.05);
        topoplotTFR_highlight(37,gravehi_audPtac_comb2,[audbasemax(ll); .55-tacbasemax(ll)],ylimhi(2,:),[-1.2 audbasemax(ll)],[],0.05);
        topoplotTFR_highlight(38,gravehi_audMSpN_comb2,[audbasemax(ll); .55-tacbasemax(ll)],ylimhi(2,:),[-1.2 audbasemax(ll)],[],0.05);
      end
      
      if 0
        topoplotTFR_highlight(45,gravelo_audPtac_fftaddbn,[audbasemax(ll); .55-tacbasemax(ll)],ylimlo(2,:),[-1.2 audbasemax(ll)],[],0.05);
        topoplotTFR_highlight(46,gravelo_audMSpN_fftaddbn,[audbasemax(ll); .55-tacbasemax(ll)],ylimlo(2,:),[-1.2 audbasemax(ll)],[],0.05);
        topoplotTFR_highlight(47,gravehi_audPtac_fftaddbn,[audbasemax(ll); .55-tacbasemax(ll)],ylimhi(2,:),[-1.2 audbasemax(ll)],[],0.05);
        topoplotTFR_highlight(48,gravehi_audMSpN_fftaddbn,[audbasemax(ll); .55-tacbasemax(ll)],ylimhi(2,:),[-1.2 audbasemax(ll)],[],0.05);
        
        topoplotTFR_highlight(55,gravelo_audPtac_comb1bn,[audbasemax(ll); .55-tacbasemax(ll)],ylimlo(2,:),[-1.2 audbasemax(ll)],[],0.05);
        topoplotTFR_highlight(56,gravelo_audMSpN_comb1bn,[audbasemax(ll); .55-tacbasemax(ll)],ylimlo(2,:),[-1.2 audbasemax(ll)],[],0.05);
        topoplotTFR_highlight(57,gravehi_audPtac_comb1bn,[audbasemax(ll); .55-tacbasemax(ll)],ylimhi(2,:),[-1.2 audbasemax(ll)],[],0.05);
        topoplotTFR_highlight(58,gravehi_audMSpN_comb1bn,[audbasemax(ll); .55-tacbasemax(ll)],ylimhi(2,:),[-1.2 audbasemax(ll)],[],0.05);
      end
    end
    
    
  end % plotflag
  
  if 0
    %       keyboard; % get multiplot going as well...want to see TFR and PLV for channel-cluster!
    figure;
    cfg=[];
    cfg.layout='elec1010.lay';
    cfg.baselinetype='relchange';
    cfg.baseline=[-1.2 tacbasemax(ll)];
    cfg.zlim='maxabs';
    cfg.ylim=[4 24];
    ft_multiplotTFR(cfg,gravelo_tacPaud_comb1);
    ft_multiplotTFR(cfg,gravelo_tacMSpN_comb1);
    cfg.parameter='plvavgabs';
    ft_multiplotTFR(cfg,gravelo_tacPaud_comb1);
    ft_multiplotTFR(cfg,gravelo_tacMSpN_comb1);
  end
  
  
  
  % % % Stats begins  % % %
  
  nsub=length(freqloall_tAudAlone_comb);
  
  if statsflag && nsub>1
  cfg=[];
  %       cfg.latency=[tacbasemax(ll); .55-audbasemax(ll)];
  cfg.latency=[-.15-audbasemax(ll); 1.05-audbasemax(ll)];
  cfg.neighbours=neighbours;
  cfg.parameter='powspctrm';
  cfg.method='montecarlo';
  cfg.numrandomization=2000;
  % cfg.correctm='holm';
  cfg.correctm='cluster';
  cfg.clusteralpha = 0.05;
  cfg.clusterstatistic = 'maxsum';
  cfg.minnbchan = 2;
  cfg.statistic='depsamplesT';
  cfg.design=zeros(2,2*nsub);
  cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
  cfg.design(2,:)=[1:nsub 1:nsub];
  cfg.ivar=1;
  cfg.uvar=2;
  cfg.randomseed=mcseed;
  if usetr==2
    grindlo_TPA_MSPN_comb1_zeros=grindlo_TPA_MSPN_comb1;
    grindlo_TPA_MSPN_comb1_zeros.powspctrm=zeros(size(grindlo_TPA_MSPN_comb1.powspctrm));
    grindlo_TPA_MSPN_comb1_zeros.plvspctrm=zeros(size(grindlo_TPA_MSPN_comb1.plvspctrm));
    grindhi_TPA_MSPN_comb1_zeros=grindhi_TPA_MSPN_comb1;
    grindhi_TPA_MSPN_comb1_zeros.powspctrm=zeros(size(grindhi_TPA_MSPN_comb1.powspctrm));
    grindhi_TPA_MSPN_comb1_zeros.plvspctrm=zeros(size(grindhi_TPA_MSPN_comb1.plvspctrm));
    stattl_mc_comb1=ft_freqstatistics(cfg, grindlo_TPA_MSPN_comb1, grindlo_TPA_MSPN_comb1_zeros);
    statth_mc_comb1=ft_freqstatistics(cfg, grindhi_TPA_MSPN_comb1, grindhi_TPA_MSPN_comb1_zeros);
    if comb2flag
      grindlo_TPA_MSPN_comb2_zeros=grindlo_TPA_MSPN_comb2;
      grindlo_TPA_MSPN_comb2_zeros.powspctrm=zeros(size(grindlo_TPA_MSPN_comb2.powspctrm));
      grindlo_TPA_MSPN_comb2_zeros.plvspctrm=zeros(size(grindlo_TPA_MSPN_comb2.plvspctrm));
      grindhi_TPA_MSPN_comb2_zeros=grindhi_TPA_MSPN_comb2;
      grindhi_TPA_MSPN_comb2_zeros.powspctrm=zeros(size(grindhi_TPA_MSPN_comb2.powspctrm));
      grindhi_TPA_MSPN_comb2_zeros.plvspctrm=zeros(size(grindhi_TPA_MSPN_comb2.plvspctrm));
      stattl_mc_comb2=ft_freqstatistics(cfg, grindlo_TPA_MSPN_comb2, grindlo_TPA_MSPN_comb2_zeros);
      statth_mc_comb2=ft_freqstatistics(cfg, grindhi_TPA_MSPN_comb2, grindhi_TPA_MSPN_comb2_zeros);
    end
  else
    if fftaddflag
      stattl_mc_fftadd=ft_freqstatistics(cfg, grindlo_tacPaud_fftadd, grindlo_tacMSpN_fftadd);
      statth_mc_fftadd=ft_freqstatistics(cfg, grindhi_tacPaud_fftadd, grindhi_tacMSpN_fftadd);
    end
    stattl_mc_comb1=ft_freqstatistics(cfg, grindlo_tacPaud_comb1, grindlo_tacMSpN_comb1);
    statth_mc_comb1=ft_freqstatistics(cfg, grindhi_tacPaud_comb1, grindhi_tacMSpN_comb1);
    if comb2flag
      stattl_mc_comb2=ft_freqstatistics(cfg, grindlo_tacPaud_comb2, grindlo_tacMSpN_comb2);
      statth_mc_comb2=ft_freqstatistics(cfg, grindhi_tacPaud_comb2, grindhi_tacMSpN_comb2);
    end
    %       stattl_mc_fftaddbn=ft_freqstatistics(cfg, grindlo_tacPaud_fftaddbn, grindlo_tacMSpN_fftaddbn);
    %       statth_mc_fftaddbn=ft_freqstatistics(cfg, grindhi_tacPaud_fftaddbn, grindhi_tacMSpN_fftaddbn);
    %       stattl_mc_comb1bn=ft_freqstatistics(cfg, grindlo_tacPaud_comb1bn, grindlo_tacMSpN_comb1bn);
    %       statth_mc_comb1bn=ft_freqstatistics(cfg, grindhi_tacPaud_comb1bn, grindhi_tacMSpN_comb1bn);
    if synchasynch && ll<5
      cfg.latency=[tacbasemax(9); .55-audbasemax(9)];
      stattl_synch_comb1=ft_freqstatistics(cfg, grindlo_tMSsynch_comb1, grindlo_tMSasynch_comb1);
      statth_synch_comb1=ft_freqstatistics(cfg, grindhi_tMSsynch_comb1, grindhi_tMSasynch_comb1);
      stattl_synch_comb2=ft_freqstatistics(cfg, grindlo_tMSsynch_comb2, grindlo_tMSasynch_comb2);
      statth_synch_comb2=ft_freqstatistics(cfg, grindhi_tMSsynch_comb2, grindhi_tMSasynch_comb2);
    end
    
    if audtacflag
      cfg.latency=[audbasemax(ll); .55-tacbasemax(ll)]; % ???
      statal_mc_fftadd=ft_freqstatistics(cfg, grindlo_audPtac_fftadd, grindlo_audMSpN_fftadd);
      statah_mc_fftadd=ft_freqstatistics(cfg, grindhi_audPtac_fftadd, grindhi_audMSpN_fftadd);
      statal_mc_comb1=ft_freqstatistics(cfg, grindlo_audPtac_comb1, grindlo_audMSpN_comb1);
      statah_mc_comb1=ft_freqstatistics(cfg, grindhi_audPtac_comb1, grindhi_audMSpN_comb1);
      if comb2flag
        statal_mc_comb2=ft_freqstatistics(cfg, grindlo_audPtac_comb2, grindlo_audMSpN_comb2);
        statah_mc_comb2=ft_freqstatistics(cfg, grindhi_audPtac_comb2, grindhi_audMSpN_comb2);
      end
      %       statal_mc_fftaddbn=ft_freqstatistics(cfg, grindlo_audPtac_fftaddbn, grindlo_audMSpN_fftaddbn);
      %       statah_mc_fftaddbn=ft_freqstatistics(cfg, grindhi_audPtac_fftaddbn, grindhi_audMSpN_fftaddbn);
      %       statal_mc_comb1bn=ft_freqstatistics(cfg, grindlo_audPtac_comb1bn, grindlo_audMSpN_comb1bn);
      %       statah_mc_comb1bn=ft_freqstatistics(cfg, grindhi_audPtac_comb1bn, grindhi_audMSpN_comb1bn);
    end
  end % usetr
  
  %       cfg.latency=[tacbasemax(ll); .55-audbasemax(ll)];
  cfg.latency=[-.15-audbasemax(ll); 1.05-audbasemax(ll)];
  cfg.statistic='diff_itc';
  cfg.parameter='plvspctrm';
  cfg.complex='diffabs'; % default; not sensitive to phase differences between conditions
  %       cfg.clusterthreshold = 'nonparametric_common'; % or 'nonparametric_individual' see clusterstat.m
  %       stattl_mc_comb1plvdac=ft_freqstatistics(cfg, grindlo_tacPaud_comb1, grindlo_tacMSpN_comb1);
  %       statth_mc_comb1plvdac=ft_freqstatistics(cfg, grindhi_tacPaud_comb1, grindhi_tacMSpN_comb1);
  %       if comb2flag
  %         stattl_mc_comb2plvdac=ft_freqstatistics(cfg, grindlo_tacPaud_comb2, grindlo_tacMSpN_comb2);
  %         statth_mc_comb2plvdac=ft_freqstatistics(cfg, grindhi_tacPaud_comb2, grindhi_tacMSpN_comb2);
  %       end
  cfg.clusterthreshold = 'nonparametric_individual'; % or 'nonparametric_individual' see clusterstat.m
  if usetr==2
    stattl_mc_comb1plvdai=ft_freqstatistics(cfg, grindlo_TPA_MSPN_comb1, grindlo_TPA_MSPN_comb1_zeros);
    statth_mc_comb1plvdai=ft_freqstatistics(cfg, grindhi_TPA_MSPN_comb1, grindhi_TPA_MSPN_comb1_zeros);
    if comb2flag
      stattl_mc_comb2plvdai=ft_freqstatistics(cfg, grindlo_TPA_MSPN_comb2, grindlo_TPA_MSPN_comb2_zeros);
      statth_mc_comb2plvdai=ft_freqstatistics(cfg, grindhi_TPA_MSPN_comb2, grindhi_TPA_MSPN_comb2_zeros);
    end
  else
  stattl_mc_comb1plvdai=ft_freqstatistics(cfg, grindlo_tacPaud_comb1, grindlo_tacMSpN_comb1);
  statth_mc_comb1plvdai=ft_freqstatistics(cfg, grindhi_tacPaud_comb1, grindhi_tacMSpN_comb1);
  if comb2flag
    stattl_mc_comb2plvdai=ft_freqstatistics(cfg, grindlo_tacPaud_comb2, grindlo_tacMSpN_comb2);
    statth_mc_comb2plvdai=ft_freqstatistics(cfg, grindhi_tacPaud_comb2, grindhi_tacMSpN_comb2);
  end
  if synchasynch && ll<5
    cfg.latency=[tacbasemax(9); .55-audbasemax(9)];
    stattl_synch_comb1plvdai=ft_freqstatistics(cfg, grindlo_tMSsynch_comb1, grindlo_tMSasynch_comb1);
    statth_synch_comb1plvdai=ft_freqstatistics(cfg, grindhi_tMSsynch_comb1, grindhi_tMSasynch_comb1);
    stattl_synch_comb2plvdai=ft_freqstatistics(cfg, grindlo_tMSsynch_comb2, grindlo_tMSasynch_comb2);
    statth_synch_comb2plvdai=ft_freqstatistics(cfg, grindhi_tMSsynch_comb2, grindhi_tMSasynch_comb2);
  end
  if audtacflag
    cfg.latency=[audbasemax(ll); .55-tacbasemax(ll)];
    %         cfg.clusterthreshold = 'nonparametric_common'; % or 'nonparametric_individual' see clusterstat.m
    %         statal_mc_comb1plvdac=ft_freqstatistics(cfg, grindlo_audPtac_comb1, grindlo_audMSpN_comb1);
    %         statah_mc_comb1plvdac=ft_freqstatistics(cfg, grindhi_audPtac_comb1, grindhi_audMSpN_comb1);
    %         if comb2flag
    %           statal_mc_comb2plvdac=ft_freqstatistics(cfg, grindlo_audPtac_comb2, grindlo_audMSpN_comb2);
    %           statah_mc_comb2plvdac=ft_freqstatistics(cfg, grindhi_audPtac_comb2, grindhi_audMSpN_comb2);
    %         end
    cfg.clusterthreshold = 'nonparametric_individual'; % or 'nonparametric_individual' see clusterstat.m
    statal_mc_comb1plvdai=ft_freqstatistics(cfg, grindlo_audPtac_comb1, grindlo_audMSpN_comb1);
    statah_mc_comb1plvdai=ft_freqstatistics(cfg, grindhi_audPtac_comb1, grindhi_audMSpN_comb1);
    if comb2flag
      statal_mc_comb2plvdai=ft_freqstatistics(cfg, grindlo_audPtac_comb2, grindlo_audMSpN_comb2);
      statah_mc_comb2plvdai=ft_freqstatistics(cfg, grindhi_audPtac_comb2, grindhi_audMSpN_comb2);
    end
  end
  
  %       cfg.latency=[tacbasemax(ll); .55-audbasemax(ll)];
  cfg.latency=[-.15-audbasemax(ll); 1.05-audbasemax(ll)];
  cfg.complex='absdiff'; % sensitive to phase difference between conditions
  %       cfg.clusterthreshold = 'nonparametric_common'; % or 'nonparametric_individual' see clusterstat.m
  %       stattl_mc_comb1plvadc=ft_freqstatistics(cfg, grindlo_tacPaud_comb1, grindlo_tacMSpN_comb1);
  %       statth_mc_comb1plvadc=ft_freqstatistics(cfg, grindhi_tacPaud_comb1, grindhi_tacMSpN_comb1);
  %       if comb2flag
  %         stattl_mc_comb2plvadc=ft_freqstatistics(cfg, grindlo_tacPaud_comb2, grindlo_tacMSpN_comb2);
  %         statth_mc_comb2plvadc=ft_freqstatistics(cfg, grindhi_tacPaud_comb2, grindhi_tacMSpN_comb2);
  %       end
  cfg.clusterthreshold = 'nonparametric_individual'; % or 'nonparametric_individual' see clusterstat.m
  stattl_mc_comb1plvadi=ft_freqstatistics(cfg, grindlo_tacPaud_comb1, grindlo_tacMSpN_comb1);
  statth_mc_comb1plvadi=ft_freqstatistics(cfg, grindhi_tacPaud_comb1, grindhi_tacMSpN_comb1);
  if comb2flag
    stattl_mc_comb2plvadi=ft_freqstatistics(cfg, grindlo_tacPaud_comb2, grindlo_tacMSpN_comb2);
    statth_mc_comb2plvadi=ft_freqstatistics(cfg, grindhi_tacPaud_comb2, grindhi_tacMSpN_comb2);
  end
  if synchasynch && ll<5
    cfg.latency=[tacbasemax(9); .55-audbasemax(9)];
    stattl_synch_comb1plvadi=ft_freqstatistics(cfg, grindlo_tMSsynch_comb1, grindlo_tMSasynch_comb1);
    statth_synch_comb1plvadi=ft_freqstatistics(cfg, grindhi_tMSsynch_comb1, grindhi_tMSasynch_comb1);
    stattl_synch_comb2plvadi=ft_freqstatistics(cfg, grindlo_tMSsynch_comb2, grindlo_tMSasynch_comb2);
    statth_synch_comb2plvadi=ft_freqstatistics(cfg, grindhi_tMSsynch_comb2, grindhi_tMSasynch_comb2);
  end
  
  if audtacflag
    cfg.latency=[audbasemax(ll); .55-tacbasemax(ll)];
    %         cfg.clusterthreshold = 'nonparametric_common'; % or 'nonparametric_individual' see clusterstat.m
    %         statal_mc_comb1plvadc=ft_freqstatistics(cfg, grindlo_audPtac_comb1, grindlo_audMSpN_comb1);
    %         statah_mc_comb1plvadc=ft_freqstatistics(cfg, grindhi_audPtac_comb1, grindhi_audMSpN_comb1);
    %         if comb2flag
    %           statal_mc_comb2plvadc=ft_freqstatistics(cfg, grindlo_audPtac_comb2, grindlo_audMSpN_comb2);
    %           statah_mc_comb2plvadc=ft_freqstatistics(cfg, grindhi_audPtac_comb2, grindhi_audMSpN_comb2);
    %         end
    cfg.clusterthreshold = 'nonparametric_individual'; % or 'nonparametric_individual' see clusterstat.m
    statal_mc_comb1plvadi=ft_freqstatistics(cfg, grindlo_audPtac_comb1, grindlo_audMSpN_comb1);
    statah_mc_comb1plvadi=ft_freqstatistics(cfg, grindhi_audPtac_comb1, grindhi_audMSpN_comb1);
    if comb2flag
      statal_mc_comb2plvadi=ft_freqstatistics(cfg, grindlo_audPtac_comb2, grindlo_audMSpN_comb2);
      statah_mc_comb2plvadi=ft_freqstatistics(cfg, grindhi_audPtac_comb2, grindhi_audMSpN_comb2);
    end
  end
  end % usetr
  end % statsflag
  
  if plotflag
    plvabsdiffmaxz=plvabsmaxz/2;
    timwin=[]; % not needed to narrow down, as stat will have already been preselected.
    base=[]; % stats weren't baseline corrected, so viewing of data shouldn't here be either.
    for yy=1:size(ylimlo,1)
      ylim=ylimlo(yy,:);
      
      try
        fig=240+yy;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb1,timwin,ylim,base,stattl_mc_comb1,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb1_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=840+yy;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb1,timwin,ylim,base,stattl_mc_comb1plvdai,.05,'plvavgabs',[0 plvabsdiffmaxz(yy)]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb1plvabsdai_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
        fig=1840+yy;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb1,timwin,ylim,base,stattl_mc_comb1plvdai,.05,'plvavgang',[-pi pi]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb1plvangdai_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=940+yy;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb1,timwin,ylim,base,stattl_mc_comb1plvadi,.05,'plvavgabs',[0 plvabsdiffmaxz(yy)]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb1plvabsadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
        fig=1940+yy;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb1,timwin,ylim,base,stattl_mc_comb1plvadi,.05,'plvavgang',[-pi pi]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb1plvangadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=340+yy;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb2,timwin,ylim,base,stattl_mc_comb2,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb2_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=2840+yy;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb2,timwin,ylim,base,stattl_mc_comb2plvdai,.05,'plvavgabs',[0 plvabsdiffmaxz(yy)]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb2plvabsdai_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
        fig=3840+yy;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb2,timwin,ylim,base,stattl_mc_comb2plvdai,.05,'plvavgang',[-pi pi]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb2plvangdai_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=2940+yy;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb2,timwin,ylim,base,stattl_mc_comb2plvadi,.05,'plvavgabs',[0 plvabsdiffmaxz(yy)]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb2plvabsadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
        fig=3940+yy;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb2,timwin,ylim,base,stattl_mc_comb2plvadi,.05,'plvavgang',[-pi pi]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb2plvangadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=1340+yy;
        topoplotTFR_highlight(fig,gravelo_TMSs_TMSa_comb1,timwin,ylim,base,stattl_synch_comb1,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_synch_ta_comb1_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=2340+yy;
        topoplotTFR_highlight(fig,gravelo_TMSs_TMSa_comb2,timwin,ylim,base,stattl_synch_comb2,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_synch_ta_comb2_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      
      
      if 0 % obsolete
        try
          fig=140+yy;
          topoplotTFR_highlight(fig,gravelo_TPA_MSPN_fftadd,timwin,ylim,base,stattl_mc_fftadd,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_fftadd_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=640+yy;
          topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb1,timwin,ylim,base,stattl_mc_comb1plvdac,.05,'plvavgabs',[0 plvabsdiffmaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb1plvdac_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=740+yy;
          topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb1,timwin,ylim,base,stattl_mc_comb1plvadc,.05,'plvavgabs',[0 plvabsdiffmaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb1plvabsadc_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
          fig=1740+yy;
          topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb1,timwin,ylim,base,stattl_mc_comb1plvadc,.05,'plvavgang',[-pi pi]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb1plvangadc_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=440+yy;
          topoplotTFR_highlight(fig,gravelo_TPA_MSPN_fftaddbn,timwin,ylim,base,stattl_mc_fftaddbn,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_fftaddbn_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=540+yy;
          topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb1bn,timwin,ylim,base,stattl_mc_comb1bn,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb1bn_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=160+yy;
          topoplotTFR_highlight(fig,gravelo_APT_MSPN_fftadd,timwin,ylim,base,statal_mc_fftadd,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_at_fftadd_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=260+yy;
          topoplotTFR_highlight(fig,gravelo_APT_MSPN_comb1,timwin,ylim,base,statal_mc_comb1,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_at_comb1_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=660+yy;
          topoplotTFR_highlight(fig,gravelo_APT_MSPN_comb1,timwin,ylim,base,statal_mc_comb1plvdac,.05,'plvavgabs',[0 plvabsdiffmaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_at_comb1plvdac_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=760+yy;
          topoplotTFR_highlight(fig,gravelo_APT_MSPN_comb1,timwin,ylim,base,statal_mc_comb1plvadc,.05,'plvavgabs',[0 plvabsdiffmaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_at_comb1plvabsadc_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
          fig=1760+yy;
          topoplotTFR_highlight(fig,gravelo_APT_MSPN_comb1,timwin,ylim,base,statal_mc_comb1plvadc,.05,'plvavgang',[-pi pi]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_at_comb1plvangadc_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=860+yy;
          topoplotTFR_highlight(fig,gravelo_APT_MSPN_comb1,timwin,ylim,base,statal_mc_comb1plvdai,.05,'plvavgabs',[0 plvabsdiffmaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_at_comb1plvdai_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=960+yy;
          topoplotTFR_highlight(fig,gravelo_APT_MSPN_comb1,timwin,ylim,base,statal_mc_comb1plvadi,.05,'plvavgabs',[0 plvabsdiffmaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_at_comb1plvabsadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
          fig=1960+yy;
          topoplotTFR_highlight(fig,gravelo_APT_MSPN_comb1,timwin,ylim,base,statal_mc_comb1plvadi,.05,'plvavgang',[-pi pi]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_at_comb1plvangadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=360+yy;
          topoplotTFR_highlight(fig,gravelo_APT_MSPN_comb2,timwin,ylim,base,statal_mc_comb2,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_at_comb2_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=460+yy;
          topoplotTFR_highlight(fig,gravelo_APT_MSPN_fftaddbn,timwin,ylim,base,statal_mc_fftaddbn,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_at_fftaddbn_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=560+yy;
          topoplotTFR_highlight(fig,gravelo_APT_MSPN_comb1bn,timwin,ylim,base,statal_mc_comb1bn,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_at_comb1bn_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
      end
      
    end
    
    for yy=1:size(ylimhi,1)
      ylim=ylimhi(yy,:);
      
      try
        fig=250+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb1,timwin,ylim,base,statth_mc_comb1,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb1_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=850+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb1,timwin,ylim,base,statth_mc_comb1plvdai,.05,'plvavgabs',[0 plvabsdiffmaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb1plvabsdai_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
        fig=1850+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb1,timwin,ylim,base,statth_mc_comb1plvdai,.05,'plvavgang',[-pi pi]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb1plvangdai_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=950+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb1,timwin,ylim,base,statth_mc_comb1plvadi,.05,'plvavgabs',[0 plvabsdiffmaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb1plvabsadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
        fig=1950+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb1,timwin,ylim,base,statth_mc_comb1plvadi,.05,'plvavgang',[-pi pi]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb1plvangadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=350+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb2,timwin,ylim,base,statth_mc_comb2,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb2_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=2850+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb2,timwin,ylim,base,statth_mc_comb2plvdai,.05,'plvavgabs',[0 plvabsdiffmaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb2plvabsdai_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
        fig=3850+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb2,timwin,ylim,base,statth_mc_comb2plvdai,.05,'plvavgang',[-pi pi]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb2plvangdai_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=2950+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb2,timwin,ylim,base,statth_mc_comb2plvadi,.05,'plvavgabs',[0 plvabsdiffmaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb2plvabsadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
        fig=3950+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb2,timwin,ylim,base,statth_mc_comb2plvadi,.05,'plvavgang',[-pi pi]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb2plvangadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=1350+yy;
        topoplotTFR_highlight(fig,gravehi_TMSs_TMSa_comb1,timwin,ylim,base,statth_synch_comb1,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_synch_ta_comb1_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=2350+yy;
        topoplotTFR_highlight(fig,gravehi_TMSs_TMSa_comb2,timwin,ylim,base,statth_synch_comb2,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_synch_ta_comb2_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      
      try
        fig=150+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_fftadd,timwin,ylim,base,statth_mc_fftadd,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_fftadd_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=650+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb1,timwin,ylim,base,statth_mc_comb1plvdac,.05,'plvavgabs',[0 plvabsdiffmaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb1plvdac_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=750+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb1,timwin,ylim,base,statth_mc_comb1plvadc,.05,'plvavgabs',[0 plvabsdiffmaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb1plvabsadc_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
        fig=1750+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb1,timwin,ylim,base,statth_mc_comb1plvadc,.05,'plvavgang',[-pi pi]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb1plvangadc_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=450+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_fftaddbn,timwin,ylim,base,statth_mc_fftaddbn,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_fftaddbn_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=550+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb1bn,timwin,ylim,base,statth_mc_comb1bn,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb1bn_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      
      if 0 % obsolete
        try
          fig=170+yy;
          topoplotTFR_highlight(fig,gravehi_APT_MSPN_fftadd,timwin,ylim,base,statah_mc_fftadd,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
          if printflag
            print(fig,[fdir 'grdiffhi_topoOverTime_mc_at_fftadd_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=270+yy;
          topoplotTFR_highlight(fig,gravehi_APT_MSPN_comb1,timwin,ylim,base,statah_mc_comb1,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
          if printflag
            print(fig,[fdir 'grdiffhi_topoOverTime_mc_at_comb1_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=670+yy;
          topoplotTFR_highlight(fig,gravehi_APT_MSPN_comb1,timwin,ylim,base,statah_mc_comb1plvdac,.05,'plvavgabs',[0 plvabsdiffmaxz(yy+3)]);
          if printflag
            print(fig,[fdir 'grdiffhi_topoOverTime_mc_at_comb1plvdac_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=770+yy;
          topoplotTFR_highlight(fig,gravehi_APT_MSPN_comb1,timwin,ylim,base,statah_mc_comb1plvadc,.05,'plvavgabs',[0 plvabsdiffmaxz(yy+3)]);
          if printflag
            print(fig,[fdir 'grdiffhi_topoOverTime_mc_at_comb1plvabsadc_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
          fig=1770+yy;
          topoplotTFR_highlight(fig,gravehi_APT_MSPN_comb1,timwin,ylim,base,statah_mc_comb1plvadc,.05,'plvavgang',[-pi pi]);
          if printflag
            print(fig,[fdir 'grdiffhi_topoOverTime_mc_at_comb1plvangadc_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=870+yy;
          topoplotTFR_highlight(fig,gravehi_APT_MSPN_comb1,timwin,ylim,base,statah_mc_comb1plvdai,.05,'plvavgabs',[0 plvabsdiffmaxz(yy+3)]);
          if printflag
            print(fig,[fdir 'grdiffhi_topoOverTime_mc_at_comb1plvdai_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=970+yy;
          topoplotTFR_highlight(fig,gravehi_APT_MSPN_comb1,timwin,ylim,base,statah_mc_comb1plvadi,.05,'plvavgabs',[0 plvabsdiffmaxz(yy+3)]);
          if printflag
            print(fig,[fdir 'grdiffhi_topoOverTime_mc_at_comb1plvabsadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
          fig=1970+yy;
          topoplotTFR_highlight(fig,gravehi_APT_MSPN_comb1,timwin,ylim,base,statah_mc_comb1plvadi,.05,'plvavgang',[-pi pi]);
          if printflag
            print(fig,[fdir 'grdiffhi_topoOverTime_mc_at_comb1plvangadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=370+yy;
          topoplotTFR_highlight(fig,gravehi_APT_MSPN_comb2,timwin,ylim,base,statah_mc_comb2,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
          if printflag
            print(fig,[fdir 'grdiffhi_topoOverTime_mc_at_comb2_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=470+yy;
          topoplotTFR_highlight(fig,gravehi_APT_MSPN_fftaddbn,timwin,ylim,base,statah_mc_fftaddbn,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
          if printflag
            print(fig,[fdir 'grdiffhi_topoOverTime_mc_at_fftaddbn_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=570+yy;
          topoplotTFR_highlight(fig,gravehi_APT_MSPN_comb1bn,timwin,ylim,base,statah_mc_comb1bn,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
          if printflag
            print(fig,[fdir 'grdiffhi_topoOverTime_mc_at_comb1bn_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
      end
      
    end % yy
  end % plotflag
  
  %   save([edir 'statsgrave_TFR_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '.mat'],'stat*','grave*');
%   save([edir 'statsgrave_TFR_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '_mcseed' num2str(mcseed) '.mat'],'stat*','grave*');
  save([edir 'statsgrave_TFR_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '_iter' num2str(iter) '_usetr' num2str(usetr) '_mcseed' num2str(mcseed) '.mat'],'stat*','grave*');
  
  clear stat*mc* gr*
  
  
  
  
  
  
end % ll
%   end % tt
% end % sleep

% keeping record of stats
tt=3;ss=10;sleep=0;usetr=3;mcseed=13;iter=31;
soalist=[1 3 4 5 6 7 9];
for ll=soalist
  load([edir 'statsgrave_TFR_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '_iter' num2str(iter) '_usetr' num2str(usetr) '_mcseed' num2str(mcseed) '.mat'],'stat*')
  comb_pvalues{ll,1,1,1}=unique(stattl_mc_comb1.prob(:));
  comb_pvalues{ll,2,1,1}=unique(stattl_mc_comb2.prob(:));
  comb_pvalues{ll,1,1,2}=unique(statth_mc_comb1.prob(:));
  comb_pvalues{ll,2,1,2}=unique(statth_mc_comb2.prob(:));
  comb_pvalues{ll,1,2,1}=unique(stattl_mc_comb1plvdai.prob(:));
  comb_pvalues{ll,2,2,1}=unique(stattl_mc_comb2plvdai.prob(:));
  comb_pvalues{ll,1,2,2}=unique(statth_mc_comb1plvdai.prob(:));
  comb_pvalues{ll,2,2,2}=unique(statth_mc_comb2plvdai.prob(:));
  clear stat*
end
save([edir 'comb_pvalues_TFR_' num2str(ss) num2str(sleep) '_iter' num2str(iter) '_usetr' num2str(usetr) '_mcseed' num2str(mcseed) '.mat'],'comb*')

tt=3;ss=10;sleep=0;
soalist=[1 3 4 5 6 7 9];
for ll=soalist
  load([edir 'statsgrave_TFR_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '.mat'],'stat*')
  comb_pvalues{ll,1,1,1}=unique(stattl_mc_comb1.prob(:));
  comb_pvalues{ll,2,1,1}=unique(stattl_mc_comb2.prob(:));
  comb_pvalues{ll,1,1,2}=unique(statth_mc_comb1.prob(:));
  comb_pvalues{ll,2,1,2}=unique(statth_mc_comb2.prob(:));
  comb_pvalues{ll,1,2,1}=unique(stattl_mc_comb1plvdai.prob(:));
  comb_pvalues{ll,2,2,1}=unique(stattl_mc_comb2plvdai.prob(:));
  comb_pvalues{ll,1,2,2}=unique(statth_mc_comb1plvdai.prob(:));
  comb_pvalues{ll,2,2,2}=unique(statth_mc_comb2plvdai.prob(:));
  clear stat*
end
save([edir 'comb_pvalues_TFR_' num2str(ss) num2str(sleep) '.mat'],'comb*')


%% Extra stats for comb2
% allcond_sameN=1;
soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];
tacbasemax=min(-.15,soades-.15);
audbasemax=min(-.15,fliplr(soades)-.15);

ylimlo=[4 7; 8 12; 14 30];
ylimhi=[30 45; 55 80];

soadesc={'Aud first by 500ms' '' 'Aud first by 70ms' 'Aud first by 20ms' 'Simultaneous' 'Tac first by 20ms' 'Tac first by 70ms' '' 'Tac first by 500ms'};
plotflag=1;
printflag=1;
statsflag=1;
audtacflag=0;
comb2flag=1;

soalist=[1 3 4 5 6 7 9];
% soalist=[3 4 5 6 7 9];



for sleep=[0]
  if sleep
    chanuse=chanuse_sleep1;
  else
    chanuse=chanuse_sleep0;
  end
  for tt=[3]
    figind=1;
    for ll=soalist
      %     for ll=[4 5]
      clearvars -except ll tt sub *dir ii*use sleep *flag figind soadesc soalist chanuse* ylim* neigh* *basemax
      
      if sleep
        subuseall=iiBuse;
      else
        subuseall=setdiff(iiSuse,[])
      end
      
      submin=subuseall(1)-1;
      subuseind=0;
      
      % Baseline correct each participant prior to entering to stats???? NO
      for ii=subuseall
        cd([edir sub{ii} ])
        %       load(['freq_diffs_averef_' sub{ii} '.mat']);
        try
          load(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_tt' num2str(tt) '.mat'])
          if audtacflag
            load(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(0) '_tt' num2str(tt) '.mat'])
          end
        catch
          if tt==2
            load(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '.mat'])
            if audtacflag
              load(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(0) '.mat'])
            end
          end
        end
        if audtacflag
          tka=load(['trialkeptTFR_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(0) '.mat']);
        end
        tkt=load(['trialkeptTFR_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '.mat']);
        
        
        if 1
          % THis is preliminary...how best to included all stages later on?
          if sleep==0
            ssuse=10; % awake
            sleepcond='Awake W';
          elseif sleep==1
            ssuse=12; % this is concatenation of N2 and N3
            sleepcond='Sleep N2';
            %             ssuse=23; % this is concatenation of N2 and N3
            %             sleepcond='Sleep (N2+N3)';
          end
        else
          if sleep
            ssuse=tkt.tr.stageuse;
          else
            ssuse=tka.tr.stageuse;
          end
        end
        
        ss=ssuse;
        %         for ss=ssuse
        subuse=subuseall; % reset to all for each sleep stage
        numtrt(ll,tt,ss,ii-submin)=numt_trials(ll,tt,ss);
        if audtacflag
          numtra(ll,tt,ss,ii-submin)=numa_trials(ll,tt,ss);
        end
        
        
        if min(numtrt(ll,tt,ss,ii-submin),[],1)<20 % what is best number to use here?
          subuse=setdiff(subuse,ii);
        else
          subuseind=subuseind+1;
          freqloall_tacPaud_comb2{subuseind}=freqlo_tacPaud_comb{ll,tt,ss,2};
          freqloall_tacMSpN_comb2{subuseind}=freqlo_tacMSpN_comb{ll,tt,ss,2};
          freqhiall_tacPaud_comb2{subuseind}=freqhi_tacPaud_comb{ll,tt,ss,2};
          freqhiall_tacMSpN_comb2{subuseind}=freqhi_tacMSpN_comb{ll,tt,ss,2};
          freqloall_tacPaud_comb2{subuseind}.dimord='chan_freq_time';
          freqloall_tacMSpN_comb2{subuseind}.dimord='chan_freq_time';
          freqhiall_tacPaud_comb2{subuseind}.dimord='chan_freq_time';
          freqhiall_tacMSpN_comb2{subuseind}.dimord='chan_freq_time';
          
          if audtacflag
            freqloall_audPtac_comb2{subuseind}=freqlo_audPtac_comb{ll,tt,ss,2};
            freqloall_audMSpN_comb2{subuseind}=freqlo_audMSpN_comb{ll,tt,ss,2};
            freqhiall_audPtac_comb2{subuseind}=freqhi_audPtac_comb{ll,tt,ss,2};
            freqhiall_audMSpN_comb2{subuseind}=freqhi_audMSpN_comb{ll,tt,ss,2};
            freqloall_audPtac_comb2{subuseind}.dimord='chan_freq_time';
            freqloall_audMSpN_comb2{subuseind}.dimord='chan_freq_time';
            freqhiall_audPtac_comb2{subuseind}.dimord='chan_freq_time';
            freqhiall_audMSpN_comb2{subuseind}.dimord='chan_freq_time';
          end
        end
        %         end % ss
        clear freqlo_* freqhi_*
      end % ii
      
      cfg=[];
      cfg.keepindividual='yes';
      cfg.parameter={'powspctrm' 'plvspctrm','stdpow','stdplv'};
      grindlo_tacPaud_comb2=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb2{:});
      grindlo_tacMSpN_comb2=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb2{:});
      grindhi_tacPaud_comb2=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb2{:});
      grindhi_tacMSpN_comb2=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb2{:});
      
      if audtacflag
        grindlo_audPtac_comb2=ft_freqgrandaverage(cfg,freqloall_audPtac_comb2{:});
        grindlo_audMSpN_comb2=ft_freqgrandaverage(cfg,freqloall_audMSpN_comb2{:});
        grindhi_audPtac_comb2=ft_freqgrandaverage(cfg,freqhiall_audPtac_comb2{:});
        grindhi_audMSpN_comb2=ft_freqgrandaverage(cfg,freqhiall_audMSpN_comb2{:});
      end
      
      
      cfg=[];
      cfg.parameter={'powspctrm' 'plvspctrm','stdpow','stdplv'};
      gravelo_tacPaud_comb2=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb2{:});
      gravelo_tacMSpN_comb2=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb2{:});
      gravehi_tacPaud_comb2=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb2{:});
      gravehi_tacMSpN_comb2=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb2{:});
      
      if audtacflag
        gravelo_audPtac_comb2=ft_freqgrandaverage(cfg,freqloall_audPtac_comb2{:});
        gravelo_audMSpN_comb2=ft_freqgrandaverage(cfg,freqloall_audMSpN_comb2{:});
        gravehi_audPtac_comb2=ft_freqgrandaverage(cfg,freqhiall_audPtac_comb2{:});
        gravehi_audMSpN_comb2=ft_freqgrandaverage(cfg,freqhiall_audMSpN_comb2{:});
      end
      
      
      if plotflag
        powalonemaxz=[0.5 0.4 0.3 0.2 0.2]; % one per frequency band
        plvabsmaxz=[.6 .4 .2 .1 .1];
        % contrast of conditions second
        % power first
        topoplotTFR_highlight(130,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(1) powalonemaxz(1)]);
        topoplotTFR_highlight(131,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(1) powalonemaxz(1)]);
        topoplotTFR_highlight(132,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(2) powalonemaxz(2)]);
        topoplotTFR_highlight(133,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(2) powalonemaxz(2)]);
        topoplotTFR_highlight(134,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(3) powalonemaxz(3)]);
        topoplotTFR_highlight(135,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(3) powalonemaxz(3)]);
        topoplotTFR_highlight(136,gravehi_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(4) powalonemaxz(4)]);
        topoplotTFR_highlight(137,gravehi_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(4) powalonemaxz(4)]);
        topoplotTFR_highlight(138,gravehi_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(5) powalonemaxz(5)]);
        topoplotTFR_highlight(139,gravehi_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(5) powalonemaxz(5)]);
        
        if printflag
          ylim=ylimlo(1,:);
          print(130,[fdir 'gravelo_tacPaud_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(131,[fdir 'gravelo_tacMSpN_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          ylim=ylimlo(2,:);
          print(132,[fdir 'gravelo_tacPaud_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(133,[fdir 'gravelo_tacMSpN_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          ylim=ylimlo(3,:);
          print(134,[fdir 'gravelo_tacPaud_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(135,[fdir 'gravelo_tacMSpN_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          ylim=ylimhi(1,:);
          print(136,[fdir 'gravehi_tacPaud_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(137,[fdir 'gravehi_tacMSpN_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          ylim=ylimhi(2,:);
          print(138,[fdir 'gravehi_tacPaud_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(139,[fdir 'gravehi_tacMSpN_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
        close all
        
        % plv second
        topoplotTFR_highlight(230,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(1)]);
        topoplotTFR_highlight(231,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(1)]);
        topoplotTFR_highlight(232,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(2)]);
        topoplotTFR_highlight(233,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(2)]);
        topoplotTFR_highlight(234,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(3)]);
        topoplotTFR_highlight(235,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(3)]);
        topoplotTFR_highlight(236,gravehi_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(4)]);
        topoplotTFR_highlight(237,gravehi_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(4)]);
        topoplotTFR_highlight(238,gravehi_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(5)]);
        topoplotTFR_highlight(239,gravehi_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(5)]);
        topoplotTFR_highlight(1230,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(1)]);
        topoplotTFR_highlight(1231,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(1)]);
        topoplotTFR_highlight(1232,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(2)]);
        topoplotTFR_highlight(1233,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(2)]);
        topoplotTFR_highlight(1234,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(3)]);
        topoplotTFR_highlight(1235,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(3)]);
        topoplotTFR_highlight(1236,gravehi_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(4)]);
        topoplotTFR_highlight(1237,gravehi_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(4)]);
        topoplotTFR_highlight(1238,gravehi_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(5)]);
        topoplotTFR_highlight(1239,gravehi_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(5)]);
        
        topoplotTFR_highlight(260,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvangavg',[-pi pi]);
        topoplotTFR_highlight(261,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvangavg',[-pi pi]);
        topoplotTFR_highlight(262,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvangavg',[-pi pi]);
        topoplotTFR_highlight(263,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvangavg',[-pi pi]);
        topoplotTFR_highlight(264,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvangavg',[-pi pi]);
        topoplotTFR_highlight(265,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvangavg',[-pi pi]);
        topoplotTFR_highlight(266,gravehi_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvangavg',[-pi pi]);
        topoplotTFR_highlight(267,gravehi_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvangavg',[-pi pi]);
        topoplotTFR_highlight(268,gravehi_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvangavg',[-pi pi]);
        topoplotTFR_highlight(269,gravehi_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvangavg',[-pi pi]);
        topoplotTFR_highlight(1260,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgang',[-pi pi]);
        topoplotTFR_highlight(1261,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgang',[-pi pi]);
        topoplotTFR_highlight(1262,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgang',[-pi pi]);
        topoplotTFR_highlight(1263,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgang',[-pi pi]);
        topoplotTFR_highlight(1264,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgang',[-pi pi]);
        topoplotTFR_highlight(1265,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgang',[-pi pi]);
        topoplotTFR_highlight(1266,gravehi_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgang',[-pi pi]);
        topoplotTFR_highlight(1267,gravehi_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgang',[-pi pi]);
        topoplotTFR_highlight(1268,gravehi_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgang',[-pi pi]);
        topoplotTFR_highlight(1269,gravehi_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgang',[-pi pi]);
        
        if printflag
          ylim=ylimlo(1,:);
          print(230,[fdir 'gravelo_tacPaud_comb2_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(231,[fdir 'gravelo_tacMSpN_comb2_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(260,[fdir 'gravelo_tacPaud_comb2_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(261,[fdir 'gravelo_tacMSpN_comb2_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          ylim=ylimlo(2,:);
          print(232,[fdir 'gravelo_tacPaud_comb2_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(233,[fdir 'gravelo_tacMSpN_comb2_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(262,[fdir 'gravelo_tacPaud_comb2_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(263,[fdir 'gravelo_tacMSpN_comb2_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          ylim=ylimlo(3,:);
          print(234,[fdir 'gravelo_tacPaud_comb2_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(235,[fdir 'gravelo_tacMSpN_comb2_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(264,[fdir 'gravelo_tacPaud_comb2_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(265,[fdir 'gravelo_tacMSpN_comb2_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          ylim=ylimhi(1,:);
          print(236,[fdir 'gravehi_tacPaud_comb2_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(237,[fdir 'gravehi_tacMSpN_comb2_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(266,[fdir 'gravehi_tacPaud_comb2_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(267,[fdir 'gravehi_tacMSpN_comb2_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          ylim=ylimhi(2,:);
          print(238,[fdir 'gravehi_tacPaud_comb2_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(239,[fdir 'gravehi_tacMSpN_comb2_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(268,[fdir 'gravehi_tacPaud_comb2_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(269,[fdir 'gravehi_tacMSpN_comb2_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          
          ylim=ylimlo(1,:);
          print(1230,[fdir 'gravelo_tacPaud_comb2_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(1231,[fdir 'gravelo_tacMSpN_comb2_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(1260,[fdir 'gravelo_tacPaud_comb2_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(1261,[fdir 'gravelo_tacMSpN_comb2_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          ylim=ylimlo(2,:);
          print(1232,[fdir 'gravelo_tacPaud_comb2_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(1233,[fdir 'gravelo_tacMSpN_comb2_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(1262,[fdir 'gravelo_tacPaud_comb2_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(1263,[fdir 'gravelo_tacMSpN_comb2_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          ylim=ylimlo(3,:);
          print(1234,[fdir 'gravelo_tacPaud_comb2_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(1235,[fdir 'gravelo_tacMSpN_comb2_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(1264,[fdir 'gravelo_tacPaud_comb2_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(1265,[fdir 'gravelo_tacMSpN_comb2_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          ylim=ylimhi(1,:);
          print(1236,[fdir 'gravehi_tacPaud_comb2_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(1237,[fdir 'gravehi_tacMSpN_comb2_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(1266,[fdir 'gravehi_tacPaud_comb2_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(1267,[fdir 'gravehi_tacMSpN_comb2_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          ylim=ylimhi(2,:);
          print(1238,[fdir 'gravehi_tacPaud_comb2_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(1239,[fdir 'gravehi_tacMSpN_comb2_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(1268,[fdir 'gravehi_tacPaud_comb2_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          print(1269,[fdir 'gravehi_tacMSpN_comb2_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
        close all
        
        
        
        if audtacflag
          topoplotTFR_highlight(35,gravelo_audPtac_comb2,[audbasemax(ll); .55-tacbasemax(ll)],ylimlo(2,:),[-1.2 audbasemax(ll)],[],0.05);
          topoplotTFR_highlight(36,gravelo_audMSpN_comb2,[audbasemax(ll); .55-tacbasemax(ll)],ylimlo(2,:),[-1.2 audbasemax(ll)],[],0.05);
          topoplotTFR_highlight(37,gravehi_audPtac_comb2,[audbasemax(ll); .55-tacbasemax(ll)],ylimhi(2,:),[-1.2 audbasemax(ll)],[],0.05);
          topoplotTFR_highlight(38,gravehi_audMSpN_comb2,[audbasemax(ll); .55-tacbasemax(ll)],ylimhi(2,:),[-1.2 audbasemax(ll)],[],0.05);
        end
      end % plotflag
      
      for ii=1:length(subuse)
        cfg=[];
        cfg.operation='subtract';
        cfg.parameter='powspctrm';
        cfg.parameter={'powspctrm' 'plvspctrm'};
        freqloall_TPA_MSPN_comb2{ii}=ft_math(cfg,freqloall_tacPaud_comb2{ii},freqloall_tacMSpN_comb2{ii});
        freqhiall_TPA_MSPN_comb2{ii}=ft_math(cfg,freqhiall_tacPaud_comb2{ii},freqhiall_tacMSpN_comb2{ii});
        
        if audtacflag
          freqloall_APT_MSPN_comb2{ii}=ft_math(cfg,freqloall_audPtac_comb2{ii},freqloall_audMSpN_comb2{ii});
          freqhiall_APT_MSPN_comb2{ii}=ft_math(cfg,freqhiall_audPtac_comb2{ii},freqhiall_audMSpN_comb2{ii});
        end
      end % end ii
      
      cfg=[];
      cfg.parameter={'powspctrm' 'plvspctrm'};
      gravelo_TPA_MSPN_comb2=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb2{:});
      gravehi_TPA_MSPN_comb2=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb2{:});
      if audtacflag
        gravelo_APT_MSPN_comb2=ft_freqgrandaverage(cfg,freqloall_APT_MSPN_comb2{:});
        gravehi_APT_MSPN_comb2=ft_freqgrandaverage(cfg,freqhiall_APT_MSPN_comb2{:});
      end
      
      % Get t-value per person, channel, time-frequency point
      % wikipedia.org/wiki/Student%27s_t-test#Equal_or_unequal_sample_sizes.2C_unequal_variances
      tvalue_lopowTPAMSPN{ll,tt,ss}=nut_ttest_meanStd(grindlo_tacPaud_comb2.powspctrm,grindlo_tacMSpN_comb2.powspctrm,grindlo_tacPaud_comb2.stdpow,grindlo_tacMSpN_comb2.stdpow,tkt.numcondtfinal,tkt.numcondtfinal);
      tvalue_loplvTPAMSPN{ll,tt,ss}=nut_ttest_meanStd(grindlo_tacPaud_comb2.plvspctrm,grindlo_tacMSpN_comb2.plvspctrm,grindlo_tacPaud_comb2.stdplv,grindlo_tacMSpN_comb2.stdplv,tkt.numcondtfinal,tkt.numcondtfinal);
      tvalue_hipowTPAMSPN{ll,tt,ss}=nut_ttest_meanStd(grindhi_tacPaud_comb2.powspctrm,grindhi_tacMSpN_comb2.powspctrm,grindhi_tacPaud_comb2.stdpow,grindhi_tacMSpN_comb2.stdpow,tkt.numcondtfinal,tkt.numcondtfinal);
      tvalue_hiplvTPAMSPN{ll,tt,ss}=nut_ttest_meanStd(grindhi_tacPaud_comb2.plvspctrm,grindhi_tacMSpN_comb2.plvspctrm,grindhi_tacPaud_comb2.stdplv,grindhi_tacMSpN_comb2.stdplv,tkt.numcondtfinal,tkt.numcondtfinal);
      % cut off t-values based on N=20(minimum) and p<0.01, means tcut=2.845
      % cut off t-values based on N=20(minimum) and p<0.05, means tcut=2.086
      tcut=2.086;
      tvmask_lopowTPAMSPN{ll,tt,ss}=tvalue_lopowTPAMSPN{ll,tt,ss}>tcut | tvalue_lopowTPAMSPN{ll,tt,ss}<-tcut;
      tvmask_loplvTPAMSPN{ll,tt,ss}=tvalue_loplvTPAMSPN{ll,tt,ss}>tcut | tvalue_loplvTPAMSPN{ll,tt,ss}<-tcut;
      tvmask_hipowTPAMSPN{ll,tt,ss}=tvalue_hipowTPAMSPN{ll,tt,ss}>tcut | tvalue_hipowTPAMSPN{ll,tt,ss}<-tcut;
      tvmask_hiplvTPAMSPN{ll,tt,ss}=tvalue_hiplvTPAMSPN{ll,tt,ss}>tcut | tvalue_hiplvTPAMSPN{ll,tt,ss}<-tcut;
      
      % Also try KL divergence between Gaussian distributions
      % KL(p,q) = log(std2/std1) + [std1.^2 + (mu1-nu2).^2]./[2*std2.^2] -1/2
      kldiv_lopowTPAMSPN{ll,tt,ss}=KLdiv_between_gaussians(grindlo_tacPaud_comb2.powspctrm,grindlo_tacMSpN_comb2.powspctrm,grindlo_tacPaud_comb2.stdpow,grindlo_tacMSpN_comb2.stdpow)
      kldiv_loplvTPAMSPN{ll,tt,ss}=KLdiv_between_gaussians(grindlo_tacPaud_comb2.plvspctrm,grindlo_tacMSpN_comb2.plvspctrm,grindlo_tacPaud_comb2.stdplv,grindlo_tacMSpN_comb2.stdplv)
      kldiv_hipowTPAMSPN{ll,tt,ss}=KLdiv_between_gaussians(grindhi_tacPaud_comb2.powspctrm,grindhi_tacMSpN_comb2.powspctrm,grindhi_tacPaud_comb2.stdpow,grindhi_tacMSpN_comb2.stdpow)
      kldiv_hiplvTPAMSPN{ll,tt,ss}=KLdiv_between_gaussians(grindhi_tacPaud_comb2.plvspctrm,grindhi_tacMSpN_comb2.plvspctrm,grindhi_tacPaud_comb2.stdplv,grindhi_tacMSpN_comb2.stdplv)
      
      % what is significantly different KLdiv?  2 std-dev away with unit-variance (0,2,1,1) gives KLdiv=2 (since .5*meandiff^2 for unitvariance)
      for cc=1:size(kldiv_lopowTPAMSPN{ll,tt,ss},2)
        for ff=1:size(kldiv_lopowTPAMSPN{ll,tt,ss},3)
          for ttime=1:size(kldiv_lopowTPAMSPN{ll,tt,ss},4)
            klmask_lowpowTPAMSPN{ll,tt,ss}(cc,ff,ttime)=signrank(kldiv_lopowTPAMSPN{ll,tt,ss}(:,cc,ff,ttime),tcut^2/2,'tail','right');
            klmask_hiwpowTPAMSPN{ll,tt,ss}(cc,ff,ttime)=signrank(kldiv_hipowTPAMSPN{ll,tt,ss}(:,cc,ff,ttime),tcut^2/2,'tail','right');
            klmask_lowplvTPAMSPN{ll,tt,ss}(cc,ff,ttime)=signrank(kldiv_loplvTPAMSPN{ll,tt,ss}(:,cc,ff,ttime),tcut^2/2,'tail','right');
            klmask_hiwplvTPAMSPN{ll,tt,ss}(cc,ff,ttime)=signrank(kldiv_hiplvTPAMSPN{ll,tt,ss}(:,cc,ff,ttime),tcut^2/2,'tail','right');
          end
        end
      end
      
      save([edir 'gravecomb2_TFR_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '.mat'],'grave*');
      clear gr*
      
    end % ll
  end % tt
  save([edir 'statcomb2_TFR_sleep' num2str(sleep) '.mat'],'tv*');
  clear tv*
end % sleep


%%  Plotting results


chanplot{1}={'Fz' 'Cz' 'F1' 'F2' 'FC1' 'FC2' 'C1' 'C2'}; % frontocentral
% chanplot{2}={'CP5' 'P5' 'PO3' 'POz' 'PO4' 'P6' 'CP6' 'Pz' 'P3' 'P4'}; % occipital
chanplot{2}={'CP5' 'POz' 'Pz' 'P3' 'P4' 'C4' 'O1' 'O2' 'P7' 'PO7'}; % occipital; % same as for ERP
% chanplot{3}={'FC6' 'FC4' 'F4' 'F6' 'F8' 'FT8' 'T8' 'C6' 'C4'}; % right frontotemporal
chanlabel{1}='Frontocentral electrodes';
chanlabel{2}='Occipital-parietal electrodes';
% chanlabel{3}='Right frontotemporal electrodes';
soalist=[1 3 4 5 6 7 9];
tt=3;
fftaddflag=0;
synchasynch=0;
tacbasemax=min(-.15,soades-.15);
% load([edir 'statsgrave_TFR_cond' num2str(ll) num2str(tt) num2str(ss) '.mat']);
stdir='/home/zumerj/audtac/eeg_data/statsgrave_oldstattimwin/';

sleep=0;
if sleep
  ss=12;
else
  ss=10;
end



% for all ll, print at least a TFR even if nothing significant.
% for ll=soalist
for ll=[5 7 9]
  close all
  clear grave* stat*
  load([edir 'statsgrave_TFR_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '.mat']);
  %   load([stdir 'statsgrave_TFR_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '.mat']);
  
  for combval=1
    close all
    clear adda
    % power
    if fftaddflag
      cfg=[];
      cfg.latency=[stattl_mc_fftadd.time(1) stattl_mc_fftadd.time(end)];
      tmp=ft_selectdata(cfg,gravelo_TPA_MSPN_fftadd);
      tmp.mask=stattl_mc_fftadd.mask;
      tmps=ft_selectdata(cfg,gravelo_tacMSpN_fftadd);
      tmps.mask=stattl_mc_fftadd.mask;
      tmpu=ft_selectdata(cfg,gravelo_tacPaud_fftadd);
      tmpu.mask=stattl_mc_fftadd.mask;
    else
      if combval==1
        cfg=[];
        cfg.latency=[stattl_mc_comb1.time(1) stattl_mc_comb1.time(end)];
        tmp=ft_selectdata(cfg,gravelo_TPA_MSPN_comb1);
        tmp.mask=stattl_mc_comb1.mask;
        tmps=ft_selectdata(cfg,gravelo_tacMSpN_comb1);
        tmps.mask=stattl_mc_comb1.mask;
        tmpu=ft_selectdata(cfg,gravelo_tacPaud_comb1);
        tmpu.mask=stattl_mc_comb1.mask;
        
        cfg=[];
        tmpA=ft_selectdata(cfg,gravelo_TPA_MSPN_comb1);
        tmpA.mask=logical(ones(size(tmpA.powspctrm)));
        tmpA.mask(:,:,dsearchn(tmpA.time',stattl_mc_comb1.time(1)):dsearchn(tmpA.time',stattl_mc_comb1.time(end)))=stattl_mc_comb1.mask;
        tmpsA=ft_selectdata(cfg,gravelo_tacMSpN_comb1);
        tmpsA.mask=logical(ones(size(tmpsA.powspctrm)));
        tmpsA.mask(:,:,dsearchn(tmpsA.time',stattl_mc_comb1.time(1)):dsearchn(tmpsA.time',stattl_mc_comb1.time(end)))=stattl_mc_comb1.mask;
        tmpuA=ft_selectdata(cfg,gravelo_tacPaud_comb1);
        tmpuA.mask=logical(ones(size(tmpuA.powspctrm)));
        tmpuA.mask(:,:,dsearchn(tmpuA.time',stattl_mc_comb1.time(1)):dsearchn(tmpuA.time',stattl_mc_comb1.time(end)))=stattl_mc_comb1.mask;
      elseif combval==2
        cfg=[];
        cfg.latency=[stattl_mc_comb2.time(1) stattl_mc_comb2.time(end)];
        tmp=ft_selectdata(cfg,gravelo_TPA_MSPN_comb2);
        tmp.mask=stattl_mc_comb2.mask;
        tmps=ft_selectdata(cfg,gravelo_tacMSpN_comb2);
        tmps.mask=stattl_mc_comb2.mask;
        tmpu=ft_selectdata(cfg,gravelo_tacPaud_comb2);
        tmpu.mask=stattl_mc_comb2.mask;
      end
      %         cfg=[];
      %         tmpn=ft_selectdata(cfg,gravelo_tNulAlone_comb );
      %         tmpt=ft_selectdata(cfg,gravelo_tTacAlone_comb );
      %         tmpa=ft_selectdata(cfg,gravelo_tAudAlone_comb );
      %         tmpm=ft_selectdata(cfg,gravelo_tMSAlone_comb );
      gravelo_tNulAlone_comb.mask=logical(ones(size(gravelo_tNulAlone_comb.powspctrm)));
      gravelo_tTacAlone_comb.mask=logical(ones(size(gravelo_tTacAlone_comb.powspctrm)));
      gravelo_tAudAlone_comb.mask=logical(ones(size(gravelo_tAudAlone_comb.powspctrm)));
      gravelo_tMSAlone_comb.mask= logical(ones(size(gravelo_tMSAlone_comb.powspctrm)));
    end
    
    %     if 0
    %       cfg=[];
    %       cfg.avgoverchan='yes';
    %       if ~isempty(tmp.label(find(ceil(mean(mean(tmp.mask(:,:,:),2),3)))))
    %         cfg.channel=tmp.label(find(ceil(mean(mean(tmp.mask(:,:,:),2),3))));
    %       else
    %         cfg.channel='all';
    %       end
    %       tmp1=ft_selectdata(cfg,tmp);
    %       tmp1.mask=logical(ceil(tmp1.mask));
    %
    %       figure(10*ll);
    %       cfg=[];
    %       cfg.parameter='powspctrm';
    %       cfg.layout='elec1010.lay';
    %       cfg.maskparameter='mask';
    %       cfg.maskalpha=0.5;
    %       cfg.zlim='maxabs';
    %       ft_singleplotTFR(cfg,tmp1);
    %       print(10*ll,[fdir 'tfrlo_final_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    %     end
    
    
    % plot TFR of powspctrm, for each suplot: u and m and difference.
    zlim=[-3 3];
    baseline1=[tacbasemax(ll) tacbasemax(ll)+.08];
    
    if 0
      pow{1}=gravelo_tMSAlone_comb;
      pow{2}=tmpu;
      pow{3}=tmps;
      pow{4}=tmp;
      
      chansel=chanplot{1};
      figinds=10*ll;
      figstrings=[];
      figstrings{1}=[fdir 'tfrlo_crop_FC_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
      tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline1)
      
      chansel=chanplot{2};
      figinds=10*ll+1;
      figstrings=[];
      figstrings{1}=[fdir 'tfrlo_crop_OP_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
      tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline1)
      
      chansel=[chanplot{1} chanplot{2}];
      figinds=10*ll+2;
      figstrings=[];
      figstrings{1}=[fdir 'tfrlo_crop_FCOP_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
      tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline1)
    end
    
    pow{1}=gravelo_tMSAlone_comb;
    pow{2}=tmpuA;
    pow{3}=tmpsA;
    pow{4}=tmpA;
    baseline2=[tacbasemax(1) tacbasemax(1)+.08];
    
    chansel=chanplot{1};
    figinds=10*ll;
    figstrings=[];
    figstrings{1}=[fdir 'tfrlo_all_FC_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
    tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline2)
    %     tfr_subchannel_3cond_plotRows_pow(pow,chansel,zlim,figinds,figstrings,baseline2)
    
    chansel=chanplot{2};
    figinds=10*ll+1;
    figstrings=[];
    figstrings{1}=[fdir 'tfrlo_all_OP_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
    tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline2)
    
    chansel=[chanplot{1} chanplot{2}];
    figinds=10*ll+2;
    figstrings=[];
    figstrings{1}=[fdir 'tfrlo_all_FCOP_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
    tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline2)
    
    masktime=find(squeeze(any(mean(tmp.mask(:,1:end,:),2),1)));
    if ~isempty(masktime)
      chansel=tmp.label(find(ceil(mean(mean(tmp.mask(:,:,dsearchn(tmp.time',tmp.time(masktime(1))):dsearchn(tmp.time',tmp.time(masktime(end)))),2),3))));
      figinds=10*ll+9;
      figstrings=[];
      figstrings{1}=[fdir 'tfrlo_all_allsig_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
      tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline2)
    end
    
    %     % tried this out but doesn't look right; baselining seems to make sense for display prior to contrast
    %     zlim=[0 10];
    %     baseline='no';
    
    
    % old way
    %     chansel=chanplot{1};
    %     figinds=10*ll;
    %     figstrings=[];
    %     figstrings{1}=[fdir 'tfrlo_final_FC_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.png'];
    %     baseline=[tmp.time(1) tmp.time(9)];
    %     tfr_subchannel_3cond_plot_pow(tmp,tmpu,tmps,chansel,zlim,figinds,figstrings,baseline)
    %
    %     chansel=chanplot{2};
    %     figinds=10*ll+1;
    %     figstrings=[];
    %     figstrings{1}=[fdir 'tfrlo_final_OP_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.png'];
    %     baseline=[tmp.time(1) tmp.time(9)];
    %     tfr_subchannel_3cond_plot_pow(tmp,tmpu,tmps,chansel,zlim,figinds,figstrings,baseline)
    
    
    % theta topo
    if ~isempty(find(squeeze(any(mean(tmp.mask(:,1:2,:),2),1))))
      cfg=[];
      cfg.parameter='powspctrm';
      cfg.layout='elec1010.lay';
      cfg.maskalpha=0.5;
      cfg.ylim=[4 6.5];
      cfg.zlim=[-1.4 1.4];
      cfg.highlight='on';
      cfg.comment='no';
      %       masktime=find(squeeze(any(mean(tmp.mask(:,1:2,:),2),1)));
      %       cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
      masktime_tmp=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1)));
      difftimes=diff(find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1))));
      clear masktime
      if any(difftimes>1)
        finddifftimes=find(difftimes>1);
        for dd=1:length(finddifftimes)+1
          if dd==1
            masktime{dd}=masktime_tmp(1):masktime_tmp(find(difftimes>1));
          elseif dd==length(finddifftimes)+1
            masktime{dd}=masktime_tmp(finddifftimes(dd-1)+1):masktime_tmp(end);
          else
            masktime{dd}=masktime_tmp(finddifftimes(dd-1)+1):masktime_tmp(finddifftimes(dd));
          end
        end
      else
        masktime{1}=masktime_tmp;
      end
      for dd=1:length(masktime)
        %     if ll==1
        %       cfg.xlim=[.1 .35];
        %     elseif ll==3
        %       cfg.xlim=[.06 .42];
        %     elseif ll==5
        %       cfg.xlim=[-.04 .36];
        %     elseif ll==7
        %       cfg.xlim=[.1 .44];
        %     end
        cfg.xlim=[tmp.time(masktime{dd}(1)) tmp.time(masktime{dd}(end))];
        cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2)) ),2),3))));
        %           cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,1:2,dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2))),2),3))));
        %       cfg.baseline=[tmp1.time(1) tmp1.time(9)];
        cfg.baseline=baseline2;
        figure(10*ll+8);
        ft_topoplotTFR(cfg,tmpuA);
        print(10*ll+8,[fdir 'tfrlo_topoU_theta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        figure(10*ll+3);
        ft_topoplotTFR(cfg,tmpsA);
        print(10*ll+3,[fdir 'tfrlo_topoM_theta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        cfg.zlim=[-1.4 1.4];
        %     cfg.zlim='maxabs';
        %     cfg.baseline='no';
        figure(10*ll+4);
        ft_topoplotTFR(cfg,tmpA);
        print(10*ll+4,[fdir 'tfrlo_topoDiff_theta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
      end
      
      
      if 0; %~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2})))
        %         chansel=setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2}));
        chansel=cfg.highlightchannel;
        zlim=[-3 3];
        figinds=10*ll+1000;
        figstrings=[];
        figstrings{1}=[fdir 'tfrlo_final_Xtheta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        baseline=[tmp.time(1) tmp.time(9)];
        tfr_subchannel_3cond_plot_pow(tmp,tmpu,tmps,chansel,zlim,figinds,figstrings,baseline)
      end
      
    end
    
    % alpha topo
    if ~isempty(find(squeeze(any(mean(tmp.mask(:,3:5,:),2),1))))
      cfg=[];
      if sleep==0 && ll==3 && combval==1  % ??
        cfg.ylim=[10 12];
        %       cfg.xlim=[0 .1];
        %       elseif ll==7
        %         cfg.ylim=[8 10];
        %       cfg.xlim=[.16 .3];
      elseif sleep==0 && ll==1 && combval==1
        cfg.ylim=[8 14];
      else
        disp('get ylim right per ll alpha')
        keyboard
      end
      masktime=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1)));
      cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
      cfg.parameter='powspctrm';
      cfg.layout='elec1010.lay';
      cfg.maskalpha=0.5;
      %       cfg.zlim=[-1.4 1.4];
      cfg.zlim=[-9 9];
      cfg.highlight='on';
      cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2)) ),2),3))));
      cfg.baseline=baseline2;
      %       cfg.baseline=[tmp1.time(1) tmp1.time(9)];
      figure(10*ll+5);
      ft_topoplotTFR(cfg,tmpuA);
      print(10*ll+5,[fdir 'tfrlo_topoU_alpha_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
      figure(10*ll+6);
      ft_topoplotTFR(cfg,tmpsA);
      print(10*ll+6,[fdir 'tfrlo_topoM_alpha_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
      %       cfg.zlim=[-1.4 1.4];
      cfg.zlim=[-9 9];
      %     cfg.zlim='maxabs';
      %     cfg.baseline='no';
      figure(10*ll+7);
      ft_topoplotTFR(cfg,tmpA);
      print(10*ll+7,[fdir 'tfrlo_topoDiff_alpha_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
      
      if 0 %~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2})))
        chansel=setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2}));
        zlim=[-3 3];
        figinds=10*ll+1000+1;
        figstrings=[];
        figstrings{1}=[fdir 'tfrlo_final_Xalpha_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        baseline=[tmp.time(1) tmp.time(9)];
        tfr_subchannel_3cond_plot_pow(tmp,tmpu,tmps,chansel,zlim,figinds,figstrings,baseline)
      end
    end
    
    % beta topo
    if ~isempty(find(squeeze(any(mean(tmp.mask(:,6:10,:),2),1))))
      cfg=[];
      if sleep==0 && ll==5
        cfg.ylim=[14 30];
        %       cfg.xlim=[.2 .28];
      elseif sleep==0 && ll==1 && combval==1
        cfg.ylim=[13 15];
      elseif sleep==0 && ll==3 && combval==1
        cfg.ylim=[14 28];
      elseif sleep==0 && ll==6 && combval==1
        cfg.ylim=[16 28];
      else
        disp('get ylim right per ll beta')
        keyboard
      end
      masktime=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1)));
      cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
      cfg.parameter='powspctrm';
      cfg.layout='elec1010.lay';
      cfg.maskalpha=0.5;
      cfg.zlim=[-1.1 1.1];
      cfg.highlight='on';
      cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2)) ),2),3))));
      %       cfg.baseline=[tmp1.time(1) tmp1.time(9)];
      cfg.baseline=baseline2;
      cfg.comment='no';
      figure(10*ll+5);
      ft_topoplotTFR(cfg,tmpuA);
      print(10*ll+5,[fdir 'tfrlo_topoU_beta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
      figure(10*ll+6);
      ft_topoplotTFR(cfg,tmpsA);
      print(10*ll+6,[fdir 'tfrlo_topoM_beta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
      cfg.zlim=[-1.1 1.1];
      %     cfg.baseline='no';
      %     cfg.zlim='maxabs';
      figure(10*ll+7);
      ft_topoplotTFR(cfg,tmpA);
      print(10*ll+7,[fdir 'tfrlo_topoDiff_beta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
      if 0 %~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2})))
        chansel=setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2}));
        zlim=[-2 2];
        figinds=10*ll+1000+2;
        figstrings=[];
        figstrings{1}=[fdir 'tfrlo_final_Xbeta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        baseline=[tmp.time(1) tmp.time(9)];
        tfr_subchannel_3cond_plot_pow(tmp,tmpu,tmps,chansel,zlim,figinds,figstrings,baseline)
      end
    end
    
    
    
    
    
    %    %  % %%%%%%% PLV   %%%%%%%%
    for adda=2
      if combval==1
        if adda==1
          cfg=[];
          cfg.latency=[stattl_mc_comb1plvadi.time(1) stattl_mc_comb1plvadi.time(end)];
          tmp=ft_selectdata(cfg,gravelo_TPA_MSPN_comb1);
          tmp.mask=stattl_mc_comb1plvadi.mask;
          tmps=ft_selectdata(cfg,gravelo_tacMSpN_comb1);
          tmps.mask=stattl_mc_comb1plvadi.mask;
          tmpu=ft_selectdata(cfg,gravelo_tacPaud_comb1);
          tmpu.mask=stattl_mc_comb1plvadi.mask;
          
          cfg=[];
          tmpA=ft_selectdata(cfg,gravelo_TPA_MSPN_comb1);
          tmpA.mask=logical(ones(size(tmpA.powspctrm)));
          tmpA.mask(:,:,dsearchn(tmpA.time',stattl_mc_comb1plvadi.time(1)):dsearchn(tmpA.time',stattl_mc_comb1plvadi.time(end)))=stattl_mc_comb1plvadi.mask;
          tmpsA=ft_selectdata(cfg,gravelo_tacMSpN_comb1);
          tmpsA.mask=logical(ones(size(tmpsA.powspctrm)));
          tmpsA.mask(:,:,dsearchn(tmpsA.time',stattl_mc_comb1plvadi.time(1)):dsearchn(tmpsA.time',stattl_mc_comb1plvadi.time(end)))=stattl_mc_comb1plvadi.mask;
          tmpuA=ft_selectdata(cfg,gravelo_tacPaud_comb1);
          tmpuA.mask=logical(ones(size(tmpuA.powspctrm)));
          tmpuA.mask(:,:,dsearchn(tmpuA.time',stattl_mc_comb1plvadi.time(1)):dsearchn(tmpuA.time',stattl_mc_comb1plvadi.time(end)))=stattl_mc_comb1plvadi.mask;
        elseif adda==2
          cfg=[];
          cfg.latency=[stattl_mc_comb1plvdai.time(1) stattl_mc_comb1plvdai.time(end)];
          tmp=ft_selectdata(cfg,gravelo_TPA_MSPN_comb1);
          tmp.mask=stattl_mc_comb1plvdai.mask;
          tmps=ft_selectdata(cfg,gravelo_tacMSpN_comb1);
          tmps.mask=stattl_mc_comb1plvdai.mask;
          tmpu=ft_selectdata(cfg,gravelo_tacPaud_comb1);
          tmpu.mask=stattl_mc_comb1plvdai.mask;
          
          cfg=[];
          tmpA=ft_selectdata(cfg,gravelo_TPA_MSPN_comb1);
          tmpA.mask=logical(ones(size(tmpA.powspctrm)));
          tmpA.mask(:,:,dsearchn(tmpA.time',stattl_mc_comb1plvdai.time(1)):dsearchn(tmpA.time',stattl_mc_comb1plvdai.time(end)))=stattl_mc_comb1plvdai.mask;
          tmpsA=ft_selectdata(cfg,gravelo_tacMSpN_comb1);
          tmpsA.mask=logical(ones(size(tmpsA.powspctrm)));
          tmpsA.mask(:,:,dsearchn(tmpsA.time',stattl_mc_comb1plvdai.time(1)):dsearchn(tmpsA.time',stattl_mc_comb1plvdai.time(end)))=stattl_mc_comb1plvdai.mask;
          tmpuA=ft_selectdata(cfg,gravelo_tacPaud_comb1);
          tmpuA.mask=logical(ones(size(tmpuA.powspctrm)));
          tmpuA.mask(:,:,dsearchn(tmpuA.time',stattl_mc_comb1plvdai.time(1)):dsearchn(tmpuA.time',stattl_mc_comb1plvdai.time(end)))=stattl_mc_comb1plvdai.mask;
        end
      elseif combval==2
        if adda==1
          cfg=[];
          cfg.latency=[stattl_mc_comb2plvadi.time(1) stattl_mc_comb2plvadi.time(end)];
          tmp=ft_selectdata(cfg,gravelo_TPA_MSPN_comb2);
          tmp.mask=stattl_mc_comb2plvadi.mask;
          tmps=ft_selectdata(cfg,gravelo_tacMSpN_comb2);
          tmps.mask=stattl_mc_comb2plvadi.mask;
          tmpu=ft_selectdata(cfg,gravelo_tacPaud_comb2);
          tmpu.mask=stattl_mc_comb2plvadi.mask;
        elseif adda==2
          cfg=[];
          cfg.latency=[stattl_mc_comb2plvdai.time(1) stattl_mc_comb2plvdai.time(end)];
          tmp=ft_selectdata(cfg,gravelo_TPA_MSPN_comb2);
          tmp.mask=stattl_mc_comb2plvdai.mask;
          tmps=ft_selectdata(cfg,gravelo_tacMSpN_comb2);
          tmps.mask=stattl_mc_comb2plvdai.mask;
          tmpu=ft_selectdata(cfg,gravelo_tacPaud_comb2);
          tmpu.mask=stattl_mc_comb2plvdai.mask;
        end
      end
      
      tmp.plvavgangwrap=wrapToPi(tmp.plvavgang);
      tmps.plvavgangwrap=wrapToPi(tmps.plvavgang);
      tmpu.plvavgangwrap=wrapToPi(tmpu.plvavgang);
      tmpA.plvavgangwrap=wrapToPi(tmpA.plvavgang);
      tmpsA.plvavgangwrap=wrapToPi(tmpsA.plvavgang);
      tmpuA.plvavgangwrap=wrapToPi(tmpuA.plvavgang);
      gravelo_tNulAlone_comb.plvavgangwrap=wrapToPi(gravelo_tNulAlone_comb.plvavgang);
      gravelo_tTacAlone_comb.plvavgangwrap=wrapToPi(gravelo_tTacAlone_comb.plvavgang);
      gravelo_tAudAlone_comb.plvavgangwrap=wrapToPi(gravelo_tAudAlone_comb.plvavgang);
      gravelo_tMSAlone_comb.plvavgangwrap= wrapToPi(gravelo_tMSAlone_comb.plvavgang);
      
      %       if 0
      %         cfg=[];
      %         cfg.avgoverchan='yes';
      %         if ~isempty(tmp.label(find(ceil(mean(mean(tmp.mask(:,:,:),2),3)))))
      %           cfg.channel=tmp.label(find(ceil(mean(mean(tmp.mask(:,:,:),2),3))));
      %         else
      %           cfg.channel='all';
      %         end
      %         tmp1=ft_selectdata(cfg,tmp);
      %         tmp1.mask=logical(ceil(tmp1.mask));
      %
      %         figure(100*ll);
      %         cfg=[];
      %         cfg.parameter='plvavgabs';
      %         cfg.layout='elec1010.lay';
      %         cfg.maskparameter='mask';
      %         cfg.maskalpha=0.5;
      %         cfg.zlim='maxabs';
      %         ft_singleplotTFR(cfg,tmp1);
      %         print(100*ll,[fdir 'plvabslo_final_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
      %         figure(100*ll+10);
      %         cfg.parameter='plvavgang';
      %         ft_singleplotTFR(cfg,tmp1);
      %         print(100*ll+10,[fdir 'plvanglo_final_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
      %       end
      
      
      % plot TFR of plv abs and plv ang, for each suplot: u and m and difference.
      zlim=[0 0.6; -.25 .25; -20 400; -10 300];
      
      if 0
        pow{1}=gravelo_tMSAlone_comb;
        pow{2}=tmpu;
        pow{3}=tmps;
        pow{4}=tmp;
        
        chansel=chanplot{1};
        figinds=[100*ll;  100*ll+10];
        figstrings{1}=[fdir 'plvabslo_crop_FC_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        figstrings{2}=[fdir 'plvanglo_crop_FC_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
        
        chansel=chanplot{2};
        figinds=[100*ll+1;  100*ll+10+1];
        figstrings{1}=[fdir 'plvabslo_crop_OP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        figstrings{2}=[fdir 'plvanglo_crop_OP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
        
        chansel=[chanplot{1} chanplot{2}];
        figinds=[100*ll+2;  100*ll+10+2];
        figstrings{1}=[fdir 'plvabslo_crop_FCOP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        figstrings{2}=[fdir 'plvanglo_crop_FCOP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
      end
      
      pow{1}=gravelo_tMSAlone_comb;
      pow{2}=tmpuA;
      pow{3}=tmpsA;
      pow{4}=tmpA;
      
      chansel=chanplot{1};
      figinds=[100*ll;  100*ll+10];
      figstrings{1}=[fdir 'plvabslo_all_FC_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
      figstrings{2}=[fdir 'plvanglo_all_FC_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
      tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
      
      chansel=chanplot{2};
      figinds=[100*ll+1;  100*ll+10+1];
      figstrings{1}=[fdir 'plvabslo_all_OP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
      figstrings{2}=[fdir 'plvanglo_all_OP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
      tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
      
      chansel=[chanplot{1} chanplot{2}];
      figinds=[100*ll+2;  100*ll+10+2];
      figstrings{1}=[fdir 'plvabslo_all_FCOP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
      figstrings{2}=[fdir 'plvanglo_all_FCOP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
      tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
      
      masktime=find(squeeze(any(mean(tmp.mask(:,1:end,:),2),1)));
      if ~isempty(masktime)
        chansel=tmp.label(find(ceil(mean(mean(tmp.mask(:,:,dsearchn(tmp.time',tmp.time(masktime(1))):dsearchn(tmp.time',tmp.time(masktime(end)))),2),3))));
        figinds=[100*ll+9;  100*ll+10+9];
        figstrings{1}=[fdir 'plvabslo_all_allsig_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        figstrings{2}=[fdir 'plvanglo_all_allsig_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
      end
      
      % theta topo
      if ~isempty(find(squeeze(any(mean(tmp.mask(:,1:2,:),2),1))))
        cfg=[];
        if sleep==0 && ll==3
          cfg.ylim=[4 6.5];
        elseif sleep==0 && ll==6 && adda==1
          cfg.ylim=[5 6.5];
        elseif sleep==0 && ll==7
          cfg.ylim=[4 6.5];
        else
          disp('get ylim right per ll theta plv')
          keyboard
        end
        masktime=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1)));
        cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
        cfg.parameter='plvavgabs';
        cfg.layout='elec1010.lay';
        cfg.maskalpha=0.5;
        cfg.zlim=[-.25 .25];
        cfg.highlight='on';
        cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2))),2),3))));
        cfg.baseline=baseline2;
        cfg.zlim=[-.25 .25];
        cfg.comment='no';
        figure(100*ll+8);
        ft_topoplotTFR(cfg,tmpuA);
        print(100*ll+8,[fdir 'plvabslo_topoU_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        figure(100*ll+3);
        ft_topoplotTFR(cfg,tmpsA);
        print(100*ll+3,[fdir 'plvabslo_topoM_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        cfg.zlim=[-.15 .15];
        figure(100*ll+4);
        ft_topoplotTFR(cfg,tmpA);
        print(100*ll+4,[fdir 'plvabslo_topoDiff_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        
        if ~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2})))
          chansel=cfg.highlightchannel;
          zlim=[-.3 .3; 0 0.35; -30 30; -60 60];
          zlim=[0 0.6; -.25 .25; -20 400; -10 300];
          figinds=[100*ll+20+1;  100*ll+30+1];
          figstrings{1}=[fdir 'plvabslo_all_Xtheta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          figstrings{2}=[fdir 'plvanglo_all_Xtheta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
        end
        
        cfg.parameter='plvavgangwrap';
        if sleep==0 && ll==3
          cfg.xlim=[0.19 0.19];
        elseif sleep==0 && ll==6 && combval==1 && adda==1
          cfg.xlim=[0.41 0.41];
        elseif sleep==0 && ll==6 && combval==2 && adda==1
          cfg.xlim=[0.32 0.32];
        elseif sleep==0 && ll==7 && combval==1
          cfg.xlim=[.21 .21];
        elseif sleep==0 && ll==7 && combval==2
          cfg.xlim=[.23 .23];
        else
          disp('get xlim right per ll theta plv')
          keyboard
        end
        cfg.zlim=[-4 4];
        figure(100*ll+8+10);
        ft_topoplotTFR(cfg,tmpuA);
        print(100*ll+8+10,[fdir 'plvanglo_topoU_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        figure(100*ll+3+10);
        ft_topoplotTFR(cfg,tmpsA);
        print(100*ll+3+10,[fdir 'plvanglo_topoM_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        figure(100*ll+4+10);
        ft_topoplotTFR(cfg,tmpA);
        print(100*ll+4+10,[fdir 'plvanglo_topoDiff_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
      end
      
      % alpha topo
      if ~isempty(find(squeeze(any(mean(tmp.mask(:,3:5,:),2),1))))
        cfg=[];
        if sleep==0 && ll==6 && [adda==1 || (combval==2 && adda==2)]
          cfg.ylim=[8 11];
        elseif sleep==0 && ll==6 && combval==1 && adda==2
          cfg.ylim=[8 14]; % slightly into beta but then no beta
        elseif sleep==0 && ll==3 && combval==1 && adda==1
          cfg.ylim=[8 12];
        elseif sleep==0 && ll==3 && adda==2
          cfg.ylim=[9 12];
        elseif sleep==0 && ll==3 && combval==2 && adda==1
          cfg.ylim=[8 9];
        else
          disp('get ylim right per ll alpha plv')
          keyboard
        end
        masktime_tmp=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1)));
        difftimes=diff(find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1))));
        clear masktime
        if any(difftimes>1)
          finddifftimes=find(difftimes>1);
          for dd=1:length(finddifftimes)+1
            if dd==1
              masktime{dd}=masktime_tmp(1):masktime_tmp(find(difftimes>1));
            elseif dd==length(finddifftimes)+1
              masktime{dd}=masktime_tmp(finddifftimes(dd-1)+1):masktime_tmp(end);
            else
              masktime{dd}=masktime_tmp(finddifftimes(dd-1)+1):masktime_tmp(finddifftimes(dd));
            end
          end
        else
          masktime{1}=masktime_tmp;
        end
        for dd=1:length(masktime)
          
          cfg.parameter='plvavgabs';
          cfg.layout='elec1010.lay';
          cfg.maskalpha=0.5;
          %     cfg.zlim='maxabs';
          cfg.highlight='on';
          cfg.xlim=[tmp.time(masktime{dd}(1)) tmp.time(masktime{dd}(end))];
          cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2)) ),2),3))));
          cfg.baseline=baseline2;
          cfg.zlim=[-.15 .15];
          cfg.comment='no';
          cfg.comment='auto';
          figure(100*ll+5);
          ft_topoplotTFR(cfg,tmpuA);
          print(100*ll+5,[fdir 'plvabslo_topoU_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '_time' num2str(dd) '.eps'],'-depsc2')
          figure(100*ll+6);
          ft_topoplotTFR(cfg,tmpsA);
          print(100*ll+6,[fdir 'plvabslo_topoM_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '_time' num2str(dd) '.eps'],'-depsc2')
          cfg.zlim=[-.1 .1];
          figure(100*ll+7);
          ft_topoplotTFR(cfg,tmpA);
          print(100*ll+7,[fdir 'plvabslo_topoDiff_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '_time' num2str(dd) '.eps'],'-depsc2')
        end
        
        
        if ~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2})))
          chansel=cfg.highlightchannel;
          zlim=[-.3 .3; 0 0.35; -30 30; -60 60];
          zlim=[0 0.6; -.25 .25; -20 400; -10 300];
          figinds=[100*ll+40+1;  100*ll+50+1];
          figstrings{1}=[fdir 'plvabslo_all_Xalpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          figstrings{2}=[fdir 'plvanglo_all_Xalpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
        end
        
        cfg.parameter='plvavgangwrap';
        if sleep==0 && ll==6
          cfg.xlim=[.34 .34];
        elseif sleep==0 && ll==3 && combval==1 && adda==1
          cfg.xlim=[0 0];
        elseif sleep==0 && ll==3 && adda==2
          cfg.xlim=[-.06 -.06];
        elseif sleep==0 && ll==3 && combval==2 && adda==1
          cfg.xlim=[.21 .21];
        else
          disp('get xlim right per ll alpha plv')
          keyboard
        end
        cfg.zlim=[-4 4];
        figure(100*ll+5+10);
        ft_topoplotTFR(cfg,tmpuA);
        print(100*ll+5+10,[fdir 'plvanglo_topoU_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        figure(100*ll+6+10);
        ft_topoplotTFR(cfg,tmpsA);
        print(100*ll+6+10,[fdir 'plvanglo_topoM_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        figure(100*ll+7+10);
        ft_topoplotTFR(cfg,tmpA);
        print(100*ll+7+10,[fdir 'plvanglo_topoDiff_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
      end
      
      % beta topo
      if ~isempty(find(squeeze(any(mean(tmp.mask(:,6:10,:),2),1))))
        cfg=[];
        if sleep==0 && ll==6 && combval==1 && adda==2
          cfg.ylim=[13 15];
        elseif sleep==0 && ll==3
          cfg.ylim=[13 15];
        elseif sleep==1 && ll==1
          cfg.ylim=[13 21];
        else
          disp('get ylim right per ll beta plv')
          tmp.freq(find(mean(mean(tmp.mask,1),3)))
          keyboard
        end
        masktime=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1)));
        cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
        cfg.parameter='plvavgabs';
        cfg.layout='elec1010.lay';
        cfg.maskalpha=0.5;
        %     cfg.zlim='maxabs';
        cfg.highlight='on';
        cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2)) ),2),3))));
        cfg.baseline=baseline2;
        cfg.zlim=[-.15 .15];
        cfg.comment='no';
        figure(100*ll+5);
        ft_topoplotTFR(cfg,tmpuA);
        print(100*ll+5,[fdir 'plvabslo_topoU_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        figure(100*ll+6);
        ft_topoplotTFR(cfg,tmpsA);
        print(100*ll+6,[fdir 'plvabslo_topoM_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        cfg.zlim=[-.1 .1];
        figure(100*ll+7);
        ft_topoplotTFR(cfg,tmpA);
        print(100*ll+7,[fdir 'plvabslo_topoDiff_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        
        if ~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2})))
          chansel=cfg.highlightchannel;
          zlim=[-.3 .3; 0 0.35; -30 30; -60 60];
          zlim=[0 0.6; -.25 .25; -20 400; -10 300];
          figinds=[100*ll+40+1;  100*ll+50+1];
          figstrings{1}=[fdir 'plvabslo_final_Xbeta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          figstrings{2}=[fdir 'plvanglo_final_Xbeta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
        end
        
        cfg.parameter='plvavgangwrap';
        if sleep==0 && ll==6 && combval==1 && adda==2
          cfg.xlim=[.32 .32];
        elseif sleep==0 && ll==3 && combval==1 && adda==1
          cfg.xlim=[0 0];
        elseif sleep==0 && ll==3 && adda==2
          cfg.xlim=[-.06 -.06];
        elseif sleep==1 && ll==1 && adda==2
          cfg.xlim=[0.05 0.61];
        else
          disp('get xlim right per ll beta plv')
          tmp.time(find(mean(mean(tmp.mask,2),3)))
          keyboard
        end
        cfg.zlim=[-4 4];
        figure(100*ll+5+10);
        ft_topoplotTFR(cfg,tmpuA);
        print(100*ll+5+10,[fdir 'plvanglo_topoU_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        figure(100*ll+6+10);
        ft_topoplotTFR(cfg,tmpsA);
        print(100*ll+6+10,[fdir 'plvanglo_topoM_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        figure(100*ll+7+10);
        ft_topoplotTFR(cfg,tmpA);
        print(100*ll+7+10,[fdir 'plvanglo_topoDiff_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
      end
      
    end % adda
    
  end % combval
  
  %%
  
  if synchasynch && ll<5
    for combval=1
      close all
      clear adda
      % power
      if combval==1
        cfg=[];
        cfg.latency=[stattl_synch_comb1.time(1) stattl_synch_comb1.time(end)];
        tmp=ft_selectdata(cfg,gravelo_TMSs_TMSa_comb1);
        tmp.mask=stattl_synch_comb1.mask;
        tmps=ft_selectdata(cfg,gravelo_tMSasynch_comb1);
        tmps.mask=stattl_synch_comb1.mask;
        tmpu=ft_selectdata(cfg,gravelo_tMSsynch_comb1);
        tmpu.mask=stattl_synch_comb1.mask;
        
        cfg=[];
        tmpA=ft_selectdata(cfg,gravelo_TMSs_TMSa_comb1);
        tmpA.mask=logical(ones(size(tmpA.powspctrm)));
        tmpA.mask(:,:,dsearchn(tmpA.time',stattl_synch_comb1.time(1)):dsearchn(tmpA.time',stattl_synch_comb1.time(end)))=stattl_synch_comb1.mask;
        tmpsA=ft_selectdata(cfg,gravelo_tMSasynch_comb1);
        tmpsA.mask=logical(ones(size(tmpsA.powspctrm)));
        tmpsA.mask(:,:,dsearchn(tmpsA.time',stattl_synch_comb1.time(1)):dsearchn(tmpsA.time',stattl_synch_comb1.time(end)))=stattl_synch_comb1.mask;
        tmpuA=ft_selectdata(cfg,gravelo_tMSsynch_comb1);
        tmpuA.mask=logical(ones(size(tmpuA.powspctrm)));
        tmpuA.mask(:,:,dsearchn(tmpuA.time',stattl_synch_comb1.time(1)):dsearchn(tmpuA.time',stattl_synch_comb1.time(end)))=stattl_synch_comb1.mask;
      elseif combval==2
        cfg=[];
        cfg.latency=[stattl_synch_comb2.time(1) stattl_synch_comb2.time(end)];
        tmp=ft_selectdata(cfg,gravelo_TMSs_TMSa_comb2);
        tmp.mask=stattl_synch_comb2.mask;
        tmps=ft_selectdata(cfg,gravelo_tMSasynch_comb2);
        tmps.mask=stattl_synch_comb2.mask;
        tmpu=ft_selectdata(cfg,gravelo_tMSsynch_comb2);
        tmpu.mask=stattl_synch_comb2.mask;
      end
      gravelo_tMSAlone_comb.mask= logical(ones(size(gravelo_tMSAlone_comb.powspctrm)));
      
      
      % plot TFR of powspctrm, for each suplot: u and m and difference.
      zlim=[-3 3];
      baseline3=[-.15 -.07];
      
      if 0
        pow{1}=gravelo_tMSAlone_comb;
        pow{2}=tmpu;
        pow{3}=tmps;
        pow{4}=tmp;
        
        chansel=chanplot{1};
        figinds=10*ll;
        figstrings=[];
        figstrings{1}=[fdir 'tfrlo_synch_crop_FC_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline3)
        
        chansel=chanplot{2};
        figinds=10*ll+1;
        figstrings=[];
        figstrings{1}=[fdir 'tfrlo_synch_crop_OP_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline3)
        
        chansel=[chanplot{1} chanplot{2}];
        figinds=10*ll+2;
        figstrings=[];
        figstrings{1}=[fdir 'tfrlo_synch_crop_FCOP_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline3)
      end
      
      pow{1}=gravelo_tMSAlone_comb;
      pow{2}=tmpuA;
      pow{3}=tmpsA;
      pow{4}=tmpA;
      
      chansel=chanplot{1};
      figinds=10*ll;
      figstrings=[];
      figstrings{1}=[fdir 'tfrlo_synch_all_FC_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
      tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline3)
      
      chansel=chanplot{2};
      figinds=10*ll+1;
      figstrings=[];
      figstrings{1}=[fdir 'tfrlo_synch_all_OP_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
      tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline3)
      
      chansel=[chanplot{1} chanplot{2}];
      figinds=10*ll+2;
      figstrings=[];
      figstrings{1}=[fdir 'tfrlo_synch_all_FCOP_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
      tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline3)
      
      masktime=find(squeeze(any(mean(tmp.mask(:,1:end,:),2),1)));
      if ~isempty(masktime)
        chansel=tmp.label(find(ceil(mean(mean(tmp.mask(:,:,dsearchn(tmp.time',tmp.time(masktime(1))):dsearchn(tmp.time',tmp.time(masktime(end)))),2),3))));
        figinds=10*ll+9;
        figstrings=[];
        figstrings{1}=[fdir 'tfrlo_synch_all_allsig_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        tfr_subchannel_3cond_plot_pow(pow,chansel,zlim,figinds,figstrings,baseline3)
      end
      
      % theta topo
      if ~isempty(find(squeeze(any(mean(tmp.mask(:,1:2,:),2),1))))
        masktime=find(squeeze(any(mean(tmp.mask(:,1:2,:),2),1)));
        cfg=[];
        cfg.parameter='powspctrm';
        cfg.layout='elec1010.lay';
        cfg.maskalpha=0.5;
        if ll==4
          cfg.ylim=[4 8.5];
        else
          cfg.ylim=[4 6.5];
        end
        cfg.zlim=[-1.4 1.4];
        cfg.highlight='on';
        cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
        %     if ll==1
        %       cfg.xlim=[.1 .35];
        %     elseif ll==3
        %       cfg.xlim=[.06 .42];
        %     elseif ll==5
        %       cfg.xlim=[-.04 .36];
        %     elseif ll==7
        %       cfg.xlim=[.1 .44];
        %     end
        cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,1:2,dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2))),2),3))));
        cfg.baseline=baseline3;
        figure(10*ll+8);
        ft_topoplotTFR(cfg,tmpuA);
        print(10*ll+8,[fdir 'tfrlo_synch_topoU_theta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        figure(10*ll+3);
        ft_topoplotTFR(cfg,tmpsA);
        print(10*ll+3,[fdir 'tfrlo_synch_topoM_theta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        cfg.zlim=[-1.4 1.4];
        %     cfg.zlim='maxabs';
        %     cfg.baseline='no';
        figure(10*ll+4);
        ft_topoplotTFR(cfg,tmpA);
        print(10*ll+4,[fdir 'tfrlo_synch_topoDiff_theta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        
        if 0 %~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2})))
          chansel=setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2}));
          zlim=[-3 3];
          figinds=10*ll+1000;
          figstrings=[];
          figstrings{1}=[fdir 'tfrlo_synch_final_Xtheta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          baseline=[tmp.time(1) tmp.time(9)];
          tfr_subchannel_3cond_plot_pow(tmp,tmpu,tmps,chansel,zlim,figinds,figstrings,baseline)
        end
        
      end
      
      % alpha topo
      if ~isempty(find(squeeze(any(mean(tmp.mask(:,3:5,:),2),1))))
        cfg=[];
        if ll==3
          cfg.ylim=[8 11];
        elseif ll==4
          cfg.ylim=[10 14];
        else
          keyboard
        end
        disp('get ylim right per ll alpha synch')
        masktime=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1)));
        if ll==3 || ll==4
          cfg.xlim=[tmp.time(masktime(1)) .5];
        else
          cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
        end
        cfg.parameter='powspctrm';
        cfg.layout='elec1010.lay';
        cfg.maskalpha=0.5;
        cfg.zlim=[-9 9];
        cfg.highlight='on';
        cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2)) ),2),3))));
        cfg.baseline=baseline3;
        figure(10*ll+5);
        ft_topoplotTFR(cfg,tmpuA);
        print(10*ll+5,[fdir 'tfrlo_synch_topoU_alpha_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        figure(10*ll+6);
        ft_topoplotTFR(cfg,tmpsA);
        print(10*ll+6,[fdir 'tfrlo_synch_topoM_alpha_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        cfg.zlim=[-9 9];
        %     cfg.zlim='maxabs';
        %     cfg.baseline='no';
        figure(10*ll+7);
        ft_topoplotTFR(cfg,tmpA);
        print(10*ll+7,[fdir 'tfrlo_synch_topoDiff_alpha_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        if ll==3 || ll==4
          cfg.xlim=[.5 tmp.time(masktime(end))];
          cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2)) ),2),3))));
          figure(10*ll+5);
          ft_topoplotTFR(cfg,tmpuA);
          print(10*ll+5,[fdir 'tfrlo_synch_topoU_alpha2_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          figure(10*ll+6);
          ft_topoplotTFR(cfg,tmpsA);
          print(10*ll+6,[fdir 'tfrlo_synch_topoM_alpha2_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          cfg.zlim=[-9 9];
          %     cfg.zlim='maxabs';
          %     cfg.baseline='no';
          figure(10*ll+7);
          ft_topoplotTFR(cfg,tmpA);
          print(10*ll+7,[fdir 'tfrlo_synch_topoDiff_alpha2_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        end
        
        if 0 %~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2})))
          chansel=setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2}));
          zlim=[-3 3];
          figinds=10*ll+1000+1;
          figstrings=[];
          figstrings{1}=[fdir 'tfrlo_synch_final_Xalpha_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          baseline=[tmp.time(1) tmp.time(9)];
          tfr_subchannel_3cond_plot_pow(tmp,tmpu,tmps,chansel,zlim,figinds,figstrings,baseline)
        end
      end
      
      % beta topo
      if ~isempty(find(squeeze(any(mean(tmp.mask(:,6:10,:),2),1))))
        cfg=[];
        if ll==3
          cfg.ylim=[14 28];
        elseif ll==4
          cfg.ylim=[16 30];
        else
          keyboard
        end
        disp('get ylim right per ll beta synch')
        masktime=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1)));
        if ll==3 || ll==4
          cfg.xlim=[tmp.time(masktime(1)) .5];
        else
          cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
        end
        cfg.parameter='powspctrm';
        cfg.layout='elec1010.lay';
        cfg.maskalpha=0.5;
        cfg.zlim=[-1.1 1.1];
        cfg.highlight='on';
        cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2)) ),2),3))));
        cfg.baseline=baseline3;
        figure(10*ll+5);
        ft_topoplotTFR(cfg,tmpuA);
        print(10*ll+5,[fdir 'tfrlo_synch_topoU_beta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        figure(10*ll+6);
        ft_topoplotTFR(cfg,tmpsA);
        print(10*ll+6,[fdir 'tfrlo_synch_topoM_beta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        cfg.zlim=[-1.1 1.1];
        %     cfg.baseline='no';
        %     cfg.zlim='maxabs';
        figure(10*ll+7);
        ft_topoplotTFR(cfg,tmpA);
        print(10*ll+7,[fdir 'tfrlo_synch_topoDiff_beta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        if ll==3 || ll==4
          cfg.xlim=[.5 tmp.time(masktime(end))];
          cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2)) ),2),3))));
          figure(10*ll+5);
          ft_topoplotTFR(cfg,tmpuA);
          print(10*ll+5,[fdir 'tfrlo_synch_topoU_beta2_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          figure(10*ll+6);
          ft_topoplotTFR(cfg,tmpsA);
          print(10*ll+6,[fdir 'tfrlo_synch_topoM_beta2_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          cfg.zlim=[-1.1 1.1];
          %     cfg.baseline='no';
          %     cfg.zlim='maxabs';
          figure(10*ll+7);
          ft_topoplotTFR(cfg,tmpA);
          print(10*ll+7,[fdir 'tfrlo_synch_topoDiff_beta2_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        end
        
        
        if 0 %~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2})))
          chansel=setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2}));
          zlim=[-2 2];
          figinds=10*ll+1000+2;
          figstrings=[];
          figstrings{1}=[fdir 'tfrlo_synch_final_Xbeta_combval' num2str(combval) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          baseline=[tmp.time(1) tmp.time(9)];
          tfr_subchannel_3cond_plot_pow(tmp,tmpu,tmps,chansel,zlim,figinds,figstrings,baseline)
        end
      end
      
      
      
      
      
      %    %  % %%%%%%% PLV   %%%%%%%%
      for adda=2
        if combval==1
          if adda==1
            cfg=[];
            cfg.latency=[stattl_synch_comb1plvadi.time(1) stattl_synch_comb1plvadi.time(end)];
            tmp=ft_selectdata(cfg,gravelo_TMSs_TMSa_comb1);
            tmp.mask=stattl_synch_comb1plvadi.mask;
            tmps=ft_selectdata(cfg,gravelo_tMSasynch_comb1);
            tmps.mask=stattl_synch_comb1plvadi.mask;
            tmpu=ft_selectdata(cfg,gravelo_tMSsynch_comb1);
            tmpu.mask=stattl_synch_comb1plvadi.mask;
            
            cfg=[];
            tmpA=ft_selectdata(cfg,gravelo_TMSs_TMSa_comb1);
            tmpA.mask=logical(ones(size(tmpA.powspctrm)));
            tmpA.mask(:,:,dsearchn(tmpA.time',stattl_synch_comb1plvadi.time(1)):dsearchn(tmpA.time',stattl_synch_comb1plvadi.time(end)))=stattl_synch_comb1plvadi.mask;
            tmpsA=ft_selectdata(cfg,gravelo_tMSasynch_comb1);
            tmpsA.mask=logical(ones(size(tmpsA.powspctrm)));
            tmpsA.mask(:,:,dsearchn(tmpsA.time',stattl_synch_comb1plvadi.time(1)):dsearchn(tmpsA.time',stattl_synch_comb1plvadi.time(end)))=stattl_synch_comb1plvadi.mask;
            tmpuA=ft_selectdata(cfg,gravelo_tMSsynch_comb1);
            tmpuA.mask=logical(ones(size(tmpuA.powspctrm)));
            tmpuA.mask(:,:,dsearchn(tmpuA.time',stattl_synch_comb1plvadi.time(1)):dsearchn(tmpuA.time',stattl_synch_comb1plvadi.time(end)))=stattl_synch_comb1plvadi.mask;
          elseif adda==2
            cfg=[];
            cfg.latency=[stattl_synch_comb1plvdai.time(1) stattl_synch_comb1plvdai.time(end)];
            tmp=ft_selectdata(cfg,gravelo_TMSs_TMSa_comb1);
            tmp.mask=stattl_synch_comb1plvdai.mask;
            tmps=ft_selectdata(cfg,gravelo_tMSasynch_comb1);
            tmps.mask=stattl_synch_comb1plvdai.mask;
            tmpu=ft_selectdata(cfg,gravelo_tMSsynch_comb1);
            tmpu.mask=stattl_synch_comb1plvdai.mask;
            
            cfg=[];
            tmpA=ft_selectdata(cfg,gravelo_TMSs_TMSa_comb1);
            tmpA.mask=logical(ones(size(tmpA.powspctrm)));
            tmpA.mask(:,:,dsearchn(tmpA.time',stattl_synch_comb1plvdai.time(1)):dsearchn(tmpA.time',stattl_synch_comb1plvdai.time(end)))=stattl_synch_comb1plvdai.mask;
            tmpsA=ft_selectdata(cfg,gravelo_tMSasynch_comb1);
            tmpsA.mask=logical(ones(size(tmpsA.powspctrm)));
            tmpsA.mask(:,:,dsearchn(tmpsA.time',stattl_synch_comb1plvdai.time(1)):dsearchn(tmpsA.time',stattl_synch_comb1plvdai.time(end)))=stattl_synch_comb1plvdai.mask;
            tmpuA=ft_selectdata(cfg,gravelo_tMSsynch_comb1);
            tmpuA.mask=logical(ones(size(tmpuA.powspctrm)));
            tmpuA.mask(:,:,dsearchn(tmpuA.time',stattl_synch_comb1plvdai.time(1)):dsearchn(tmpuA.time',stattl_synch_comb1plvdai.time(end)))=stattl_synch_comb1plvdai.mask;
          end
        elseif combval==2
          if adda==1
            cfg=[];
            cfg.latency=[stattl_synch_comb2plvadi.time(1) stattl_synch_comb2plvadi.time(end)];
            tmp=ft_selectdata(cfg,gravelo_TMSs_TMSa_comb2);
            tmp.mask=stattl_synch_comb2plvadi.mask;
            tmps=ft_selectdata(cfg,gravelo_tMSasynch_comb2);
            tmps.mask=stattl_synch_comb2plvadi.mask;
            tmpu=ft_selectdata(cfg,gravelo_tMSsynch_comb2);
            tmpu.mask=stattl_synch_comb2plvadi.mask;
          elseif adda==2
            cfg=[];
            cfg.latency=[stattl_synch_comb2plvdai.time(1) stattl_synch_comb2plvdai.time(end)];
            tmp=ft_selectdata(cfg,gravelo_TMSs_TMSa_comb2);
            tmp.mask=stattl_synch_comb2plvdai.mask;
            tmps=ft_selectdata(cfg,gravelo_tMSasynch_comb2);
            tmps.mask=stattl_synch_comb2plvdai.mask;
            tmpu=ft_selectdata(cfg,gravelo_tMSsynch_comb2);
            tmpu.mask=stattl_synch_comb2plvdai.mask;
          end
        end
        
        tmp.plvavgangwrap=wrapToPi(tmp.plvavgang);
        tmps.plvavgangwrap=wrapToPi(tmps.plvavgang);
        tmpu.plvavgangwrap=wrapToPi(tmpu.plvavgang);
        tmpA.plvavgangwrap=wrapToPi(tmpA.plvavgang);
        tmpsA.plvavgangwrap=wrapToPi(tmpsA.plvavgang);
        tmpuA.plvavgangwrap=wrapToPi(tmpuA.plvavgang);
        
        % plot TFR of plv abs and plv ang, for each suplot: u and m and difference.
        zlim=[0 0.6; -.25 .25; -20 400; -10 300];
        
        
        if 0
          pow{1}=gravelo_tMSAlone_comb;
          pow{2}=tmpu;
          pow{3}=tmps;
          pow{4}=tmp;
          
          chansel=chanplot{1};
          figinds=[100*ll;  100*ll+10];
          figstrings{1}=[fdir 'plvabslo_synch_crop_FC_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          figstrings{2}=[fdir 'plvanglo_synch_crop_FC_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
          
          chansel=chanplot{2};
          figinds=[100*ll+1;  100*ll+10+1];
          figstrings{1}=[fdir 'plvabslo_synch_crop_OP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          figstrings{2}=[fdir 'plvanglo_synch_crop_OP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
          
          chansel=[chanplot{1} chanplot{2}];
          figinds=[100*ll+2;  100*ll+10+2];
          figstrings{1}=[fdir 'plvabslo_synch_crop_FCOP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          figstrings{2}=[fdir 'plvanglo_synch_crop_FCOP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
        end
        
        pow{1}=gravelo_tMSAlone_comb;
        pow{2}=tmpuA;
        pow{3}=tmpsA;
        pow{4}=tmpA;
        
        chansel=chanplot{1};
        figinds=[100*ll;  100*ll+10];
        figstrings{1}=[fdir 'plvabslo_synch_all_FC_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        figstrings{2}=[fdir 'plvanglo_synch_all_FC_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
        
        chansel=chanplot{2};
        figinds=[100*ll+1;  100*ll+10+1];
        figstrings{1}=[fdir 'plvabslo_synch_all_OP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        figstrings{2}=[fdir 'plvanglo_synch_all_OP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
        
        chansel=[chanplot{1} chanplot{2}];
        figinds=[100*ll+2;  100*ll+10+2];
        figstrings{1}=[fdir 'plvabslo_synch_all_FCOP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        figstrings{2}=[fdir 'plvanglo_synch_all_FCOP_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
        tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
        
        masktime=find(squeeze(any(mean(tmp.mask(:,1:end,:),2),1)));
        if ~isempty(masktime)
          chansel=tmp.label(find(ceil(mean(mean(tmp.mask(:,:,dsearchn(tmp.time',tmp.time(masktime(1))):dsearchn(tmp.time',tmp.time(masktime(end)))),2),3))));
          figinds=[100*ll+9;  100*ll+10+9];
          figstrings{1}=[fdir 'plvabslo_synch_all_allsig_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          figstrings{2}=[fdir 'plvanglo_synch_all_allsig_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
          tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
        end
        
        % theta topo
        if ~isempty(find(squeeze(any(mean(tmp.mask(:,1:2,:),2),1))))
          cfg=[];
          if ll==3 && adda==1
            cfg.ylim=[4 6.5];
          elseif ll==3 && adda==2
            cfg.ylim=[5.5 6.5];
          else
            keyboard
          end
          disp('get ylim right per ll theta synch plv')
          masktime=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1)));
          cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
          cfg.parameter='plvavgabs';
          cfg.layout='elec1010.lay';
          cfg.maskalpha=0.5;
          cfg.zlim=[-.25 .25];
          cfg.highlight='on';
          cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2))),2),3))));
          cfg.baseline=baseline3;
          cfg.zlim=[-.25 .25];
          figure(100*ll+8);
          ft_topoplotTFR(cfg,tmpuA);
          print(100*ll+8,[fdir 'plvabslo_synch_topoU_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          figure(100*ll+3);
          ft_topoplotTFR(cfg,tmpsA);
          print(100*ll+3,[fdir 'plvabslo_synch_topoM_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          cfg.zlim=[-.15 .15];
          figure(100*ll+4);
          ft_topoplotTFR(cfg,tmpA);
          print(100*ll+4,[fdir 'plvabslo_synch_topoDiff_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          
          if ~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2})))
            chansel=cfg.highlightchannel;
            zlim=[0 0.6; -.25 .25; -20 400; -10 300];
            figinds=[100*ll+20+1;  100*ll+30+1];
            figstrings{1}=[fdir 'plvabslo_synch_final_Xtheta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
            figstrings{2}=[fdir 'plvanglo_synch_final_Xtheta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
            tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
          end
          
          cfg.parameter='plvavgangwrap';
          if ll==3 && adda==1
            cfg.xlim=[0.19 0.19];
          elseif ll==3 && adda==2
            cfg.xlim=[0 0];
          else
            keyboard
          end
          disp('get xlim right per ll theta synch plv')
          cfg.zlim=[-4 4];
          figure(100*ll+8+10);
          ft_topoplotTFR(cfg,tmpuA);
          print(100*ll+8+10,[fdir 'plvanglo_synch_topoU_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          figure(100*ll+3+10);
          ft_topoplotTFR(cfg,tmpsA);
          print(100*ll+3+10,[fdir 'plvanglo_synch_topoM_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          figure(100*ll+4+10);
          ft_topoplotTFR(cfg,tmpA);
          print(100*ll+4+10,[fdir 'plvanglo_synch_topoDiff_theta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        end
        
        % alpha topo
        if ~isempty(find(squeeze(any(mean(tmp.mask(:,3:5,:),2),1))))
          cfg=[];
          if ll==3 && adda==1
            cfg.ylim=[8 9];
          else
            keyboard
          end
          disp('get ylim right per ll alpha synch plv')
          masktime=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1)));
          cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
          cfg.parameter='plvavgabs';
          cfg.layout='elec1010.lay';
          cfg.maskalpha=0.5;
          %     cfg.zlim='maxabs';
          cfg.highlight='on';
          cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2)) ),2),3))));
          cfg.baseline=baseline3;
          cfg.zlim=[-.15 .15];
          figure(100*ll+5);
          ft_topoplotTFR(cfg,tmpuA);
          print(100*ll+5,[fdir 'plvabslo_synch_topoU_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          figure(100*ll+6);
          ft_topoplotTFR(cfg,tmpsA);
          print(100*ll+6,[fdir 'plvabslo_synch_topoM_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          cfg.zlim=[-.1 .1];
          figure(100*ll+7);
          ft_topoplotTFR(cfg,tmpA);
          print(100*ll+7,[fdir 'plvabslo_synch_topoDiff_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          
          if ~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2})))
            chansel=cfg.highlightchannel;
            %             zlim=[-.3 .3; 0 0.35; -30 30; -60 60];
            zlim=[0 0.6; -.25 .25; -20 400; -10 300];
            figinds=[100*ll+40+1;  100*ll+50+1];
            figstrings{1}=[fdir 'plvabslo_synch_final_Xalpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
            figstrings{2}=[fdir 'plvanglo_synch_final_Xalpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
            tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
          end
          
          cfg.parameter='plvavgangwrap';
          if ll==3 && adda==1
            cfg.xlim=[0.19 0.19];
          else
            keyboard
          end
          disp('get xlim right per ll alpha synch plv')
          cfg.zlim=[-4 4];
          figure(100*ll+5+10);
          ft_topoplotTFR(cfg,tmpuA);
          print(100*ll+5+10,[fdir 'plvanglo_synch_topoU_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          figure(100*ll+6+10);
          ft_topoplotTFR(cfg,tmpsA);
          print(100*ll+6+10,[fdir 'plvanglo_synch_topoM_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          figure(100*ll+7+10);
          ft_topoplotTFR(cfg,tmpA);
          print(100*ll+7+10,[fdir 'plvanglo_synch_topoDiff_alpha_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        end
        
        % beta topo
        if ~isempty(find(squeeze(any(mean(tmp.mask(:,6:10,:),2),1))))
          cfg=[];
          if ll==6 && combval==1 && adda==2
            cfg.ylim=[13 15];
            %         elseif ll==3
            %           cfg.ylim=[13 20];
            %           keyboard
          else
            keyboard
          end
          disp('get ylim right per ll beta synch plv')
          masktime=find(squeeze(any(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),:),2),1)));
          cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
          cfg.parameter='plvavgabs';
          cfg.layout='elec1010.lay';
          cfg.maskalpha=0.5;
          %     cfg.zlim='maxabs';
          cfg.highlight='on';
          cfg.highlightchannel=tmp.label(find(ceil(mean(mean(tmp.mask(:,dsearchn(tmp.freq',cfg.ylim(1)):dsearchn(tmp.freq',cfg.ylim(end)),dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2)) ),2),3))));
          cfg.baseline=baseline;
          cfg.zlim=[-.15 .15];
          figure(100*ll+5);
          ft_topoplotTFR(cfg,tmpuA);
          print(100*ll+5,[fdir 'plvabslo_synch_topoU_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          figure(100*ll+6);
          ft_topoplotTFR(cfg,tmpsA);
          print(100*ll+6,[fdir 'plvabslo_synch_topoM_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          cfg.zlim=[-.1 .1];
          figure(100*ll+7);
          ft_topoplotTFR(cfg,tmpA);
          print(100*ll+7,[fdir 'plvabslo_synch_topoDiff_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          
          if ~isempty(setdiff(cfg.highlightchannel,union(chanplot{1},chanplot{2})))
            chansel=cfg.highlightchannel;
            zlim=[-.3 .3; 0 0.35; -30 30; -60 60];
            zlim=[0 0.6; -.25 .25; -20 400; -10 300];
            figinds=[100*ll+40+1;  100*ll+50+1];
            figstrings{1}=[fdir 'plvabslo_synch_final_Xbeta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
            figstrings{2}=[fdir 'plvanglo_synch_final_Xbeta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'];
            tfr_subchannel_3cond_plot(pow,chansel,zlim,figinds,figstrings)
          end
          
          cfg.parameter='plvavgangwrap';
          if ll==6 && combval==1 && adda==2
            cfg.xlim=[.32 .32];
          elseif ll==3 && combval==1 && adda==1
            cfg.xlim=[0 0];
          elseif ll==3 && adda==2
            cfg.xlim=[-.06 -.06];
          else
            keyboard
          end
          disp('get xlim right per ll beta synch plv')
          cfg.zlim=[-4 4];
          figure(100*ll+5+10);
          ft_topoplotTFR(cfg,tmpuA);
          print(100*ll+5+10,[fdir 'plvanglo_synch_topoU_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          figure(100*ll+6+10);
          ft_topoplotTFR(cfg,tmpsA);
          print(100*ll+6+10,[fdir 'plvanglo_synch_topoM_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
          figure(100*ll+7+10);
          ft_topoplotTFR(cfg,tmpA);
          print(100*ll+7+10,[fdir 'plvanglo_synch_topoDiff_beta_combval' num2str(combval) '_adda' num2str(adda) '_' num2str(ll) num2str(tt) num2str(ss) '.eps'],'-depsc2')
        end
        
      end % adda
      
    end % combval
  end  % if ll<5
  
end % ll

%% phaset0

% Q1 uniform distribution for unisensory at time 0?
%    Which method / channel selection / time point prior to 0 to use?
%
%  Answer:  FFT PCA at -100ms for 10Hz

% load elec1010_neighb.mat

% allcond_sameN=1;
% soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];
% tacbasemax=min(-.15,soades-.15);
% audbasemax=min(-.15,fliplr(soades)-.15);

% soadesc={'Aud first by 500ms' '' 'Aud first by 70ms' 'Aud first by 20ms' 'Simultaneous' 'Tac first by 20ms' 'Tac first by 70ms' '' 'Tac first by 500ms'};
plotflag=1;
printflag=1;
statsflag=1;
audtacflag=0;

soalist=[1 3 4 5 6 7 9];

% chanuse_sleep0={'all' '-F4'};
% chanuse_sleep1={'all' '-AF7' '-AF3' '-Fp1'};


for sleep=[0]
  %   if sleep
  %     chanuse=chanuse_sleep1;
  %   else
  %     chanuse=chanuse_sleep0;
  %   end
  for tt=[3]
    clearvars -except ll tt sub *dir ii*use sleep *flag figind soadesc soalist chanuse* ylim* neigh* *basemax
    
    if sleep
      subuseall=iiBuse;
      iter=11;
      ss=12;
    else
      %       subuseall=setdiff(iiSuse,[16:32])
      subuseall=iiSuse;
      iter=27;
      ss=10;
    end
    usetr=1;
    submin=subuseall(1)-1;
    
    
    phasedist_tac=nan(13,6,9,11,length(subuseall)); % freq x type x soalist x time x subject
    phasedist_aud=nan(13,6,9,11,length(subuseall));
    phasedist_nul=nan(13,6,9,11,length(subuseall));
    phasedist_ms1=nan(13,6,9,11,length(subuseall));
    phasedist_ms2=nan(13,6,9,11,length(subuseall));
    
    subuseind=0;
    for ii=subuseall
      subuseind=subuseind+1;
      cd([edir sub{ii} ])
      load(['freqtime0_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_tt' num2str(tt) '.mat'])
      load(['tlock_unimultsensHilbert_' sub{ii} '_sleep0_tt' num2str(tt) '_tacaud1_iter' num2str(iter) '_trialkc-1.mat'])
      for ll=soalist
        % question 1: uniform distribution at time 0?
        phasedist_tac(:,1,ll,:,subuseind)=squeeze(mean(exp(i*angle(freqlo_tTacAlone_time0{ll,tt,ss}.fourierspctrm(:,18,:,:))),1));
        phasedist_tac(:,2,ll,:,subuseind)=squeeze(mean(exp(i*angle(freqlo_tTacAlone_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,:,:))),1));
        phasedist_tac(:,3,ll,:,subuseind)=squeeze(mean(exp(i*angle(freqlo_tTacAlone_chanPCA_time0{ll,tt,ss}.fourierspctrm(:,1,:,:))),1));
        phasedist_aud(:,1,ll,:,subuseind)=squeeze(mean(exp(i*angle(freqlo_tAudAlone_time0{ll,tt,ss}.fourierspctrm(:,18,:,:))),1));
        phasedist_aud(:,2,ll,:,subuseind)=squeeze(mean(exp(i*angle(freqlo_tAudAlone_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,:,:))),1));
        phasedist_aud(:,3,ll,:,subuseind)=squeeze(mean(exp(i*angle(freqlo_tAudAlone_chanPCA_time0{ll,tt,ss}.fourierspctrm(:,1,:,:))),1));
        phasedist_nul(:,1,ll,:,subuseind)=squeeze(mean(exp(i*angle(freqlo_tNulAlone_time0{ll,tt,ss}.fourierspctrm(:,18,:,:))),1));
        phasedist_nul(:,2,ll,:,subuseind)=squeeze(mean(exp(i*angle(freqlo_tNulAlone_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,:,:))),1));
        phasedist_nul(:,3,ll,:,subuseind)=squeeze(mean(exp(i*angle(freqlo_tNulAlone_chanPCA_time0{ll,tt,ss}.fourierspctrm(:,1,:,:))),1));
        phasedist_ms1(:,1,ll,:,subuseind)=squeeze(mean(exp(i*angle(freqlo_tMSAlone_first_time0{ll,tt,ss}.fourierspctrm(:,18,:,:))),1));
        phasedist_ms1(:,2,ll,:,subuseind)=squeeze(mean(exp(i*angle(freqlo_tMSAlone_first_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,:,:))),1));
        phasedist_ms1(:,3,ll,:,subuseind)=squeeze(mean(exp(i*angle(freqlo_tMSAlone_first_chanPCA_time0{ll,tt,ss}.fourierspctrm(:,1,:,:))),1));
        phasedist_ms2(:,1,ll,:,subuseind)=squeeze(mean(exp(i*angle(freqlo_tMSAlone_second_time0{ll,tt,ss}.fourierspctrm(:,18,:,:))),1));
        phasedist_ms2(:,2,ll,:,subuseind)=squeeze(mean(exp(i*angle(freqlo_tMSAlone_second_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,:,:))),1));
        phasedist_ms2(:,3,ll,:,subuseind)=squeeze(mean(exp(i*angle(freqlo_tMSAlone_second_chanPCA_time0{ll,tt,ss}.fourierspctrm(:,1,:,:))),1));
        
        phasedist_tac(1:4,4,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_Tac{ll,tt,ss}(:,18,1,:)))))',[4 1]);
        phasedist_tac(5:8,4,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_Tac{ll,tt,ss}(:,18,2,:)))))',[4 1]);
        phasedist_tac(9:12,4,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_Tac{ll,tt,ss}(:,18,3,:)))))',[4 1]);
        phasedist_aud(1:4,4,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_Aud{ll,tt,ss}(:,18,1,:)))))',[4 1]);
        phasedist_aud(5:8,4,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_Aud{ll,tt,ss}(:,18,2,:)))))',[4 1]);
        phasedist_aud(9:12,4,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_Aud{ll,tt,ss}(:,18,3,:)))))',[4 1]);
        phasedist_nul(1:4,4,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_Nul{ll,tt,ss}(:,18,1,:)))))',[4 1]);
        phasedist_nul(5:8,4,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_Nul{ll,tt,ss}(:,18,2,:)))))',[4 1]);
        phasedist_nul(9:12,4,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_Nul{ll,tt,ss}(:,18,3,:)))))',[4 1]);
        phasedist_ms1(1:4,4,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_MS{ll,tt,ss}(:,18,1,:)))))',[4 1]);
        phasedist_ms1(5:8,4,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_MS{ll,tt,ss}(:,18,2,:)))))',[4 1]);
        phasedist_ms1(9:12,4,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_MS{ll,tt,ss}(:,18,3,:)))))',[4 1]);
        phasedist_ms2(1:4,4,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_secondstim_MS{ll,tt,ss}(:,18,1,:)))))',[4 1]);
        phasedist_ms2(5:8,4,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_secondstim_MS{ll,tt,ss}(:,18,2,:)))))',[4 1]);
        phasedist_ms2(9:12,4,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_secondstim_MS{ll,tt,ss}(:,18,3,:)))))',[4 1]);
        
        phasedist_tac(1:4,5,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanFC_Tac{ll,tt,ss}(:,1,1,:)))))',[4 1]);
        phasedist_tac(5:8,5,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanFC_Tac{ll,tt,ss}(:,1,2,:)))))',[4 1]);
        phasedist_tac(9:12,5,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanFC_Tac{ll,tt,ss}(:,1,3,:)))))',[4 1]);
        phasedist_aud(1:4,5,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanFC_Aud{ll,tt,ss}(:,1,1,:)))))',[4 1]);
        phasedist_aud(5:8,5,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanFC_Aud{ll,tt,ss}(:,1,2,:)))))',[4 1]);
        phasedist_aud(9:12,5,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanFC_Aud{ll,tt,ss}(:,1,3,:)))))',[4 1]);
        phasedist_nul(1:4,5,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanFC_Nul{ll,tt,ss}(:,1,1,:)))))',[4 1]);
        phasedist_nul(5:8,5,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanFC_Nul{ll,tt,ss}(:,1,2,:)))))',[4 1]);
        phasedist_nul(9:12,5,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanFC_Nul{ll,tt,ss}(:,1,3,:)))))',[4 1]);
        phasedist_ms1(1:4,5,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanFC_MS{ll,tt,ss}(:,1,1,:)))))',[4 1]);
        phasedist_ms1(5:8,5,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanFC_MS{ll,tt,ss}(:,1,2,:)))))',[4 1]);
        phasedist_ms1(9:12,5,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanFC_MS{ll,tt,ss}(:,1,3,:)))))',[4 1]);
        phasedist_ms2(1:4,5,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_secondstim_chanFC_MS{ll,tt,ss}(:,1,1,:)))))',[4 1]);
        phasedist_ms2(5:8,5,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_secondstim_chanFC_MS{ll,tt,ss}(:,1,2,:)))))',[4 1]);
        phasedist_ms2(9:12,5,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_secondstim_chanFC_MS{ll,tt,ss}(:,1,3,:)))))',[4 1]);
        
        phasedist_tac(1:4,6,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanPCA_Tac{ll,tt,ss}(:,1,1,:)))))',[4 1]);
        phasedist_tac(5:8,6,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanPCA_Tac{ll,tt,ss}(:,1,2,:)))))',[4 1]);
        phasedist_tac(9:12,6,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanPCA_Tac{ll,tt,ss}(:,1,3,:)))))',[4 1]);
        phasedist_aud(1:4,6,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanPCA_Aud{ll,tt,ss}(:,1,1,:)))))',[4 1]);
        phasedist_aud(5:8,6,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanPCA_Aud{ll,tt,ss}(:,1,2,:)))))',[4 1]);
        phasedist_aud(9:12,6,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanPCA_Aud{ll,tt,ss}(:,1,3,:)))))',[4 1]);
        phasedist_nul(1:4,6,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanPCA_Nul{ll,tt,ss}(:,1,1,:)))))',[4 1]);
        phasedist_nul(5:8,6,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanPCA_Nul{ll,tt,ss}(:,1,2,:)))))',[4 1]);
        phasedist_nul(9:12,6,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanPCA_Nul{ll,tt,ss}(:,1,3,:)))))',[4 1]);
        phasedist_ms1(1:4,6,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanPCA_MS{ll,tt,ss}(:,1,1,:)))))',[4 1]);
        phasedist_ms1(5:8,6,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanPCA_MS{ll,tt,ss}(:,1,2,:)))))',[4 1]);
        phasedist_ms1(9:12,6,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_firststim_chanPCA_MS{ll,tt,ss}(:,1,3,:)))))',[4 1]);
        phasedist_ms2(1:4,6,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_secondstim_chanPCA_MS{ll,tt,ss}(:,1,1,:)))))',[4 1]);
        phasedist_ms2(5:8,6,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_secondstim_chanPCA_MS{ll,tt,ss}(:,1,2,:)))))',[4 1]);
        phasedist_ms2(9:12,6,ll,:,subuseind)=repmat(squeeze(mean(exp(i*angle(hilb_time0_secondstim_chanPCA_MS{ll,tt,ss}(:,1,3,:)))))',[4 1]);
        
        
        
        
        %         plv_persub_tac.label=freqlo_tAudAlone_time0{ll,tt,ss}.label;
        %         plv_persub_tac.freq=freqlo_tAudAlone_time0{ll,tt,ss}.freq;
        %         plv_persub_tac.dimord=freqlo_tAudAlone_time0{ll,tt,ss}.dimord;
        %         plv_persub_tac.time=1:9;
        %         plv_persub_aud=plv_persub_tac;
        %         plv_persub_nul=plv_persub_tac;
        %         plv_persub_ms1=plv_persub_tac;
        %         plv_persub_ms2=plv_persub_tac;
        %           freqloall_tAudAlone_time0{ll,subuseind}      =freqlo_tAudAlone_time0{ll,tt,ss}.fourierspctrm;
        %           freqloall_tNulAlone_time0{ll,subuseind}      =freqlo_tNulAlone_time0{ll,tt,ss}.fourierspctrm;
        %           freqloall_tTacAlone_time0{ll,subuseind}      =freqlo_tTacAlone_time0{ll,tt,ss}.fourierspctrm;
        %           freqloall_tMSAlone_first_time0{ll,subuseind} =freqlo_tMSAlone_first_time0{ll,tt,ss}.fourierspctrm;
        %           freqloall_tMSAlone_second_time0{ll,subuseind}=freqlo_tMSAlone_second_time0{ll,tt,ss}.fourierspctrm;
        %
        %           plv_persub_tac.fourierspctrm(subuseind,:,:,ll)=squeeze(mean(exp(i*angle(freqloall_tTacAlone_time0{ll,subuseind})),1));
        %           plv_persub_aud.fourierspctrm(subuseind,:,:,ll)=squeeze(mean(exp(i*angle(freqloall_tAudAlone_time0{ll,subuseind})),1));
        %           plv_persub_nul.fourierspctrm(subuseind,:,:,ll)=squeeze(mean(exp(i*angle(freqloall_tNulAlone_time0{ll,subuseind})),1));
        %           plv_persub_ms1.fourierspctrm(subuseind,:,:,ll) =squeeze(mean(exp(i*angle(freqloall_tMSAlone_first_time0{ll,subuseind})),1));
        %           plv_persub_ms2.fourierspctrm(subuseind,:,:,ll) =squeeze(mean(exp(i*angle(freqloall_tMSAlone_second_time0{ll,subuseind})),1));
        %
        %
        %           cfg=[];
        %           cfg.ylim=[10 10];
        %           cfg.layout='eeg1010.lay';
        %           cfg.parameter='fourierspctrm';
        %           ft_topoplotTFR(cfg,freqlo_tAudAlone_time0{ll,tt,ss});
        
      end % ll
      if plotflag
        
        if 0
          for yy=1:size(phasedist_ms2,3) % index of time ind relative to time=0 [..... -100ms -50ms 0ms]
            figure(ii+100*yy);
            for ll=1:length(soalist)
              subplot(5,7,ll);plot(abs(phasedist_nul(:,:,soalist(ll),yy,subuseind)));axis([0 15 0 0.8])
              switch soalist(ll)
                case 1
                  title('AT500')
                case 3
                  title('AT70')
                case 4
                  title('AT20')
                case 5
                  title('AT0')
                case 6
                  title('TA20')
                case 7
                  title('TA70')
                case 9
                  title('TA500')
              end
              if ll==1,ylabel('Null');end
              subplot(5,7,ll+7);plot(abs(phasedist_tac(:,:,soalist(ll),yy,subuseind)));axis([0 15 0 0.8])
              if ll==1,ylabel('TacAlone');end
              subplot(5,7,ll+14);plot(abs(phasedist_aud(:,:,soalist(ll),yy,subuseind)));axis([0 15 0 0.8])
              if ll==1,ylabel('AudAlone');end
              subplot(5,7,ll+21);plot(abs(phasedist_ms1(:,:,soalist(ll),yy,subuseind)));axis([0 15 0 0.8])
              if ll==1,ylabel('MS first stim');end
              subplot(5,7,ll+28);plot(abs(phasedist_ms2(:,:,soalist(ll),yy,subuseind)));axis([0 15 0 0.8])
              if ll==1,ylabel('MS second stim');end
            end
            legend({'FFT Cz' 'FFT FrontCent' 'FFT PCA' 'Hilb Cz'})
            
            figure(2000+ii+100*yy);
            for ll=1:length(soalist)
              subplot(5,7,ll);plot(angle(phasedist_nul(:,:,soalist(ll),yy,subuseind)),'o');axis([0 15 -pi pi])
              switch soalist(ll)
                case 1
                  title('AT500')
                case 3
                  title('AT70')
                case 4
                  title('AT20')
                case 5
                  title('AT0')
                case 6
                  title('TA20')
                case 7
                  title('TA70')
                case 9
                  title('TA500')
              end
              if ll==1,ylabel('Null');end
              subplot(5,7,ll+7);plot(angle(phasedist_tac(:,:,soalist(ll),yy,subuseind)),'o');axis([0 15 -pi pi])
              if ll==1,ylabel('TacAlone');end
              subplot(5,7,ll+14);plot(angle(phasedist_aud(:,:,soalist(ll),yy,subuseind)),'o');axis([0 15 -pi pi])
              if ll==1,ylabel('AudAlone');end
              subplot(5,7,ll+21);plot(angle(phasedist_ms1(:,:,soalist(ll),yy,subuseind)),'o');axis([0 15 -pi pi])
              if ll==1,ylabel('MS first stim');end
              subplot(5,7,ll+28);plot(angle(phasedist_ms2(:,:,soalist(ll),yy,subuseind)),'o');axis([0 15 -pi pi])
              if ll==1,ylabel('MS second stim');end
            end
            legend({'FFT Cz' 'FFT FrontCent' 'FFT PCA' 'Hilb Cz'})
          end % yy
        end
        
        for yy=1:size(phasedist_ms2,4) % index of time ind relative to time=0 [..... -100ms -50ms 0ms]
          figure(ii);
          for ll=4 % soalist(ll)
            subplot(5,11,yy);plot(abs(phasedist_nul(:,:,soalist(ll),yy,subuseind)));axis([0 15 0 0.8])
            title([num2str((yy-11)*50) ' ms'])
            if yy==1,ylabel('Null');end
            subplot(5,11,yy+11);plot(abs(phasedist_tac(:,:,soalist(ll),yy,subuseind)));axis([0 15 0 0.8])
            if yy==1,ylabel('TacAlone');end
            subplot(5,11,yy+22);plot(abs(phasedist_aud(:,:,soalist(ll),yy,subuseind)));axis([0 15 0 0.8])
            if yy==1,ylabel('AudAlone');end
            subplot(5,11,yy+33);plot(abs(phasedist_ms1(:,:,soalist(ll),yy,subuseind)));axis([0 15 0 0.8])
            if yy==1,ylabel('MS first stim');end
            subplot(5,11,yy+44);plot(abs(phasedist_ms2(:,:,soalist(ll),yy,subuseind)));axis([0 15 0 0.8])
            if yy==1,ylabel('MS second stim');end
          end
        end % yy
        legend({'FFT Cz' 'FFT FrontCent' 'FFT PCA' 'Hilb Cz' 'Hilb FC' 'Hilb PCA'})
        
        for yy=1:size(phasedist_ms2,4) % index of time ind relative to time=0 [..... -100ms -50ms 0ms]
          figure(ii+100);
          for ll=4 % soalist(ll)
            subplot(5,11,yy);plot(angle(phasedist_nul(:,:,soalist(ll),yy,subuseind)),'o');axis([0 15 -pi pi])
            title([num2str((yy-11)*50) ' ms'])
            if yy==1,ylabel('Null');end
            subplot(5,11,yy+11);plot(angle(phasedist_tac(:,:,soalist(ll),yy,subuseind)),'o');axis([0 15 -pi pi])
            if yy==1,ylabel('TacAlone');end
            subplot(5,11,yy+22);plot(angle(phasedist_aud(:,:,soalist(ll),yy,subuseind)),'o');axis([0 15 -pi pi])
            if yy==1,ylabel('AudAlone');end
            subplot(5,11,yy+33);plot(angle(phasedist_ms1(:,:,soalist(ll),yy,subuseind)),'o');axis([0 15 -pi pi])
            if yy==1,ylabel('MS first stim');end
            subplot(5,11,yy+44);plot(angle(phasedist_ms2(:,:,soalist(ll),yy,subuseind)),'o');axis([0 15 -pi pi])
            if yy==1,ylabel('MS second stim');end
          end
        end % yy
        legend({'FFT Cz' 'FFT FrontCent' 'FFT PCA' 'Hilb Cz' 'Hilb FC' 'Hilb PCA'})
        
      end % plotflag
      
    end % ii
    
    %     htac=nan(9,12,4);ptac=htac;for ll=soalist,for ff=1:12,for yy=1:4,[htac(ll,ff,yy),ptac(ll,ff,yy)]=ttest(phasedist_tac(ff,yy,ll,:),phasedist_nul(ff,yy,ll,:));end;end;end
    %     haud=nan(9,12,4);paud=haud;for ll=soalist,for ff=1:12,for yy=1:4,[haud(ll,ff,yy),paud(ll,ff,yy)]=ttest(phasedist_aud(ff,yy,ll,:),phasedist_nul(ff,yy,ll,:));end;end;end
    %     hms1=nan(9,12,4);pms1=hms1;for ll=soalist,for ff=1:12,for yy=1:4,[hms1(ll,ff,yy),pms1(ll,ff,yy)]=ttest(phasedist_ms1(ff,yy,ll,:),phasedist_nul(ff,yy,ll,:));end;end;end
    %     hms2=nan(9,12,4);pms2=hms2;for ll=soalist,for ff=1:12,for yy=1:4,[hms2(ll,ff,yy),pms2(ll,ff,yy)]=ttest(phasedist_ms2(ff,yy,ll,:),phasedist_nul(ff,yy,ll,:));end;end;end
    
    % Test if different from Nul; if so, then biased, but if not, then good.
    %  Thus, which is most similar to (i.e. most not-significantly-different from) nul, especially in alpha?
    htac=nan(9,12,6,11);ptac=htac;for ll=soalist,for ff=1:12,for yy=1:6,for tm=1:11,try [htac(ll,ff,yy,tm),ptac(ll,ff,yy,tm)]=ttest(abs(phasedist_tac(ff,yy,ll,tm,:)),abs(phasedist_nul(ff,yy,ll,tm,:)),'alpha',.01);end;end;end;end;end
    haud=nan(9,12,6,11);paud=haud;for ll=soalist,for ff=1:12,for yy=1:6,for tm=1:11,try [haud(ll,ff,yy,tm),paud(ll,ff,yy,tm)]=ttest(abs(phasedist_aud(ff,yy,ll,tm,:)),abs(phasedist_nul(ff,yy,ll,tm,:)),'alpha',.01);end;end;end;end;end
    hms1=nan(9,12,6,11);pms1=hms1;for ll=soalist,for ff=1:12,for yy=1:6,for tm=1:11,try [hms1(ll,ff,yy,tm),pms1(ll,ff,yy,tm)]=ttest(abs(phasedist_ms1(ff,yy,ll,tm,:)),abs(phasedist_nul(ff,yy,ll,tm,:)),'alpha',.01);end;end;end;end;end
    hms2=nan(9,12,6,11);pms2=hms2;for ll=soalist,for ff=1:12,for yy=1:6,for tm=1:11,try [hms2(ll,ff,yy,tm),pms2(ll,ff,yy,tm)]=ttest(abs(phasedist_ms2(ff,yy,ll,tm,:)),abs(phasedist_nul(ff,yy,ll,tm,:)),'alpha',.01);end;end;end;end;end
    
    if plotflag
      figure;bar(squeeze(any(ptac(soalist,2,:,:)<.05,1))+squeeze(any(paud(soalist,2,:,:)<.05,1))+squeeze(any(pms1(soalist,2,:,:)<.05,1)))
      legend({'-500ms' '-450ms' '-400ms' '-350ms' '-300ms' '-250ms' '-200ms' '-150ms' '-100ms' '-50ms' '-0ms'}); title('2 Hz')
      figure;bar(squeeze(any(ptac(soalist,6,:,:)<.05,1))+squeeze(any(paud(soalist,6,:,:)<.05,1))+squeeze(any(pms1(soalist,6,:,:)<.05,1)))
      legend({'-500ms' '-450ms' '-400ms' '-350ms' '-300ms' '-250ms' '-200ms' '-150ms' '-100ms' '-50ms' '-0ms'}) ; title('6 Hz')
      figure;bar(squeeze(any(ptac(soalist,10,:,:)<.05,1))+squeeze(any(paud(soalist,10,:,:)<.05,1))+squeeze(any(pms1(soalist,10,:,:)<.05,1)))
      legend({'-500ms' '-450ms' '-400ms' '-350ms' '-300ms' '-250ms' '-200ms' '-150ms' '-100ms' '-50ms' '-0ms'}) ; title('10 Hz')
    end
    
    
    
  end % tt
end % sleep


%%  Sort ERP based on phase from above analysis

plotflag=1;
printflag=1;
statsflag=1;
audtacflag=0;

soalist=[1 3 4 5 6 7 9];

% chanuse_sleep0={'all' '-F4'};
% chanuse_sleep1={'all' '-AF7' '-AF3' '-Fp1'};


for sleep=[0]
  %   if sleep
  %     chanuse=chanuse_sleep1;
  %   else
  %     chanuse=chanuse_sleep0;
  %   end
  for tt=[3]
    clearvars -except ll tt sub *dir ii*use sleep *flag figind soadesc soalist chanuse* ylim* neigh* *basemax
    
    if sleep
      subuseall=iiBuse;
      iter=11;
      ss=12;
    else
      %       subuseall=setdiff(iiSuse,[16:32])
      subuseall=iiSuse;
      iter=27;
      ss=10;
    end
    usetr=1;
    submin=subuseall(1)-1;
    
    
    subuseind=0;
    for ii=subuseall
      subuseind=subuseind+1;
      cd([edir sub{ii} ])
      load(['freqtime0_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_tt' num2str(tt) '.mat'])
      load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud1_iter' num2str(iter) '.mat']);
      
      for ll=soalist
        angle(freqlo_tAudAlone_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,10,9))
      end % ll
    end % ii
  end % tt
end % sleep