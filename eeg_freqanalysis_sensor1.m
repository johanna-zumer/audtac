% preprocessing of EEG data from Hills, 64ch MR-cap

eeg_legomagic_preamble

%%

minnumcomb=1; % in effect just the normal combination of trials
minnumcomb=2; % don't let minnumcomb go higher than 1000.
% minnumcomb=100; % don't let minnumcomb go higher than 1000.

timestepfft=0.02;
timestepplv=0.01;
toibeg=min(soades,0)-.7;
toiend=max(soades,0)+1.3;

dofftadd=0;
usetr=1; % see inside eeg_legomagic_trialSelection2_wakeSleep_sepTacAud for what this does/means
synchasynch=0;
use23=0;
phaset0=0;

% for sleep=0
sleep=1;
if sleep
  iiuse=iiBuse;
  %     iiuse=[32];
  iteruse=11;
  trialkc=-1;  % obsolete (used to vary this from -1, 0, and 1) but now use trlfeat/trlkeep instead
else
  iiuse=iiSuse;
  %         iiuse=setdiff(iiSuse,1:30);
  iteruse=31;  % iteruse=31 with usetr3 for final iter31; iteruse=31 with usetr2 for final iter32
  % 1 june 2016, re-run iter 31/32 with updated cfg.toi
  iteruse=27; % final use; keep it simple
  trialkc=-1;
end
for ii=iiuse;
  % for ii=setdiff(iiuse,1:8)
%   for ii=8
  cd([edir sub{ii} ])
  clearvars -except ii sub *dir ii*use sleep minnumcomb hostname timestep* soades dofftadd statst use* synchasynch iteruse trialkc phaset0 toibeg toiend  subuseall trlkeep*
  [raw_tac, raw_aud, raw_nul, artflag]=eeg_legomagic_epoching2(ii,sleep,1,0); % featfull=1, saveflag =0;
  
  %% Do TacPlusAud first
  tacaud=1;
  %   for tt=[3]
  tt=3;
  for iter=iteruse
    %       [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection2_wakeSleep_sepTacAud(ii,sleep,tt,tacaud);
    [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection2_wakeSleep_sepTacAud(ii,sleep,tt,tacaud,raw_tac, raw_aud, raw_nul, iter, usetr, trialkc);
    
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
%         ssuse=[tr.stageuse+10];
        ssuse=12;
      end
      load(['trlfeat_' sub{ii} '_iter' num2str(iteruse) '.mat']);      
    else
      tlock_aud{10,tt,12}=[];
      tlock_aud{10,tt,13}=[];
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
              cfg.toi=toibeg(ll):timestepfft:toiend(ll);
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
              cfg.toi=toibeg(ll):timestepfft:toiend(ll);
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
              cfg.toi=toibeg(ll):timestepfft:toiend(ll);
              cfg.t_ftimwin=4./cfg.foi;
              cfg.output='powandcsd';
              freqlo_aud{ll,tt,ss}=ft_freqanalysis(cfg,tlock_aud{ll,tt,ss});
              cfg=[];
              cfg.method='mtmconvol';
              cfg.pad=4;
              cfg.foi=[30:5:45 55:5:80];
              cfg.taper='dpss';
              cfg.tapsmofrq=7*ones(1,length(cfg.foi));
              cfg.toi=toibeg(ll):timestepfft:toiend(ll);
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
              cfg.toi=toibeg(ll):timestepfft:toiend(ll);
              cfg.t_ftimwin=4./cfg.foi;
              cfg.output='powandcsd';
              freqlo_nul{ll,tt,ss}=ft_freqanalysis(cfg,tlock_nul{ll,tt,ss});
              cfg=[];
              cfg.method='mtmconvol';
              cfg.pad=4;
              cfg.foi=[30:5:45 55:5:80];
              cfg.taper='dpss';
              cfg.tapsmofrq=7*ones(1,length(cfg.foi));
              cfg.toi=toibeg(ll):timestepfft:toiend(ll);
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
                cfg.toi=toibeg(ll):timestepfft:toiend(ll);
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
                cfg.toi=toibeg(ll):timestepfft:toiend(ll);
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
            if ll==1 % time of first stimulus
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
            if ll==1 || ll==3 || ll==4 || ll==5  % time of second stimulus
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
          cfg.toi=toibeg(ll):timestepplv:toiend(ll);
          cfg.t_ftimwin=4./cfg.foi;
          %         cfg.keeptrials='yes';  % not necessary at this stage for keeping trials (yes later for source analysis and phase resetting)
          cfg.output='fourier';
          freqlo_tTacAlone_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tac{ll+20,tt,ss});
          freqlo_tAudAlone_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_aud{ll+40,tt,ss});
          freqlo_tMSAlone_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tac{ll,tt,ss});
          freqlo_tNulAlone_comb{ll,tt,ss}=ft_freqanalysis(cfg,tlock_nul{ll+50,tt,ss});
          
          %           disp('study getplv'), keyboard
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
          cfg.toi=toibeg(ll):timestepplv:toiend(ll);
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
          keyboard
          cfg=[];cfg.randomseed=13;ft_preamble randomseed
          for cc=1:numtests
            %               tmp=Shuffle(combindex(:));
            %               combuse=tmp(1:numt_trials(ll,tt,ss));
            allnumtr=1:numt_trials(ll,tt,ss);
            combuse=Shuffle(allnumtr);
            
            if cc==1
              tlock_fake=addbeforeFFT(tlock_tac{ll+20,tt,ss},tlock_aud{ll+40,tt,ss},allnumtr,allnumtr);
              tlock_fake_nKD=addbeforeFFT(tlock_tac{ll+20,tt,ss},tlock_aud{ll+40,tt,ss},trlkeep_tac_nKD{ll+20,ii},trlkeep_aud_nKD{ll+40,ii});
              tlock_fake_nSD=addbeforeFFT(tlock_tac{ll+20,tt,ss},tlock_aud{ll+40,tt,ss},trlkeep_tac_nSD{ll+20,ii},trlkeep_aud_nSD{ll+40,ii});
              tlock_fake_nKD_nSD=addbeforeFFT(tlock_tac{ll+20,tt,ss},tlock_aud{ll+40,tt,ss},trlkeep_tac_nKD_nSD{ll+20,ii},trlkeep_aud_nKD_nSD{ll+40,ii});
            else
              tlock_fake=addbeforeFFT(tlock_tac{ll+20,tt,ss},tlock_aud{ll+40,tt,ss},allnumtr,combuse);
            end
%             for at=1:numt_trials(ll,tt,ss)
%               if cc==1  % do it as 'normal'
%                 tind=at;aind=at;
%               else
%                 %                   [tind,aind]=find(combindex==combuse(at));
%                 tind=at;
%                 aind=combuse(at);
%                 %                   [tind,aind]=find(combindex==combuse(at));
%               end
%               cfg=[];
%               cfg.trials=tind;
%               tmpt=ft_selectdata(cfg,tlock_tac{ll+20,tt,ss});
%               cfg.trials=aind;
%               tmpa=ft_selectdata(cfg,tlock_aud{ll+40,tt,ss});
%               cfg=[];
%               cfg.operation='add';
%               cfg.parameter='trial';
%               tmpsum=ft_math(cfg,tmpt,tmpa);
%               if at==1
%                 tlock_fake=tmpsum;
%               end
%               tlock_fake.trial(at,:,:)=tmpsum.trial(1,:,:);
%             end
            
            cfg=[];
            cfg.method='mtmconvol';
            cfg.pad=4;
            cfg.foi=4:2:30;
            cfg.taper='hanning';
            cfg.toi=toibeg(ll):timestepplv:toiend(ll);
            cfg.t_ftimwin=4./cfg.foi;
            %         cfg.keeptrials='yes';  % not necessary at this stage for keeping trials (yes later for source analysis and phase resetting)
            cfg.output='fourier';
            freqlo_tacPaud_tmp{cc}= powNplv_freqanalysis(cfg,tlock_fake);
            if cc==1
              freqlo_tacPaud_nKD{ll,tt,ss}= powNplv_freqanalysis(cfg,tlock_fake_nKD);
              freqlo_tacPaud_nSD{ll,tt,ss}= powNplv_freqanalysis(cfg,tlock_fake_nSD);
              freqlo_tacPaud_nKD_nSD{ll,tt,ss}= powNplv_freqanalysis(cfg,tlock_fake_nKD_nSD);
            end            
            
            cfg=[];
            cfg.method='mtmconvol';
            cfg.pad=4;
            cfg.foi=[30:5:45 55:5:80];
            cfg.taper='dpss';
            cfg.tapsmofrq=7*ones(1,length(cfg.foi));
            cfg.toi=toibeg(ll):timestepplv:toiend(ll);
            cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
            cfg.output='fourier';
            freqhi_tacPaud_tmp{cc}= powNplv_freqanalysis(cfg,tlock_fake);
            if cc==1
              freqhi_tacPaud_nKD{ll,tt,ss}= powNplv_freqanalysis(cfg,tlock_fake_nKD);
              freqhi_tacPaud_nSD{ll,tt,ss}= powNplv_freqanalysis(cfg,tlock_fake_nSD);
              freqhi_tacPaud_nKD_nSD{ll,tt,ss}= powNplv_freqanalysis(cfg,tlock_fake_nKD_nSD);
            end            

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
          freqlo_tacMSpN_tmp{1}=[];
          freqhi_tacMSpN_tmp{1}=[];
          for cc=1:numtests
            allnumtr=1:numt_trials(ll,tt,ss);
            combuse=Shuffle(allnumtr);
            
            if cc==1
              tlock_fake=addbeforeFFT(tlock_tac{ll,tt,ss},tlock_nul{ll+50,tt,ss},allnumtr,allnumtr);
              tlock_fake_nKD=addbeforeFFT(tlock_tac{ll,tt,ss},tlock_nul{ll+50,tt,ss},trlkeep_tac_nKD{ll,ii},trlkeep_nul_nKD{ll+50,ii});
              tlock_fake_nSD=addbeforeFFT(tlock_tac{ll,tt,ss},tlock_nul{ll+50,tt,ss},trlkeep_tac_nSD{ll,ii},trlkeep_nul_nSD{ll+50,ii});
              tlock_fake_nKD_nSD=addbeforeFFT(tlock_tac{ll,tt,ss},tlock_nul{ll+50,tt,ss},trlkeep_tac_nKD_nSD{ll,ii},trlkeep_nul_nKD_nSD{ll+50,ii});
            else
              tlock_fake=addbeforeFFT(tlock_tac{ll,tt,ss},tlock_nul{ll+50,tt,ss},allnumtr,combuse);
            end
            
%             for at=1:numt_trials(ll,tt,ss)
%               if cc==1  % do it as 'normal'
%                 tind=at;aind=at;
%               else
%                 %                   [tind,aind]=find(combindex==combuse(at));
%                 tind=at;
%                 aind=combuse(at);
%               end
%               cfg=[];
%               cfg.trials=tind;
%               tmpt=ft_selectdata(cfg,tlock_tac{ll,tt,ss});
%               cfg.trials=aind;
%               tmpa=ft_selectdata(cfg,tlock_nul{ll+50,tt,ss});
%               cfg=[];
%               cfg.operation='add';
%               cfg.parameter='trial';
%               tmpsum=ft_math(cfg,tmpt,tmpa);
%               if at==1
%                 tlock_fake=tmpsum;
%               end
%               tlock_fake.trial(at,:,:)=tmpsum.trial(1,:,:);
%               %           tlock_fake.trial(at,:,:)=tlock_tac{ll,tt,ss}.trial(tind,:,:)+tlock_nul{ll+50,tt,ss}.trial(aind,:,:);
%             end
            
            cfg=[];
            cfg.method='mtmconvol';
            cfg.pad=4;
            cfg.foi=4:2:30;
            cfg.taper='hanning';
            cfg.toi=toibeg(ll):timestepplv:toiend(ll);
            cfg.t_ftimwin=4./cfg.foi;
            %         cfg.keeptrials='yes';  % not necessary at this stage for keeping trials (yes later for source analysis and phase resetting)
            cfg.output='fourier'; % this forces keeptrials to be 'yes'
            
            freqlo_tacMSpN_tmp{cc}= powNplv_freqanalysis(cfg,tlock_fake);
            if cc==1
              freqlo_tacMSpN_nKD{ll,tt,ss}= powNplv_freqanalysis(cfg,tlock_fake_nKD);
              freqlo_tacMSpN_nSD{ll,tt,ss}= powNplv_freqanalysis(cfg,tlock_fake_nSD);
              freqlo_tacMSpN_nKD_nSD{ll,tt,ss}= powNplv_freqanalysis(cfg,tlock_fake_nKD_nSD);
            end
            
            cfg=[];
            cfg.method='mtmconvol';
            cfg.pad=4;
            cfg.foi=[30:5:45 55:5:80];
            cfg.taper='dpss';
            cfg.tapsmofrq=7*ones(1,length(cfg.foi));
            cfg.toi=toibeg(ll):timestepplv:toiend(ll);
            cfg.t_ftimwin=.2*ones(1,length(cfg.foi));
            cfg.output='fourier';
            freqhi_tacMSpN_tmp{cc}= powNplv_freqanalysis(cfg,tlock_fake);
            if cc==1
              freqhi_tacMSpN_nKD{ll,tt,ss}= powNplv_freqanalysis(cfg,tlock_fake_nKD);
              freqhi_tacMSpN_nSD{ll,tt,ss}= powNplv_freqanalysis(cfg,tlock_fake_nSD);
              freqhi_tacMSpN_nKD_nSD{ll,tt,ss}= powNplv_freqanalysis(cfg,tlock_fake_nKD_nSD);
            end
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
              cfg.toi=toibeg(ll):timestepplv:toiend(ll);
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
              cfg.toi=toibeg(ll):timestepplv:toiend(ll);
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
              cfg.toi=toibeg(ll):timestepplv:toiend(ll);
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
              cfg.toi=toibeg(ll):timestepplv:toiend(ll);
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
                save(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_tt' num2str(tt) '_trialkc' num2str(trialkc) '.mat'],'freq*comb','stats*','num*trials','freq*_n*','-v7.3')
              else
                %                 save(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_tt' num2str(tt) '_trialkc' num2str(trialkc) '.mat'],'freq*comb','num*trials','-v7.3')
                save(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'freq*comb','num*trials','freq*_n*','-v7.3')
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
          save(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_tt' num2str(tt) '_trialkc' num2str(trialkc) '.mat'],'freq*comb','stats*','num*trials','freq*_n*','-v7.3')
        else
          %           save(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_tt' num2str(tt) '_trialkc' num2str(trialkc) '.mat'],'freq*comb','num*trials','-v7.3')
          save(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'freq*comb','num*trials','freq*_n*','-v7.3')
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
                cfg.toi=toibeg(ll):timestepfft:toiend(ll);
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
                cfg.toi=toibeg(ll):timestepfft:toiend(ll);
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
                cfg.toi=toibeg(ll):timestepfft:toiend(ll);
                cfg.t_ftimwin=4./cfg.foi;
                cfg.output='powandcsd';
                freqlo_tac{ll,tt,ss}=ft_freqanalysis(cfg,tlock_tac{ll,tt,ss});
                cfg=[];
                cfg.method='mtmconvol';
                cfg.pad=4;
                cfg.foi=[30:5:45 55:5:80];
                cfg.taper='dpss';
                cfg.tapsmofrq=7*ones(1,length(cfg.foi));
                cfg.toi=toibeg(ll):timestepfft:toiend(ll);
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
                cfg.toi=toibeg(ll):timestepfft:toiend(ll);
                cfg.t_ftimwin=4./cfg.foi;
                cfg.output='powandcsd';
                freqlo_nul{ll,tt,ss}=ft_freqanalysis(cfg,tlock_nul{ll,tt,ss});
                cfg=[];
                cfg.method='mtmconvol';
                cfg.pad=4;
                cfg.foi=[30:5:45 55:5:80];
                cfg.taper='dpss';
                cfg.tapsmofrq=7*ones(1,length(cfg.foi));
                cfg.toi=toibeg(ll):timestepfft:toiend(ll);
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
                cfg.toi=toibeg(ll):timestepfft:toiend(ll);
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
                cfg.toi=toibeg(ll):timestepfft:toiend(ll);
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
            cfg.toi=toibeg(ll):timestepplv:toiend(ll);
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
            cfg.toi=toibeg(ll):timestepplv:toiend(ll);
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
              cfg.toi=toibeg(ll):timestepplv:toiend(ll);
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
              cfg.toi=toibeg(ll):timestepplv:toiend(ll);
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
              cfg.toi=toibeg(ll):timestepplv:toiend(ll);
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
              cfg.toi=toibeg(ll):timestepplv:toiend(ll);
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
                cfg.toi=toibeg(ll):timestepplv:toiend(ll);
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
                cfg.toi=toibeg(ll):timestepplv:toiend(ll);
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
                cfg.toi=toibeg(ll):timestepplv:toiend(ll);
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
                cfg.toi=toibeg(ll):timestepplv:toiend(ll);
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






