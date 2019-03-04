% statistics of EEG awake data

eeg_legomagic_preamble

%%
usetr=1; % = 0 for use new random sampling for that iter;  = 1 use existing trial sampling if that iter has already been run;  = 2 if 'smart sampling' (including non-included trials of that existing iter)
ttuse=[3];
tacaloneproc=0;
phaset0=0; % compute phase at time 0
phaset0use=2; % 0 = compute Hilbert; 1 = FFT Cz, 2 = FFT FC, 3 = FFT PCA, 4 = Hilb Cz, 5 = Hilb FC, 6 = Hilb PCA
synchasynch=0;
binonly=0;
pcmflag=0;
plothist=0;


% for sleep=[1]
sleep=1;
if sleep
  iiuse=setdiff(iiBuse,3:7);
  %     iiuse=setdiff(iiBuse,[3:7 8]);
    trialkc=-1;  % CHANGE ME: -1 (use all), 0 (only non-Kc), or 1 (only Kc)
    iteruse=11;
else
  iiuse=iiSuse;
  %         iiuse=setdiff(iiSuse,1:9);
  %             iiuse=8;
    trialkc=-1; % always -1
    iteruse=27; % for main analysis
    %       iteruse=30; % for tacalone analysis
    % for new 'smart sampling' use 31 and 33
    %       iteruse=32; % smart sampling;  set to value of usetr=0 in combination with usetr=2 to get new paired sample (iteruse+1)
end
freqsub=nan(1,max(iiuse));
for ii=iiuse
%     for ii=setdiff(iiuse,[8:16])
%   for ii=16
  cd([edir sub{ii} ])
  clearvars -except ii sub edir ddir ii*use sleep featurestats* ttuse soades usetr tacaloneproc synchasynch phaset0* binonly pcmflag plothist iter* trialkc*
  [raw_tac, raw_aud, raw_nul, artflag]=eeg_legomagic_epoching2(ii,sleep,1,0); % featfull=1, saveflag =0;
  
  %   try
  %     load(['tlock_trialSel_' sub{ii} '.mat']);
  %   catch
  
  %     [tlock_tac_s0,tlock_aud_s0,tlock_nul_s0,tr]=eeg_legomagic_trialSelection1(ii);
  %     save(['tlock_trialSel_' sub{ii} '.mat'],'tlock*s0','tr','-v7.3');
  %     [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection1_wakeSleep(ii,sleep);

  % For memory reasons, do TacPlusAud separately from AudPlusTac
  %% Do TacPlusAud first
  tacaud=1; % tacaud=1 means triggered on tactile; tacaud=0 means triggered on auditory
  %     for tt=ttuse % refers to lim=[no-limit .01 .005 .002];
  tt=ttuse;
  %   profile on
  %       [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection2_wakeSleep_sepTacAud(ii,sleep,tt,tacaud);
  
  % % iter 1:10 same but with trialSelection eog of prestim (13) only
  % % iter 11:20 same but with trialSelection eog of all trial (13 and 23)
  % % iter 21:30 now 'epoching' uses  'artfctdef_sit_automan' rather than 'auto';
  %       for iter=iteruse
  iter=iteruse;
  if ~isempty(whos('raw*'))
    [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection2_wakeSleep_sepTacAud(ii,sleep,tt,tacaud,raw_tac, raw_aud, raw_nul, iter, usetr, trialkc);
  else
    error('with sleep=1 cant run more than one iter or one tt at once for memory reasons')
  end
  if usetr==2 && isfield(tr,'iter') && tr.iter~=iter
    if length(iteruse)>1
      error('cannot have iteruse more than scalar and usetr=2')
    end
    iter=tr.iter;
  end
  
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
  
  if plothist
    sl=1;
    for ll=[soalist 10]
      [tmp,mind]=min(squeeze(mean(tlock_aud{ll,3,12}.trial(:,17,:),1)))
      figure(ii);
      subplot(8,2,sl*2);hist(tlock_aud{ll,3,12}.trial(:,17,mind),20);xlim([-150 60]);ylim([0 100])
      sl=sl+1;
    end
    sl=1;
    for ll=[soalist 10]
      [tmp,mind]=min(squeeze(mean(tlock_tac{ll,3,12}.trial(:,17,:),1)))
      figure(ii);
      subplot(8,2,sl*2-1);hist(tlock_tac{ll,3,12}.trial(:,17,mind),20);xlim([-150 60]);ylim([0 100])
      sl=sl+1;
    end
    continue
  end
  
  if sleep
    ssuse=[tr.stageuse+10 23]; % not enough memory for 23, plus is it useful?
    %           ssuse=[tr.stageuse+10];
    ssuse=[12];
    clear raw*
  else
    ssuse=tr.stageuse+10;
    if pcmflag
      ssuse=10;
    end
  end
  
  
  for ss=ssuse
    try
      fsample=1/diff(tlock_tac{10,tt,ss}.time(1:2)); % sampling rate
    catch
      try
        fsample=1/diff(tlock_aud{10,tt,ss}.time(1:2)); % sampling rate
      catch
        try
          fsample=1/diff(tlock_nul{10,tt,ss}.time(1:2)); % sampling rate
        catch
          fsample=1/diff(tlock_nul{10,tt,12}.time(1:2)); % sampling rate
        end
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
            keyboard
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
    
    if phaset0 && sleep==0 && ss~=10
      continue
    elseif phaset0 && sleep==1 && ss~=12
      continue
    end
    
    if phaset0 && phaset0use==0
      % For PCA of data, hopefully picking up channels greatly contribution to evoked response
      % This is not dependent on 'll' but use all data
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
      if mxc<.6 || mind==1
        if ii==16
        else
          keyboard
        end
        %                   if     [ii==11 && ll==9 && ss==10 && sleep==0]
        %                   elseif [ii==26 && ll==9 && ss==10 && sleep==0]
        %                   elseif [ii==31 && ll==9 && ss==10 && sleep==0]
        %                     mind=2;
        %                   elseif [ii==8 && ll==3 && ss==10 && sleep==0]
        %                     mind=2;
        %                   else
        %                   end
      end
      montage     = [];
      montage.tra = pcatacaud_tlock.topo(:, mind)';
      montage.labelorg = pcatacaud.topolabel;
      montage.labelnew = pcatacaud.label(mind);
    end
    
    if phaset0 && phaset0use && [(sleep==0 && ss==10) || (sleep==1 && ss==12)]
      [angbin,absbin,freqtr]=eeg_legomagic_power_phase_bins(ii,sub,sleep,ss,trialkc,phaset0use);
      
      tacbin = angbin.tacbin;
      audbin = angbin.audbin;
      nulbin = angbin.nulbin;
      ms1bin = angbin.ms1bin;
      ms2bin = angbin.ms2bin;
      tac4mscon1bin = angbin.tac4mscon1bin;
      tac4mscon2bin = angbin.tac4mscon2bin;
      aud4mscon1bin = angbin.aud4mscon1bin;
      aud4mscon2bin = angbin.aud4mscon2bin;
      nul4mscon1bin = angbin.nul4mscon1bin;
      nul4mscon2bin = angbin.nul4mscon2bin;
      ms4mscon1bin = angbin.ms4mscon1bin;
      ms4mscon2bin = angbin.ms4mscon2bin;
      
      tacabsbin = absbin.tacabsbin;
      audabsbin = absbin.audabsbin;
      nulabsbin = absbin.nulabsbin;
      ms1absbin = absbin.ms1absbin;
      ms2absbin = absbin.ms2absbin;
      tac4mscon1absbin = absbin.tac4mscon1absbin;
      aud4mscon1absbin = absbin.aud4mscon1absbin;
      nul4mscon1absbin = absbin.nul4mscon1absbin;
      ms4mscon1absbin = absbin.ms4mscon1absbin;
      tac4mscon2absbin = absbin.tac4mscon2absbin;
      aud4mscon2absbin = absbin.aud4mscon2absbin;
      nul4mscon2absbin = absbin.nul4mscon2absbin;
      ms4mscon2absbin = absbin.ms4mscon2absbin;
      
      if binonly
        %               save(['absbin_' num2str(sub{ii}) '_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'*abs*bin','*abs','*absIAF');
        %               save(['absbin_' num2str(sub{ii}) '_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'*abs*bin','*abs');
        save(['absbin_' num2str(sub{ii}) '_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'angbin','absbin','freqtr');
        continue
      end
      
    end
    
    for ll=[soalist soalist+20 soalist+40]
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
            keyboard
            continue
          end
          if ll==max(soalist) && ~tacaloneproc
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
          if ll==max(soalist) && ~tacaloneproc
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
            keyboard
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
          if ll==max(soalist)  && ~tacaloneproc
            tlock_nul{10,tt,12}=[];
            tlock_nul{10,tt,13}=[];
          end
          
        else % ss==23
          if ~tr.t10trialkept{ll,tt,ss}
            numt_trials(ll,tt,ss)=0;
            keyboard
            continue
          else
            cfg=[];
            cfg.trials=tr.t10trialkept{ll,tt,ss};
            if length(cfg.trials)<2
              tlock_tac{ll+20,tt,ss}=[]; % create +20 here as 10 will be different for every ll,tt
            else
              tlock_tac{ll+20,tt,ss}=ft_selectdata(cfg,tlock_tac{10,tt,ss}); % create +20 here as 10 will be different for every ll,tt
            end
            if ~tacaloneproc
              if ll==max(soalist) && ss<12 && ~phaset0
                tlock_tac{10,tt,ss}=[];
              end
              
              %               cfg=[];
              %               cfg.trials=tr.a10trialkept{ll,tt,ss};
              %               tlock_aud{ll+20,tt,ss}=ft_selectdata(cfg,tlock_aud{10,tt,ss}); % create +20 here as 10 will be different for every ll,tt
              if ll==max(soalist) && ss<12 && ~phaset0
                tlock_aud{10,tt,ss}=[];
              end
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
            if ~tacaloneproc
              if ll==max(soalist) && ss<12
                tlock_nul{10,tt,ss}=[];
              end
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
          keyboard
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
          keyboard
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
      %               keyboard; % CHECK baseline correct for EVERYWHERE
      if ll<40
        if isfield(tlock_tac{ll,tt,ss},'trial') && size(tlock_tac{ll,tt,ss}.trial,1)
          if phaset0 && phaset0use==0
            
            if ll<10
              tlock_ms_pca=ft_apply_montage(tlock_tac{ll,tt,ss},montage);
              tlock_ms_pca.label{1}='pca_maxcorrerp';
            else
              tlock_tac_pca=ft_apply_montage(tlock_tac{ll,tt,ss},montage);
              tlock_tac_pca.label{1}='pca_maxcorrerp';
            end
            if sign(corr_erppca(mind))==-1
              cfg=[];
              cfg.operation='-1*x1';
              cfg.parameter='trial';
              if ll<10
                tlock_ms_pca =ft_math(cfg,tlock_ms_pca);
              else
                tlock_tac_pca=ft_math(cfg,tlock_tac_pca);
              end
            end
            cfg=[];
            cfg.keeptrials='yes';
            if ll<10
              tlock_ms_pca =ft_timelockanalysis(cfg,tlock_ms_pca);
            else
              tlock_tac_pca=ft_timelockanalysis(cfg,tlock_tac_pca);
            end
            
            ccfg=[];  % for channel subselection to get phase for channel grouping
            ccfg.channel={'Fz' 'F1' 'F2' 'FC1' 'FC2' 'C1' 'Cz' 'C2'};
            ccfg.avgoverchan='yes';
            
            timetest=-0.5:0.05:0;
            freqtest=[1 4; 4 8; 8 12];
            for ff=1:size(freqtest,1)
              for yy=1:length(timetest)
                scfg=[];
                if ll==1
                  scfg.latency=[-0.25 0.25]-0.5+timetest(yy);  % choose asymmetric window as we want phase in range of -0.1 to 0
                elseif ll==3
                  scfg.latency=[-0.25 0.25]-0.07+timetest(yy);
                elseif ll==4
                  scfg.latency=[-0.25 0.25]-0.02+timetest(yy);
                elseif ll==5 || ll==6 || ll==7 || ll==9 || ll>20
                  scfg.latency=[-0.25 0.25]-0+timetest(yy);
                end
                cfg=[];
                cfg.bpfilter='yes';
                cfg.bpfreq=freqtest(ff,:);
                cfg.bpfilttype='fir';
                cfg.bpfiltord=[3*fsample/cfg.bpfreq(1)]; % 3* is default
                cfg.bpfiltdir='onepass-zerophase';
                cfg.demean='yes';
                cfg.hilbert='complex';
                tlocktmp=ft_preprocessing(cfg,ft_selectdata(scfg,tlock_tac{ll,tt,ss}));
                tlocktmpchan=ft_preprocessing(cfg,ft_selectdata(ccfg,ft_selectdata(scfg,tlock_tac{ll,tt,ss}))); % avg over chan then take filter/hilbert
                indsamp=round(length(tlocktmp.time)/2);
                if ll<10
                  tlocktmppca=ft_preprocessing(cfg,ft_selectdata(scfg,tlock_ms_pca)); % avg over chan then take filter/hilbert
                  hilb_time0_firststim_MS{ll,tt,ss}(:,:,ff,yy)       =tlocktmp.trial(:,:,indsamp);
                  hilb_time0_firststim_chanFC_MS{ll,tt,ss}(:,:,ff,yy)=tlocktmpchan.trial(:,:,indsamp);
                  hilb_time0_firststim_chanPCA_MS{ll,tt,ss}(:,:,ff,yy)=tlocktmppca.trial(:,:,indsamp);
                else
                  tlocktmppca=ft_preprocessing(cfg,ft_selectdata(scfg,tlock_tac_pca)); % avg over chan then take filter/hilbert
                  hilb_time0_firststim_Tac{ll-20,tt,ss}(:,:,ff,yy)       =tlocktmp.trial(:,:,indsamp);
                  hilb_time0_firststim_chanFC_Tac{ll-20,tt,ss}(:,:,ff,yy)=tlocktmpchan.trial(:,:,indsamp);
                  hilb_time0_firststim_chanPCA_Tac{ll-20,tt,ss}(:,:,ff,yy)=tlocktmppca.trial(:,:,indsamp);
                end
              end %yy
              %                     if ll==1
              %                       scfg.latency=[-0.75 0.25]-0.5;  % choose asymmetric window as we want phase in range of -0.1 to 0
              %                     elseif ll==3
              %                       scfg.latency=[-0.75 0.25]-0.07;
              %                     elseif ll==4
              %                       scfg.latency=[-0.75 0.25]-0.02;
              %                     elseif ll==5 || ll==6 || ll==7 || ll==9 || ll>20
              %                       scfg.latency=[-0.75 0.25]-0;
              %                     end
              %                     cfg=[];
              %                     cfg.bpfilter='yes';
              %                     cfg.bpfreq=freqtest(ff,:);
              %                     cfg.bpfilttype='fir';
              %                     cfg.bpfiltord=[3*fsample/cfg.bpfreq(1)]; % 3* is default
              %                     cfg.bpfiltdir='onepass-zerophase';
              %                     cfg.demean='yes';
              %                     cfg.hilbert='complex';
              %                     tlocktmp=ft_preprocessing(cfg,ft_selectdata(scfg,tlock_tac{ll,tt,ss}));
              %                     tlocktmpchan=ft_preprocessing(cfg,ft_selectdata(ccfg,ft_selectdata(scfg,tlock_tac{ll,tt,ss}))); % avg over chan then take filter/hilbert
              %                     indsamp=round(length(tlocktmp.time)/2)+[-250:50:250];
              %                     if ll<10
              %                       tlocktmppca=ft_preprocessing(cfg,ft_selectdata(scfg,tlock_ms_pca)); % avg over chan then take filter/hilbert
              %                       hilb_time0_firststim_MS{ll,tt,ss}(:,:,ff,:)       =tlocktmp.trial(:,:,indsamp);
              %                       hilb_time0_firststim_chanFC_MS{ll,tt,ss}(:,:,ff,:)=tlocktmpchan.trial(:,:,indsamp);
              %                       hilb_time0_firststim_chanPCA_MS{ll,tt,ss}(:,:,ff,:)=tlocktmppca.trial(:,:,indsamp);
              %                     else
              %                       tlocktmppca=ft_preprocessing(cfg,ft_selectdata(scfg,tlock_tac_pca)); % avg over chan then take filter/hilbert
              %                       hilb_time0_firststim_Tac{ll-20,tt,ss}(:,:,ff,:)       =tlocktmp.trial(:,:,indsamp);
              %                       hilb_time0_firststim_chanFC_Tac{ll-20,tt,ss}(:,:,ff,:)=tlocktmpchan.trial(:,:,indsamp);
              %                       hilb_time0_firststim_chanPCA_Tac{ll-20,tt,ss}(:,:,ff,:)=tlocktmppca.trial(:,:,indsamp);
              %                     end
              
              if ll<10  % also get at time0 of second stim within multisensory condition
                for yy=1:length(timetest)
                  scfg=[];
                  if ll==1 || ll==3 || ll==4 || ll==5
                    scfg.latency=[-0.25 0.25]+0+timetest(yy);
                  elseif ll==6
                    scfg.latency=[-0.25 0.25]+0.02+timetest(yy);
                  elseif ll==7
                    scfg.latency=[-0.25 0.25]+0.07+timetest(yy);
                  elseif ll==9
                    scfg.latency=[-0.25 0.25]+0.5+timetest(yy);
                  end
                  cfg=[];
                  cfg.bpfilter='yes';
                  cfg.bpfreq=freqtest(ff,:);
                  cfg.bpfilttype='fir';
                  cfg.bpfiltord=[3*fsample/cfg.bpfreq(1)]; % 3* is default
                  cfg.bpfiltdir='onepass-zerophase';
                  cfg.demean='yes';
                  cfg.hilbert='complex';
                  tlocktmp=ft_preprocessing(cfg,ft_selectdata(scfg,tlock_tac{ll,tt,ss}));
                  tlocktmpchan=ft_preprocessing(cfg,ft_selectdata(ccfg,ft_selectdata(scfg,tlock_tac{ll,tt,ss})));
                  tlocktmppca=ft_preprocessing(cfg,ft_selectdata(scfg,tlock_ms_pca)); % avg over chan then take filter/hilbert
                  hilb_time0_secondstim_MS{ll,tt,ss}(:,:,ff,yy)=tlocktmp.trial(:,:,indsamp);
                  hilb_time0_secondstim_chanFC_MS{ll,tt,ss}(:,:,ff,yy)=tlocktmpchan.trial(:,:,indsamp);
                  hilb_time0_secondstim_chanPCA_MS{ll,tt,ss}(:,:,ff,yy)=tlocktmppca.trial(:,:,indsamp);
                end %yy
                %                       scfg=[];
                %                       if ll==1 || ll==3 || ll==4 || ll==5
                %                         scfg.latency=[-0.75 0.25]+0;
                %                       elseif ll==6
                %                         scfg.latency=[-0.75 0.25]+0.02;
                %                       elseif ll==7
                %                         scfg.latency=[-0.75 0.25]+0.07;
                %                       elseif ll==9
                %                         scfg.latency=[-0.75 0.25]+0.5;
                %                       end
                %                       cfg=[];
                %                       cfg.bpfilter='yes';
                %                       cfg.bpfreq=freqtest(ff,:);
                %                       cfg.bpfilttype='fir';
                %                       cfg.bpfiltord=[3*fsample/cfg.bpfreq(1)]; % 3* is default
                %                       cfg.bpfiltdir='onepass-zerophase';
                %                       cfg.demean='yes';
                %                       cfg.hilbert='complex';
                %                       tlocktmp=ft_preprocessing(cfg,ft_selectdata(scfg,tlock_tac{ll,tt,ss}));
                %                       tlocktmpchan=ft_preprocessing(cfg,ft_selectdata(ccfg,ft_selectdata(scfg,tlock_tac{ll,tt,ss})));
                %                       tlocktmppca=ft_preprocessing(cfg,ft_selectdata(scfg,tlock_ms_pca)); % avg over chan then take filter/hilbert
                %                       hilb_time0_secondstim_MS{ll,tt,ss}(:,:,ff,:)=tlocktmp.trial(:,:,indsamp);
                %                       hilb_time0_secondstim_chanFC_MS{ll,tt,ss}(:,:,ff,:)=tlocktmpchan.trial(:,:,indsamp);
                %                       hilb_time0_secondstim_chanPCA_MS{ll,tt,ss}(:,:,ff,:)=tlocktmppca.trial(:,:,indsamp);
              end
            end
          else
            cfg=[];
            cfg.lpfilter='yes';
            cfg.lpfreq=40;
            cfg.demean='yes';
            %cfg.baselinewindow=[-1.7 -0.6];
            if ll==1 || ll==21
              cfg.baselinewindow=[-.4 0]-.6;
            elseif ll==3 || ll==23
              cfg.baselinewindow=[-.4 0]-.17;
            elseif ll==4 || ll==24
              cfg.baselinewindow=[-.4 0]-.12;
            elseif ll==5 || ll==6 || ll==7 || ll==9 || ll>24
              cfg.baselinewindow=[-.4 0]-.1;
            end
            
            disp('ft_preprocessing tac')
            tlock_tac{ll,tt,ss}=ft_preprocessing(cfg,tlock_tac{ll,tt,ss});
            featurestats_tac(:,ll,tt,ss,ii)=mean(tlock_tac{ll,tt,ss}.trialinfo(:,12:19));
            cfg=[];
            cfg.covariance='yes';
            cfg.covariancewindow=[-0.6 0.7]; % a full window valid for all ll
            tlock_tactlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
            if pcmflag
              if rem(ll,10)==1 || rem(ll,10)==3 || rem(ll,10)==4 || rem(ll,10)==5
                starttime=0;     endtime=.5;
              elseif rem(ll,10)==6
                starttime=0.02;  endtime=.52;
              elseif rem(ll,10)==7
                starttime=0.07;  endtime=.57;
              elseif rem(ll,10)==9
                starttime=0.5;   endtime=1;
              end
              cfg.keeptrials='yes';
              cfg.covariance='no';
              tmp_tactlock=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
              cfg=[];
              cfg.latency=[starttime endtime];
              pcm_tactlock{ll,tt,ss}=ft_selectdata(cfg,tmp_tactlock);
            end
            if phaset0 && phaset0use
              %                       if ll<10
              %                         for at=1:numt_trials(ll,tt,ss)
              %                           tind=at;aind=at;
              %                           cfg=[];
              %                           cfg.trials=tind;
              %                           tmpt=ft_selectdata(cfg,tlock_tac{ll+20,tt,ss});
              %                           tmpms=ft_selectdata(cfg,tlock_tac{ll,tt,ss});
              %                           cfg.trials=aind;
              %                           tmpa=ft_selectdata(cfg,tlock_aud{ll+40,tt,ss});
              %                           tmpn=ft_selectdata(cfg,tlock_nul{ll+50,tt,ss});
              %                           cfg=[];
              %                           cfg.operation='add';
              %                           cfg.parameter='trial';
              %                           tmpsumU=ft_math(cfg,tmpt,tmpa);
              %                           tmpsumM=ft_math(cfg,tmpms,tmpn);
              %                           cfg=[];
              %                           cfg.operation='subtract';
              %                           cfg.parameter='trial';
              %                           tmpconMT=ft_math(cfg,tmpms,tmpt); % these are different contrasts, since add/subtact is done before FFT (nonlinear) step
              %                           tmpconMA=ft_math(cfg,tmpms,tmpa);
              %                           tmpconAN=ft_math(cfg,tmpa,tmpn);
              %                           tmpconTN=ft_math(cfg,tmpt,tmpn);
              %                           if at==1
              %                             tlock_fakeU=tmpsumU;
              %                             tlock_fakeM=tmpsumM;
              %                             tlock_fakeMT=tmpconMT;
              %                             tlock_fakeMA=tmpconMA;
              %                             tlock_fakeTN=tmpconTN;
              %                             tlock_fakeAN=tmpconAN;
              %                           end
              %                           tlock_fakeU.trial(at,:,:)=tmpsumU.trial(1,:,:);
              %                           tlock_fakeM.trial(at,:,:)=tmpsumM.trial(1,:,:);
              %                           tlock_fakeMT.trial(at,:,:)=tmpconMT.trial(1,:,:);
              %                           tlock_fakeMA.trial(at,:,:)=tmpconMA.trial(1,:,:);
              %                           tlock_fakeTN.trial(at,:,:)=tmpconTN.trial(1,:,:);
              %                           tlock_fakeAN.trial(at,:,:)=tmpconAN.trial(1,:,:);
              %                         end % at
              %                         clear tmp*
              %                       end
              
              for bb=1:4
                if ll<10
                  cfg=[];
                  %                         cfg.trials=find(ms1bin{ll,tt,ss}(:,bb));
                  %                         tlock_ms1tlock_phasebin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
                  %                         cfg.trials=find(ms2bin{ll,tt,ss}(:,bb));
                  % %                         if ~isempty(cfg.trials)
                  %                           tlock_ms2tlock_phasebin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
                  % %                         else
                  % %                           tlock_ms2tlock_phasebin{ll,tt,ss,bb}=[];
                  % %                         end
                  
                  cfg.trials=find(ms4mscon1bin{ll,tt,ss}(:,bb));
                  tlock_ms4mscon1_phasebin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
                  cfg.trials=find(ms4mscon2bin{ll,tt,ss}(:,bb));
                  tlock_ms4mscon2_phasebin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
                  
                  %                         cfg.trials=find(ms4mscon1IAFbin{ll,tt,ss}(:,bb));
                  %                         tlock_ms4mscon1IAF_phasebin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
                  %                         cfg.trials=find(ms4mscon2IAFbin{ll,tt,ss}(:,bb));
                  % %                         if ~isempty(cfg.trials)
                  %                           tlock_ms4mscon2IAF_phasebin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
                  % %                         else
                  % %                           tlock_ms4mscon2IAF_phasebin{ll,tt,ss,bb}=[];
                  % %                         end
                  
                  
                  cfg=[];
                  %                         cfg.trials=find(ms1absbin{ll,tt,ss}(:,bb));
                  %                         tlock_ms1tlock_absbin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
                  %                         cfg.trials=find(ms2absbin{ll,tt,ss}(:,bb));
                  %                         tlock_ms2tlock_absbin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
                  
                  cfg.trials=find(ms4mscon1absbin{ll,tt,ss}(:,bb));
                  tlock_ms4mscon1_absbin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
                  cfg.trials=find(ms4mscon2absbin{ll,tt,ss}(:,bb));
                  %                         if ~isempty(cfg.trials)
                  tlock_ms4mscon2_absbin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
                  %                         else
                  %                           tlock_ms4mscon2_absbin{ll,tt,ss,bb}=[];
                  %                         end
                  
                  %                         cfg.trials=find(ms4mscon1absIAFbin{ll,tt,ss}(:,bb));
                  %                         tlock_ms4mscon1IAF_absbin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
                  %                         cfg.trials=find(ms4mscon2absIAFbin{ll,tt,ss}(:,bb));
                  % %                         if ~isempty(cfg.trials)
                  %                           tlock_ms4mscon2IAF_absbin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
                  % %                         else
                  % %                           tlock_ms4mscon2IAF_absbin{ll,tt,ss,bb}=[];
                  % %                         end
                  
                else
                  cfg=[];
                  %                         cfg.trials=find(tacbin{ll-20,tt,ss}(:,bb));
                  %                         tlock_tactlock_phasebin{ll-20,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
                  cfg.trials=find(tac4mscon1bin{ll-20,tt,ss}(:,bb));
                  tlock_tac4mscon1_phasebin{ll-20,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
                  cfg.trials=find(tac4mscon2bin{ll-20,tt,ss}(:,bb));
                  tlock_tac4mscon2_phasebin{ll-20,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
                  
                  %                         cfg.trials=find(tac4mscon1IAFbin{ll-20,tt,ss}(:,bb));
                  %                         tlock_tac4mscon1IAF_phasebin{ll-20,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
                  %                         cfg.trials=find(tac4mscon2IAFbin{ll-20,tt,ss}(:,bb));
                  %                         tlock_tac4mscon2IAF_phasebin{ll-20,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
                  
                  cfg=[];
                  %                         cfg.trials=find(tacabsbin{ll-20,tt,ss}(:,bb));
                  % %                         if ~isempty(cfg.trials)
                  %                           tlock_tactlock_absbin{ll-20,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
                  % %                         else
                  % %                           tlock_tactlock_absbin{ll-20,tt,ss,bb}=[];
                  % %                         end
                  cfg.trials=find(tac4mscon1absbin{ll-20,tt,ss}(:,bb));
                  tlock_tac4mscon1_absbin{ll-20,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
                  cfg.trials=find(tac4mscon2absbin{ll-20,tt,ss}(:,bb));
                  tlock_tac4mscon2_absbin{ll-20,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
                  
                  %                         cfg.trials=find(tac4mscon1absIAFbin{ll-20,tt,ss}(:,bb));
                  %                         tlock_tac4mscon1IAF_absbin{ll-20,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
                  %                         cfg.trials=find(tac4mscon2absIAFbin{ll-20,tt,ss}(:,bb));
                  %                         tlock_tac4mscon2IAF_absbin{ll-20,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
                end
                %                       % we need to do all analysis here before tac{ll,tt,ss} gets deleted
                %                       if ll<10
                %                         cfg=[];
                %                         cfg.trials=find(tacPaud1bin{ll,tt,ss}(:,bb));
                %                         tlock_tacPaud1_phasebin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_fakeU);
                %                         cfg.trials=find(tacPaud2bin{ll,tt,ss}(:,bb));
                %                         tlock_tacPaud2_phasebin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_fakeU);
                %                         cfg.trials=find(tacMSpN1bin{ll,tt,ss}(:,bb));
                %                         tlock_tacMSpN1_phasebin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_fakeM);
                %                         cfg.trials=find(tacMSpN2bin{ll,tt,ss}(:,bb));
                %                         tlock_tacMSpN2_phasebin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_fakeM);
                %                         cfg.trials=find(tacAmN1bin{ll,tt,ss}(:,bb));
                %                         tlock_tacAmN1_phasebin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_fakeAN);
                %                         cfg.trials=find(tacAmN2bin{ll,tt,ss}(:,bb));
                %                         tlock_tacAmN2_phasebin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_fakeAN);
                %                         cfg.trials=find(tacTmN1bin{ll,tt,ss}(:,bb));
                %                         tlock_tactmN1_phasebin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_fakeTN);
                %                         cfg.trials=find(tactmN2bin{ll,tt,ss}(:,bb));
                %                         tlock_tacTmN2_phasebin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_fakeTN);
                %                         cfg.trials=find(tacMSmA1bin{ll,tt,ss}(:,bb));
                %                         tlock_tacMSmA1_phasebin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_fakeMA);
                %                         cfg.trials=find(tacMSmA2bin{ll,tt,ss}(:,bb));
                %                         tlock_tacMSmA2_phasebin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_fakeMA);
                %                         cfg.trials=find(tacMSmT1bin{ll,tt,ss}(:,bb));
                %                         tlock_tacMSmT1_phasebin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_fakeMT);
                %                         cfg.trials=find(tacMSmT2bin{ll,tt,ss}(:,bb));
                %                         tlock_tacMSmT2_phasebin{ll,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_fakeMT);
                %                       end
                
              end % bb
            end
          end
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
          if phaset0 && phaset0use==0
            tlock_aud_pca=ft_apply_montage(tlock_aud{ll,tt,ss},montage);
            tlock_aud_pca.label{1}='pca_maxcorrerp';
            if sign(corr_erppca(mind))==-1
              cfg=[];
              cfg.operation='-1*x1';
              cfg.parameter='trial';
              tlock_aud_pca=ft_math(cfg,tlock_aud_pca);
            end
            cfg=[];
            cfg.keeptrials='yes';
            tlock_aud_pca=ft_timelockanalysis(cfg,tlock_aud_pca);
            
            freqtest=[1 4; 4 8; 8 12];
            for ff=1:size(freqtest,1)
              for yy=1:length(timetest)
                scfg=[];
                if ll==41
                  scfg.latency=[-0.25 0.25]-0.5+timetest(yy);
                elseif ll==43
                  scfg.latency=[-0.25 0.25]-0.07+timetest(yy);
                elseif ll==44
                  scfg.latency=[-0.25 0.25]-0.02+timetest(yy);
                elseif ll==45
                  scfg.latency=[-0.25 0.25]-0+timetest(yy);
                elseif ll==46
                  scfg.latency=[-0.25 0.25]+0.02+timetest(yy);
                elseif ll==47
                  scfg.latency=[-0.25 0.25]+0.07+timetest(yy);
                elseif ll==49
                  scfg.latency=[-0.25 0.25]+0.5+timetest(yy);
                end
                cfg=[];
                cfg.bpfilter='yes';
                cfg.bpfreq=freqtest(ff,:);
                cfg.bpfilttype='fir';
                cfg.bpfiltord=[3*fsample/cfg.bpfreq(1)]; % 3* is default
                cfg.bpfiltdir='onepass-zerophase';
                cfg.demean='yes';
                cfg.hilbert='complex';
                tlocktmp=ft_preprocessing(cfg,ft_selectdata(scfg,tlock_aud{ll,tt,ss}));
                tlocktmpchan=ft_preprocessing(cfg,ft_selectdata(ccfg,ft_selectdata(scfg,tlock_aud{ll,tt,ss})));
                tlocktmppca=ft_preprocessing(cfg,ft_selectdata(scfg,tlock_aud_pca));
                hilb_time0_firststim_Aud{ll-40,tt,ss}(:,:,ff,yy)=tlocktmp.trial(:,:,indsamp);
                hilb_time0_firststim_chanFC_Aud{ll-40,tt,ss}(:,:,ff,yy)=tlocktmpchan.trial(:,:,indsamp);
                hilb_time0_firststim_chanPCA_Aud{ll-40,tt,ss}(:,:,ff,yy)=tlocktmppca.trial(:,:,indsamp);
              end
              %                     scfg=[];
              %                     if ll==41
              %                       scfg.latency=[-0.75 0.25]-0.5;
              %                     elseif ll==43
              %                       scfg.latency=[-0.75 0.25]-0.07;
              %                     elseif ll==44
              %                       scfg.latency=[-0.75 0.25]-0.02;
              %                     elseif ll==45
              %                       scfg.latency=[-0.75 0.25]-0;
              %                     elseif ll==46
              %                       scfg.latency=[-0.75 0.25]+0.02;
              %                     elseif ll==47
              %                       scfg.latency=[-0.75 0.25]+0.07;
              %                     elseif ll==49
              %                       scfg.latency=[-0.75 0.25]+0.5;
              %                     end
              %                     cfg=[];
              %                     cfg.bpfilter='yes';
              %                     cfg.bpfreq=freqtest(ff,:);
              %                     cfg.bpfilttype='fir';
              %                     cfg.bpfiltord=[3*fsample/cfg.bpfreq(1)]; % 3* is default
              %                     cfg.bpfiltdir='onepass-zerophase';
              %                     cfg.demean='yes';
              %                     cfg.hilbert='complex';
              %                     tlocktmp=ft_preprocessing(cfg,ft_selectdata(scfg,tlock_aud{ll,tt,ss}));
              %                     tlocktmpchan=ft_preprocessing(cfg,ft_selectdata(ccfg,ft_selectdata(scfg,tlock_aud{ll,tt,ss})));
              %                     tlocktmppca=ft_preprocessing(cfg,ft_selectdata(scfg,tlock_aud_pca));
              %                     indsamp=round(length(tlocktmp.time)/2)+[-250:50:250];
              %                     hilb_time0_firststim_Aud{ll-40,tt,ss}(:,:,ff,:)=tlocktmp.trial(:,:,indsamp);
              %                     hilb_time0_firststim_chanFC_Aud{ll-40,tt,ss}(:,:,ff,:)=tlocktmpchan.trial(:,:,indsamp);
              %                     hilb_time0_firststim_chanPCA_Aud{ll-40,tt,ss}(:,:,ff,:)=tlocktmppca.trial(:,:,indsamp);
            end
          else
            cfg=[];
            cfg.lpfilter='yes';
            cfg.lpfreq=40;
            %         cfg.bpfilter='yes';
            %         cfg.bpfreq=[1 40];
            cfg.demean='yes';
            %               cfg.baselinewindow=[-1.7 -0.6];
            if ll==41
              cfg.baselinewindow=[-.4 0]-.6;
            elseif ll==43
              cfg.baselinewindow=[-.4 0]-.17;
            elseif ll==44
              cfg.baselinewindow=[-.4 0]-.12;
            elseif ll==45 || ll==46 || ll==47 || ll==49
              cfg.baselinewindow=[-.4 0]-.1;
            end
            disp('ft_preprocessing aud')
            tlock_aud{ll,tt,ss}=ft_preprocessing(cfg,tlock_aud{ll,tt,ss});
            featurestats_aud(:,ll,tt,ss,ii)=mean(tlock_aud{ll,tt,ss}.trialinfo(:,12:19));
            cfg=[];
            cfg.covariance='yes';
            cfg.covariancewindow=[-0.6 0.7]; % a full window valid for all ll
            tlock_audtlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_aud{ll,tt,ss});
            if pcmflag
              if rem(ll,10)==1 || rem(ll,10)==3 || rem(ll,10)==4 || rem(ll,10)==5
                starttime=0;     endtime=.5;
              elseif rem(ll,10)==6
                starttime=0.02;  endtime=.52;
              elseif rem(ll,10)==7
                starttime=0.07;  endtime=.57;
              elseif rem(ll,10)==9
                starttime=0.5;   endtime=1;
              end
              cfg.keeptrials='yes';
              cfg.covariance='no';
              tmp_audtlock=ft_timelockanalysis(cfg,tlock_aud{ll,tt,ss});
              cfg=[];
              cfg.latency=[starttime endtime];
              pcm_audtlock{ll,tt,ss}=ft_selectdata(cfg,tmp_audtlock);
            end
            if phaset0 && phaset0use
              for bb=1:4
                cfg=[];
                %                       cfg.trials=find(audbin{ll-40,tt,ss}(:,bb));
                % %                       if ~isempty(cfg.trials)
                %                         tlock_audtlock_phasebin{ll-40,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_aud{ll,tt,ss});
                % %                       else
                % %                         tlock_audtlock_phasebin{ll-40,tt,ss,bb}=[];
                % %                       end
                cfg.trials=find(aud4mscon1bin{ll-40,tt,ss}(:,bb));
                tlock_aud4mscon1_phasebin{ll-40,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_aud{ll,tt,ss});
                cfg.trials=find(aud4mscon2bin{ll-40,tt,ss}(:,bb));
                %                       if ~isempty(cfg.trials)
                tlock_aud4mscon2_phasebin{ll-40,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_aud{ll,tt,ss});
                %                       else
                %                         tlock_aud4mscon2_phasebin{ll-40,tt,ss,bb}=[];
                %                       end
                %                       cfg.trials=find(aud4mscon1IAFbin{ll-40,tt,ss}(:,bb));
                %                       tlock_aud4mscon1IAF_phasebin{ll-40,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_aud{ll,tt,ss});
                %                       cfg.trials=find(aud4mscon2IAFbin{ll-40,tt,ss}(:,bb));
                % %                       if ~isempty(cfg.trials)
                %                         tlock_aud4mscon2IAF_phasebin{ll-40,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_aud{ll,tt,ss});
                % %                       else
                % %                         tlock_aud4mscon2IAF_phasebin{ll-40,tt,ss,bb}=[];
                % %                       end
                
                cfg=[];
                %                       cfg.trials=find(audabsbin{ll-40,tt,ss}(:,bb));
                %                       tlock_audtlock_absbin{ll-40,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_aud{ll,tt,ss});
                cfg.trials=find(aud4mscon1absbin{ll-40,tt,ss}(:,bb));
                tlock_aud4mscon1_absbin{ll-40,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_aud{ll,tt,ss});
                cfg.trials=find(aud4mscon2absbin{ll-40,tt,ss}(:,bb));
                tlock_aud4mscon2_absbin{ll-40,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_aud{ll,tt,ss});
                %                       cfg.trials=find(aud4mscon1absIAFbin{ll-40,tt,ss}(:,bb));
                %                       tlock_aud4mscon1IAF_absbin{ll-40,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_aud{ll,tt,ss});
                %                       cfg.trials=find(aud4mscon2absIAFbin{ll-40,tt,ss}(:,bb));
                %                       tlock_aud4mscon2IAF_absbin{ll-40,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_aud{ll,tt,ss});
              end
            end
          end
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
      if length(tr.nllttrialkept{ll-50,tt,ss})>=2 && isfield(tlock_nul{ll,tt,ss},'trial') && size(tlock_nul{ll,tt,ss}.trial,1)
        if phaset0 && phaset0use==0
          tlock_nul_pca=ft_apply_montage(tlock_nul{ll,tt,ss},montage);
          tlock_nul_pca.label{1}='pca_maxcorrerp';
          if sign(corr_erppca(mind))==-1
            cfg=[];
            cfg.operation='-1*x1';
            cfg.parameter='trial';
            tlock_nul_pca=ft_math(cfg,tlock_nul_pca);
          end
          cfg=[];
          cfg.keeptrials='yes';
          tlock_nul_pca=ft_timelockanalysis(cfg,tlock_nul_pca);
          
          freqtest=[1 4; 4 8; 8 12];
          for ff=1:size(freqtest,1)
            for yy=1:length(timetest)
              scfg=[];
              scfg.latency=[-0.25 0.25]-0+timetest(yy);
              cfg=[];
              cfg.bpfilter='yes';
              cfg.bpfreq=freqtest(ff,:);
              cfg.bpfilttype='fir';
              cfg.bpfiltord=[3*fsample/cfg.bpfreq(1)]; % 3* is default
              cfg.bpfiltdir='onepass-zerophase';
              cfg.demean='yes';
              cfg.hilbert='complex';
              tlocktmp=ft_preprocessing(cfg,ft_selectdata(scfg,tlock_nul{ll,tt,ss}));
              tlocktmpchan=ft_preprocessing(cfg,ft_selectdata(ccfg,ft_selectdata(scfg,tlock_nul{ll,tt,ss})));
              tlocktmppca=ft_preprocessing(cfg,ft_selectdata(scfg,tlock_nul_pca));
              hilb_time0_firststim_Nul{ll-50,tt,ss}(:,:,ff,yy)=tlocktmp.trial(:,:,indsamp);
              hilb_time0_firststim_chanFC_Nul{ll-50,tt,ss}(:,:,ff,yy)=tlocktmpchan.trial(:,:,indsamp);
              hilb_time0_firststim_chanPCA_Nul{ll-50,tt,ss}(:,:,ff,yy)=tlocktmppca.trial(:,:,indsamp);
            end
            %                   scfg=[];
            %                   scfg.latency=[-0.75 0.25]-0;
            %                   cfg=[];
            %                   cfg.bpfilter='yes';
            %                   cfg.bpfreq=freqtest(ff,:);
            %                   cfg.bpfilttype='fir';
            %                   cfg.bpfiltord=[3*fsample/cfg.bpfreq(1)]; % 3* is default
            %                   cfg.bpfiltdir='onepass-zerophase';
            %                   cfg.demean='yes';
            %                   cfg.hilbert='complex';
            %                   tlocktmp=ft_preprocessing(cfg,ft_selectdata(scfg,tlock_nul{ll,tt,ss}));
            %                   tlocktmpchan=ft_preprocessing(cfg,ft_selectdata(ccfg,ft_selectdata(scfg,tlock_nul{ll,tt,ss})));
            %                   tlocktmppca=ft_preprocessing(cfg,ft_selectdata(scfg,tlock_nul_pca));
            %                   indsamp=round(length(tlocktmp.time)/2)+[-250:50:250];
            %                   hilb_time0_firststim_Nul{ll-50,tt,ss}(:,:,ff,:)=tlocktmp.trial(:,:,indsamp);
            %                   hilb_time0_firststim_chanFC_Nul{ll-50,tt,ss}(:,:,ff,:)=tlocktmpchan.trial(:,:,indsamp);
            %                   hilb_time0_firststim_chanPCA_Nul{ll-50,tt,ss}(:,:,ff,:)=tlocktmppca.trial(:,:,indsamp);
          end
        else
          cfg=[];
          cfg.lpfilter='yes';
          cfg.lpfreq=40;
          cfg.demean='yes';
          %             cfg.baselinewindow=[-1.7 -0.6];
          if ll==51
            cfg.baselinewindow=[-.4 0]-.6;
          elseif ll==53
            cfg.baselinewindow=[-.4 0]-.17;
          elseif ll==54
            cfg.baselinewindow=[-.4 0]-.12;
          elseif ll==55 || ll==56 || ll==57 || ll==59
            cfg.baselinewindow=[-.4 0]-.1;
          end
          %           cfg.bpfilter='yes';
          %           cfg.bpfreq=[1 40];
          disp('ft_preprocessing nul')
          tlock_nul{ll,tt,ss}=ft_preprocessing(cfg,tlock_nul{ll,tt,ss});
          featurestats_nul(:,ll,tt,ss,ii)=mean(tlock_nul{ll,tt,ss}.trialinfo(:,12:19));
          cfg=[];
          cfg.covariance='yes';
          cfg.covariancewindow=[-0.6 0.7]; % a full window valid for all ll
          tlock_nultlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_nul{ll,tt,ss});
          if pcmflag
            if rem(ll,10)==1 || rem(ll,10)==3 || rem(ll,10)==4 || rem(ll,10)==5
              starttime=0;     endtime=.5;
            elseif rem(ll,10)==6
              starttime=0.02;  endtime=.52;
            elseif rem(ll,10)==7
              starttime=0.07;  endtime=.57;
            elseif rem(ll,10)==9
              starttime=0.5;   endtime=1;
            end
            cfg.keeptrials='yes';
            cfg.covariance='no';
            tmp_nultlock=ft_timelockanalysis(cfg,tlock_nul{ll,tt,ss});
            cfg=[];
            cfg.latency=[starttime endtime];
            pcm_nultlock{ll,tt,ss}=ft_selectdata(cfg,tmp_nultlock);
          end
          if phaset0 && phaset0use
            for bb=1:4
              cfg=[];
              %                     cfg.trials=find(nulbin{ll-50,tt,ss}(:,bb));
              % %                     if ~isempty(cfg.trials)
              %                       tlock_nultlock_phasebin{ll-50,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_nul{ll,tt,ss});
              % %                     else
              % %                       tlock_nultlock_phasebin{ll-50,tt,ss,bb}=[];
              % %                     end
              cfg.trials=find(nul4mscon1bin{ll-50,tt,ss}(:,bb));
              %                     if ~isempty(cfg.trials)
              tlock_nul4mscon1_phasebin{ll-50,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_nul{ll,tt,ss});
              %                     else
              %                       tlock_nul4mscon1_phasebin{ll-50,tt,ss,bb}=[];
              %                     end
              cfg.trials=find(nul4mscon2bin{ll-50,tt,ss}(:,bb));
              %                     if ~isempty(cfg.trials)
              tlock_nul4mscon2_phasebin{ll-50,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_nul{ll,tt,ss});
              %                     else
              %                       tlock_nul4mscon2_phasebin{ll-50,tt,ss,bb}=[];
              %                     end
              %                     cfg.trials=find(nul4mscon1IAFbin{ll-50,tt,ss}(:,bb));
              %                     tlock_nul4mscon1IAF_phasebin{ll-50,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_nul{ll,tt,ss});
              %                     cfg.trials=find(nul4mscon2IAFbin{ll-50,tt,ss}(:,bb));
              %                     tlock_nul4mscon2IAF_phasebin{ll-50,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_nul{ll,tt,ss});
              
              cfg=[];
              %                     cfg.trials=find(nulabsbin{ll-50,tt,ss}(:,bb));
              %                     tlock_nultlock_absbin{ll-50,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_nul{ll,tt,ss});
              cfg.trials=find(nul4mscon1absbin{ll-50,tt,ss}(:,bb));
              tlock_nul4mscon1_absbin{ll-50,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_nul{ll,tt,ss});
              cfg.trials=find(nul4mscon2absbin{ll-50,tt,ss}(:,bb));
              tlock_nul4mscon2_absbin{ll-50,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_nul{ll,tt,ss});
              %                     cfg.trials=find(nul4mscon1absIAFbin{ll-50,tt,ss}(:,bb));
              %                     tlock_nul4mscon1IAF_absbin{ll-50,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_nul{ll,tt,ss});
              %                     cfg.trials=find(nul4mscon2absIAFbin{ll-50,tt,ss}(:,bb));
              % %                     if ~isempty(cfg.trials)
              %                       tlock_nul4mscon2IAF_absbin{ll-50,tt,ss,bb}=ft_timelockanalysis(cfg,tlock_nul{ll,tt,ss});
              % %                     else
              % %                       tlock_nul4mscon2IAF_absbin{ll-50,tt,ss,bb}=[];
              % %                     end
            end
          end
        end
      else
        tlock_nultlock{ll,tt,ss}=[];
      end
      tlock_nul{ll,tt,ss}=[];
    end
    clear tmp_*
    
    if phaset0 && phaset0use==0
      save(['tlock_unimultsensHilbert_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'hilb*','tr')
      continue  % continue to next ss
    elseif phaset0 && phaset0use>0
      %             save(['tlock_ERPphasebin_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'tlock_*phasebin','tr','phaset0use');
      %             continue  % continue to next ss
    end
    
    
    if tacaloneproc
      if isfield(tlock_tac{102,tt,ss},'trial') && size(tlock_tac{102,tt,ss}.trial,1)>2 && isfield(tlock_tac{10,tt,ss},'trial') && size(tlock_tac{10,tt,ss}.trial,1)>2 && isfield(tlock_aud{102,tt,ss},'trial') && size(tlock_aud{102,tt,ss}.trial,1)>2 && isfield(tlock_aud{10,tt,ss},'trial') && size(tlock_aud{10,tt,ss}.trial,1)>2 && isfield(tlock_nul{10,tt,ss},'trial') && size(tlock_nul{10,tt,ss}.trial,1)>2
        
        
        if ~isempty(tlock_tac{102,tt,ss}), tlock_tac{102,tt,ss}=trim_nans(tlock_tac{102,tt,ss}); end
        if ~isempty(tlock_tac{10,tt,ss}),  tlock_tac{10,tt,ss}=trim_nans(tlock_tac{10,tt,ss}); end
        if ~isempty(tlock_aud{102,tt,ss}), tlock_aud{102,tt,ss}=trim_nans(tlock_aud{102,tt,ss}); end
        if ~isempty(tlock_aud{10,tt,ss}),  tlock_aud{10,tt,ss}=trim_nans(tlock_aud{10,tt,ss}); end
        
        cfg=[];
        try
          tlock_tac9T_tlock{tt,ss}=ft_timelockanalysis(cfg,tlock_tac{102,tt,ss});
        catch
          tlock_tac9T_tlock{tt,ss}=[];
        end
        try tlock_aud1A_tlock{tt,ss}=ft_timelockanalysis(cfg,tlock_aud{102,tt,ss}); catch
          tlock_aud1A_tlock{tt,ss}=[];end
        try tlock_tacT_tlock{tt,ss}=ft_timelockanalysis(cfg,tlock_tac{10,tt,ss});   catch
          tlock_tacT_tlock{tt,ss}=[];end
        try tlock_audA_tlock{tt,ss}=ft_timelockanalysis(cfg,tlock_aud{10,tt,ss});   catch
          tlock_audA_tlock{tt,ss}=[];end
        try tlock_nul_tlock{tt,ss}=ft_timelockanalysis(cfg,tlock_nul{10,tt,ss});    catch
          tlock_nul_tlock{tt,ss}=[];end
        
        cfg=[];
        cfg.lpfilter='yes';
        cfg.lpfreq=40;
        cfg.demean='yes';
        cfg.baselinewindow=[-.4 0]-.6;
        tlock_tac9T{tt,ss}=ft_preprocessing(cfg,tlock_tac{102,tt,ss});
        tlock_tacT{tt,ss}=ft_preprocessing(cfg,tlock_tac{10,tt,ss});
        tlock_aud1A{tt,ss}=ft_preprocessing(cfg,tlock_aud{102,tt,ss});
        tlock_audA{tt,ss}=ft_preprocessing(cfg,tlock_aud{10,tt,ss});
        tlock_nul{tt,ss}=ft_preprocessing(cfg,tlock_nul{10,tt,ss});
        
        cfg=[];
        cfg.method = 'mtmconvol';
        cfg.keeptrials = 'yes';
        cfg.taper = 'hanning';
        cfg.foi = 1:13;
        cfg.toi = -0.5./cfg.foi;
        cfg.t_ftimwin = 1./cfg.foi; % use e.g. 4 cycles to estimate the phase
        cfg.output = 'fourier';
        cfg.pad = 4;
        freq_tac9T{tt,ss} = ft_freqanalysis(cfg, tlock_tac{102,tt,ss});
        freq_tacT{tt,ss} = ft_freqanalysis(cfg, tlock_tac{10,tt,ss});
        freq_aud1A{tt,ss} = ft_freqanalysis(cfg, tlock_aud{102,tt,ss});
        freq_audA{tt,ss} = ft_freqanalysis(cfg, tlock_aud{10,tt,ss});
        freq_nul{tt,ss} = ft_freqanalysis(cfg, tlock_nul{10,tt,ss});
        
        
        freqhil=[1 3; 4 7; 9 12];
        cfg=[];
        cfg.hilbert='complex';
        cfg.bpfilter='yes';
        cfg.bpfilttype='firws';
        cfg.bpfiltdir='onepass-zerophase'; % default for firws
        for ff=1:size(freqhil,1)
          cfg.bpfreq=freqhil(ff,:);
          tlock_hilbert_tac9T{tt,ss}{ff}=ft_preprocessing(cfg,tlock_tac{102,tt,ss});
          tlock_hilbert_tacT{tt,ss}{ff} =ft_preprocessing(cfg,tlock_tac{10,tt,ss});
          tlock_hilbert_aud1A{tt,ss}{ff}=ft_preprocessing(cfg,tlock_aud{102,tt,ss});
          tlock_hilbert_audA{tt,ss}{ff} =ft_preprocessing(cfg,tlock_aud{10,tt,ss});
          tlock_hilbert_nul{tt,ss}{ff}  =ft_preprocessing(cfg,tlock_nul{10,tt,ss});
        end
        
        % 'angle' gives angle of cosine.  0 and 2pi are peaks, pi and
        % -pi are troughs; pi/2 and -pi/2 are zero-crossings.
        % subdivide into 4 quadrants based on peak and trough
        
        % use shifting time for phase peak location for each frequency
        erpamp_phasefft_tac9T{tt,ss} = erpamp_phasesorted_fft(freq_tac9T{tt,ss}, tlock_tac9T{tt,ss}, 'Fz', pi);
        erpamp_phasefft_tacT{tt,ss}  = erpamp_phasesorted_fft(freq_tacT{tt,ss},  tlock_tacT{tt,ss}, 'Fz', pi);
        erpamp_phasefft_aud1A{tt,ss} = erpamp_phasesorted_fft(freq_aud1A{tt,ss}, tlock_aud1A{tt,ss}, 'Fz', pi);
        erpamp_phasefft_audA{tt,ss}  = erpamp_phasesorted_fft(freq_audA{tt,ss},  tlock_audA{tt,ss}, 'Fz', pi);
        erpamp_phasefft_nul{tt,ss}   = erpamp_phasesorted_fft(freq_nul{tt,ss},   tlock_nul{tt,ss}, 'Fz', pi);
        
        % picking channel Fz and time point 0 after Hilbert
        erpamp_phasehilb_tac9T{tt,ss} = erpamp_phasesorted_hilbert(tlock_hilbert_tac9T{tt,ss}, tlock_tac9T{tt,ss}, 'Fz', 0 );
        erpamp_phasehilb_tacT{tt,ss} =  erpamp_phasesorted_hilbert(tlock_hilbert_tacT{tt,ss},  tlock_tacT{tt,ss}, 'Fz', 0 );
        erpamp_phasehilb_aud1A{tt,ss} = erpamp_phasesorted_hilbert(tlock_hilbert_aud1A{tt,ss}, tlock_aud1A{tt,ss}, 'Fz', 0 );
        erpamp_phasehilb_audA{tt,ss} =  erpamp_phasesorted_hilbert(tlock_hilbert_audA{tt,ss},  tlock_audA{tt,ss}, 'Fz', 0 );
        erpamp_phasehilb_nul{tt,ss} =   erpamp_phasesorted_hilbert(tlock_hilbert_nul{tt,ss},   tlock_nul{tt,ss}, 'Fz', 0 );
        
        if ss==12 % we want all data available stages W, N1, N2
          erpamp_powrat_tac9T = erpamp_powratsorted( freq_tac9T, tlock_tac9T, 'Fz', [-1 -.5]);
          erpamp_powrat_tacT =  erpamp_powratsorted( freq_tacT,  tlock_tacT, 'Fz', [-1 -.5]);
          erpamp_powrat_aud1A = erpamp_powratsorted( freq_aud1A, tlock_aud1A, 'Fz', [-1 -.5]);
          erpamp_powrat_audA =  erpamp_powratsorted( freq_audA,  tlock_audA, 'Fz', [-1 -.5]);
          erpamp_powrat_nul =   erpamp_powratsorted( freq_nul,   tlock_nul, 'Fz', [-1 -.5]);
        end
        
        % get amplitude sorted by sleep stage
        erpamp_tac9T(ss,1) = median(tlock_tac9T{tt,ss}.trial(:,match_str(tlock_tac9T{tt,ss}.label,'Fz'),dsearchn(tlock_tac9T{tt,ss}.time',0.1)));
        erpamp_tac9T(ss,2) = median(tlock_tac9T{tt,ss}.trial(:,match_str(tlock_tac9T{tt,ss}.label,'Fz'),dsearchn(tlock_tac9T{tt,ss}.time',0.2)));
        erpamp_tacT(ss,1)  = median(tlock_tacT{tt,ss}.trial(:, match_str(tlock_tacT{tt,ss}.label,'Fz'), dsearchn(tlock_tacT{tt,ss}.time',0.1)));
        erpamp_tacT(ss,2)  = median(tlock_tacT{tt,ss}.trial(:, match_str(tlock_tacT{tt,ss}.label,'Fz'), dsearchn(tlock_tacT{tt,ss}.time',0.2)));
        erpamp_aud1A(ss,1) = median(tlock_aud1A{tt,ss}.trial(:,match_str(tlock_aud1A{tt,ss}.label,'Fz'),dsearchn(tlock_aud1A{tt,ss}.time',0.1)));
        erpamp_aud1A(ss,2) = median(tlock_aud1A{tt,ss}.trial(:,match_str(tlock_aud1A{tt,ss}.label,'Fz'),dsearchn(tlock_aud1A{tt,ss}.time',0.2)));
        erpamp_audA(ss,1)  = median(tlock_audA{tt,ss}.trial(:, match_str(tlock_audA{tt,ss}.label,'Fz'), dsearchn(tlock_audA{tt,ss}.time',0.1)));
        erpamp_audA(ss,2)  = median(tlock_audA{tt,ss}.trial(:, match_str(tlock_audA{tt,ss}.label,'Fz'), dsearchn(tlock_audA{tt,ss}.time',0.2)));
        erpamp_nul(ss,1)   = median(tlock_nul{tt,ss}.trial(:,  match_str(tlock_nul{tt,ss}.label,'Fz'),  dsearchn(tlock_nul{tt,ss}.time',0.1)));
        erpamp_nul(ss,2)   = median(tlock_nul{tt,ss}.trial(:,  match_str(tlock_nul{tt,ss}.label,'Fz'),  dsearchn(tlock_nul{tt,ss}.time',0.2)));
        
        
      else
        tlock_tac9T_tlock{tt,ss}=[];
        erpamp_phasefft_tac9T=[];
      end
      save(['tlock_tacaloneERPamp_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'tlock_*_tlock','erpamp*')
      continue  % continue to next ss
    end
    
    for ll=soalist
      % create sum of unisensory conditions
      cfg=[];
      cfg.operation='add';
      cfg.parameter='avg';
      %           cfg.parameter={'avg' 'cov'};
      % the call to ft_timelockanalysis is b/c .avg field not present after ft_redefinetrial
      if numt_trials(ll,tt,ss)>=2
        tlock_tacPaud{ll,tt,ss}=ft_math(cfg,tlock_tactlock{ll+20,tt,ss},tlock_audtlock{ll+40,tt,ss}); % keeping tac centred but jitter aud (to match tlock_tac{X,tt})
        cfg.parameter='cov';
        tmp=ft_math(cfg,tlock_tactlock{ll+20,tt,ss},tlock_audtlock{ll+40,tt,ss}); % keeping tac centred but jitter aud (to match tlock_tac{X,tt})
        tlock_tacPaud{ll,tt,ss}.cov=tmp.cov;
        featurestats_tacPaud(:,ll,tt,ss,ii)=mean([featurestats_tac(:,ll+20,tt,ss,ii), featurestats_aud(:,ll+40,tt,ss,ii)]');
        tlock_tTacAlone{ll,tt,ss}=tlock_tactlock{ll+20,tt,ss}; % keeping tac centred but jitter aud (to match tlock_tac{X,tt})
        tlock_tAudAlone{ll,tt,ss}=tlock_audtlock{ll+40,tt,ss}; % keeping tac centred but jitter aud (to match tlock_tac{X,tt})
        if phaset0
          cfg.parameter='avg';
          for bb=1:4
            tlock_tacPaud_4mscon1_phasebin{ll,tt,ss,bb}=ft_math(cfg,tlock_tac4mscon1_phasebin{ll,tt,ss,bb},tlock_aud4mscon1_phasebin{ll,tt,ss,bb});
            tlock_tacPaud_4mscon2_phasebin{ll,tt,ss,bb}=ft_math(cfg,tlock_tac4mscon2_phasebin{ll,tt,ss,bb},tlock_aud4mscon2_phasebin{ll,tt,ss,bb});
            %                   tlock_tacPaud_4mscon1IAF_phasebin{ll,tt,ss,bb}=ft_math(cfg,tlock_tac4mscon1IAF_phasebin{ll,tt,ss,bb},tlock_aud4mscon1IAF_phasebin{ll,tt,ss,bb});
            %                   try
            %                     tlock_tacPaud_4mscon2IAF_phasebin{ll,tt,ss,bb}=ft_math(cfg,tlock_tac4mscon2IAF_phasebin{ll,tt,ss,bb},tlock_aud4mscon2IAF_phasebin{ll,tt,ss,bb});
            %                   catch
            %                     tlock_tacPaud_4mscon2IAF_phasebin{ll,tt,ss,bb}=[];
            %                   end
            
            tlock_tac4mscon1_phasebin{ll,tt,ss,bb}=rmfield(tlock_tac4mscon1_phasebin{ll,tt,ss,bb},{'var' 'dof'});
            tlock_tac4mscon2_phasebin{ll,tt,ss,bb}=rmfield(tlock_tac4mscon2_phasebin{ll,tt,ss,bb},{'var' 'dof'});
            tlock_aud4mscon1_phasebin{ll,tt,ss,bb}=rmfield(tlock_aud4mscon1_phasebin{ll,tt,ss,bb},{'var' 'dof'});
            tlock_aud4mscon2_phasebin{ll,tt,ss,bb}=rmfield(tlock_aud4mscon2_phasebin{ll,tt,ss,bb},{'var' 'dof'});
            %                   tlock_tac4mscon1IAF_phasebin{ll,tt,ss,bb}=rmfield(tlock_tac4mscon1IAF_phasebin{ll,tt,ss,bb},{'var' 'dof'});
            %                   tlock_tac4mscon2IAF_phasebin{ll,tt,ss,bb}=rmfield(tlock_tac4mscon2IAF_phasebin{ll,tt,ss,bb},{'var' 'dof'});
            %                   tlock_aud4mscon1IAF_phasebin{ll,tt,ss,bb}=rmfield(tlock_aud4mscon1IAF_phasebin{ll,tt,ss,bb},{'var' 'dof'});
            %                   try tlock_aud4mscon2IAF_phasebin{ll,tt,ss,bb}=rmfield(tlock_aud4mscon2IAF_phasebin{ll,tt,ss,bb},{'var' 'dof'}); catch end
            
            tlock_tacPaud_4mscon1_absbin{ll,tt,ss,bb}=ft_math(cfg,tlock_tac4mscon1_absbin{ll,tt,ss,bb},tlock_aud4mscon1_absbin{ll,tt,ss,bb});
            tlock_tacPaud_4mscon2_absbin{ll,tt,ss,bb}=ft_math(cfg,tlock_tac4mscon2_absbin{ll,tt,ss,bb},tlock_aud4mscon2_absbin{ll,tt,ss,bb});
            %                   tlock_tacPaud_4mscon1IAF_absbin{ll,tt,ss,bb}=ft_math(cfg,tlock_tac4mscon1IAF_absbin{ll,tt,ss,bb},tlock_aud4mscon1IAF_absbin{ll,tt,ss,bb});
            %                   tlock_tacPaud_4mscon2IAF_absbin{ll,tt,ss,bb}=ft_math(cfg,tlock_tac4mscon2IAF_absbin{ll,tt,ss,bb},tlock_aud4mscon2IAF_absbin{ll,tt,ss,bb});
            
            tlock_tac4mscon1_absbin{ll,tt,ss,bb}=rmfield(tlock_tac4mscon1_absbin{ll,tt,ss,bb},{'var' 'dof'});
            tlock_tac4mscon2_absbin{ll,tt,ss,bb}=rmfield(tlock_tac4mscon2_absbin{ll,tt,ss,bb},{'var' 'dof'});
            tlock_aud4mscon1_absbin{ll,tt,ss,bb}=rmfield(tlock_aud4mscon1_absbin{ll,tt,ss,bb},{'var' 'dof'});
            tlock_aud4mscon2_absbin{ll,tt,ss,bb}=rmfield(tlock_aud4mscon2_absbin{ll,tt,ss,bb},{'var' 'dof'});
            %                   tlock_tac4mscon1IAF_absbin{ll,tt,ss,bb}=rmfield(tlock_tac4mscon1IAF_absbin{ll,tt,ss,bb},{'var' 'dof'});
            %                   tlock_tac4mscon2IAF_absbin{ll,tt,ss,bb}=rmfield(tlock_tac4mscon2IAF_absbin{ll,tt,ss,bb},{'var' 'dof'});
            %                   tlock_aud4mscon1IAF_absbin{ll,tt,ss,bb}=rmfield(tlock_aud4mscon1IAF_absbin{ll,tt,ss,bb},{'var' 'dof'});
            %                   tlock_aud4mscon2IAF_absbin{ll,tt,ss,bb}=rmfield(tlock_aud4mscon2IAF_absbin{ll,tt,ss,bb},{'var' 'dof'});
          end
        end
      else
        tlock_tacPaud{ll,tt,ss}=[];
        tlock_tTacAlone{ll,tt,ss}=[]; % keeping tac centred but jitter aud (to match tlock_tac{X,tt})
        tlock_tAudAlone{ll,tt,ss}=[]; % keeping tac centred but jitter aud (to match tlock_tac{X,tt})
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
        cfg.parameter='cov';
        tmp=ft_math(cfg,tlock_tactlock{ll,tt,ss},tlock_nultlock{ll+50,tt,ss}); % TA+N
        tlock_tacMSpN{ll,tt,ss}.cov=tmp.cov;
        featurestats_tacMSpN(:,ll,tt,ss,ii)=mean([featurestats_tac(:,ll,tt,ss,ii), featurestats_nul(:,ll+50,tt,ss,ii)]');
        tlock_tMSAlone{ll,tt,ss}=tlock_tactlock{ll,tt,ss}; % TA+N
        tlock_tNulAlone{ll,tt,ss}=tlock_nultlock{ll+50,tt,ss}; % TA+N
        if phaset0
          cfg.parameter='avg';
          for bb=1:4
            tlock_tacMSpN_4mscon1_phasebin{ll,tt,ss,bb}=ft_math(cfg,tlock_ms4mscon1_phasebin{ll,tt,ss,bb},tlock_nul4mscon1_phasebin{ll,tt,ss,bb});
            tlock_tacMSpN_4mscon2_phasebin{ll,tt,ss,bb}=ft_math(cfg,tlock_ms4mscon2_phasebin{ll,tt,ss,bb},tlock_nul4mscon2_phasebin{ll,tt,ss,bb});
            %                   tlock_tacMSpN_4mscon1IAF_phasebin{ll,tt,ss,bb}=ft_math(cfg,tlock_ms4mscon1IAF_phasebin{ll,tt,ss,bb},tlock_nul4mscon1IAF_phasebin{ll,tt,ss,bb});
            %                   tlock_tacMSpN_4mscon2IAF_phasebin{ll,tt,ss,bb}=ft_math(cfg,tlock_ms4mscon2IAF_phasebin{ll,tt,ss,bb},tlock_nul4mscon2IAF_phasebin{ll,tt,ss,bb});
            
            tlock_ms4mscon1_phasebin{ll,tt,ss,bb}=rmfield(tlock_ms4mscon1_phasebin{ll,tt,ss,bb},{'var' 'dof'});
            tlock_ms4mscon2_phasebin{ll,tt,ss,bb}=rmfield(tlock_ms4mscon2_phasebin{ll,tt,ss,bb},{'var' 'dof'});
            tlock_nul4mscon1_phasebin{ll,tt,ss,bb}=rmfield(tlock_nul4mscon1_phasebin{ll,tt,ss,bb},{'var' 'dof'});
            tlock_nul4mscon2_phasebin{ll,tt,ss,bb}=rmfield(tlock_nul4mscon2_phasebin{ll,tt,ss,bb},{'var' 'dof'});
            %                   tlock_ms4mscon1IAF_phasebin{ll,tt,ss,bb}=rmfield(tlock_ms4mscon1IAF_phasebin{ll,tt,ss,bb},{'var' 'dof'});
            %                   tlock_ms4mscon2IAF_phasebin{ll,tt,ss,bb}=rmfield(tlock_ms4mscon2IAF_phasebin{ll,tt,ss,bb},{'var' 'dof'});
            %                   tlock_nul4mscon1IAF_phasebin{ll,tt,ss,bb}=rmfield(tlock_nul4mscon1IAF_phasebin{ll,tt,ss,bb},{'var' 'dof'});
            %                   tlock_nul4mscon2IAF_phasebin{ll,tt,ss,bb}=rmfield(tlock_nul4mscon2IAF_phasebin{ll,tt,ss,bb},{'var' 'dof'});
            
            tlock_tacMSpN_4mscon1_absbin{ll,tt,ss,bb}=ft_math(cfg,tlock_ms4mscon1_absbin{ll,tt,ss,bb},tlock_nul4mscon1_absbin{ll,tt,ss,bb});
            tlock_tacMSpN_4mscon2_absbin{ll,tt,ss,bb}=ft_math(cfg,tlock_ms4mscon2_absbin{ll,tt,ss,bb},tlock_nul4mscon2_absbin{ll,tt,ss,bb});
            %                   tlock_tacMSpN_4mscon1IAF_absbin{ll,tt,ss,bb}=ft_math(cfg,tlock_ms4mscon1IAF_absbin{ll,tt,ss,bb},tlock_nul4mscon1IAF_absbin{ll,tt,ss,bb});
            %                   tlock_tacMSpN_4mscon2IAF_absbin{ll,tt,ss,bb}=ft_math(cfg,tlock_ms4mscon2IAF_absbin{ll,tt,ss,bb},tlock_nul4mscon2IAF_absbin{ll,tt,ss,bb});
            
            tlock_ms4mscon1_absbin{ll,tt,ss,bb}=rmfield(tlock_ms4mscon1_absbin{ll,tt,ss,bb},{'var' 'dof'});
            tlock_ms4mscon2_absbin{ll,tt,ss,bb}=rmfield(tlock_ms4mscon2_absbin{ll,tt,ss,bb},{'var' 'dof'});
            tlock_nul4mscon1_absbin{ll,tt,ss,bb}=rmfield(tlock_nul4mscon1_absbin{ll,tt,ss,bb},{'var' 'dof'});
            tlock_nul4mscon2_absbin{ll,tt,ss,bb}=rmfield(tlock_nul4mscon2_absbin{ll,tt,ss,bb},{'var' 'dof'});
            %                   tlock_ms4mscon1IAF_absbin{ll,tt,ss,bb}=rmfield(tlock_ms4mscon1IAF_absbin{ll,tt,ss,bb},{'var' 'dof'});
            %                   tlock_ms4mscon2IAF_absbin{ll,tt,ss,bb}=rmfield(tlock_ms4mscon2IAF_absbin{ll,tt,ss,bb},{'var' 'dof'});
            %                   tlock_nul4mscon1IAF_absbin{ll,tt,ss,bb}=rmfield(tlock_nul4mscon1IAF_absbin{ll,tt,ss,bb},{'var' 'dof'});
            %                   tlock_nul4mscon2IAF_absbin{ll,tt,ss,bb}=rmfield(tlock_nul4mscon2IAF_absbin{ll,tt,ss,bb},{'var' 'dof'});
          end
        end
      else
        tlock_tacMSpN{ll,tt,ss}=[];
        tlock_tMSAlone{ll,tt,ss}=[]; % TA+N
        tlock_tNulAlone{ll,tt,ss}=[]; % TA+N
      end
      %           if numa_trials(ll,tt,ss)
      %             tlock_audMSpN{ll,tt,ss}=ft_math(cfg,tlock_audtlock{ll,tt,ss},tlock_nultlock{ll+60,tt,ss}); % AT+N
      %             featurestats_audMSpN(:,ll,tt,ss,ii)=mean([featurestats_aud(:,ll,tt,ss,ii), featurestats_nul(:,ll+60,tt,ss,ii)]');
      %           else
      %             tlock_audMSpN{ll,tt,ss}=[];
      %           end
      
      
      if pcmflag
        disp('pcmflag')
        mintrl=min([size(pcm_tactlock{ll,tt,ss}.trial,1) size(pcm_nultlock{ll+50,tt,ss}.trial,1) size(pcm_tactlock{ll+20,tt,ss}.trial,1) size(pcm_audtlock{ll+40,tt,ss}.trial,1)]);
        cfg=[];
        cfg.trials=1:mintrl;
        pcm_tactlock{ll,tt,ss}=ft_selectdata(cfg,pcm_tactlock{ll,tt,ss});
        pcm_tactlock{ll+20,tt,ss}=ft_selectdata(cfg,pcm_tactlock{ll+20,tt,ss});
        pcm_nultlock{ll+50,tt,ss}=ft_selectdata(cfg,pcm_nultlock{ll+50,tt,ss});
        pcm_audtlock{ll+40,tt,ss}=ft_selectdata(cfg,pcm_audtlock{ll+40,tt,ss});
        cfg=[];
        cfg.parameter='trial';
        cfg.operation='add';
        pcm_tacMSpN=ft_math(cfg,pcm_tactlock{ll,tt,ss},   pcm_nultlock{ll+50,tt,ss}); % AT + N
        pcm_tacPaud=ft_math(cfg,pcm_tactlock{ll+20,tt,ss},pcm_audtlock{ll+40,tt,ss}); % T + A
        cfg.operation='subtract';
        pcm_TPAmMSPN{ll,tt,ss}=ft_math(cfg,pcm_tacPaud,pcm_tacMSpN);
        %               pcm_TPAmMSPN{ll,tt,ss}.dimord='rpt_chan_time';
        %               pcm_TPAmMSPN{ll,tt,ss}.label=pcm_tacPaud.label;
        %               for jj=1:25 % every 20 ms
        %                 pcm_TPAmMSPN{ll,tt,ss}.trial(:,:,jj)=mean(TPAmMSPN.trial(:,:,jj*20-19:jj*20+1),3);
        %                 pcm_TPAmMSPN{ll,tt,ss}.time(jj)=TPAmMSPN.time(jj*20-10);
        %               end
        clear TPAmMSPN pcm_tacMSpN pcm_tacPaud
      end
      
    end % end ll
    
    if synchasynch
      for ll=[1 3 4]
        if ~isempty(tlock_tacMSshift1{ll,tt,ss}) && size(tlock_tacMSshift1{ll,tt,ss}.trial,1)
          cfg=[];
          cfg.lpfilter='yes';
          cfg.lpfreq=40;
          cfg.demean='yes';
          cfg.baselinewindow=[-.4 0]-.1; % can use same baseline for all, as shifting has already occured
          %               if ll==1
          %                 cfg.baselinewindow=[-.4 0]-.6;
          %               elseif ll==3
          %                 cfg.baselinewindow=[-.4 0]-.17;
          %               elseif ll==4
          %                 cfg.baselinewindow=[-.4 0]-.12;
          %               elseif ll==5
          %                 cfg.baselinewindow=[-.4 0]-.1;
          %               end
          
          disp('ft_preprocessing tac')
          tlock_tacMSshift1{ll,tt,ss}=ft_preprocessing(cfg,tlock_tacMSshift1{ll,tt,ss});
          tlock_tacMSshift2{ll,tt,ss}=ft_preprocessing(cfg,tlock_tacMSshift2{ll,tt,ss});
          tlock_tacMSshift3{ll,tt,ss}=ft_preprocessing(cfg,tlock_tacMSshift3{ll,tt,ss});
          tlock_tacMSshift4{ll,tt,ss}=ft_preprocessing(cfg,tlock_tacMSshift4{ll,tt,ss});
          featurestats_tacMSshift1(:,ll,tt,ss,ii)=mean(tlock_tacMSshift1{ll,tt,ss}.trialinfo(:,12:19));
          featurestats_tacMSshift2(:,ll,tt,ss,ii)=mean(tlock_tacMSshift2{ll,tt,ss}.trialinfo(:,12:19));
          featurestats_tacMSshift3(:,ll,tt,ss,ii)=mean(tlock_tacMSshift3{ll,tt,ss}.trialinfo(:,12:19));
          featurestats_tacMSshift4(:,ll,tt,ss,ii)=mean(tlock_tacMSshift4{ll,tt,ss}.trialinfo(:,12:19));
          cfg=[];
          cfg.covariance='yes';
          cfg.covariancewindow=[-0.1 0.7]; % a full window valid for all ll
          tlock_tacMSshift1tlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_tacMSshift1{ll,tt,ss});
          tlock_tacMSshift2tlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_tacMSshift2{ll,tt,ss});
          tlock_tacMSshift3tlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_tacMSshift3{ll,tt,ss});
          tlock_tacMSshift4tlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_tacMSshift4{ll,tt,ss});
        else
          tlock_tacMSshift1tlock{ll,tt,ss}=[];
          tlock_tacMSshift2tlock{ll,tt,ss}=[];
          tlock_tacMSshift3tlock{ll,tt,ss}=[];
          tlock_tacMSshift4tlock{ll,tt,ss}=[];
        end
        
        if ~isempty(tlock_tacMSshift1tlock{ll,tt,ss})
          % create sum of asynchronous multisensory conditions (e.g. AT70_shifted plus TA70)
          cfg=[];
          cfg.operation='add';
          cfg.parameter='avg';
          tlock_tMSasynch{ll,tt,ss}=ft_math(cfg,tlock_tacMSshift1tlock{ll,tt,ss},tlock_tacMSshift2tlock{ll,tt,ss});
          cfg.parameter='cov';
          tmp=ft_math(cfg,tlock_tacMSshift1tlock{ll,tt,ss},tlock_tacMSshift2tlock{ll,tt,ss});
          tlock_tMSasynch{ll,tt,ss}.cov=tmp.cov;
          featurestats_tMSasynch(:,ll,tt,ss,ii)=mean([featurestats_tacMSshift1(:,ll,tt,ss,ii), featurestats_tacMSshift2(:,ll,tt,ss,ii)]');
          % create sum of synchronous multisensory conditions (e.g. TA0_shifted plus TA0)
          cfg=[];
          cfg.operation='add';
          cfg.parameter='avg';
          tlock_tMSsynch{ll,tt,ss}=ft_math(cfg,tlock_tacMSshift3tlock{ll,tt,ss},tlock_tacMSshift4tlock{ll,tt,ss});
          cfg.parameter='cov';
          tmp=ft_math(cfg,tlock_tacMSshift3tlock{ll,tt,ss},tlock_tacMSshift4tlock{ll,tt,ss});
          tlock_tMSsynch{ll,tt,ss}.cov=tmp.cov;
          featurestats_tMSsynch(:,ll,tt,ss,ii)=mean([featurestats_tacMSshift3(:,ll,tt,ss,ii), featurestats_tacMSshift4(:,ll,tt,ss,ii)]');
        else
          tlock_tMSasynch{ll,tt,ss}=[];
          tlock_tMSsynch{ll,tt,ss}=[];
        end
      end % ll
    end % synchasynch
    
    
  end  % end ss
  
  %         if binonly
  %           continue
  %         end
  
  
  %clear tlock_tac tlock_aud tlock_nul tr tlock*tlock
  clear tlock_tac tlock_aud tlock_nul
  clear tlock_tacMSshift1 tlock_tacMSshift2 tlock_tacMSshift3 tlock_tacMSshift4
  if ~tacaloneproc && ~phaset0 && ~pcmflag
    if isempty(tr.stageuse)
      save(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'tr')
    else
      save(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'tlock_*P*','tlock_*N*','tlock_t*','tlock*tlock','num*trials','featurestats_*')
    end
  elseif phaset0 && phaset0use>0
    if binonly
    else
      save(['tlock_ERPphasebin_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'tlock_*con*phasebin','tlock_*con*absbin','tr','phaset0use','freqs*');
    end
  elseif pcmflag
    save(['tlock_pcm_' sub{ii} '_sleep' num2str(sleep) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'pcm_TPAmMSPN','tr');
  end
  %       end % end iter
  %     end % end tt
  
  continue
  
  %% Do AudPlusTac second
  tacaud=0;
  for tt=ttuse % refers to lim=[no-limit .01 .005 .002];
    %   profile on
    %       [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection2_wakeSleep_sepTacAud(ii,sleep,tt,tacaud);
    [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection2_wakeSleep_sepTacAud(ii,sleep,tt,tacaud,raw_tac, raw_aud, raw_nul);
    if tt==ttuse(end)
      clear raw_* % not needed anymore
    end
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
            %               cfg.baselinewindow=[-1.7 -0.6];
            if ll==49
              cfg.baselinewindow=[-.4 0]-.6;
            elseif ll==47
              cfg.baselinewindow=[-.4 0]-.17;
            elseif ll==46
              cfg.baselinewindow=[-.4 0]-.12;
            elseif ll==45 || ll==44 || ll==43 || ll==41
              cfg.baselinewindow=[-.4 0]-.1;
            end
            disp('ft_preprocessing tac')
            tlock_tac{ll,tt,ss}=ft_preprocessing(cfg,tlock_tac{ll,tt,ss});
            featurestats_tac(:,ll,tt,ss,ii)=mean(tlock_tac{ll,tt,ss}.trialinfo(:,12:19));
            cfg=[];
            cfg.covariance='yes';
            cfg.covariancewindow=[-0.6 0.7]; % a full window valid for all ll
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
            %               cfg.baselinewindow=[-1.7 -0.6];
            if ll==9 || ll==29
              cfg.baselinewindow=[-.4 0]-.6;
            elseif ll==7 || ll==27
              cfg.baselinewindow=[-.4 0]-.17;
            elseif ll==6 || ll==26
              cfg.baselinewindow=[-.4 0]-.12;
            elseif ll==5 || ll==4 || ll==3 || ll==1 || ll==25 || ll==24 || ll==23 || ll==21
              cfg.baselinewindow=[-.4 0]-.1;
            end
            disp('ft_preprocessing aud')
            tlock_aud{ll,tt,ss}=ft_preprocessing(cfg,tlock_aud{ll,tt,ss});
            featurestats_aud(:,ll,tt,ss,ii)=mean(tlock_aud{ll,tt,ss}.trialinfo(:,12:19));
            cfg=[];
            cfg.covariance='yes';
            cfg.covariancewindow=[-0.6 0.7]; % a full window valid for all ll
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
          %             cfg.baselinewindow=[-1.7 -0.6];
          if ll==69
            cfg.baselinewindow=[-.4 0]-.6;
          elseif ll==67
            cfg.baselinewindow=[-.4 0]-.17;
          elseif ll==66
            cfg.baselinewindow=[-.4 0]-.12;
          elseif ll==65 || ll==64 || ll==63 || ll==61
            cfg.baselinewindow=[-.4 0]-.1;
          end
          %           cfg.bpfilter='yes';
          %           cfg.bpfreq=[1 40];
          disp('ft_preprocessing nul')
          tlock_nul{ll,tt,ss}=ft_preprocessing(cfg,tlock_nul{ll,tt,ss});
          featurestats_nul(:,ll,tt,ss,ii)=mean(tlock_nul{ll,tt,ss}.trialinfo(:,12:19));
          cfg=[];
          cfg.covariance='yes';
          cfg.covariancewindow=[-0.6 0.7]; % a full window valid for all ll
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
          cfg.parameter='cov';
          tmp=ft_math(cfg,tlock_audtlock{ll+20,tt,ss},tlock_tactlock{ll+40,tt,ss}); % keeping aud centred but jitter tac (to match tlock_aud{X,tt})
          tlock_audPtac{ll,tt,ss}.cov=tmp.cov;
          featurestats_audPtac(:,ll,tt,ss,ii)=mean([featurestats_aud(:,ll+20,tt,ss,ii), featurestats_tac(:,ll+40,tt,ss,ii)]');
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
          cfg.parameter='cov';
          tmp=ft_math(cfg,tlock_audtlock{ll,tt,ss},tlock_nultlock{ll+60,tt,ss}); % AT+N
          tlock_audMSpN{ll,tt,ss}.cov=tmp.cov;
          featurestats_audMSpN(:,ll,tt,ss,ii)=mean([featurestats_aud(:,ll,tt,ss,ii), featurestats_nul(:,ll+60,tt,ss,ii)]');
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
      
      for ll=[1 3 4]
        if ~isempty(tlock_audMSshift1{ll,tt,ss})&& size(tlock_audMSshift1{ll,tt,ss}.trial,1)
          cfg=[];
          cfg.lpfilter='yes';
          cfg.lpfreq=40;
          cfg.demean='yes';
          cfg.baselinewindow=[-.4 0]-.1; % can use same baseline for all, as shifting has already occured
          %               if ll==1
          %                 cfg.baselinewindow=[-.4 0]-.6;
          %               elseif ll==3
          %                 cfg.baselinewindow=[-.4 0]-.17;
          %               elseif ll==4
          %                 cfg.baselinewindow=[-.4 0]-.12;
          %               elseif ll==5
          %                 cfg.baselinewindow=[-.4 0]-.1;
          %               end
          
          disp('ft_preprocessing aud')
          tlock_audMSshift1{ll,tt,ss}=ft_preprocessing(cfg,tlock_audMSshift1{ll,tt,ss});
          tlock_audMSshift2{ll,tt,ss}=ft_preprocessing(cfg,tlock_audMSshift2{ll,tt,ss});
          tlock_audMSshift3{ll,tt,ss}=ft_preprocessing(cfg,tlock_audMSshift3{ll,tt,ss});
          tlock_audMSshift4{ll,tt,ss}=ft_preprocessing(cfg,tlock_audMSshift4{ll,tt,ss});
          featurestats_audMSshift1(:,ll,tt,ss,ii)=mean(tlock_audMSshift1{ll,tt,ss}.trialinfo(:,12:19));
          featurestats_audMSshift2(:,ll,tt,ss,ii)=mean(tlock_audMSshift2{ll,tt,ss}.trialinfo(:,12:19));
          featurestats_audMSshift3(:,ll,tt,ss,ii)=mean(tlock_audMSshift3{ll,tt,ss}.trialinfo(:,12:19));
          featurestats_audMSshift4(:,ll,tt,ss,ii)=mean(tlock_audMSshift4{ll,tt,ss}.trialinfo(:,12:19));
          cfg=[];
          cfg.covariance='yes';
          cfg.covariancewindow=[-0.1 0.7]; % a full window valid for all ll
          tlock_audMSshift1tlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_audMSshift1{ll,tt,ss});
          tlock_audMSshift2tlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_audMSshift2{ll,tt,ss});
          tlock_audMSshift3tlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_audMSshift3{ll,tt,ss});
          tlock_audMSshift4tlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_audMSshift4{ll,tt,ss});
        else
          tlock_audMSshift1tlock{ll,tt,ss}=[];
          tlock_audMSshift2tlock{ll,tt,ss}=[];
          tlock_audMSshift3tlock{ll,tt,ss}=[];
          tlock_audMSshift4tlock{ll,tt,ss}=[];
        end
        
        if ~isempty(tlock_audMSshift1tlock{ll,tt,ss})
          % create sum of asynchronous multisensory conditions (e.g. AT70_shifted plus TA70)
          cfg=[];
          cfg.operation='add';
          cfg.parameter='avg';
          tlock_aMSasynch{ll,tt,ss}=ft_math(cfg,tlock_audMSshift1tlock{ll,tt,ss},tlock_audMSshift2tlock{ll,tt,ss});
          cfg.parameter='cov';
          tmp=ft_math(cfg,tlock_audMSshift1tlock{ll,tt,ss},tlock_audMSshift2tlock{ll,tt,ss});
          tlock_aMSasynch{ll,tt,ss}.cov=tmp.cov;
          featurestats_aMSasynch(:,ll,tt,ss,ii)=mean([featurestats_audMSshift1(:,ll,tt,ss,ii), featurestats_audMSshift2(:,ll,tt,ss,ii)]');
          % create sum of synchronous multisensory conditions (e.g. TA0_shifted plus TA0)
          cfg=[];
          cfg.operation='add';
          cfg.parameter='avg';
          tlock_aMSsynch{ll,tt,ss}=ft_math(cfg,tlock_audMSshift3tlock{ll,tt,ss},tlock_audMSshift4tlock{ll,tt,ss});
          cfg.parameter='cov';
          tmp=ft_math(cfg,tlock_audMSshift3tlock{ll,tt,ss},tlock_audMSshift4tlock{ll,tt,ss});
          tlock_aMSsynch{ll,tt,ss}.cov=tmp.cov;
          featurestats_aMSsynch(:,ll,tt,ss,ii)=mean([featurestats_audMSshift3(:,ll,tt,ss,ii), featurestats_audMSshift4(:,ll,tt,ss,ii)]');
        else
          tlock_aMSasynch{ll,tt,ss}=[];
          tlock_aMSsynch{ll,tt,ss}=[];
        end
      end % ll
      
    end  % end ss
    
    %clear tlock_tac tlock_aud tlock_nul tr tlock*tlock
    clear tlock_tac tlock_aud tlock_nul
    save(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_trialkc' num2str(trialkc) '.mat'],'tlock_*P*','tlock_*N*','tlock_a*','tlock*tlock','num*trials','featurestats_*')
  end % end tt
  
  
  
  %%
  
  
  
  %   save(['tlock_diffs_averef_' sub{ii} '.mat'],'*trialkept','tlock_*P*','tlock_*N*')
  % save(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '.mat'],'tlock_*P*','tlock_*N*','tlock_t*','tlock_a*','num*trials','featurestats_*')
  if ~exist(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'file')
    save(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'tlock_*P*','tlock_*N*','tlock_t*','tlock_a*','tlock*tlock','num*trials','featurestats_*')
  end
  
  
end % ii
% end % sleep

return





