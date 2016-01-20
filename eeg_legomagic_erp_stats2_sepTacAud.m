% statistics of EEG awake data

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
  rmpath(genpath('/mnt/hgfs/D/matlab/spm8/external/fieldtrip/'))
  rmpath(genpath('/mnt/hgfs/D/fieldtrip_svn/'))
  addpath('/mnt/hgfs/D/fieldtrip_svn/')
end
which ft_defaults.m
ft_defaults;

load([edir 'iikeep.mat'])
soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];

%%
usetr=1;
ttuse=[3];
tacaloneproc=0;
phaset0=1;
phaset0use=2; % 0 = compute Hilbert; 1 = FFT Cz, 2 = FFT FC, 3 = FFT PCA, 4 = Hilb Cz, 5 = Hilb FC, 6 = Hilb PCA
synchasynch=0;
binonly=1;

for sleep=[1]
  if sleep
    iiuse=setdiff(iiBuse,3:7);
    %     iiuse=setdiff(iiBuse,[3:7 8]);
  else
    iiuse=iiSuse;
    %         iiuse=setdiff(iiSuse,1:9);
    %         iiuse=8;
  end
  freqsub=nan(1,max(iiuse));
  for ii=iiuse
    cd([edir sub{ii} ])
    clearvars -except ii sub edir ddir ii*use sleep featurestats* ttuse soades usetr tacaloneproc synchasynch phaset0* binonly
    [raw_tac, raw_aud, raw_nul, artflag]=eeg_legomagic_epoching2(ii,sleep,1,0); % featfull=1, saveflag =0;
    
    %   try
    %     load(['tlock_trialSel_' sub{ii} '.mat']);
    %   catch
    
    %     [tlock_tac_s0,tlock_aud_s0,tlock_nul_s0,tr]=eeg_legomagic_trialSelection1(ii);
    %     save(['tlock_trialSel_' sub{ii} '.mat'],'tlock*s0','tr','-v7.3');
    %     [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection1_wakeSleep(ii,sleep);
    
    if sleep
      trialkc=-1;  % -1 (use all), 0 (only non-Kc), or 1 (only Kc)
      iteruse=11;
    elseif sleep==0
      trialkc=-1; % always -1
      iteruse=27; % for main analysis
      %       iteruse=30; % for tacalone analysis
    end
    
    % For memory reasons, do TacPlusAud separately from AudPlusTac
    %% Do TacPlusAud first
    tacaud=1; % tacaud=1 means triggered on tactile; tacaud=0 means triggered on auditory
    for tt=ttuse % refers to lim=[no-limit .01 .005 .002];
      %   profile on
      %       [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection2_wakeSleep_sepTacAud(ii,sleep,tt,tacaud);
      
      % % iter 1:10 same but with trialSelection eog of prestim (13) only
      % % iter 11:20 same but with trialSelection eog of all trial (13 and 23)
      % % iter 21:30 now 'epoching' uses  'artfctdef_sit_automan' rather than 'auto';
      for iter=iteruse
        if ~isempty(whos('raw*'))
          [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection2_wakeSleep_sepTacAud(ii,sleep,tt,tacaud,raw_tac, raw_aud, raw_nul, iter, usetr, trialkc);
        else
          error('with sleep=1 cant run more than one iter or one tt at once for memory reasons')
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
          ssuse=[tr.stageuse+10 23]; % not enough memory for 23, plus is it useful?
          %           ssuse=[tr.stageuse+10];
          ssuse=[12];
          clear raw*
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
        
        if binonly
          continue
        end
        
        %clear tlock_tac tlock_aud tlock_nul tr tlock*tlock
        clear tlock_tac tlock_aud tlock_nul
        clear tlock_tacMSshift1 tlock_tacMSshift2 tlock_tacMSshift3 tlock_tacMSshift4
        if ~tacaloneproc && ~phaset0
          if isempty(tr.stageuse)
            save(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'tr')
          else
            save(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'tlock_*P*','tlock_*N*','tlock_t*','tlock*tlock','num*trials','featurestats_*')
          end
        elseif phaset0 && phaset0use>0
          save(['tlock_ERPphasebin_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'tlock_*con*phasebin','tlock_*con*absbin','tr','phaset0use','freqs*');
        end
      end % end iter
    end % end tt
    
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
end % sleep

return


%%  Testing the power distribution of delta prestimulus

sleep=1;
if sleep==1
  subuse=setdiff(iiBuse,3:7);
  iter=11;
  trialkc=-1; % -1 all, 0 no Kc, 1 only Kc
elseif sleep==0
  subuse=iiSuse;
  iter=27;
end
audbineach=nan(4,9,19);
tacbineach=nan(4,9,19);
nulbineach=nan(4,9,19);
ms1bineach=nan(4,9,19);
soalist=[1 3 4 5 6 7 9];
subuseind=0;
for ii=subuse
  subuseind=subuseind+1;
  cd([edir sub{ii} ])
  %   try
  load(['absbin_' sub{ii} '_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat']);
  %   catch
  %     if trialkc==0
  %       load(['absbin_' sub{ii} '_sleep' num2str(sleep) '.mat']);
  %     else
  %       error('where is the file?')
  %     end
  %   end
  
  poweach4{subuseind}=[aud4mscon1abs{5,3,12}; tac4mscon1abs{5,3,12}; nul4mscon1abs{5,3,12}; ms4mscon1abs{5,3,12}]';
  poweach10{subuseind}=[aud4mscon1abs{5,3,12}; tac4mscon1abs{5,3,12}; nul4mscon1abs{5,3,12}; ms4mscon1abs{1,3,12}; ms4mscon1abs{3,3,12}; ms4mscon1abs{4,3,12}; ms4mscon1abs{5,3,12}; ms4mscon1abs{6,3,12}; ms4mscon1abs{7,3,12}; ms4mscon1abs{9,3,12}]';
  poweachIAF4{subuseind}=[aud4mscon1absIAF{5,3,12}; tac4mscon1absIAF{5,3,12}; nul4mscon1absIAF{5,3,12}; ms4mscon1absIAF{5,3,12}]';
  poweachIAF10{subuseind}=[aud4mscon1absIAF{5,3,12}; tac4mscon1absIAF{5,3,12}; nul4mscon1absIAF{5,3,12}; ms4mscon1absIAF{1,3,12}; ms4mscon1absIAF{3,3,12}; ms4mscon1absIAF{4,3,12}; ms4mscon1absIAF{5,3,12}; ms4mscon1absIAF{6,3,12}; ms4mscon1absIAF{7,3,12}; ms4mscon1absIAF{9,3,12}]';
  
  
  audeach{subuseind}=aud4mscon1abs{5,3,12}';
  taceach{subuseind}=tac4mscon1abs{5,3,12}';
  nuleach{subuseind}=nul4mscon1abs{5,3,12}';
  ms1each{subuseind}=ms4mscon1abs{1,3,12}';
  ms3each{subuseind}=ms4mscon1abs{3,3,12}';
  ms4each{subuseind}=ms4mscon1abs{4,3,12}';
  ms5each{subuseind}=ms4mscon1abs{5,3,12}';
  ms6each{subuseind}=ms4mscon1abs{6,3,12}';
  ms7each{subuseind}=ms4mscon1abs{7,3,12}';
  ms9each{subuseind}=ms4mscon1abs{9,3,12}';
  
  for ll=soalist
    audbineach(:,ll,subuseind)=mean(aud4mscon1absbin{ll,3,12});
    tacbineach(:,ll,subuseind)=mean(tac4mscon1absbin{ll,3,12});
    nulbineach(:,ll,subuseind)=mean(nul4mscon1absbin{ll,3,12});
    ms1bineach(:,ll,subuseind)=mean(ms4mscon1absbin{ll,3,12});
  end
  
end % ii

thresh=.05;
mean(audbineach<thresh,3)
mean(tacbineach<thresh,3)
mean(nulbineach<thresh,3)
mean(ms1bineach<thresh,3)


for ii=1:length(poweach4)
  figure(1);
  subplot(4,5,ii);hist(poweach4{ii},[0:10 20]);
  figure(2);
  subplot(4,5,ii);hist(poweach10{ii},[0:10 20]);
  figure(3);
  subplot(4,5,ii);hist(poweachIAF4{ii},[0:10 20]);
  figure(4);
  subplot(4,5,ii);hist(poweachIAF10{ii},[0:10 20]);
end
figure;hist([poweach10{:}],[0:10 20])
figure;hist([poweachIAF10{:}],[0:10 20])

sortpow=sort([poweach10{:}]);
sortpow(round(1*length(sortpow)/4))
%     0.9973
sortpow(round(2*length(sortpow)/4))
%     1.6479
sortpow(round(3*length(sortpow)/4))
%     2.5349

% % with Kcomplexes (trialkc = -1)
%     1.1966
%     2.0076
%     3.3088

figure;
subplot(2,5,1);hist([nuleach{:}],[0:10 20]);
subplot(2,5,2);hist([taceach{:}],[0:10 20]);
subplot(2,5,3);hist([audeach{:}],[0:10 20]);
subplot(2,5,4);hist([ms1each{:}],[0:10 20]);
subplot(2,5,5);hist([ms3each{:}],[0:10 20]);
subplot(2,5,6);hist([ms4each{:}],[0:10 20]);
subplot(2,5,7);hist([ms5each{:}],[0:10 20]);
subplot(2,5,8);hist([ms6each{:}],[0:10 20]);
subplot(2,5,9);hist([ms7each{:}],[0:10 20]);
subplot(2,5,10);hist([ms9each{:}],[0:10 20]);





%% Group level: awake and asleep testing tactile alone response
% This is not for the main finding, but rather testing strange
% finding of no seeming Tactile response for N23 bed data.

plotflag=1;
printflag=1;
statsflag=0;
tacaud=1; % tacaud=1 means triggered on tactile; tacaud=0 means triggered on auditory

tt=3;
sleepss_tacAll=zeros(13,11,max([iiSuse iiBuse]));
sleepss_tac19T=zeros(13,11,max([iiSuse iiBuse]));
sleepss_tacAlone=zeros(13,11,max([iiSuse iiBuse]));
for sleep=[1]
  clearvars -except tt sub edir ddir ii*use *flag sleep* tacaud
  if sleep==1
    subuse=setdiff(iiBuse,3:7);
    iter=11;
    trialkc=-1; % -1 all, 0 no Kc, 1 only Kc
  elseif sleep==0
    subuse=iiSuse;
    iter=27;
  end
  subind=1;
  for ii=subuse
    cd([edir sub{ii} ])
    %   load(['tlock_diffs_' sub{ii} '.mat']);
    %       load(['tlock_diffs_averef_' sub{ii} '.mat']);
    %     load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '.mat'])
    if sleep
      load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_trialkc' num2str(trialkc) '.mat'])
      tk0=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(0) '_trialkc' num2str(trialkc) '.mat']);
      tk1=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_trialkc' num2str(trialkc) '.mat']);
      ssuse=tk1.tr.stageuse;
    else
      load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud1' '_iter' num2str(iter) '.mat'])
      tk1=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter) '.mat']);
      ssuse=tk1.tr.stageuse;
    end
    %         load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '.mat'],'tr')
    
    for ss=ssuse+10
      try
        sleepss_tacAll(ss,sleep+10,subind)=tlock_tacAll{tt,ss}.dof(1,dsearchn(tlock_tacAll{tt,ss}.time',0));
        %         if sleepss_tacAll(ss,sleep+10,subind)>=20
        tlock_tacAll_each{ss,sleep+10,subind}=tlock_tacAll{tt,ss};
        %         end
      end
      try
        sleepss_tac19T(ss,sleep+10,subind)=tlock_tac19T{tt,ss}.dof(1,dsearchn(tlock_tac19T{tt,ss}.time',0));
        %         if sleepss_tac19T(ss,sleep+10,subind)>=20
        tlock_tac19T_each{ss,sleep+10,subind}=tlock_tac19T{tt,ss};
        %         end
      end
      try
        sleepss_tacAlone(ss,sleep+10,subind)=tlock_tTacAlone{5,tt,ss}.dof(1,dsearchn(tlock_tTacAlone{5,tt,ss}.time',0));
        %         if sleepss_tacAlone(ss,sleep+10,subind)>=20
        tlock_tacAlone_each{ss,sleep+10,subind}=tlock_tTacAlone{5,tt,ss}; % the 'll' is arbitrary here
        %         end
      end
    end % ss
    
    clear tlock_tacAll tlock_tac19T tlock_tTacAlone tlock_a* tlock_tacMSpN tlock_tacPaud tlock_t*Alone
  end % sleep
  % if any(any(sleepss_tacAll(:,:,subind))) || any(any(sleepss_tac19T(:,:,subind))) || any(any(sleepss_Alone(:,:,subind)))
  subind=subind+1;
end % ii
subind=subind-1;

thresh=10;
for ii=1:subind
  sleep=0;
  ss=10;
  if sleepss_tacAll(ss,sleep+10,ii)>=thresh
    tacAll_s0{ii}=tlock_tacAll_each{ss,sleep+10,ii};
  else
    tacAll_s0{ii}=[];
  end
  if sleepss_tac19T(ss,sleep+10,ii)>=thresh
    tac19T_s0{ii}=tlock_tac19T_each{ss,sleep+10,ii};
  else
    tac19T_s0{ii}=[];
  end
  if sleepss_tacAlone(ss,sleep+10,ii)>=thresh
    tacAlone_s0{ii}=tlock_tacAlone_each{ss,sleep+10,ii};
  else
    tacAlone_s0{ii}=[]
  end
  
  sleep=0;
  ss=11;
  if sleepss_tacAll(ss,sleep+10,ii)>=thresh
    tacAll_s1{ii}=tlock_tacAll_each{ss,sleep+10,ii};
  else
    tacAll_s1{ii}=[];
  end
  if sleepss_tac19T(ss,sleep+10,ii)>=thresh
    tac19T_s1{ii}=tlock_tac19T_each{ss,sleep+10,ii};
  else
    tac19T_s1{ii}=[];
  end
  if sleepss_tacAlone(ss,sleep+10,ii)>=thresh
    tacAlone_s1{ii}=tlock_tacAlone_each{ss,sleep+10,ii};
  else
    tacAlone_s1{ii}=[];
  end
  
  sleep=0;
  ss=12;
  if sleepss_tacAll(ss,sleep+10,ii)>=thresh
    tacAll_s2{ii}=tlock_tacAll_each{ss,sleep+10,ii};
  else
    tacAll_s2{ii}=[];
  end
  if sleepss_tac19T(ss,sleep+10,ii)>=thresh
    tac19T_s2{ii}=tlock_tac19T_each{ss,sleep+10,ii};
  else
    tac19T_s2{ii}=[];
  end
  if sleepss_tacAlone(ss,sleep+10,ii)>=thresh
    tacAlone_s2{ii}=tlock_tacAlone_each{ss,sleep+10,ii};
  else
    tacAlone_s2{ii}=[];
  end
  
  sleep=1;
  ss=10;
  if sleepss_tacAll(ss,sleep+10,ii)>=thresh
    tacAll_b0{ii}=tlock_tacAll_each{ss,sleep+10,ii};
  else
    tacAll_b0{ii}=[];
  end
  if sleepss_tac19T(ss,sleep+10,ii)>=thresh
    tac19T_b0{ii}=tlock_tac19T_each{ss,sleep+10,ii};
  else
    tac19T_b0{ii}=[];
  end
  if sleepss_tacAlone(ss,sleep+10,ii)>=thresh
    tacAlone_b0{ii}=tlock_tacAlone_each{ss,sleep+10,ii};
  else
    tacAlone_b0{ii}=[];
  end
  
  sleep=1;
  ss=11;
  if sleepss_tacAll(ss,sleep+10,ii)>=thresh
    tacAll_b1{ii}=tlock_tacAll_each{ss,sleep+10,ii};
  else
    tacAll_b1{ii}=[];
  end
  if sleepss_tac19T(ss,sleep+10,ii)>=thresh
    tac19T_b1{ii}=tlock_tac19T_each{ss,sleep+10,ii};
  else
    tac19T_b1{ii}=[];
  end
  if sleepss_tacAlone(ss,sleep+10,ii)>=thresh
    tacAlone_b1{ii}=tlock_tacAlone_each{ss,sleep+10,ii};
  else
    tacAlone_b1{ii}=[];
  end
  
  sleep=1;
  ss=12;
  if sleepss_tacAll(ss,sleep+10,ii)>=thresh
    tacAll_b2{ii}=tlock_tacAll_each{ss,sleep+10,ii};
  else
    tacAll_b2{ii}=[];
  end
  if sleepss_tac19T(ss,sleep+10,ii)>=thresh
    tac19T_b2{ii}=tlock_tac19T_each{ss,sleep+10,ii};
  else
    tac19T_b2{ii}=[];
  end
  if sleepss_tacAlone(ss,sleep+10,ii)>=thresh
    tacAlone_b2{ii}=tlock_tacAlone_each{ss,sleep+10,ii};
  else
    tacAlone_b2{ii}=[];
  end
  
  sleep=1;
  ss=13;
  if sleepss_tacAll(ss,sleep+10,ii)>=thresh
    tacAll_b3{ii}=tlock_tacAll_each{ss,sleep+10,ii};
  else
    tacAll_b3{ii}=[];
  end
  if sleepss_tac19T(ss,sleep+10,ii)>=thresh
    tac19T_b3{ii}=tlock_tac19T_each{ss,sleep+10,ii};
  else
    tac19T_b3{ii}=[];
  end
  if sleepss_tacAlone(ss,sleep+10,ii)>=thresh
    tacAlone_b3{ii}=tlock_tacAlone_each{ss,sleep+10,ii};
  else
    tacAlone_b3{ii}=[];
  end
end % ii


% insert these into lines below.
figure(100);
figcfg=[];
figcfg.xlim=[-0.7 1.1];
figcfg.channel={'Fz' 'Cz' 'F1' 'F2' 'FC1' 'FC2' 'C1' 'C2'};
figcfg.ylim=[-7 7];

cfg=[];
tmp={tacAll_s0{squeeze(sleepss_tacAll(10,10,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAll_s0=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,1);ft_singleplotER(figcfg,grave_tacAll_s0);
  ylabel('TacAll')
  title(num2str(sum(squeeze(sleepss_tacAll(10,10,:)))))
else
end
cfg=[];
tmp={tacAll_s1{squeeze(sleepss_tacAll(11,10,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAll_s1=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,3);ft_singleplotER(figcfg,grave_tacAll_s1);
  title(num2str(sum(squeeze(sleepss_tacAll(11,10,:)))))
else
end
cfg=[];
tmp={tacAll_s2{squeeze(sleepss_tacAll(12,10,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAll_s2=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,5);ft_singleplotER(figcfg,grave_tacAll_s2);
  title(num2str(sum(squeeze(sleepss_tacAll(12,10,:)))))
else
end
% cfg=[];
% tmp={tacAll_s3{squeeze(sleepss_tacAll(13,10,:))>=thresh}};
% grave_tacAll_s3=ft_timelockgrandaverage(cfg,tmp{:});
% subplot(3,8,7);ft_singleplotER(figcfg,grave_tacAll_s3);

cfg=[];
tmp={tacAll_b0{squeeze(sleepss_tacAll(10,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAll_b0=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,2);ft_singleplotER(figcfg,grave_tacAll_b0);
  title(num2str(sum(squeeze(sleepss_tacAll(10,11,:)))))
else
end
cfg=[];
tmp={tacAll_b1{squeeze(sleepss_tacAll(11,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAll_b1=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,4);ft_singleplotER(figcfg,grave_tacAll_b1);
  title(num2str(sum(squeeze(sleepss_tacAll(11,11,:)))))
else
end
cfg=[];
tmp={tacAll_b2{squeeze(sleepss_tacAll(12,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAll_b2=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,6);ft_singleplotER(figcfg,grave_tacAll_b2);
  title(num2str(sum(squeeze(sleepss_tacAll(12,11,:)))))
else
end
cfg=[];
tmp={tacAll_b3{squeeze(sleepss_tacAll(13,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAll_b3=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,8);ft_singleplotER(figcfg,grave_tacAll_b3);
  title(num2str(sum(squeeze(sleepss_tacAll(13,11,:)))))
else
end

cfg=[];
tmp={tac19T_s0{squeeze(sleepss_tac19T(10,10,:))>=thresh}};
if ~isempty(tmp)
  grave_tac19T_s0=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,9); ft_singleplotER(figcfg,grave_tac19T_s0);
  title(num2str(sum(squeeze(sleepss_tac19T(10,10,:)))))
  ylabel('Tac PlusMinus500 and Alone')
else
end
cfg=[];
tmp={tac19T_s1{squeeze(sleepss_tac19T(11,10,:))>=thresh}};
if ~isempty(tmp)
  grave_tac19T_s1=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,11);ft_singleplotER(figcfg,grave_tac19T_s1);
  title(num2str(sum(squeeze(sleepss_tac19T(11,10,:)))))
else
end
cfg=[];
tmp={tac19T_s2{squeeze(sleepss_tac19T(12,10,:))>=thresh}};
if ~isempty(tmp)
  grave_tac19T_s2=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,13);ft_singleplotER(figcfg,grave_tac19T_s2);
  title(num2str(sum(squeeze(sleepss_tac19T(12,10,:)))))
else
end
% cfg=[];
% tmp={tac19T_s3{squeeze(sleepss_tac19T(13,10,:))>=thresh}};
% grave_tac19T_s3=ft_timelockgrandaverage(cfg,tmp{:});
% subplot(3,8,15);ft_singleplotER(figcfg,grave_tac19T_s3);

cfg=[];
tmp={tac19T_b0{squeeze(sleepss_tac19T(10,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tac19T_b0=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,10);ft_singleplotER(figcfg,grave_tac19T_b0);
  title(num2str(sum(squeeze(sleepss_tac19T(10,11,:)))))
else
end
cfg=[];
tmp={tac19T_b1{squeeze(sleepss_tac19T(11,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tac19T_b1=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,12);ft_singleplotER(figcfg,grave_tac19T_b1);
  title(num2str(sum(squeeze(sleepss_tac19T(11,11,:)))))
else
end
cfg=[];
tmp={tac19T_b2{squeeze(sleepss_tac19T(12,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tac19T_b2=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,14);ft_singleplotER(figcfg,grave_tac19T_b2);
  title(num2str(sum(squeeze(sleepss_tac19T(12,11,:)))))
else
end
cfg=[];
tmp={tac19T_b3{squeeze(sleepss_tac19T(13,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tac19T_b3=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,16);ft_singleplotER(figcfg,grave_tac19T_b3);
  title(num2str(sum(squeeze(sleepss_tac19T(13,11,:)))))
else
end

cfg=[];
tmp={tacAlone_s0{squeeze(sleepss_tacAlone(10,10,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAlone_s0=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,17);ft_singleplotER(figcfg,grave_tacAlone_s0);
  title(num2str(sum(squeeze(sleepss_tacAlone(10,10,:)))))
  ylabel('Tac Alone')
else
end
cfg=[];
tmp={tacAlone_s1{squeeze(sleepss_tacAlone(11,10,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAlone_s1=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,19);ft_singleplotER(figcfg,grave_tacAlone_s1);
  title(num2str(sum(squeeze(sleepss_tacAlone(11,10,:)))))
else
end
cfg=[];
tmp={tacAlone_s2{squeeze(sleepss_tacAlone(12,10,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAlone_s2=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,21);ft_singleplotER(figcfg,grave_tacAlone_s2);
  title(num2str(sum(squeeze(sleepss_tacAlone(12,10,:)))))
else
end
% cfg=[];
% tmp={tacAlone_s3{squeeze(sleepss_tacAlone(13,10,:))>=thresh}};
% grave_tacAlone_s3=ft_timelockgrandaverage(cfg,tmp{:});
% subplot(3,8,23);ft_singleplotER(figcfg,grave_tacAlone_s3);

cfg=[];
tmp={tacAlone_b0{squeeze(sleepss_tacAlone(10,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAlone_b0=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,18);ft_singleplotER(figcfg,grave_tacAlone_b0);
  title(num2str(sum(squeeze(sleepss_tacAlone(10,11,:)))))
else
end
cfg=[];
tmp={tacAlone_b1{squeeze(sleepss_tacAlone(11,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAlone_b1=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,20);ft_singleplotER(figcfg,grave_tacAlone_b1);
  title(num2str(sum(squeeze(sleepss_tacAlone(11,11,:)))))
else
end
cfg=[];
tmp={tacAlone_b2{squeeze(sleepss_tacAlone(12,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAlone_b2=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,22);ft_singleplotER(figcfg,grave_tacAlone_b2);
  title(num2str(sum(squeeze(sleepss_tacAlone(12,11,:)))))
else
end
cfg=[];
tmp={tacAlone_b3{squeeze(sleepss_tacAlone(13,11,:))>=thresh}};
if ~isempty(tmp)
  grave_tacAlone_b3=ft_timelockgrandaverage(cfg,tmp{:});
  subplot(3,8,24);ft_singleplotER(figcfg,grave_tacAlone_b3);
  title(num2str(sum(squeeze(sleepss_tacAlone(13,11,:)))))
else
end

%% More testing of tactile-alone response, for sleep=1 only, using tacaloneproc=1 flag
% sorting ERP amplitudes

tt=3;
sleep=1;
iiuse=setdiff(iiBuse,[3:7]);
tacaud=1;
iter=11;
trialkc=-1;  % redo for -1 also

clearvars -except tt sub* *dir ii*use *flag sleep* tacaud iter trialkc

submin=iiuse(1)-1;
subuseind=0;
subuse=iiuse;

fsample=1000;
latuse=[-0.5 1];
time=latuse(1):1/fsample:latuse(end);
coloruse=varycolor(length(subuse));
numtr=10;
for ii=iiuse
  cd([edir sub{ii} ])
  subuseind=subuseind+1;
  for tk=1:3
    try
      load(['tlock_tacaloneERPamp_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(tk-2) '.mat'])
    catch
      tlock_tacT_tlockall(:,:,10:12,tk,subuseind)= nan(63,length(time),3);
      tlock_tac9T_tlockall(:,:,10:12,tk,subuseind)= nan(63,length(time),3);
      tlock_audA_tlockall(:,:,10:12,tk,subuseind)= nan(63,length(time),3);
      tlock_aud1A_tlockall(:,:,10:12,tk,subuseind)= nan(63,length(time),3);
      tlock_nul_tlockall(:,:,10:12,tk,subuseind)= nan(63,length(time),3);
      continue
    end
    for ss=10:12
      clear *tmp
      cfg=[];
      cfg.latency=latuse;
      if ~isempty(tlock_tacT_tlock{tt,ss}) && tlock_tacT_tlock{tt,ss}.dof(17,dsearchn(tlock_tacT_tlock{tt,ss}.time',0))>numtr
        tlock_tacT_tlocktmp=ft_selectdata(cfg,tlock_tacT_tlock{tt,ss});
        numtr_tacT(ss,tk,subuseind)=tlock_tacT_tlock{tt,ss}.dof(17,dsearchn(tlock_tacT_tlock{tt,ss}.time',0));
      elseif ~isempty(tlock_tacT_tlock{tt,ss})
        tlock_tacT_tlocktmp.avg=nan(63,length(time));
        numtr_tacT(ss,tk,subuseind)=tlock_tacT_tlock{tt,ss}.dof(17,dsearchn(tlock_tacT_tlock{tt,ss}.time',0));
      else
        tlock_tacT_tlocktmp.avg=nan(63,length(time));
        numtr_tacT(ss,tk,subuseind)=0;
      end
      if ~isempty(tlock_tac9T_tlock{tt,ss}) && tlock_tac9T_tlock{tt,ss}.dof(17,dsearchn(tlock_tac9T_tlock{tt,ss}.time',0))>numtr
        tlock_tac9T_tlocktmp=ft_selectdata(cfg,tlock_tac9T_tlock{tt,ss});
        numtr_tac9T(ss,tk,subuseind)=tlock_tac9T_tlock{tt,ss}.dof(17,dsearchn(tlock_tac9T_tlock{tt,ss}.time',0));
      elseif ~isempty(tlock_tac9T_tlock{tt,ss})
        tlock_tac9T_tlocktmp.avg=nan(63,length(time));
        numtr_tac9T(ss,tk,subuseind)=tlock_tac9T_tlock{tt,ss}.dof(17,dsearchn(tlock_tac9T_tlock{tt,ss}.time',0));
      else
        tlock_tac9T_tlocktmp.avg=nan(63,length(time));
        numtr_tac9T(ss,tk,subuseind)=0;
      end
      if ~isempty(tlock_audA_tlock{tt,ss}) && tlock_audA_tlock{tt,ss}.dof(17,dsearchn(tlock_audA_tlock{tt,ss}.time',0))>numtr
        tlock_audA_tlocktmp=ft_selectdata(cfg,tlock_audA_tlock{tt,ss});
        numtr_audA(ss,tk,subuseind)=tlock_audA_tlock{tt,ss}.dof(17,dsearchn(tlock_audA_tlock{tt,ss}.time',0));
      elseif ~isempty(tlock_audA_tlock{tt,ss})
        tlock_audA_tlocktmp.avg=nan(63,length(time));
        numtr_audA(ss,tk,subuseind)=tlock_audA_tlock{tt,ss}.dof(17,dsearchn(tlock_audA_tlock{tt,ss}.time',0));
      else
        tlock_audA_tlocktmp.avg=nan(63,length(time));
        numtr_audA(ss,tk,subuseind)=0;
      end
      if ~isempty(tlock_aud1A_tlock{tt,ss}) && tlock_aud1A_tlock{tt,ss}.dof(17,dsearchn(tlock_aud1A_tlock{tt,ss}.time',0))>numtr
        tlock_aud1A_tlocktmp=ft_selectdata(cfg,tlock_aud1A_tlock{tt,ss});
        numtr_aud1A(ss,tk,subuseind)=tlock_aud1A_tlock{tt,ss}.dof(17,dsearchn(tlock_aud1A_tlock{tt,ss}.time',0));
      elseif ~isempty(tlock_aud1A_tlock{tt,ss})
        tlock_aud1A_tlocktmp.avg=nan(63,length(time));
        numtr_aud1A(ss,tk,subuseind)=tlock_aud1A_tlock{tt,ss}.dof(17,dsearchn(tlock_aud1A_tlock{tt,ss}.time',0));
      else
        tlock_aud1A_tlocktmp.avg=nan(63,length(time));
        numtr_aud1A(ss,tk,subuseind)=0;
      end
      if ~isempty(tlock_nul_tlock{tt,ss}) && tlock_nul_tlock{tt,ss}.dof(17,dsearchn(tlock_nul_tlock{tt,ss}.time',0))>numtr
        tlock_nul_tlocktmp=ft_selectdata(cfg,tlock_nul_tlock{tt,ss});
        numtr_nul(ss,tk,subuseind)=tlock_nul_tlock{tt,ss}.dof(17,dsearchn(tlock_nul_tlock{tt,ss}.time',0));
      elseif ~isempty(tlock_nul_tlock{tt,ss})
        tlock_nul_tlocktmp.avg=nan(63,length(time));
        numtr_nul(ss,tk,subuseind)=tlock_nul_tlock{tt,ss}.dof(17,dsearchn(tlock_nul_tlock{tt,ss}.time',0));
      else
        tlock_nul_tlocktmp.avg=nan(63,length(time));
        numtr_nul(ss,tk,subuseind)=0;
      end
      
      tlock_tacT_tlockall(:,:,ss,tk,subuseind)= tlock_tacT_tlocktmp.avg;
      tlock_tac9T_tlockall(:,:,ss,tk,subuseind)= tlock_tac9T_tlocktmp.avg;
      tlock_audA_tlockall(:,:,ss,tk,subuseind)= tlock_audA_tlocktmp.avg;
      tlock_aud1A_tlockall(:,:,ss,tk,subuseind)= tlock_aud1A_tlocktmp.avg;
      tlock_nul_tlockall(:,:,ss,tk,subuseind)= tlock_nul_tlocktmp.avg;
    end % ss
  end % tk
end % ii

% figure;plot(time,squeeze(tlock_tacT_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,:)))
% figure;plot(time,squeeze(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,:)))
% figure;plot(time,squeeze(tlock_audA_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,:)))
% figure;plot(time,squeeze(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,:)))
% figure;plot(time,squeeze(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,:)))
% figure;imagesc(time,1:subuseind,squeeze(tlock_tacT_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,:))');caxis([-5 5])
% figure;imagesc(time,1:subuseind,squeeze(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,:))');caxis([-5 5])
% figure;imagesc(time,1:subuseind,squeeze(tlock_audA_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,:))');caxis([-5 5])
% figure;imagesc(time,1:subuseind,squeeze(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,:))');caxis([-5 5])
% figure;imagesc(time,1:subuseind,squeeze(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,:))');caxis([-5 5])

% Get max/min of peaks from awake data
[maxval,maxind]=max(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,1,:),5));
[minval,minind]=min(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,10,1,:),5));

for ss=10:12
  [sortval(:,ss),sortind(:,ss)]=sort(squeeze(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),maxind,ss,1,:)-tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),minind,ss,1,:)));
end

% wideband P200 versus prestim
for ss=12
  [sortwideval(:,ss),sortwideind(:,ss)]=sort(squeeze(mean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),maxind-30:maxind+30,ss,1,:),2)-mean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),1:500,ss,1,:),2)));
end

% wideband P200 versus wideband N100
for ss=12
  [sortwideP2N1val(:,ss),sortwideP2N1ind(:,ss)]=sort(squeeze(mean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),maxind-30:maxind+30,ss,1,:),2)-mean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),minind-30:minind+30,ss,1,:),2)));
end

% for tk=1:3
%   figure(tk)
%   for ss=10:12
%     subplot(3,3,1+(ss-10)*3);plot(time,squeeze(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:)));axis([-inf inf -5 5])
%     hold on;plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k');
%     subplot(3,3,2+(ss-10)*3);plot(time,squeeze(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:)));axis([-inf inf -5 5])
%     hold on;plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k');
%     subplot(3,3,3+(ss-10)*3);plot(time,squeeze(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:)));axis([-inf inf -5 5])
%     hold on;plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k');
%   end
% end

sortTacN2=subuse(sortind(:,12))';
save([edir 'sortTacN2.mat'],'sort*');

% tophalf=sortind(end-8:end,12);
% bothalf=sortind(1:9,12);
% figure(10);
% for ss=10:12
%   subplot(3,3,1+(ss-10)*3);axis([-inf inf -5 5])
%   hold on;plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,tophalf),5)),'Color',cmapuse(end,:));
%   hold on;plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,bothalf),5)),'Color',cmapuse(1,:));
%   subplot(3,3,2+(ss-10)*3);axis([-inf inf -5 5])
%   hold on;plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,tophalf),5)),'Color',cmapuse(end,:));
%   hold on;plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,bothalf),5)),'Color',cmapuse(1,:));
%   subplot(3,3,3+(ss-10)*3);axis([-inf inf -5 5])
%   hold on;plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,tophalf),5)),'Color',cmapuse(end,:));
%   hold on;plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,bothalf),5)),'Color',cmapuse(1,:));
% end
%
% tophalfwide=sortwideind(end-8:end,12);
% bothalfwide=sortwideind(1:9,12);
% figure(11);
% for ss=10:12
%   subplot(3,3,1+(ss-10)*3);axis([-inf inf -5 5])
%   hold on;plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,tophalfwide),5)),'Color',cmapuse(end,:));
%   hold on;plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,bothalfwide),5)),'Color',cmapuse(1,:));
%   subplot(3,3,2+(ss-10)*3);axis([-inf inf -5 5])
%   hold on;plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,tophalfwide),5)),'Color',cmapuse(end,:));
%   hold on;plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,bothalfwide),5)),'Color',cmapuse(1,:));
%   subplot(3,3,3+(ss-10)*3);axis([-inf inf -5 5])
%   hold on;plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,tophalfwide),5)),'Color',cmapuse(end,:));
%   hold on;plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,bothalfwide),5)),'Color',cmapuse(1,:));
% end
%
% tophalfwideP2N1=sortwideP2N1ind(end-8:end,12);
% bothalfwideP2N1=sortwideP2N1ind(1:9,12);
% figure(12);
% for ss=10:12
%   subplot(3,3,1+(ss-10)*3);axis([-inf inf -5 5])
%   hold on;plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,tophalfwideP2N1),5)),'Color',cmapuse(end,:));
%   hold on;plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,bothalfwideP2N1),5)),'Color',cmapuse(1,:));
%   subplot(3,3,2+(ss-10)*3);axis([-inf inf -5 5])
%   hold on;plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,tophalfwideP2N1),5)),'Color',cmapuse(end,:));
%   hold on;plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,bothalfwideP2N1),5)),'Color',cmapuse(1,:));
%   subplot(3,3,3+(ss-10)*3);axis([-inf inf -5 5])
%   hold on;plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,tophalfwideP2N1),5)),'Color',cmapuse(end,:));
%   hold on;plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,bothalfwideP2N1),5)),'Color',cmapuse(1,:));
% end

tophalfwide=sortwideind(end-8:end,12);
bothalfwide=sortwideind(1:9,12);
figure(10);
for ss=10:12
  subplot(3,3,1+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,tophalfwide),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,bothalfwide),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
  subplot(3,3,2+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,tophalfwide),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,bothalfwide),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
  subplot(3,3,3+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,tophalfwide),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,bothalfwide),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
end

tophalfwideP2N1=sortwideP2N1ind(end-8:end,12);
bothalfwideP2N1=sortwideP2N1ind(1:9,12);
figure(11);
for ss=10:12
  subplot(3,3,1+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,tophalfwideP2N1),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,bothalfwideP2N1),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
  subplot(3,3,2+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,tophalfwideP2N1),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,bothalfwideP2N1),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
  subplot(3,3,3+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,tophalfwideP2N1),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,bothalfwideP2N1),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
end

tophalf=sortind(end-8:end,12);
bothalf=sortind(1:9,12);
figure(12);
for ss=10:12
  subplot(3,3,1+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,tophalf),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,bothalf),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
  subplot(3,3,2+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,tophalf),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,bothalf),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
  subplot(3,3,3+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,tophalf),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,bothalf),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
end




% contrast tophalf bothalf for conditions not used for sorting.
for ss=10:12
  for tk=1:2
    tlock_tac9T_gndavgtop{ss,tk}=[];
    tlock_tac9T_gndavgtop{ss,tk}.label=tlock_nul_tlocktmp.label;
    tlock_tac9T_gndavgtop{ss,tk}.dimord='subj_chan_time';
    tlock_tac9T_gndavgtop{ss,tk}.time=tlock_nul_tlocktmp.time;
    
    tlock_tac9T_gndavgbot{ss,tk}=tlock_tac9T_gndavgtop{ss,tk};
    tlock_aud1A_gndavgtop{ss,tk}=tlock_tac9T_gndavgtop{ss,tk};
    tlock_aud1A_gndavgbot{ss,tk}=tlock_tac9T_gndavgtop{ss,tk};
    tlock_nul_gndavgtop{ss,tk}=tlock_tac9T_gndavgtop{ss,tk};
    tlock_nul_gndavgbot{ss,tk}=tlock_tac9T_gndavgtop{ss,tk};
    
    tlock_tac9T_gndavgtop{ss,tk}.individual=permute(squeeze(tlock_tac9T_tlockall(:,:,10,1,tophalf)),[3 1 2]);
    tlock_tac9T_gndavgbot{ss,tk}.individual=permute(squeeze(tlock_tac9T_tlockall(:,:,10,1,bothalf)),[3 1 2]);
    tlock_aud1A_gndavgtop{ss,tk}.individual=permute(squeeze(tlock_aud1A_tlockall(:,:,10,1,tophalf)),[3 1 2]);
    tlock_aud1A_gndavgbot{ss,tk}.individual=permute(squeeze(tlock_aud1A_tlockall(:,:,10,1,bothalf)),[3 1 2]);
    tlock_nul_gndavgtop{ss,tk}.individual  =permute(squeeze(tlock_nul_tlockall(:,:,10,1,tophalf)),[3 1 2]);
    tlock_nul_gndavgbot{ss,tk}.individual  =permute(squeeze(tlock_nul_tlockall(:,:,10,1,bothalf)),[3 1 2]);
  end
end

load eeg1010_neighb
nsub=length(tophalf);

cfg=[];
cfg.latency=[.03 .43];
cfg.neighbours=neighbours;
cfg.parameter='individual';
cfg.method='montecarlo';
cfg.numrandomization=2000;
cfg.correctm='cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.statistic='indepsamplesT';
cfg.design=zeros(2,2*nsub);
cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
% cfg.design(2,:)=[1:nsub 1:nsub];
cfg.ivar=1;
% cfg.uvar=2;
for ss=10:12
  for tk=1:2
    statt_tactopbot{ss,tk}=ft_timelockstatistics(cfg, tlock_tac9T_gndavgtop{ss,tk}, tlock_tac9T_gndavgbot{ss,tk});
    statt_audtopbot{ss,tk}=ft_timelockstatistics(cfg, tlock_aud1A_gndavgtop{ss,tk}, tlock_aud1A_gndavgbot{ss,tk});
    statt_nultopbot{ss,tk}=ft_timelockstatistics(cfg, tlock_nul_gndavgtop{ss,tk},   tlock_nul_gndavgbot{ss,tk});
  end
end
% nothing signif above.

% try for just one channel of main effect.
cfg=[];
cfg.channel='Fz';
cfg.latency=[.03 .43];
cfg.parameter='individual';
cfg.method='montecarlo';
cfg.numrandomization=2000;
cfg.correctm='cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.statistic='indepsamplesT';
cfg.design=zeros(2,2*nsub);
cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
cfg.ivar=1;
for ss=10:12
  for tk=1:2
    statt_tactopbotF{ss,tk}=ft_timelockstatistics(cfg, tlock_tac9T_gndavgtop{ss,tk}, tlock_tac9T_gndavgbot{ss,tk});
    statt_audtopbotF{ss,tk}=ft_timelockstatistics(cfg, tlock_aud1A_gndavgtop{ss,tk}, tlock_aud1A_gndavgbot{ss,tk});
    statt_nultopbotF{ss,tk}=ft_timelockstatistics(cfg, tlock_nul_gndavgtop{ss,tk},   tlock_nul_gndavgbot{ss,tk});
  end
end
% things only trend signif at p=0.07 here

% try for just one channel of main effect and around 100ms and 200ms peaks
cfg=[];
cfg.channel='Fz';
cfg.latency=[.08 .23];
cfg.parameter='individual';
cfg.method='montecarlo';
cfg.numrandomization=2000;
cfg.correctm='cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.statistic='indepsamplesT';
cfg.design=zeros(2,2*nsub);
cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
cfg.ivar=1;
for ss=10:12
  for tk=1:2
    statt_tactopbotF12{ss,tk}=ft_timelockstatistics(cfg, tlock_tac9T_gndavgtop{ss,tk}, tlock_tac9T_gndavgbot{ss,tk});
    statt_audtopbotF12{ss,tk}=ft_timelockstatistics(cfg, tlock_aud1A_gndavgtop{ss,tk}, tlock_aud1A_gndavgbot{ss,tk});
    statt_nultopbotF12{ss,tk}=ft_timelockstatistics(cfg, tlock_nul_gndavgtop{ss,tk},   tlock_nul_gndavgbot{ss,tk});
  end
end
% tac is now signif, aud trend, but nul highly signif!  ???



% %%%%%%%%%%%%%%%%% this code was run with 'sitting' data loaded on one machine and 'lying' data loaded on other machine.

% Uta suggested sorting based on P200 tactile awake data peak versus prestim
tt=3;
sleep=0;
iiuse=intersect(setdiff(iiBuse,[3:7]),iiSuse);
tacaud=1;
iter=27;
trialkc=-1;  % redo for -1 also

submin=iiuse(1)-1;
subuseind=0;
subuse=iiuse;

clearvars -except tt sub* edir ddir ii*use *flag sleep* tacaud iter trialkc

coloruse=varycolor(length(subuse));
cmap=colormap('parula');
sortplot=[iiuse_bothalfWwideP2N1  iiuse_tophalfWwideP2N1]
cmapuse=cmap(1:4:4*length(sortplot),:);
% numtr=10;
for ii=iiuse
  cd([edir sub{ii} ])
  subuseind=subuseind+1;
  try
    load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'])
  catch
    load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '.mat'])
  end
  cfg=[];
  cfg.latency=[-.5 0.5];
  tlock_tacall{subuseind}=ft_selectdata(cfg,tlock_tTacAlone{5,3,10});
  tlock_tacavg(:,:,subuseind)=tlock_tacall{subuseind}.avg;
end % ii
time=tlock_tacall{subuseind}.time;

% Get max/min of peaks from awake data
[maxval,maxind]=max(nanmean(tlock_tacavg(match_str(tlock_tacall{subuseind}.label,'Fz'),:,:),3));
[minval,minind]=min(nanmean(tlock_tacavg(match_str(tlock_tacall{subuseind}.label,'Fz'),:,:),3));

[sortwideWval,sortwideWind]=sort(squeeze(mean(tlock_tacavg(match_str(tlock_tacall{subuseind}.label,'Fz'),maxind-30:maxind+30,:),2)-mean(tlock_tacavg(match_str(tlock_tacall{subuseind}.label,'Fz'),1:500,:),2)));

[sortP2N1Wval,sortP2N1Wind]=sort(squeeze(tlock_tacavg(match_str(tlock_tacall{subuseind}.label,'Fz'),maxind,:)-tlock_tacavg(match_str(tlock_tacall{subuseind}.label,'Fz'),minind,:)));

[sortwideP2N1Wval,sortwideP2N1Wind]=sort(squeeze(mean(tlock_tacavg(match_str(tlock_tacall{subuseind}.label,'Fz'),maxind-30:maxind+30,:),2)-mean(tlock_tacavg(match_str(tlock_tacall{subuseind}.label,'Fz'),minind-30:minind+30,:),2)));

tophalfWwide=sortwideWind(end-7:end);
bothalfWwide=sortwideWind(1:8);

tophalfWP2N1=sortP2N1Wind(end-7:end);
bothalfWP2N1=sortP2N1Wind(1:8);

tophalfWwideP2N1=sortwideP2N1Wind(end-7:end);
bothalfWwideP2N1=sortwideP2N1Wind(1:8);

iiuse_tophalfWwide=iiuse(tophalfWwide);
iiuse_bothalfWwide=iiuse(bothalfWwide);
iiuse_tophalfWP2N1=iiuse(tophalfWP2N1);
iiuse_bothalfWP2N1=iiuse(bothalfWP2N1);
iiuse_tophalfWwideP2N1=iiuse(tophalfWwideP2N1);
iiuse_bothalfWwideP2N1=iiuse(bothalfWwideP2N1);
save([edir 'sortTacWP200.mat'],'iiuse_*');



load([edir 'sortTacWP200.mat']);
iiBfinal=iiBuse(5:end);

figure(20);
for ss=10:12
  subplot(3,3,1+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWwide')),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWwide')),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
  subplot(3,3,2+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWwide')),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWwide')),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
  subplot(3,3,3+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWwide')),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWwide')),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
end
% Conclusion:  P200 in W sitting data partly explains P200 in W lying data, but not N2 lying data


% FINAL in BACN 2015 poster.
tk=1;
figure(42);
for ss=10:12
  subplot(3,3,1+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_tacT_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWwideP2N1')),5),1)),'Color',cmapuse(end,:));set(lh(1),'linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_tacT_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWwideP2N1')),5),1)),'Color',cmapuse(1,:));set(lh(1),'linewidth',3);
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Tactile alone');end
  if ss==10,ylabel('Awake'),elseif ss==11,ylabel('N1');elseif ss==12,ylabel('N2');end
  set(gca,'FontSize',20)
  subplot(3,3,2+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_audA_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWwideP2N1')),5),1)),'Color',cmapuse(end,:));set(lh(1),'linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_audA_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWwideP2N1')),5),1)),'Color',cmapuse(1,:));set(lh(1),'linewidth',3);
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Auditory alone');end
  set(gca,'FontSize',20)
  subplot(3,3,3+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWwideP2N1')),5),1)),'Color',cmapuse(end,:));set(lh(1),'linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWwideP2N1')),5),1)),'Color',cmapuse(1,:));set(lh(1),'linewidth',3);
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Null');end
  if ss==12,lg=legend({'Top half (Part 2, P200-N100)' 'Bottom half (Part 2, P200-N100)'});set(lg,'fontsize',14);end
  set(gca,'FontSize',20)
end
print(42,[fdir 'unisensory_mediansplitWideP2N1_sleep1.png'],'-dpng');
print(42,[fdir 'unisensory_mediansplitWideP2N1_sleep1.eps'],'-depsc');


tk=1;
figure(22);
for ss=10:12
  subplot(3,3,1+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_tacT_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWP2N1')),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_tacT_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWP2N1')),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
  subplot(3,3,2+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_audA_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWP2N1')),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_audA_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWP2N1')),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
  subplot(3,3,3+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWP2N1')),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWP2N1')),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
end

tk=1;
figure(23);
for ss=10:12
  subplot(3,3,1+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_tacT_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWP2N1')),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_tacT_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWP2N1')),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
  subplot(3,3,2+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_audA_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWP2N1')),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_audA_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWP2N1')),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
  subplot(3,3,3+(ss-10)*3);axis([-inf inf -5 5])
  hold on;lh=plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_tophalfWP2N1')),5)),'Color',cmapuse(end,:));set(lh(2),'linestyle',':','linewidth',3);
  hold on;lh=plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,{'Fz' 'Cz'}),:,ss,tk,dsearchn(iiBfinal',iiuse_bothalfWP2N1')),5)),'Color',cmapuse(1,:));set(lh(2),'linestyle',':','linewidth',3);
end

% choose which 'iiuse_*' to use
% decision: iiuse_tophalfWP2N1 iiuse_bothalfWP2N1

% test if difference can be explained in trial numbers
for ss=10:12
  [hh,pp]=ttest2(squeeze(numtr_tac9T(ss,1,dsearchn(iiuse',iiuse_bothalfWwideP2N1'))),squeeze(numtr_tac9T(ss,1,dsearchn(iiuse',iiuse_tophalfWwideP2N1'))))
end

sortplot=[iiuse_bothalfWwideP2N1  iiuse_tophalfWwideP2N1]
cmap=colormap('parula');
cmapuse=cmap(1:4:4*length(sortplot),:);

% this also in final BACN poster
tk=1;figure(1)
for ss=10:12
  subplot(3,3,1+(ss-10)*3);
  for ii=1:length(sortplot),plot(time,squeeze(tlock_tacT_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,dsearchn(subuse',sortplot(ii)))),'Color',cmapuse(ii,:));hold on;end;axis([-inf inf -7 7])
  hold on;
  plot(time,squeeze(nanmean(tlock_tacT_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k','LineWidth',3');
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Tactile alone');end
  if ss==10,ylabel('Awake'),elseif ss==11,ylabel('N1');elseif ss==12,ylabel('N2');end
  set(gca,'FontSize',20)
  subplot(3,3,2+(ss-10)*3);
  for ii=1:length(sortplot),plot(time,squeeze(tlock_audA_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,dsearchn(subuse',sortplot(ii)))),'Color',cmapuse(ii,:));hold on;end;axis([-inf inf -7 7])
  hold on;
  plot(time,squeeze(nanmean(tlock_audA_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k','LineWidth',3');
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Auditory alone');end
  set(gca,'FontSize',20)
  subplot(3,3,3+(ss-10)*3);
  for ii=1:length(sortplot),plot(time,squeeze(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,dsearchn(subuse',sortplot(ii)))),'Color',cmapuse(ii,:));hold on;end;axis([-inf inf -7 7])
  hold on;
  plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k','LineWidth',3');
  if ss==12,xlabel('Time (s)');end
  if ss==10,title('Null');end
  set(gca,'FontSize',20)
end
print(1,[fdir 'unisensory_sortedWideP2N1_sleep1.png'],'-dpng');
print(1,[fdir 'unisensory_sortedWideP2N1_sleep1.eps'],'-depsc');

% tk=1;figure(2)
% for ss=10:12
%   subplot(3,3,1+(ss-10)*3);
%   for ii=1:length(subuse),plot(time,squeeze(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,sortwideind(ii,12))),'Color',cmapuse(ii,:));hold on;end;axis([-inf inf -7 7])
%   hold on;
%   plot(time,squeeze(nanmean(tlock_tac9T_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k','LineWidth',3');
%   if ss==12,xlabel('Time (s)');end
%   if ss==10,title('Tactile alone');end
%   if ss==10,ylabel('Awake'),elseif ss==11,ylabel('N1');elseif ss==12,ylabel('N2');end
%   subplot(3,3,2+(ss-10)*3);
%   for ii=1:length(subuse),plot(time,squeeze(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,sortwideind(ii,12))),'Color',cmapuse(ii,:));hold on;end;axis([-inf inf -7 7])
%   hold on;
%   plot(time,squeeze(nanmean(tlock_aud1A_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k','LineWidth',3');
%   if ss==12,xlabel('Time (s)');end
%   if ss==10,title('Auditory alone');end
%   subplot(3,3,3+(ss-10)*3);
%   for ii=1:length(subuse),plot(time,squeeze(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,sortwideind(ii,12))),'Color',cmapuse(ii,:));hold on;end;axis([-inf inf -7 7])
%   hold on;
%   plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k','LineWidth',3');
%   if ss==12,xlabel('Time (s)');end
%   if ss==10,title('Null');end
% end

% contrast tophalf bothalf for conditions not used for sorting.

for ss=10:12
  for tk=1:2
    tlock_tac9T_gndavgtopWwideP2N1{ss,tk}=[];
    tlock_tac9T_gndavgtopWwideP2N1{ss,tk}.label=tlock_nul_tlocktmp.label;
    tlock_tac9T_gndavgtopWwideP2N1{ss,tk}.dimord='subj_chan_time';
    tlock_tac9T_gndavgtopWwideP2N1{ss,tk}.time=tlock_nul_tlocktmp.time;
    
    tlock_tac9T_gndavgbotWwideP2N1{ss,tk}=tlock_tac9T_gndavgtopWwideP2N1{ss,tk};
    tlock_aud1A_gndavgtopWwideP2N1{ss,tk}=tlock_tac9T_gndavgtopWwideP2N1{ss,tk};
    tlock_aud1A_gndavgbotWwideP2N1{ss,tk}=tlock_tac9T_gndavgtopWwideP2N1{ss,tk};
    tlock_nul_gndavgtopWwideP2N1{ss,tk}=tlock_tac9T_gndavgtopWwideP2N1{ss,tk};
    tlock_nul_gndavgbotWwideP2N1{ss,tk}=tlock_tac9T_gndavgtopWwideP2N1{ss,tk};
    
    tlock_tac9T_gndavgtopWwideP2N1{ss,tk}.individual=permute(squeeze(tlock_tac9T_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_tophalfWwideP2N1'))),[3 1 2]);
    tlock_tac9T_gndavgbotWwideP2N1{ss,tk}.individual=permute(squeeze(tlock_tac9T_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_bothalfWwideP2N1'))),[3 1 2]);
    tlock_aud1A_gndavgtopWwideP2N1{ss,tk}.individual=permute(squeeze(tlock_aud1A_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_tophalfWwideP2N1'))),[3 1 2]);
    tlock_aud1A_gndavgbotWwideP2N1{ss,tk}.individual=permute(squeeze(tlock_aud1A_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_bothalfWwideP2N1'))),[3 1 2]);
    tlock_nul_gndavgtopWwideP2N1{ss,tk}.individual  =permute(squeeze(tlock_nul_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_tophalfWwideP2N1'))),[3 1 2]);
    tlock_nul_gndavgbotWwideP2N1{ss,tk}.individual  =permute(squeeze(tlock_nul_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_bothalfWwideP2N1'))),[3 1 2]);
  end
end

% % contrast tophalf bothalf for conditions not used for sorting.
% for ss=10:12
%   for tk=1:2
%     tlock_tac9T_gndavgtopWwide{ss,tk}=[];
%     tlock_tac9T_gndavgtopWwide{ss,tk}.label=tlock_nul_tlocktmp.label;
%     tlock_tac9T_gndavgtopWwide{ss,tk}.dimord='subj_chan_time';
%     tlock_tac9T_gndavgtopWwide{ss,tk}.time=tlock_nul_tlocktmp.time;
%
%     tlock_tac9T_gndavgbotWwide{ss,tk}=tlock_tac9T_gndavgtopWwide{ss,tk};
%     tlock_aud1A_gndavgtopWwide{ss,tk}=tlock_tac9T_gndavgtopWwide{ss,tk};
%     tlock_aud1A_gndavgbotWwide{ss,tk}=tlock_tac9T_gndavgtopWwide{ss,tk};
%     tlock_nul_gndavgtopWwide{ss,tk}=tlock_tac9T_gndavgtopWwide{ss,tk};
%     tlock_nul_gndavgbotWwide{ss,tk}=tlock_tac9T_gndavgtopWwide{ss,tk};
%
%     tlock_tac9T_gndavgtopWwide{ss,tk}.individual=permute(squeeze(tlock_tac9T_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_tophalfWwide'))),[3 1 2]);
%     tlock_tac9T_gndavgbotWwide{ss,tk}.individual=permute(squeeze(tlock_tac9T_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_bothalfWwide'))),[3 1 2]);
%     tlock_aud1A_gndavgtopWwide{ss,tk}.individual=permute(squeeze(tlock_aud1A_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_tophalfWwide'))),[3 1 2]);
%     tlock_aud1A_gndavgbotWwide{ss,tk}.individual=permute(squeeze(tlock_aud1A_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_bothalfWwide'))),[3 1 2]);
%     tlock_nul_gndavgtopWwide{ss,tk}.individual  =permute(squeeze(tlock_nul_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_tophalfWwide'))),[3 1 2]);
%     tlock_nul_gndavgbotWwide{ss,tk}.individual  =permute(squeeze(tlock_nul_tlockall(:,:,ss,tk,dsearchn(subuse',iiuse_bothalfWwide'))),[3 1 2]);
%   end
% end

nsub=size(tlock_tac9T_gndavgtopWwideP2N1{ss,tk}.individual,1);
% try for just one channel of main effect.
cfg=[];
cfg.channel={'Fz' 'C1' 'C2' 'FC1' 'FC2'};
cfg.avgoverchan='yes';
cfg.latency=[.05 .25];
cfg.parameter='individual';
cfg.method='montecarlo';
cfg.numrandomization=2000;
% cfg.correctm='fdr';
cfg.correctm='cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.statistic='indepsamplesT';
cfg.design=zeros(2,2*nsub);
cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
cfg.ivar=1;
% for ss=10:12
%   for tk=1:2
%     statt_tactopbotWP2N1F{ss,tk}=ft_timelockstatistics(cfg, tlock_tac9T_gndavgtopWP2N1{ss,tk}, tlock_tac9T_gndavgbotWP2N1{ss,tk});
%     statt_audtopbotWP2N1F{ss,tk}=ft_timelockstatistics(cfg, tlock_aud1A_gndavgtopWP2N1{ss,tk}, tlock_aud1A_gndavgbotWP2N1{ss,tk});
%     statt_nultopbotWP2N1F{ss,tk}=ft_timelockstatistics(cfg, tlock_nul_gndavgtopWP2N1{ss,tk},   tlock_nul_gndavgbotWP2N1{ss,tk});
%   end
% end

for ss=10:12
  for tk=1:2
    statt_tactopbotWwideP2N1F{ss,tk}=ft_timelockstatistics(cfg, tlock_tac9T_gndavgtopWwideP2N1{ss,tk}, tlock_tac9T_gndavgbotWwideP2N1{ss,tk});
    statt_audtopbotWwideP2N1F{ss,tk}=ft_timelockstatistics(cfg, tlock_aud1A_gndavgtopWwideP2N1{ss,tk}, tlock_aud1A_gndavgbotWwideP2N1{ss,tk});
    statt_nultopbotWwideP2N1F{ss,tk}=ft_timelockstatistics(cfg, tlock_nul_gndavgtopWwideP2N1{ss,tk},   tlock_nul_gndavgbotWwideP2N1{ss,tk});
  end
end

% for ss=10:12
%   for tk=1:2
%     statt_tactopbotWwideF{ss,tk}=ft_timelockstatistics(cfg, tlock_tac9T_gndavgtopWwide{ss,tk}, tlock_tac9T_gndavgbotWwide{ss,tk});
%     statt_audtopbotWwideF{ss,tk}=ft_timelockstatistics(cfg, tlock_aud1A_gndavgtopWwide{ss,tk}, tlock_aud1A_gndavgbotWwide{ss,tk});
%     statt_nultopbotWwideF{ss,tk}=ft_timelockstatistics(cfg, tlock_nul_gndavgtopWwide{ss,tk},   tlock_nul_gndavgbotWwide{ss,tk});
%   end
% end

% examine peaks across participants at times when statistically significant
ss=10;
timeind=dsearchn(tlock_tac9T_gndavgbotWwideP2N1{ss,1}.time',statt_tactopbotWwideP2N1F{ss,1}.time((statt_tactopbotWwideP2N1F{ss,1}.prob<.05))');
chanind=match_str(tlock_tac9T_gndavgbotWwideP2N1{ss,1}.label,cfg.channel);

[mx,mxind]=max(nanmean(nanmean(tlock_aud1A_tlockall(chanind,:,ss,1,:),1),5))
[mn,mnind]=min(nanmean(nanmean(tlock_aud1A_tlockall(chanind,:,ss,1,:),1),5))

sub_aud1A_tactime=squeeze(mean(mean(tlock_aud1A_tlockall(chanind,timeind,ss,1,:),2),1));
sub_aud1A_amintime=squeeze(mean(mean(tlock_aud1A_tlockall(chanind,mnind-length(timeind)/2:mnind+length(timeind)/2,ss,1,:),2),1));
sub_aud1A_amaxtime=squeeze(mean(mean(tlock_aud1A_tlockall(chanind,mxind-length(timeind)/2:mxind+length(timeind)/2,ss,1,:),2),1));

[cc,pp]=corr(sub_aud1A_tactime(~isnan(sub_aud1A_amaxtime)),sub_tac9T(~isnan(sub_aud1A_amaxtime)))
[cc,pp]=corr(sub_aud1A_amintime(~isnan(sub_aud1A_amaxtime)),sub_tac9T(~isnan(sub_aud1A_amaxtime)))
[cc,pp]=corr(sub_aud1A_amaxtime(~isnan(sub_aud1A_amaxtime)),sub_tac9T(~isnan(sub_aud1A_amaxtime)))

sub_aud1A_tactime=squeeze(mean(mean(tlock_aud1A_tlockall(chanind,timeind,ss,1,:),2),1));
sub_aud1A_amintime=squeeze(mean(mean(tlock_aud1A_tlockall(chanind,mnind-length(timeind)/2:mnind+length(timeind)/2,ss,1,:),2),1));
sub_aud1A_amaxtime=squeeze(mean(mean(tlock_aud1A_tlockall(chanind,mxind-length(timeind)/2:mxind+length(timeind)/2,ss,1,:),2),1));

[mx,mxind]=max(nanmean(nanmean(tlock_tac9T_tlockall(chanind,501:1000,ss,1,:),1),5))
[mn,mnind]=min(nanmean(nanmean(tlock_tac9T_tlockall(chanind,501:1000,ss,1,:),1),5))

sub_tac9T_tactime=squeeze(mean(mean(tlock_tac9T_tlockall(chanind,timeind,ss,1,:),2),1));
sub_tac9T_amintime=squeeze(mean(mean(tlock_tac9T_tlockall(chanind,mnind-length(timeind)/2:mnind+length(timeind)/2,ss,1,:),2),1));
sub_tac9T_amaxtime=squeeze(mean(mean(tlock_tac9T_tlockall(chanind,mxind-length(timeind)/2:mxind+length(timeind)/2,ss,1,:),2),1));

[cc,pp]=corr(sub_aud1A_tactime(~isnan(sub_aud1A_amaxtime)),sub_tac9T_tactime(~isnan(sub_aud1A_amaxtime)))
[cc,pp]=corr(sub_aud1A_amintime(~isnan(sub_aud1A_amaxtime)),sub_tac9T_tactime(~isnan(sub_aud1A_amaxtime)))
[cc,pp]=corr(sub_aud1A_amaxtime(~isnan(sub_aud1A_amaxtime)),sub_tac9T_tactime(~isnan(sub_aud1A_amaxtime)))

[cc,pp]=corr(sub_aud1A_tactime(~isnan(sub_aud1A_amaxtime)),sub_tac9T_amintime(~isnan(sub_aud1A_amaxtime)))
[cc,pp]=corr(sub_aud1A_amintime(~isnan(sub_aud1A_amaxtime)),sub_tac9T_amintime(~isnan(sub_aud1A_amaxtime)))
[cc,pp]=corr(sub_aud1A_amaxtime(~isnan(sub_aud1A_amaxtime)),sub_tac9T_amintime(~isnan(sub_aud1A_amaxtime)))

[cc,pp]=corr(sub_aud1A_tactime(~isnan(sub_aud1A_amaxtime)),sub_tac9T_amaxtime(~isnan(sub_aud1A_amaxtime))) % significant
[cc,pp]=corr(sub_aud1A_amintime(~isnan(sub_aud1A_amaxtime)),sub_tac9T_amaxtime(~isnan(sub_aud1A_amaxtime))) % significant
[cc,pp]=corr(sub_aud1A_amaxtime(~isnan(sub_aud1A_amaxtime)),sub_tac9T_amaxtime(~isnan(sub_aud1A_amaxtime)))

figure;plot(sub_aud1A_amintime(~isnan(sub_aud1A_amintime)),sub_tac9T_amaxtime(~isnan(sub_aud1A_amaxtime)),'o')

% compare across stages


%% More testing of tactile-alone response, using tacaloneproc=1 flag, phase relationship

tt=3;
sleep=1;
iiuse=setdiff(iiBuse,[3:7]);
tacaud=1;
iter=11;
trialkc=-1;  % redo for -1 also

submin=iiuse(1)-1;
subuseind=0;
subuse=iiuse;
erpamp_phasefft_aud1A_all=nan(4,13,2,12,length(iiuse));
erpamp_phasefft_audA_all=nan(4,13,2,12,length(iiuse));
erpamp_phasefft_tac9T_all=nan(4,13,2,12,length(iiuse));
erpamp_phasefft_tacT_all=nan(4,13,2,12,length(iiuse));
erpamp_phasefft_nul_all=nan(4,13,2,12,length(iiuse));
erpamp_phasehilb_aud1A_all=nan(4,3,2,12,length(iiuse));
erpamp_phasehilb_audA_all=nan(4,3,2,12,length(iiuse));
erpamp_phasehilb_tac9T_all=nan(4,3,2,12,length(iiuse));
erpamp_phasehilb_tacT_all=nan(4,3,2,12,length(iiuse));
erpamp_phasehilb_nul_all=nan(4,3,2,12,length(iiuse));
erpamp_powrat_audA_all=nan(26,3,2,length(iiuse));
erpamp_powrat_aud1A_all=nan(26,3,2,length(iiuse));
erpamp_powrat_tac9T_all=nan(26,3,2,length(iiuse));
erpamp_powrat_tacT_all=nan(26,3,2,length(iiuse));
erpamp_powrat_nul_all=nan(26,3,2,length(iiuse));
for ii=iiuse
  cd([edir sub{ii} ])
  load(['tlock_tacaloneERPamp_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'])
  subuseind=subuseind+1;
  for ss=10:12
    if ~isempty(erpamp_phasefft_aud1A{tt,ss}), erpamp_phasefft_aud1A_all(:,:,:,ss,subuseind)=erpamp_phasefft_aud1A{tt,ss}; end
    if ~isempty(erpamp_phasefft_audA{tt,ss}),  erpamp_phasefft_audA_all(:,:,:,ss,subuseind) =erpamp_phasefft_audA{tt,ss};  end
    if ~isempty(erpamp_phasefft_tac9T{tt,ss}),  erpamp_phasefft_tac9T_all(:,:,:,ss,subuseind)=erpamp_phasefft_tac9T{tt,ss};   end
    if ~isempty(erpamp_phasefft_tacT{tt,ss}),  erpamp_phasefft_tacT_all(:,:,:,ss,subuseind) =erpamp_phasefft_tacT{tt,ss};  end
    if ~isempty(erpamp_phasefft_nul{tt,ss}),   erpamp_phasefft_nul_all(:,:,:,ss,subuseind)  =erpamp_phasefft_nul{tt,ss};   end
    if ~isempty(erpamp_phasehilb_aud1A{tt,ss}), erpamp_phasehilb_aud1A_all(:,:,:,ss,subuseind)=erpamp_phasehilb_aud1A{tt,ss}; end
    if ~isempty(erpamp_phasehilb_audA{tt,ss}),  erpamp_phasehilb_audA_all(:,:,:,ss,subuseind) =erpamp_phasehilb_audA{tt,ss};  end
    if ~isempty(erpamp_phasehilb_tac9T{tt,ss}), erpamp_phasehilb_tac9T_all(:,:,:,ss,subuseind)=erpamp_phasehilb_tac9T{tt,ss};   end
    if ~isempty(erpamp_phasehilb_tacT{tt,ss}),  erpamp_phasehilb_tacT_all(:,:,:,ss,subuseind) =erpamp_phasehilb_tacT{tt,ss};  end
    if ~isempty(erpamp_phasehilb_nul{tt,ss}),   erpamp_phasehilb_nul_all(:,:,:,ss,subuseind)  =erpamp_phasehilb_nul{tt,ss};   end
  end
  erpamp_powrat_audA_all(:,:,:,subuseind) =erpamp_powrat_audA;
  erpamp_powrat_aud1A_all(:,:,:,subuseind)=erpamp_powrat_aud1A;
  erpamp_powrat_tac9T_all(:,:,:,subuseind)=erpamp_powrat_tac9T;
  erpamp_powrat_tacT_all(:,:,:,subuseind) =erpamp_powrat_tacT;
  erpamp_powrat_nul_all(:,:,:,subuseind)  =erpamp_powrat_nul;
end

for ss=[10:12 24:26]
  ss
  for pk=1:2
    pk
    [hh,pp]=ttest(squeeze(erpamp_powrat_tac9T_all(ss,1,pk,:)),squeeze(erpamp_powrat_tac9T_all(ss,3,pk,:)))
    [hh,pp]=ttest(squeeze(erpamp_powrat_tacT_all(ss,1,pk,:)), squeeze(erpamp_powrat_tacT_all(ss,3,pk,:)))
    [hh,pp]=ttest(squeeze(erpamp_powrat_aud1A_all(ss,1,pk,:)),squeeze(erpamp_powrat_aud1A_all(ss,3,pk,:)))
    [hh,pp]=ttest(squeeze(erpamp_powrat_audA_all(ss,1,pk,:)), squeeze(erpamp_powrat_audA_all(ss,3,pk,:)))
    [hh,pp]=ttest(squeeze(erpamp_powrat_nul_all(ss,1,pk,:)),  squeeze(erpamp_powrat_nul_all(ss,3,pk,:)))
  end
end
% results:
% for trialkc=0
% for ss=24, pk=1, then audA is sign
% for ss=26, pk=2, then tac9T is sign (tacT borders sign):  Does this say anything more than for N2 the P2 goes away?
%
% for trialkc=-1
% for ss=10, pk=1, tac9T
% for ss=24, pk=1, then tacT, and pk=2 tacT and aud1A
% for ss=26, pk=1, then tac9T (tacT borders) and aud1A
%
% for trialkc=-1 (only tac9T and aud1A and nul)
% for ss=10, pk=1, tac9T
% for ss=24, pk=2, aud1A
% for ss=26, pk=1, then tac9T (tacT borders) and aud1A


close all
alphaval=.005;
for ss=[10:12]
  for ff=1:13
    for pk=1:2
      ppfa1(ff,pk,ss)=anova1(squeeze(erpamp_phasefft_aud1A_all(:,ff,pk,ss,:))',[],'off');
      %       ppfa(ff,pk,ss)=anova1(squeeze(erpamp_phasefft_audA_all(:,ff,pk,ss,:))',[],'off');
      ppft9(ff,pk,ss)=anova1(squeeze(erpamp_phasefft_tac9T_all(:,ff,pk,ss,:))',[],'off');
      %       ppft(ff,pk,ss)=anova1(squeeze(erpamp_phasefft_tacT_all(:,ff,pk,ss,:))',[],'off');
      ppfn(ff,pk,ss)=anova1(squeeze(erpamp_phasefft_nul_all(:,ff,pk,ss,:))',[],'off');
      if ff<4
        ppha1(ff,pk,ss)=anova1(squeeze(erpamp_phasehilb_aud1A_all(:,ff,pk,ss,:))',[],'off');
        %         ppha(ff,pk,ss)=anova1(squeeze(erpamp_phasehilb_audA_all(:,ff,pk,ss,:))',[],'off');
        ppht9(ff,pk,ss)=anova1(squeeze(erpamp_phasehilb_tac9T_all(:,ff,pk,ss,:))',[],'off');
        %         ppht(ff,pk,ss)=anova1(squeeze(erpamp_phasehilb_tacT_all(:,ff,pk,ss,:))',[],'off');
        pphn(ff,pk,ss)=anova1(squeeze(erpamp_phasehilb_nul_all(:,ff,pk,ss,:))',[],'off');
      end
    end
    figure(1);subplot(3,13,(ss-10)*13+ff);bar(squeeze(nanmean(erpamp_phasefft_tac9T_all(:,ff,:,ss,:),5)));axis([0 5 -3 3])
    if ppft9(ff,1,ss)<alphaval, hold on; subplot(3,13,(ss-10)*13+ff);plot(2,-2.5,'b*');end
    if ppft9(ff,2,ss)<alphaval, hold on; subplot(3,13,(ss-10)*13+ff);plot(3,-2.5,'r*');end
    figure(2);subplot(3,13,(ss-10)*13+ff);bar(squeeze(nanmean(erpamp_phasefft_aud1A_all(:,ff,:,ss,:),5)));axis([0 5 -3 3])
    if ppfa1(ff,1,ss)<alphaval, hold on; subplot(3,13,(ss-10)*13+ff);plot(2,-2.5,'b*');end
    if ppfa1(ff,2,ss)<alphaval, hold on; subplot(3,13,(ss-10)*13+ff);plot(3,-2.5,'r*');end
    figure(3);subplot(3,13,(ss-10)*13+ff);bar(squeeze(nanmean(erpamp_phasefft_nul_all(:,ff,:,ss,:),5)));axis([0 5 -3 3])
    if ppfn(ff,1,ss)<alphaval, hold on; subplot(3,13,(ss-10)*13+ff);plot(2,-2.5,'b*');end
    if ppfn(ff,2,ss)<alphaval, hold on; subplot(3,13,(ss-10)*13+ff);plot(3,-2.5,'r*');end
    if ff<4
      figure(4);subplot(3,3,(ss-10)*3+ff);bar(squeeze(nanmean(erpamp_phasehilb_tac9T_all(:,ff,:,ss,:),5)));axis([0 5 -3 3])
      if ppht9(ff,1,ss)<alphaval, hold on; subplot(3,3,(ss-10)*3+ff);plot(2,-2.5,'b*');end
      if ppht9(ff,2,ss)<alphaval, hold on; subplot(3,3,(ss-10)*3+ff);plot(3,-2.5,'r*');end
      figure(5);subplot(3,3,(ss-10)*3+ff);bar(squeeze(nanmean(erpamp_phasehilb_aud1A_all(:,ff,:,ss,:),5)));axis([0 5 -3 3])
      if ppha1(ff,1,ss)<alphaval, hold on; subplot(3,3,(ss-10)*3+ff);plot(2,-2.5,'b*');end
      if ppha1(ff,2,ss)<alphaval, hold on; subplot(3,3,(ss-10)*3+ff);plot(3,-2.5,'r*');end
      figure(6);subplot(3,3,(ss-10)*3+ff);bar(squeeze(nanmean(erpamp_phasehilb_nul_all(:,ff,:,ss,:),5)));axis([0 5 -3 3])
      if pphn(ff,1,ss)<alphaval, hold on; subplot(3,3,(ss-10)*3+ff);plot(2,-2.5,'b*');end
      if pphn(ff,2,ss)<alphaval, hold on; subplot(3,3,(ss-10)*3+ff);plot(3,-2.5,'r*');end
    end
  end
end
% interesting!  some N2 delta phase dependence and W alpha phase dependence;
% however, how much of it is above and beyond what is seen in nul trials?

%% Group level: awake and asleep separately, for each asynchrony separately

soadesc={'Aud first by 500ms' '' 'Aud first by 70ms' 'Aud first by 20ms' 'Simultaneous' 'Tac first by 20ms' 'Tac first by 70ms' '' 'Tac first by 500ms'};
plotflag=0;
printflag=0;
statsflag=1;
loadprevstatsflag=0;
tacnulmsaudstatsflag=0;
% tacaudmsnul_earlyflag=1;
savegrindflag=1;
iterflag=1; % multiple iterations of random trial assignment
audtacflag=0;
runagainflag=0; % 2nd run of erp results (resampling of trials into ERP)
usetr=0;
trialkcflag=0; % 0 for ignore trialkc (old results; all trials), or 1 for use trialkc value
tophalfflag=0; % only relevant for sleep=1; tophalf of participants only (see sortTacN2.mat)
synchasynch=0;
soalist=[1 3 4 5 6 7 9];
statwinorig=0; % =1 means the latency range of .1 to .45 (first old way);  =0 means 0 to .5 (final /better way)

ftver=0;  % 0 means the svn version on the local machine that is already on path
if ftver
  if ispc
    rmpath(genpath('D:\fieldtrip_svn\'))
    addpath(['D:\fieldtrip-' num2str(ftver)]);
  else
    rmpath(genpath('/mnt/hgfs/D/fieldtrip_svn/'))
    addpath(['/mnt/hgfs/D/fieldtrip-' num2str(ftver)]);
  end
end

mcseed=17;

chanuse_sleep0={'all' '-F4'};
chanuse_sleep1={'all' '-AF7' '-AF3' '-Fp1'};

% iiSuse=setdiff(iiSuse, 11);

% allcond_sameN=1; % means using ii>=8, but we doing that anyway with iiuse
% for tacaud=[1 0]
% for sleep=[1]
sleep=0;
if sleep
  chanuse=chanuse_sleep1;
  iteruse=11;
  if trialkcflag
    trialkc=0;  % 0 (no Kc) or 1 (with Kc)
  else
    trialkc=-1;
  end
else
  chanuse=chanuse_sleep0;
%   warning('change me back to 27!!')
  iteruse=27;
  trialkc=-1;
end
iter=iteruse;
tt=3;

% for tt=[3]
close all
figind=1;

for ll=soalist
  %       for ll=[5]
  %   for tt=1:4
  clearvars -except ll tt sub edir ddir ii* sleep *flag figind soa* chanuse* stat* grave*T* grind_*save plv iter* usetr trial* synch* ftver mcseed
  
  %     if ll==1 | ll==9
  %       subuse=8:32;
  %     else
  %       subuse=5:32;
  %     end
  %     if allcond_sameN
  %       subuse=8:32;
  %     end
  if sleep
    if tophalfflag
      load sortTacN2.mat
      subuseall=sort(sortTacN2(end-8:end)');
    else
      subuseall=setdiff(iiBuse,[3:7]);
    end
  else
    subuseall=iiSuse;
  end
  %       subuseall=setdiff(iiuse,[ 27    28    30    31]);
  
  submin=subuseall(1)-1;
  subuseind=0;
  subuse=subuseall; % reset to all for each sleep stage
  
  for ii=subuseall
    %       for ii=[8 9 10 12 14 15  16 17 18]
    cd([edir sub{ii} ])
    %   load(['tlock_diffs_' sub{ii} '.mat']);
    %       load(['tlock_diffs_averef_' sub{ii} '.mat']);
    
    if audtacflag==0
      tacaud=1;
    end
    
    if ~iterflag
      if runagainflag==1
        load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '.mat'])
      else
        load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '.mat'])
        tk0=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(0) '.mat']);
        tk1=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '.mat']);
      end
    else
      %             while ~exist(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '.mat'],'file')
      %               pause(600)
      %             end
      if trialkcflag
        load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'])
        tk1=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter)  '_usetr' num2str(usetr) '_trialkc' num2str(trialkc) '.mat']);
      else % either trialkc=-1 or not exist
        try
          load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'])
        catch
          load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '.mat'])
        end
        try
          tk1=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter)  '_usetr' num2str(usetr) '_trialkc' num2str(trialkc) '.mat']);
        catch
          try
            tk1=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter)  '_usetr' num2str(usetr) '.mat']);
          catch
            tk1=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter) '.mat']);
          end
        end
      end
    end
    
    
    if 1
      % THis is preliminary...how best to included all stages later on?
      if sleep==0
        ss=10; % awake
        sleepcond='Awake';
      elseif sleep==1
        % % % Which Stage??  % %
        ss=12; % N2
        ss=11; % N1
        ss=10; % W
        sleepcond='Sleep W';
        
        %               ss=11; % N1
        %               sleepcond='Sleep N1';
        %             ssuse=23; % this is concatenation of N2 and N3
        %             sleepcond='Sleep (N2+N3)';
      end
    else
      if sleep
        ss=tk1.tr.stageuse;
      else
        ss=tk0.tr.stageuse;
      end
    end
    
    [tk1.tr.stageuse ii]
    
    %                 for ss=ssuse
    if isempty(tk1.tr.stageuse)
      subuse=setdiff(subuse,ii);
      continue
    else
      numtrt(ll,tt,ss,ii-submin)=numt_trials(ll,tt,ss); % does this need to be per 'sleep01' as well?
      if audtacflag
        numtra(ll,tt,ss,ii-submin)=numa_trials(ll,tt,ss);
      end
      
      if numtrt(ll,tt,ss,ii-submin)<20 % what is best number to use here?
        subuse=setdiff(subuse,ii);
        %         tlock_tacPaud{ii-submin}=[];
        %         tlock_audPtac{ii-submin}=[];
        %         tlock_tacMSpN{ii-submin}=[];
        %         tlock_audMSpN{ii-submin}=[];
      else
        subuseind=subuseind+1;
        %         tlock_tacPaud{subuseind}=tlock_tacPaud_s0{ll,tt,ss};
        %         tlock_audPtac{subuseind}=tlock_audPtac_s0{ll,tt,ss};
        %         tlock_tacMSpN{subuseind}=tlock_tacMSpN_s0{ll,tt,ss};
        %         tlock_audMSpN{subuseind}=tlock_audMSpN_s0{ll,tt,ss};
        
        %           if length(ssuse)==1
        %           for ll=soalist
        %             tlock_tacPaud_each{ll,subuseind}=tlock_tacPaud{ll,tt,ss};
        %             tlock_audPtac_each{ll,subuseind}=tlock_audPtac{ll,tt,ss};
        %             tlock_tacMSpN_each{ll,subuseind}=tlock_tacMSpN{ll,tt,ss};
        %             tlock_audMSpN_each{ll,subuseind}=tlock_audMSpN{ll,tt,ss};
        %           end
        tlock_tacPaud_each{subuseind}=tlock_tacPaud{ll,tt,ss};
        tlock_tacMSpN_each{subuseind}=tlock_tacMSpN{ll,tt,ss};
        tlock_MStlock_each{subuseind}=tlock_tMSAlone{ll,tt,ss}; % ll is tac plus auditory multisensory, tac-locked
        tlock_tactlock_each{subuseind}=tlock_tTacAlone{ll,tt,ss}; % ll+20 is tac alone, tac-locked
        tlock_audtlock_each{subuseind}=tlock_tAudAlone{ll,tt,ss}; % ll+20 is aud alone, aud-locked
        tlock_nulttlock_each{subuseind}=tlock_tNulAlone{ll,tt,ss}; % ll+50 is nul, tac-locked
        if audtacflag
          tlock_nulatlock_each{subuseind}=tlock_aNulAlone{ll,tt,ss}; % ll+60 is nul, aud-locked
          tlock_audPtac_each{subuseind}=tlock_audPtac{ll,tt,ss};
          tlock_audMSpN_each{subuseind}=tlock_audMSpN{ll,tt,ss};
        end
        
        cfg=[];
        cfg.operation='subtract';
        cfg.parameter='avg';
        tlock_tacVSnul_each{subuseind}=ft_math(cfg,tlock_tTacAlone{ll,tt,ss},tlock_tNulAlone{ll,tt,ss});
        tlock_audVSnul_each{subuseind}=ft_math(cfg,tlock_tAudAlone{ll,tt,ss},tlock_tNulAlone{ll,tt,ss});
        tlock_msVSnul_each{subuseind}=ft_math(cfg,tlock_tMSAlone{ll,tt,ss},tlock_tNulAlone{ll,tt,ss});
        if audtacflag
          tlock_tacVSnul_each{subuseind}=ft_math(cfg,tlock_aTacAlone{ll,tt,ss},tlock_aNulAlone{ll,tt,ss});
          tlock_audVSnul_each{subuseind}=ft_math(cfg,tlock_aAudAlone{ll,tt,ss},tlock_aNulAlone{ll,tt,ss});
          tlock_msVSnul_each{subuseind}=ft_math(cfg,tlock_aMSAlone{ll,tt,ss},tlock_aNulAlone{ll,tt,ss});
        end
        %           elseif length(ssuse)==2
        %             cfg=[];
        %             cfg.operation='(x1+x2)/s';
        %             cfg.parameter='avg';
        %             cfg.scalar=2;
        %             tlock_tacPaud_each{subuseind}=ft_math(cfg,tlock_tacPaud{ll,tt,ssuse(1)},tlock_tacPaud{ll,tt,ssuse(2)})
        %             tlock_audPtac_each{subuseind}=ft_math(cfg,tlock_audPtac{ll,tt,ssuse(1)},tlock_audPtac{ll,tt,ssuse(2)});
        %             tlock_tacMSpN_each{subuseind}=ft_math(cfg,tlock_tacMSpN{ll,tt,ssuse(1)},tlock_tacMSpN{ll,tt,ssuse(2)});
        %             tlock_audMSpN_each{subuseind}=ft_math(cfg,tlock_audMSpN{ll,tt,ssuse(1)},tlock_audMSpN{ll,tt,ssuse(2)});
        %           end
        
        % this makes up for an error in main code.  once error fixed,
        % then this shouldn't be necessary but leaving it in for now.
        if synchasynch && ll<5
          if ~isfield(tlock_tMSsynch{ll,tt,ss},'avg')
            cfg=[];
            cfg.operation='add';
            cfg.parameter='avg';
            tlock_tMSsynch{ll,tt,ss}=ft_math(cfg,tlock_tacMSshift3tlock{ll,tt,ss},tlock_tacMSshift4tlock{ll,tt,ss});
            cfg.parameter='cov';
            tmp=ft_math(cfg,tlock_tacMSshift3tlock{ll,tt,ss},tlock_tacMSshift4tlock{ll,tt,ss});
            tlock_tMSsynch{ll,tt,ss}.cov=tmp.cov;
          end
          if audtacflag
            if ~isfield(tlock_aMSsynch{ll,tt,ss},'avg')
              cfg=[];
              cfg.operation='add';
              cfg.parameter='avg';
              tlock_aMSsynch{ll,tt,ss}=ft_math(cfg,tlock_audMSshift3tlock{ll,tt,ss},tlock_audMSshift4tlock{ll,tt,ss});
              cfg.parameter='cov';
              tmp=ft_math(cfg,tlock_audMSshift3tlock{ll,tt,ss},tlock_audMSshift4tlock{ll,tt,ss});
              tlock_aMSsynch{ll,tt,ss}.cov=tmp.cov;
            end
          end
          
          % save out temporal-shift-comparison
          tlock_tMSsynch_each{subuseind}=tlock_tMSsynch{ll,tt,ss};
          tlock_tMSasynch_each{subuseind}=tlock_tMSasynch{ll,tt,ss};
          if audtacflag
            tlock_aMSsynch_each{subuseind}=tlock_aMSsynch{ll,tt,ss};
            tlock_aMSasynch_each{subuseind}=tlock_aMSasynch{ll,tt,ss};
          end
        end
        
      end
    end
    %       clear *_s0
    clear tlock*N tlock*tac tlock*aud
  end % ii
  subuseindfinal=subuseind
  %     end
  
  
  
  
  %       if plotflag
  if 0
    %         figure(20);
    %         %     plotind=1;
    %         %     for ii=1:length(subuse)
    %         for ii=1:subuseindfinal
    %           cfg=[];
    %           cfg.parameter='avg';
    %           cfg.operation='subtract';
    %           diff{ii}=ft_math(cfg,tlock_tacPaud_each{ii},tlock_tacMSpN_each{ii});
    %           subplot(3,subuseindfinal,ii);imagesc(tlock_tacPaud_each{ii}.time,1:63,tlock_tacPaud_each{ii}.avg);caxis([-6 6])
    %           subplot(3,subuseindfinal,subuseindfinal+ii);imagesc(tlock_tacMSpN_each{ii}.time,1:63,tlock_tacMSpN_each{ii}.avg);caxis([-6 6])
    %           subplot(3,subuseindfinal,2*subuseindfinal+ii);imagesc(diff{ii}.time,1:63,diff{ii}.avg);caxis([-6 6])
    %           %       plotind=plotind+1;
    %         end
    %         %       if allcond_sameN
    %         if printflag
    %           print(20,['D:\audtac\figs\sumuni_mspn_diff_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
    %         end
    %         %       else
    %         %         print(20,['D:\audtac\figs\sumuni_mspn_diff_ta_cond' num2str(ll) num2str(tt) 'a.png'],'-dpng')
    %         %       end
    %         if audtacflag
    %           figure(21);
    %           %     plotind=1;
    %           for ii=1:subuseindfinal
    %             cfg=[];
    %             cfg.parameter='avg';
    %             cfg.operation='subtract';
    %             diff{ii}=ft_math(cfg,tlock_audPtac_each{ii},tlock_audMSpN_each{ii});
    %             subplot(3,subuseindfinal,ii);imagesc(tlock_audPtac_each{ii}.time,1:63,tlock_audPtac_each{ii}.avg);caxis([-6 6])
    %             subplot(3,subuseindfinal,subuseindfinal+ii);imagesc(tlock_audMSpN_each{ii}.time,1:63,tlock_audMSpN_each{ii}.avg);caxis([-6 6])
    %             subplot(3,subuseindfinal,2*subuseindfinal+ii);imagesc(diff{ii}.time,1:63,diff{ii}.avg);caxis([-6 6])
    %             %       plotind=plotind+1;
    %           end
    %           %       if allcond_sameN
    %           if printflag
    %             print(21,['D:\audtac\figs\sumuni_mspn_diff_at_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
    %           end
    %           %       else
    %           %         print(21,['D:\audtac\figs\sumuni_mspn_diff_at_cond' num2str(ll) num2str(tt) 'a.png'],'-dpng')
    %           %       end
    %         end
  end
  
  cfg=[];
  cfg.keepindividual='yes';
  cfg.channel=chanuse;
  grind_tacPaud=ft_timelockgrandaverage(cfg,tlock_tacPaud_each{:});
  grind_tacMSpN=ft_timelockgrandaverage(cfg,tlock_tacMSpN_each{:});
  grind_tactlock_save{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_tactlock_each{:});
  grind_audtlock_save{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_audtlock_each{:});
  grind_nultlock_save{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_nulttlock_each{:});
  grind_MStlock_save{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_MStlock_each{:});
  if synchasynch && ll<5
    grind_tMSsynch=ft_timelockgrandaverage(cfg,tlock_tMSsynch_each{:});
    grind_tMSasynch=ft_timelockgrandaverage(cfg,tlock_tMSasynch_each{:});
  end
  if audtacflag
    grind_audPtac=ft_timelockgrandaverage(cfg,tlock_audPtac_each{:});
    grind_audMSpN=ft_timelockgrandaverage(cfg,tlock_audMSpN_each{:});
    if ll<5
      grind_aMSsynch=ft_timelockgrandaverage(cfg,tlock_aMSsynch_each{:});
      grind_aMSasynch=ft_timelockgrandaverage(cfg,tlock_aMSasynch_each{:});
    end
  end
  
  
  cfg=[];
  cfg.channel=chanuse;
  grave_tacPaud=ft_timelockgrandaverage(cfg,tlock_tacPaud_each{:});
  grave_tacMSpN=ft_timelockgrandaverage(cfg,tlock_tacMSpN_each{:});
  grave_tactlock=ft_timelockgrandaverage(cfg,tlock_tactlock_each{:});
  grave_audtlock=ft_timelockgrandaverage(cfg,tlock_audtlock_each{:});
  grave_nultlock=ft_timelockgrandaverage(cfg,tlock_nulttlock_each{:});
  grave_MStlock=ft_timelockgrandaverage(cfg,tlock_MStlock_each{:});
  if synchasynch && ll<5
    grave_tMSsynch=ft_timelockgrandaverage(cfg,tlock_tMSsynch_each{:});
    grave_tMSasynch=ft_timelockgrandaverage(cfg,tlock_tMSasynch_each{:});
  end
  if audtacflag
    grave_audPtac=ft_timelockgrandaverage(cfg,tlock_audPtac_each{:});
    grave_audMSpN=ft_timelockgrandaverage(cfg,tlock_audMSpN_each{:});
    if synchasynch && ll<5
      grave_aMSsynch=ft_timelockgrandaverage(cfg,tlock_aMSsynch_each{:});
      grave_aMSasynch=ft_timelockgrandaverage(cfg,tlock_aMSasynch_each{:});
    end
  end
  grave_tacVSnul=ft_timelockgrandaverage(cfg,tlock_tacVSnul_each{:});
  grave_audVSnul=ft_timelockgrandaverage(cfg,tlock_audVSnul_each{:});
  grave_msVSnul=ft_timelockgrandaverage(cfg,tlock_msVSnul_each{:});
  
  
  % unisensory Figure 1 sleep paper
  if plotflag
    figure(52);
    subplot(3,3,1+(ss-10)*3);axis([-inf inf -5 5])
    chanplothere=match_str(grave_tactlock.label,{'Fz' 'C1' 'C2' 'FC1' 'FC2'})
    timeindhere=dsearchn(grave_tactlock.time',-0.5):dsearchn(grave_tactlock.time',0.5);
    
    hold on;plot(grave_tactlock.time(timeindhere),mean(grave_tactlock.avg(chanplothere,timeindhere),1),'k','LineWidth',3');
    
    for ii=1:length(sortplot),plot(grave_tactlock.time(timeindhere),squeeze(mean(grind_tactlock_save{ll,tt,ss}.individual(ii,chanplothere,timeindhere),2)))
      squeeze(tlock_tacT_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,dsearchn(subuse',sortplot(ii)))),'Color',cmapuse(ii,:);hold on;end;axis([-inf inf -7 7])
    
    
    hold on;
    plot(time,squeeze(nanmean(tlock_tacT_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k','LineWidth',3');
    if ss==12,xlabel('Time (s)');end
    if ss==10,title('Tactile alone');end
    if ss==10,ylabel('Awake'),elseif ss==11,ylabel('N1');elseif ss==12,ylabel('N2');end
    set(gca,'FontSize',20)
    subplot(3,3,2+(ss-10)*3);
    for ii=1:length(sortplot),plot(time,squeeze(tlock_audA_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,dsearchn(subuse',sortplot(ii)))),'Color',cmapuse(ii,:));hold on;end;axis([-inf inf -7 7])
    hold on;
    plot(time,squeeze(nanmean(tlock_audA_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k','LineWidth',3');
    if ss==12,xlabel('Time (s)');end
    if ss==10,title('Auditory alone');end
    set(gca,'FontSize',20)
    subplot(3,3,3+(ss-10)*3);
    for ii=1:length(sortplot),plot(time,squeeze(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,dsearchn(subuse',sortplot(ii)))),'Color',cmapuse(ii,:));hold on;end;axis([-inf inf -7 7])
    hold on;
    plot(time,squeeze(nanmean(tlock_nul_tlockall(match_str(tlock_nul_tlocktmp.label,'Fz'),:,ss,tk,:),5)),'k','LineWidth',3');
    if ss==12,xlabel('Time (s)');end
    if ss==10,title('Null');end
    set(gca,'FontSize',20)
    print(1,[fdir 'unisensory_sortedWideP2N1_sleep1.png'],'-dpng');
    print(1,[fdir 'unisensory_sortedWideP2N1_sleep1.eps'],'-depsc');
  end
  
  
  if plotflag
    if sleep
      topoplot_highlight(111,grave_tacPaud,[.07 .77],[]);
      topoplot_highlight(113,grave_tacMSpN,[.07 .77],[]);
      if printflag
        print(111,['D:\audtac\figs\gravetacPaud_topoOverSTime_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
        print(113,['D:\audtac\figs\gravetacMSpN_topoOverSTime_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      end
    else
      topoplot_highlight(111,grave_tacPaud,[.07 .37],[]);
      topoplot_highlight(113,grave_tacMSpN,[.07 .37],[]);
      if printflag
        print(111,['D:\audtac\figs\gravetacPaud_topoOverSTime_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
        print(113,['D:\audtac\figs\gravetacMSpN_topoOverSTime_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      end
    end
  end
  
  
  if 0
    %         if plotflag
    %           %       if 0
    %           topoplot_highlight(11,grave_tacPaud,[-0.5 0.8],[]);
    %           topoplot_highlight(13,grave_tacMSpN,[-0.5 0.8],[]);
    %           if printflag
    %             %       if allcond_sameN
    %             print(11,['D:\audtac\figs\gravetacPaud_topoOverTime_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
    %             %       else
    %             %         print(11,['D:\audtac\figs\gravetacPaud_topoOverTime_cond' num2str(ll) num2str(tt)  '_28.png'],'-dpng')
    %             %       end
    %             %       if allcond_sameN
    %             print(13,['D:\audtac\figs\gravetacMSpN_topoOverTime_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
    %             %       else
    %             %         print(13,['D:\audtac\figs\gravetacMSpN_topoOverTime_cond' num2str(ll) num2str(tt)  '_28.png'],'-dpng')
    %             %       end
    %           end
    %
    %
    %           if audtacflag
    %             topoplot_highlight(12,grave_audPtac,[-0.1 1.1],[]);
    %             topoplot_highlight(14,grave_audMSpN,[-0.1 1.1],[]);
    %             if printflag
    %               %       if allcond_sameN
    %               print(12,['D:\audtac\figs\graveaudPtac_topoOverTime_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
    %               %       else
    %               %         print(12,['D:\audtac\figs\graveaudPtac_topoOverTime_cond' num2str(ll) num2str(tt)  '_28.png'],'-dpng')
    %               %       end
    %               %       if allcond_sameN
    %               print(14,['D:\audtac\figs\graveaudMSpN_topoOverTime_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
    %               %       else
    %               %         print(14,['D:\audtac\figs\graveaudMSpN_topoOverTime_cond' num2str(ll) num2str(tt)  '_28.png'],'-dpng')
    %               %       end
    %             end
    %           end
    %         end
  end
  
  
  % butteryfly plots (as in MEG-UK poster)
  if 0
    %         if sleep
    %           figure(140);
    %         else
    %           figure(130);
    %         end
    %         plot(grave_tacPaud.time,grave_tacPaud.avg,'r'); hold on;
    %         plot(grave_tacMSpN.time,grave_tacMSpN.avg,'b'); hold on;
    %         axis([-0.5 1.1 -8 8])
    %         if audtacflag
    %           if sleep
    %             figure(141);
    %           else
    %             figure(131);
    %           end
    %           plot(grave_audPtac.time,grave_audPtac.avg,'r'); hold on;
    %           plot(grave_audMSpN.time,grave_audMSpN.avg,'b'); hold on;
    %           axis([-0.5 1.1 -8 8])
    %         end
    %
    %         if sleep
    %           figure(150);
    %         else
    %           figure(151);
    %         end
    %         plot(grave_tactlock.time,grave_tactlock.avg,'r'); hold on;
    %         plot(grave_audtlock.time,grave_audtlock.avg,'g'); hold on;
    %         plot(grave_nultlock.time,grave_nultlock.avg,'k'); hold on;
    %         plot(grave_MStlock.time,grave_MStlock.avg,'b'); hold on;
    %         axis([-0.5 1.1 -8 8])
    %
    %
    %         if printflag
    %           if sleep
    %             print(140,['D:\audtac\figs\gravetacMSpN_erpbutterfly_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
    %             if audtacflag
    %               print(141,['D:\audtac\figs\graveaudMSpN_erpbutterfly_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
    %             end
    %           else
    %             print(130,['D:\audtac\figs\gravetacMSpN_erpbutterfly_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
    %             if audtacflag
    %               print(131,['D:\audtac\figs\graveaudMSpN_erpbutterfly_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
    %             end
    %           end
    %         end
  end
  
  % singleplotER averaged over channel groups
  %       if 0
  chanplot{1}={'Fz' 'Cz' 'F1' 'F2' 'FC1' 'FC2' 'C1' 'C2'}; % frontocentral
  chanplot{2}={'CP5' 'P5' 'PO3' 'POz' 'PO4' 'P6' 'CP6' 'Pz' 'P3' 'P4'}; % occipital
  chanplot{3}={'FC6' 'FC4' 'F4' 'F6' 'F8' 'FT8' 'T8' 'C6' 'C4'}; % right frontotemporal
  chanlabel{1}='Frontocentral electrodes';
  chanlabel{2}='Occipital-parietal electrodes';
  chanlabel{3}='Right frontotemporal electrodes';
  if ~tacnulmsaudstatsflag && plotflag
    for cc=1:length(chanplot)
      if cc==2 % posterior, compute phase of alpha at time where signif findings
        %           keyboard
        grindnames=whos('grind*');
        cfg=[];
        cfg.hilbert='complex';
        cfg.bpfilter='yes';
        cfg.bpfreq=[8 13];
        for gg=1:length(grindnames)
          if length(findstr(grindnames(gg).name,'_'))==1
            hilgrind.(grindnames(gg).name)=ft_preprocessing(cfg,eval(grindnames(gg).name));
            tuse=dsearchn(hilgrind.(grindnames(gg).name).time',0):dsearchn(hilgrind.(grindnames(gg).name).time',.7);
            ccuse=match_str(hilgrind.(grindnames(gg).name).label,chanplot{cc});
            plv{ll,tt,ss}.(grindnames(gg).name)=squeeze(mean(hilgrind.(grindnames(gg).name).trial(:,ccuse,tuse)./abs(hilgrind.(grindnames(gg).name).trial(:,ccuse,tuse)),2));
          end
        end
        plv{ll,tt,ss}.time=0:.001:.7;
        if plotflag
          if sleep
            figure(270);
          else
            figure(280);
          end
          subplot(2,7,figind);
          plot(plv{ll,tt,ss}.time,[abs(mean(plv{ll,tt,ss}.grind_tacPaud,1)); abs(mean(plv{ll,tt,ss}.grind_tacMSpN,1))]);
          axis([0 .7 -inf inf]);
          legend({'T Plus A' 'TA Plus Nul'})
          subplot(2,7,figind+7);
          plot(plv{ll,tt,ss}.time,[angle(mean(plv{ll,tt,ss}.grind_tacPaud,1)); angle(mean(plv{ll,tt,ss}.grind_tacMSpN,1))])
          axis([0 .7 -pi pi]);
          legend({'T Plus A' 'TA Plus Nul'})
        end
      end
      
      %             if plotflag
      if sleep
        figure(50+cc*10); % 60, 70, 80
      else
        figure(20+cc*10); % 30, 40, 50
      end
      subplot(1,7,figind);
      cfg=[];
      cfg.channel=chanplot{cc}
      cfg.xlim=[-0.5 1.1];
      cfg.ylim=[-7 7];
      ft_singleplotER(cfg, grave_tacMSpN,grave_tacPaud)
      legend({'TA Plus Nul' 'T Plus A'})
      title([soadesc{ll}])
      xlabel(['Tactile at time 0, ' sleepcond])
      ylabel(chanlabel{cc})
      
      if sleep
        figure(80+cc*10); % 90, 100, 110
      else
        figure(110+cc*10); % 120, 130, 140
      end
      subplot(1,7,figind);
      cfg=[];
      cfg.channel=chanplot{cc}
      cfg.xlim=[-0.5 1.1];
      cfg.ylim=[-7 7];
      ft_singleplotER(cfg, grave_MStlock,grave_tactlock, grave_audtlock, grave_nultlock)
      legend({'TA' 'T' 'A' 'N'})
      title([soadesc{ll}])
      xlabel(['Tactile at time 0, ' sleepcond])
      ylabel(chanlabel{cc})
      
      if synchasynch && ll<5
        if sleep
          figure(140+cc*10); % 150, 160, 170
        else
          figure(170+cc*10); % 180, 190, 200
        end
        subplot(1,7,figind);
        cfg=[];
        cfg.channel=chanplot{cc}
        cfg.xlim=[-0.5 1.1];
        if cc==1
          cfg.ylim=[-10 10];
        else
          cfg.ylim=[-7 7];
        end
        ft_singleplotER(cfg, grave_tMSsynch,grave_tMSasynch)
        legend({'MS synch + shifted' 'MS asynch + shifted'})
        title([soadesc{ll}])
        xlabel(['Tactile at time 0, ' sleepcond])
        ylabel(chanlabel{cc})
      end
      
      
      if audtacflag
        if sleep
          figure(41);
        else
          figure(31);
        end
        subplot(1,7,figind);
        cfg=[];
        cfg.channel=chanplot{cc}
        cfg.xlim=[-0.5 1.1];
        cfg.ylim=[-7 7];
        ft_singleplotER(cfg, grave_audMSpN,grave_audPtac)
        legend({'TA Plus Nul' 'T Plus A'})
        title([soadesc{ll}])
        xlabel(['Auditory at time 0, ' sleepcond])
        ylabel('Frontocentral electrodes')
        
        if synchasynch && ll<5
          if sleep
            figure(200+cc*10); % 210, 220, 230
          else
            figure(230+cc*10); % 240, 250, 260
          end
          subplot(1,7,figind);
          cfg=[];
          cfg.channel=chanplot{cc}
          cfg.xlim=[-0.5 1.1];
          cfg.ylim=[-10 10];
          ft_singleplotER(cfg, grave_tMSsynch,grave_tMSasynch)
          legend({'MS synch + shifted' 'MS asynch + shifted'})
          title([soadesc{ll}])
          xlabel(['Tactile at time 0, ' sleepcond])
          ylabel(chanlabel{cc})
        end
      end
      
      
      
      
      
      %             end % plotflag
    end % cc
  else
    plv=[];
  end
  
  if plotflag && printflag
    if sleep
      print(60,['D:\audtac\figs\gravetacMSpN_FC_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      print(70,['D:\audtac\figs\gravetacMSpN_OP_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      print(80,['D:\audtac\figs\gravetacMSpN_RT_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      print(90,['D:\audtac\figs\graveUniSens_FC_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      print(100,['D:\audtac\figs\graveUniSens_OP_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      print(110,['D:\audtac\figs\graveUniSens_RT_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      print(150,['D:\audtac\figs\gravetacMSsynch_FC_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      print(160,['D:\audtac\figs\gravetacMSsynch_OP_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      print(170,['D:\audtac\figs\gravetacMSsynch_RT_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      if audtacflag
        print(41,['D:\audtac\figs\graveaudMSpN_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
        print(210,['D:\audtac\figs\graveaudMSsynch_FC_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
        print(220,['D:\audtac\figs\graveaudMSsynch_OP_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
        print(230,['D:\audtac\figs\graveaudMSsynch_RT_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      end
    else
      print(30,['D:\audtac\figs\gravetacMSpN_FC_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      print(40,['D:\audtac\figs\gravetacMSpN_OP_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      print(50,['D:\audtac\figs\gravetacMSpN_RT_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      print(120,['D:\audtac\figs\graveUniSens_FC_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      print(130,['D:\audtac\figs\graveUniSens_OP_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      print(140,['D:\audtac\figs\graveUniSens_RT_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      print(180,['D:\audtac\figs\gravetacMSsynch_FC_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      print(190,['D:\audtac\figs\gravetacMSsynch_OP_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      print(200,['D:\audtac\figs\gravetacMSsynch_RT_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      if audtacflag
        print(31,['D:\audtac\figs\graveaudMSpN_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
        print(240,['D:\audtac\figs\graveaudMSsynch_FC_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
        print(250,['D:\audtac\figs\graveaudMSsynch_OP_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
        print(260,['D:\audtac\figs\graveaudMSsynch_RT_erp_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      end
    end
  end
  
  figind=figind+1;
  %       end
  
  
  for ii=1:subuseindfinal
    cfg=[];
    cfg.operation='subtract';
    cfg.parameter='avg';
    cfg.channel=chanuse;
    tlock_TPA_MSPN{ii}=ft_math(cfg,tlock_tacPaud_each{ii},tlock_tacMSpN_each{ii});
    if synchasynch && ll<5
      tlock_TMSs_TMSa{ii}=ft_math(cfg,tlock_tMSsynch_each{ii},tlock_tMSasynch_each{ii});
    end
    if audtacflag
      tlock_APT_MSPN{ii}=ft_math(cfg,tlock_audPtac_each{ii},tlock_audMSpN_each{ii});
      if ll<5
        tlock_AMSs_AMSa{ii}=ft_math(cfg,tlock_aMSsynch_each{ii},tlock_aMSasynch_each{ii});
      end
    end
  end
  cfg=[];
  cfg.channel=chanuse;
  grave_TPA_MSPN{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_TPA_MSPN{:});
  if synchasynch && ll<5
    grave_TMSs_TMSa{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_TMSs_TMSa{:});
  end
  if audtacflag
    grave_APT_MSPN{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_APT_MSPN{:});
    if ll<5
      grave_AMSs_AMSa{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_AMSs_AMSa{:});
    end
  end
  
  
  % figure;
  % % define parameters for plotting
  % timestep = 0.025;      %(in seconds)
  % fsample=1/(tlock_tacMSpN{1}.time(2)-tlock_tacMSpN{1}.time(1));
  % sample_count = length(grave_TPA_MSPN.time);
  % j = [-0.1:timestep:0.5];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
  % m = [1:timestep*fsample:sample_count];  % temporal endpoints in MEEG samples
  % % plot
  % for k = 1:24;
  %      subplot(4,6,k);
  %      cfg = [];
  %      cfg.xlim=[j(k) j(k+1)];
  %      cfg.zlim = [-1.5 1.5];
  % %      pos_int = all(pos(:, m(k):m(k+1)), 2);
  % %      cfg.highlight = 'on';
  % %      cfg.highlightchannel = find(pos_int);
  %      cfg.comment = 'xlim';
  %      cfg.commentpos = 'title';
  %      cfg.layout = 'elec1010.lay';
  %      ft_topoplotER(cfg, grave_TPA_MSPN);
  % end
  %
  % figure;
  % % define parameters for plotting
  % timestep = 0.025;      %(in seconds)
  % fsample=1/(tlock_tacMSpN{1}.time(2)-tlock_tacMSpN{1}.time(1));
  % sample_count = length(grave_TAmMS.time);
  % j = [-0.1:timestep:0.5];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
  % m = [1:timestep*fsample:sample_count];  % temporal endpoints in MEEG samples
  % % plot
  % for k = 1:24;
  %      subplot(4,6,k);
  %      cfg = [];
  %      cfg.xlim=[j(k) j(k+1)];
  %      cfg.zlim = [-1.5 1.5];
  % %      pos_int = all(pos(:, m(k):m(k+1)), 2);
  % %      cfg.highlight = 'on';
  % %      cfg.highlightchannel = find(pos_int);
  %      cfg.comment = 'xlim';
  %      cfg.commentpos = 'title';
  %      cfg.layout = 'elec1010.lay';
  %      ft_topoplotER(cfg, grave_APT_MSPN);
  % end
  %
  
  
  % save out for later
  cfg=[];
  %         cfg.latency=[statt_mc{ll,tt,ss}.time(1) statt_mc{ll,tt,ss}.time(end)];
  cfg.latency=[-.5 1];
  grind_tacPaud_save{ll,tt,ss}=ft_selectdata(cfg,grind_tacPaud);
  grind_tacMSpN_save{ll,tt,ss}=ft_selectdata(cfg,grind_tacMSpN);
  if synchasynch && ll<5
    cfg=[];
    %         cfg.latency=[statt_synch{ll,tt,ss}.time(1) statt_synch{ll,tt,ss}.time(end)];
    cfg.latency=[-.5 1];
    grind_tMSsynch_save{ll,tt,ss}=ft_selectdata(cfg,grind_tMSsynch);
    grind_tMSasynch_save{ll,tt,ss}=ft_selectdata(cfg,grind_tMSasynch);
  end
  
  
  if tacnulmsaudstatsflag
    load eeg1010_neighb
    nsub=length(tlock_tacMSpN_each);
    cfg=[];
    cfg.latency=[.05 .5]+soades(ll);
    cfg.channel=chanuse;
    cfg.neighbours=neighbours;
    cfg.parameter='individual';
    cfg.method='montecarlo';
    cfg.numrandomization=2000;
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
    statt_msaudmc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_MStlock_save{ll,tt,ss},  grind_audtlock_save{ll,tt,ss});
    
    if sleep
      save([edir 'tlock_statmsaud_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '_ftver' num2str(ftver) '.mat'],'statt_msaudmc');
    else
      save([edir 'tlock_statmsaud_sleep' num2str(sleep) '_iter' num2str(iter) '_ftver' num2str(ftver) '.mat'],'statt_msaudmc');
    end
    %
    disp('still add something for tac alone with tac9T')
  end
  
  
  nsub=length(tlock_tacMSpN_each);
  if statsflag && nsub>1
    load eeg1010_neighb
    
    
    
    cfg=[];
    %           if sleep
    %             %           cfg.latency=[.1 .8]; % longer to allow for Kc
    %             if ll==1 || ll==3 || ll==4 || ll==5
    %               cfg.latency=[.1 .8];
    %             elseif ll==6
    %               cfg.latency=[.12 .82];
    %             elseif ll==7
    %               cfg.latency=[.17 .87];
    %             elseif ll==9
    %               cfg.latency=[.6 .13];
    %             end
    %           else
    %           cfg.latency=[.1 .5]; % previously [-.1 .5]
    if statwinorig
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[.1 .45];
      elseif ll==6
        cfg.latency=[.12 .47];
      elseif ll==7
        cfg.latency=[.17 .52];
      elseif ll==9
        cfg.latency=[.6 .95];
      end
    else
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[0 .5];
      elseif ll==6
        cfg.latency=[0 .5]+.02;
      elseif ll==7
        cfg.latency=[0 .5]+.07;
      elseif ll==9
        cfg.latency=[0 .5]+.5;
      end
    end
    
    %           end
    cfg.channel=chanuse;
    cfg.neighbours=neighbours;
    % cfg.parameter='avg';
    cfg.parameter='individual';
    cfg.method='montecarlo';
    % cfg.method='analytic';
    cfg.numrandomization=2000;
    % cfg.correctm='holm';
    cfg.correctm='cluster';
    cfg.clusteralpha = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan = 2;
    cfg.statistic='depsamplesT';
    % cfg.statistic='indepsamplesregrT';
    % cfg.statistic='indepsamplesT';
    cfg.design=zeros(2,2*nsub);
    cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
    cfg.design(2,:)=[1:nsub 1:nsub];
    cfg.ivar=1;
    cfg.uvar=2;
    cfg.randomseed=mcseed;
    disp('test1')
    statt_mc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tacPaud, grind_tacMSpN);
    if audtacflag
      stata_mc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_audPtac, grind_audMSpN);
    end
    %           if tacaudmsnul_earlyflag
    statt_tacmc_early{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tactlock_save{ll,tt,ss}, grind_nultlock_save{ll,tt,ss});
    statt_audmc_early{ll,tt,ss}=ft_timelockstatistics(cfg, grind_audtlock_save{ll,tt,ss}, grind_nultlock_save{ll,tt,ss});
    statt_msmc_early{ll,tt,ss}=ft_timelockstatistics(cfg, grind_MStlock_save{ll,tt,ss},  grind_nultlock_save{ll,tt,ss});
    save([edir 'tlock_statmc_early_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '_statwinorig' num2str(statwinorig) '_ftver' num2str(ftver) 'mcseed' num2str(mcseed) '.mat'],'stat*early');
    %             continue
    %           end
    
    if ~length(statt_mc{ll,tt,ss}.time)
      keyboard
    end
    
    % late component
    cfg=[];
    if audtacflag==0
      if sleep
        %           cfg.latency=[.1 .8]; % longer to allow for Kc
        if ll==1 || ll==3 || ll==4 || ll==5
          cfg.latency=[.1 .8];
        elseif ll==6
          cfg.latency=[.12 .82];
        elseif ll==7
          cfg.latency=[.17 .87];
        elseif ll==9
          cfg.latency=[.6 1.3];
        end
      else
        %           cfg.latency=[.1 .5]; % previously [-.1 .5]
        if ll==1 || ll==3 || ll==4 || ll==5
          cfg.latency=[.45 .9];
        elseif ll==6
          cfg.latency=[.47 .92];
        elseif ll==7
          cfg.latency=[.47 .97];
        elseif ll==9
          cfg.latency=[.95 1.4];
        end
      end
    else
      disp('need to make appropriate latencies here')
      keyboard
    end
    cfg.channel=chanuse;
    cfg.neighbours=neighbours;
    % cfg.parameter='avg';
    cfg.parameter='individual';
    cfg.method='montecarlo';
    % cfg.method='analytic';
    cfg.numrandomization=2000;
    % cfg.correctm='holm';
    cfg.correctm='cluster';
    cfg.clusteralpha = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan = 2;
    cfg.statistic='depsamplesT';
    % cfg.statistic='indepsamplesregrT';
    % cfg.statistic='indepsamplesT';
    cfg.design=zeros(2,2*nsub);
    cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
    cfg.design(2,:)=[1:nsub 1:nsub];
    cfg.ivar=1;
    cfg.uvar=2;
    cfg.randomseed=mcseed;
    disp('test2')
    statt_latemc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tacPaud, grind_tacMSpN);
    if audtacflag
      stata_latemc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_audPtac, grind_audMSpN);
    end
    
    cfg.avgovertime='yes';
    disp('test3')
    statt_late{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tacPaud, grind_tacMSpN);
    if audtacflag
      stata_late{ll,tt,ss}=ft_timelockstatistics(cfg, grind_audPtac, grind_audMSpN);
    end
    
    % early and late all together
    cfg=[];
    if audtacflag==0
      if sleep
        %           cfg.latency=[.1 .8]; % longer to allow for Kc
        if ll==1 || ll==3 || ll==4 || ll==5
          cfg.latency=[.05 1];
        elseif ll==6
          cfg.latency=[.05 1]+.02;
        elseif ll==7
          cfg.latency=[.05 1]+.07;
        elseif ll==9
          cfg.latency=[.05 1]+.5;
        end
      else
        %           cfg.latency=[.1 .5]; % previously [-.1 .5]
        if ll==1 || ll==3 || ll==4 || ll==5
          cfg.latency=[.05 .8];
        elseif ll==6
          cfg.latency=[.07 .82];
        elseif ll==7
          cfg.latency=[.12 .87];
        elseif ll==9
          cfg.latency=[.55 1.3];
        end
      end
    else
      disp('need to make appropriate latencies here')
      keyboard
    end
    cfg.channel=chanuse;
    cfg.neighbours=neighbours;
    % cfg.parameter='avg';
    cfg.parameter='individual';
    cfg.method='montecarlo';
    % cfg.method='analytic';
    cfg.numrandomization=2000;
    % cfg.correctm='holm';
    cfg.correctm='cluster';
    cfg.clusteralpha = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan = 2;
    cfg.statistic='depsamplesT';
    % cfg.statistic='indepsamplesregrT';
    % cfg.statistic='indepsamplesT';
    cfg.design=zeros(2,2*nsub);
    cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
    cfg.design(2,:)=[1:nsub 1:nsub];
    cfg.ivar=1;
    cfg.uvar=2;
    cfg.randomseed=mcseed;
    disp('test4')
    statt_allmc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tacPaud, grind_tacMSpN);
    disp('test5')
    statt_tacmc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tactlock_save{ll,tt,ss}, grind_nultlock_save{ll,tt,ss});
    disp('test6')
    statt_audmc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_audtlock_save{ll,tt,ss}, grind_nultlock_save{ll,tt,ss});
    disp('test7')
    statt_msmc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_MStlock_save{ll,tt,ss},  grind_nultlock_save{ll,tt,ss});
    if audtacflag
      stata_allmc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_audPtac, grind_audMSpN);
      stata_tacmc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tactlock_save{ll,tt,ss}, grind_nultlock_save{ll,tt,ss});
      stata_audmc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_audtlock_save{ll,tt,ss}, grind_nultlock_save{ll,tt,ss});
      stata_msmc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_MStlock_save{ll,tt,ss},  grind_nultlock_save{ll,tt,ss});
    end
    
    % Funny/Uta temporal (in)congruent contrast: (AT70 + TA70) vs (Simult + Simult_shift70)
    if synchasynch && ll<5
      cfg=[];
      if audtacflag==0
        if sleep
          if ll==4
            cfg.latency=[.12 .82];
          elseif ll==3
            cfg.latency=[.17 .87];
          elseif ll==1
            cfg.latency=[.6 1.3];
          end
        else
          if ll==4
            cfg.latency=[.12 .47];
          elseif ll==3
            cfg.latency=[.17 .52];
          elseif ll==1
            cfg.latency=[.6 .95];
          end
        end
      else
        disp('need to make appropriate latencies here')
        keyboard
      end
      cfg.channel=chanuse;
      cfg.neighbours=neighbours;
      % cfg.parameter='avg';
      cfg.parameter='individual';
      cfg.method='montecarlo';
      % cfg.method='analytic';
      cfg.numrandomization=2000;
      % cfg.correctm='holm';
      cfg.correctm='cluster';
      cfg.clusteralpha = 0.05;
      cfg.clusterstatistic = 'maxsum';
      cfg.minnbchan = 2;
      cfg.statistic='depsamplesT';
      % cfg.statistic='indepsamplesregrT';
      % cfg.statistic='indepsamplesT';
      cfg.design=zeros(2,2*nsub);
      cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
      cfg.design(2,:)=[1:nsub 1:nsub];
      cfg.ivar=1;
      cfg.uvar=2;
      disp('test8')
      statt_synch{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tMSsynch, grind_tMSasynch);
      if audtacflag
        stata_synch{ll,tt,ss}=ft_timelockstatistics(cfg, grind_aMSsynch, grind_aMSasynch);
      end
    end
    %   save([edir 'tlock_statmc.mat'],'stat*','grave*'); % no point, as grave* isn't ll,tt,ss dependent
    if iterflag
      if sleep
        %               if trialkcflag
        save([edir 'tlock_statmc_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '_statwinorig' num2str(statwinorig) '_ftver' num2str(ftver) 'mcseed' num2str(mcseed) '.mat'],'stat*');
        %               else
        %                 save([edir 'tlock_statmc_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '.mat'],'stat*');
        %               end
      else
        save([edir 'tlock_statmc_sleep' num2str(sleep) '_iter' num2str(iter) '_statwinorig' num2str(statwinorig) '_ftver' num2str(ftver) 'mcseed' num2str(mcseed) '.mat'],'stat*');
      end
    else
      if runagain==1
        if exist([edir 'tlock_statmc_sleep' num2str(sleep) '_runagain' num2str(runagain) '.mat'],'file')
          if exist([edir 'tlock_statmc_sleep' num2str(sleep) '_runagain' num2str(runagain) '_rerun.mat'],'file')
            save([edir 'tlock_statmc_sleep' num2str(sleep) '_runagain' num2str(runagain) '_rerun2.mat'],'stat*');
          else
            save([edir 'tlock_statmc_sleep' num2str(sleep) '_runagain' num2str(runagain) '_rerun.mat'],'stat*');
          end
        else
          save([edir 'tlock_statmc_sleep' num2str(sleep) '_runagain' num2str(runagain) '.mat'],'stat*');
        end
      else
        if exist([edir 'tlock_statmc_sleep' num2str(sleep) '.mat'],'file')
          if exist([edir 'tlock_statmc_sleep' num2str(sleep) '_rerun.mat'],'file')
            save([edir 'tlock_statmc_sleep' num2str(sleep) '_rerun2.mat'],'stat*');
          else
            save([edir 'tlock_statmc_sleep' num2str(sleep) '_rerun.mat'],'stat*');
          end
        else
          save([edir 'tlock_statmc_sleep' num2str(sleep) '.mat'],'stat*');
        end
      end
    end
    %   if audtacflag
    %     save([edir 'tlock_statmc_sleep' num2str(sleep) '.mat'],'stat*');
    %   end
  end % statsflag
  
  if loadprevstatsflag
    if iterflag
      if sleep
        if trialkcflag
          load([edir 'tlock_statmc_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'statt_*');
        else
          try
            load([edir 'tlock_statmc_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'statt_*');
          catch
            load([edir 'tlock_statmc_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '.mat'],'statt_*');
          end
        end
      else
        load([edir 'tlock_statmc_sleep' num2str(sleep) '_iter' num2str(iter) '.mat'],'statt_*');
      end
    else
      error('please type out other options')
    end
    
    if plotflag
      try close 22;end
      topoplot_highlight(22,grave_TPA_MSPN{ll,tt,ss},[statt_mc{ll,tt,ss}.time(1) statt_mc{ll,tt,ss}.time(end)],statt_mc{ll,tt,ss});
      try close 24;end
      scfg=[];
      scfg.avgovertime='yes';
      scfg.latency=cfg.latency;
      topoplot_highlight(24,ft_selectdata(scfg,grave_TPA_MSPN{ll,tt,ss}),[statt_late{ll,tt,ss}.time statt_late{ll,tt,ss}.time],statt_late{ll,tt,ss});
      try close 26;end
      topoplot_highlight(26,grave_TPA_MSPN{ll,tt,ss},[statt_latemc{ll,tt,ss}.time(1) statt_latemc{ll,tt,ss}.time(end)],statt_latemc{ll,tt,ss});
      try close 28;end
      topoplot_highlight(28,grave_TPA_MSPN{ll,tt,ss},[statt_allmc{ll,tt,ss}.time(1) statt_allmc{ll,tt,ss}.time(end)],statt_allmc{ll,tt,ss});
      try close 10;end
      topoplot_highlight(10,grave_tacVSnul,[statt_tacmc{ll,tt,ss}.time(1) statt_tacmc{ll,tt,ss}.time(end)],statt_tacmc{ll,tt,ss});
      try close 12;end
      topoplot_highlight(12,grave_audVSnul,[statt_audmc{ll,tt,ss}.time(1) statt_audmc{ll,tt,ss}.time(end)],statt_audmc{ll,tt,ss});
      try close 14;end
      topoplot_highlight(14,grave_msVSnul,[statt_msmc{ll,tt,ss}.time(1) statt_msmc{ll,tt,ss}.time(end)],statt_msmc{ll,tt,ss});
      if ll<5
        try close 32;end
        topoplot_highlight(32,grave_TMSs_TMSa{ll,tt,ss},[statt_synch{ll,tt,ss}.time(1) statt_synch{ll,tt,ss}.time(end)],statt_synch{ll,tt,ss});
      end
      
      %     % get relevant (significant) values
      %     pos_cluster_pvals = [statt_mc.posclusters(:).prob];
      %     % pos_signif_clust = find(pos_cluster_pvals < statt.cfg.alpha);
      %     pos_signif_clust = find(pos_cluster_pvals < 0.05);
      %     pos = ismember(statt_mc.posclusterslabelmat, pos_signif_clust);
      %     % define parameters for plotting
      %     try close(22)
      %     catch
      %     end
      %     figure(22);
      %     timestep = 0.025;      %(in seconds)
      %     fsample=1/(tlock_tacMSpN{1}.time(2)-tlock_tacMSpN{1}.time(1));
      %     sample_count = length(grave_TPA_MSPN.time);
      %     j = [-0.1:timestep:0.5];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
      %     m = [1:timestep*fsample:sample_count];  % temporal endpoints in MEEG samples
      %     % plot
      %     for k = 1:24;
      %       subplot(4,6,k);
      %       cfg = [];
      %       cfg.xlim=[j(k) j(k+1)];
      %       cfg.zlim = [-1.5 1.5];
      %       %      pos_int = all(pos(:, m(k):m(k+1)), 2);
      %       pos_int = any(pos(:, m(k):m(k+1)), 2);
      %       cfg.highlight = 'on';
      %       cfg.highlightchannel = find(pos_int)
      %       cfg.comment = 'xlim';
      %       cfg.commentpos = 'title';
      %       cfg.layout = 'elec1010.lay';
      %       ft_topoplotER(cfg, grave_TPA_MSPN);
      %     end
      %       if allcond_sameN
      if printflag
        print(22,['D:\audtac\figs\grdiff_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '.png'],'-dpng')
        print(24,['D:\audtac\figs\grdiff_topoOverTime_late_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '.png'],'-dpng')
        print(26,['D:\audtac\figs\grdiff_topoOverTime_latemc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '.png'],'-dpng')
        print(28,['D:\audtac\figs\grdiff_topoOverTime_allmc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '.png'],'-dpng')
        print(10,['D:\audtac\figs\grdiff_topoOverTime_tacmc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '.png'],'-dpng')
        print(12,['D:\audtac\figs\grdiff_topoOverTime_audmc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '.png'],'-dpng')
        print(14,['D:\audtac\figs\grdiff_topoOverTime_msmc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '.png'],'-dpng')
        if ll<5, print(32,['D:\audtac\figs\grdiff_topoOverTime_synch_ta_cond' num2str(ll) num2str(ss) num2str(tt) num2str(sleep) '.png'],'-dpng');end
      end
      %       else
      %         print(22,['D:\audtac\figs\grdiff_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) 'a.png'],'-dpng')
      %       end
      
      if audtacflag
        topoplot_highlight(23,grave_APT_MSPN{ll,tt,ss},[stata_mc{ll,tt,ss}.time(1) stata_mc{ll,tt,ss}.time(end)],stata_mc{ll,tt,ss});
        scfg=[];
        scfg.avgovertime='yes';
        scfg.latency=cfg.latency;
        topoplot_highlight(25,ft_selectdata(scfg,grave_APT_MSPN{ll,tt,ss}),[stata_late{ll,tt,ss}.time stata_late{ll,tt,ss}.time],stata_late{ll,tt,ss});
        topoplot_highlight(27,grave_APT_MSPN{ll,tt,ss},[stata_latemc{ll,tt,ss}.time(1) stata_latemc{ll,tt,ss}.time(end)],stata_latemc{ll,tt,ss});
        topoplot_highlight(29,grave_APT_MSPN{ll,tt,ss},[stata_allmc{ll,tt,ss}.time(1) stata_allmc{ll,tt,ss}.time(end)],stata_allmc{ll,tt,ss});
        topoplot_highlight(11,grave_tacVSnul,[stata_tacmc{ll,tt,ss}.time(1) stata_tacmc{ll,tt,ss}.time(end)],stata_tacmc{ll,tt,ss});
        topoplot_highlight(13,grave_audVSnul,[stata_audmc{ll,tt,ss}.time(1) stata_audmc{ll,tt,ss}.time(end)],stata_audmc{ll,tt,ss});
        topoplot_highlight(15,grave_msVSnul,[stata_msmc{ll,tt,ss}.time(1) stata_msmc{ll,tt,ss}.time(end)],stata_msmc{ll,tt,ss});
        if ll<5,topoplot_highlight(33,grave_AMSs_AMSa{ll,tt,ss},[stata_synch{ll,tt,ss}.time(1) stata_synch{ll,tt,ss}.time(end)],stata_synch{ll,tt,ss});end
        %     % get relevant (significant) values
        %     pos_cluster_pvals = [stata_mc.posclusters(:).prob];
        %     % pos_signif_clust = find(pos_cluster_pvals < statt.cfg.alpha);
        %     pos_signif_clust = find(pos_cluster_pvals < 0.05);
        %     pos = ismember(stata_mc.posclusterslabelmat, pos_signif_clust);
        %     % define parameters for plotting
        %     try close(23)
        %     catch
        %     end
        %     figure(23);
        %     timestep = 0.025;      %(in seconds)
        %     fsample=1/(tlock_tacMSpN{1}.time(2)-tlock_tacMSpN{1}.time(1));
        %     sample_count = length(grave_APT_MSPN.time);
        %     j = [-0.1:timestep:0.5];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
        %     m = [1:timestep*fsample:sample_count];  % temporal endpoints in MEEG samples
        %     % plot
        %     for k = 1:24;
        %       subplot(4,6,k);
        %       cfg = [];
        %       cfg.xlim=[j(k) j(k+1)];
        %       cfg.zlim = [-1.5 1.5];
        %       %      pos_int = all(pos(:, m(k):m(k+1)), 2);
        %       pos_int = any(pos(:, m(k):m(k+1)), 2);
        %       cfg.highlight = 'on';
        %       cfg.highlightchannel = find(pos_int)
        %       cfg.comment = 'xlim';
        %       cfg.commentpos = 'title';
        %       cfg.layout = 'elec1010.lay';
        %       ft_topoplotER(cfg, grave_APT_MSPN);
        %     end
        %       if allcond_sameN
        if printflag
          print(23,['D:\audtac\figs\grdiff_topoOverTime_mc_at_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '.png'],'-dpng')
          print(25,['D:\audtac\figs\grdiff_topoOverTime_late_at_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '.png'],'-dpng')
          print(27,['D:\audtac\figs\grdiff_topoOverTime_latemc_at_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '.png'],'-dpng')
          print(29,['D:\audtac\figs\grdiff_topoOverTime_allmc_at_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '.png'],'-dpng')
          print(11,['D:\audtac\figs\grdiff_topoOverTime_tacmc_at_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '.png'],'-dpng')
          print(13,['D:\audtac\figs\grdiff_topoOverTime_audmc_at_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '.png'],'-dpng')
          print(15,['D:\audtac\figs\grdiff_topoOverTime_msmc_at_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '.png'],'-dpng')
          if ll<5, print(33,['D:\audtac\figs\grdiff_topoOverTime_synch_at_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '.png'],'-dpng');end
        end
      end
      
    end % if plotflag
    
  end % if loadprevstatsflag
  %       else
  %         print(23,['D:\audtac\figs\grdiff_topoOverTime_mc_at_cond' num2str(ll) num2str(tt) 'a.png'],'-dpng')
  %       end
  
  
  %Ignore Analytic for now
  %     cfg=[];
  %     cfg.latency=[-.1 .5];
  %     % cfg.neighbours=neighbours;
  %     % cfg.parameter='avg';
  %     cfg.parameter='individual';
  %     % cfg.method='montecarlo';
  %     cfg.method='analytic';
  %     cfg.alpha=0.05;
  %     % cfg.numrandomization=200;
  %     cfg.correctm='fdr';
  %     % cfg.correctm='cluster';
  %     % cfg.clusteralpha = 0.05;
  %     % cfg.clusterstatistic = 'maxsum';
  %     % cfg.minnbchan = 2;
  %     cfg.statistic='depsamplesT';
  %     % cfg.statistic='indepsamplesregrT';
  %     % cfg.statistic='indepsamplesT';
  %     cfg.design=zeros(2,2*nsub);
  %     cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
  %     cfg.design(2,:)=[1:nsub 1:nsub];
  %     cfg.ivar=1;
  %     cfg.uvar=2;
  %     statt_an=ft_timelockstatistics(cfg, grind_tacPaud, grind_tacMSpN);
  %     stata_an=ft_timelockstatistics(cfg, grind_audPtac, grind_audMSpN);
  %
  %     % get relevant (significant) values
  %     pos = statt_an.mask;
  %     % define parameters for plotting
  %     figure(24);
  %     timestep = 0.025;      %(in seconds)
  %     fsample=1/(tlock_tacMSpN{1}.time(2)-tlock_tacMSpN{1}.time(1));
  %     sample_count = length(grave_TPA_MSPN.time);
  %     j = [-0.1:timestep:0.5];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
  %     m = [1:timestep*fsample:sample_count];  % temporal endpoints in MEEG samples
  %     % plot
  %     for k = 1:24;
  %       subplot(4,6,k);
  %       cfg = [];
  %       cfg.xlim=[j(k) j(k+1)];
  %       cfg.zlim = [-1.5 1.5];
  %       %      pos_int = all(pos(:, m(k):m(k+1)), 2);
  %       pos_int = any(pos(:, m(k):m(k+1)), 2);
  %       cfg.highlight = 'on';
  %       cfg.highlightchannel = find(pos_int);
  %       cfg.comment = 'xlim';
  %       cfg.commentpos = 'title';
  %       cfg.layout = 'elec1010.lay';
  %       ft_topoplotER(cfg, grave_TPA_MSPN);
  %     end
  %     print(24,['D:\audtac\figs\grdiff_topoOverTime_an_ta_cond' num2str(ll) num2str(tt) 'a.png'],'-dpng')
  %
  %     % get relevant (significant) values
  %     pos = stata_an.mask;
  %     % define parameters for plotting
  %     figure(25);
  %     timestep = 0.025;      %(in seconds)
  %     fsample=1/(tlock_tacMSpN{1}.time(2)-tlock_tacMSpN{1}.time(1));
  %     sample_count = length(grave_APT_MSPN.time);
  %     j = [-0.1:timestep:0.5];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
  %     m = [1:timestep*fsample:sample_count];  % temporal endpoints in MEEG samples
  %     % plot
  %     for k = 1:24;
  %       subplot(4,6,k);
  %       cfg = [];
  %       cfg.xlim=[j(k) j(k+1)];
  %       cfg.zlim = [-1.5 1.5];
  %       %      pos_int = all(pos(:, m(k):m(k+1)), 2);
  %       pos_int = any(pos(:, m(k):m(k+1)), 2);
  %       cfg.highlight = 'on';
  %       cfg.highlightchannel = find(pos_int);
  %       cfg.comment = 'xlim';
  %       cfg.commentpos = 'title';
  %       cfg.layout = 'elec1010.lay';
  %       ft_topoplotER(cfg, grave_APT_MSPN);
  %     end
  %     print(25,['D:\audtac\figs\grdiff_topoOverTime_an_at_cond' num2str(ll) num2str(tt) 'a.png'],'-dpng')
  
  
  if 0 % sometimes useful to pause here
    keyboard
  end
  try
    close(111);
    close(113);
  end
  %       if sleep
  %         try
  %           close(140);
  %           close(141)
  %         end
  %       else
  %         try
  %           close(130);
  %           close(131);
  %         end
  %       end
end % ll
%     keyboard
%     end % tt

if savegrindflag
  if ~iterflag
    save([edir 'tlock_grind_sleep' num2str(sleep) '_runagain' num2str(runagain) '.mat'],'grave*T*','grind_*save','plv');
    break
  else
    if sleep
      save([edir 'tlock_grind_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '_statwinorig' num2str(statwinorig) '_ftver' num2str(ftver) '.mat'],'grave*T*','grind_*save','plv');
    else
      save([edir 'tlock_grind_sleep' num2str(sleep) '_iter' num2str(iter) '_statwinorig' num2str(statwinorig) '_ftver' num2str(ftver) '.mat'],'grave*T*','grind_*save','plv');
    end
  end
end
%   end % iter
%   save([edir 'tlock_grind_sleep' num2str(sleep) '.mat'],'grave*T*','grind_*save','plv');

%   if audtacflag
%     save([edir 'tlock_grind_sleep' num2str(sleep) '.mat'],'grave*T*','grind_*save','plv');
%   end
% end % sleep
% end
% save([edir 'tlock_numtrlltt.mat'],'numtr*','grind*'); % no point as grind* isn't ll,tt,ss dependent
save([edir 'tlock_numtrlltt.mat'],'numtr*');

close all


%% Compare Stages during 'sleep' data (W vs N1 vs N2)
% This had been written and then lost during a file transfer.
% Thus only bits of original code remain or are gussed at as to how it was set up.
% At least this was discovered only a few days after it was written, so still very fresh in my mind!

% To generate
% tlock_statmc_STAGECOMP_sleep1_ss12_iter11_trialkc-1.mat
% tlock_grind_STAGECOMP_sleep1_iter11.mat

plotflag=0;
printflag=0;
statsflag=1;
soalist=[1 3 4 5 6 7 9];
tophalfflag=0;


chanuse_sleep0={'all' '-F4'};
chanuse_sleep1={'all' '-AF7' '-AF3' '-Fp1'};
tt=3;

sleep=1;
chanuse=chanuse_sleep1;
iteruse=11;
iter=iteruse;
trialkc=-1;
usetr=1;

tacaud=1;
if sleep
  if tophalfflag
    load sortTacN2.mat
    subuseall=sort(sortTacN2(end-8:end)');
  else
    subuseall=setdiff(iiBuse,[3:7]);
  end
else
  subuseall=iiSuse;
end



% for ll=soalist
for ll=[7 9]
  clearvars -except ll tt sub edir ddir ii* sleep *flag soa* chanuse* iter* usetr trial* synch* tacaud  subuseall stat* grind*
  submin=subuseall(1)-1;
  subuseind=0;
  %         subuse=subuseall;
  subuse=nan(12,length(subuseall));
  subuse(10:12,:)=repmat(subuseall,[3 1]);
  for ii=subuseall
    subuseind=subuseind+1;
    cd([edir sub{ii} ])
    
    load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'])
    tk1=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter)  '_usetr' num2str(usetr) '_trialkc' num2str(trialkc) '.mat']);
    
    
    for ss=10:12
      
      numtrt(ll,tt,ss,subuseind)=numt_trials(ll,tt,ss); % does this need to be per 'sleep01' as well?
      if numtrt(ll,tt,ss,subuseind)<20
        subuse(ss,subuseind)=nan;
      end
      
      if ~isnan(subuse(ss,subuseind))
        
        tlock_tacPaud_each{ss}{subuseind}=tlock_tacPaud{ll,tt,ss};
        tlock_tacMSpN_each{ss}{subuseind}=tlock_tacMSpN{ll,tt,ss};
        tlock_MStlock_each{ss}{subuseind}=tlock_tMSAlone{ll,tt,ss}; % ll is tac plus auditory multisensory, tac-locked
        tlock_tactlock_each{ss}{subuseind}=tlock_tTacAlone{ll,tt,ss}; % ll+20 is tac alone, tac-locked
        tlock_audtlock_each{ss}{subuseind}=tlock_tAudAlone{ll,tt,ss}; % ll+20 is aud alone, aud-locked
        tlock_nulttlock_each{ss}{subuseind}=tlock_tNulAlone{ll,tt,ss}; % ll+50 is nul, tac-locked
        
        cfg=[];
        cfg.operation='subtract';
        cfg.parameter='avg';
        tlock_tacVSnul_each{ss}{subuseind}=ft_math(cfg,tlock_tTacAlone{ll,tt,ss},tlock_tNulAlone{ll,tt,ss});
        tlock_audVSnul_each{ss}{subuseind}=ft_math(cfg,tlock_tAudAlone{ll,tt,ss},tlock_tNulAlone{ll,tt,ss});
        tlock_msVSnul_each{ss}{subuseind}=ft_math(cfg,tlock_tMSAlone{ll,tt,ss},tlock_tNulAlone{ll,tt,ss});
      end
    end % ss
    
    clear tlock*N tlock*tac tlock*aud
  end % ii
  subuseindfinal=subuseind
  
  
  useAll=~isnan(subuse(10,:)) & ~isnan(subuse(11,:));
  useN1N2=~isnan(subuse(11,:));
  
  for ss=10:12
    
    cfg=[];
    cfg.keepindividual='yes';
    cfg.channel=chanuse;
    %     grind_tacPaud_stageAll{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_tacPaud_each{ss}{useAll});
    %     grind_tacMSpN_stageAll{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_tacMSpN_each{ss}{useAll});
    grind_tactlock_stageAll_save{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_tactlock_each{ss}{useAll});
    grind_audtlock_stageAll_save{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_audtlock_each{ss}{useAll});
    grind_nultlock_stageAll_save{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_nulttlock_each{ss}{useAll});
    grind_MStlock_stageAll_save{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_MStlock_each{ss}{useAll});
    
    if ss>10
      %       grind_tacPaud_stageN1N2{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_tacPaud_each{ss}{useN1N2});
      %       grind_tacMSpN_stageN1N2{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_tacMSpN_each{ss}{useN1N2});
      grind_tactlock_stageN1N2_save{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_tactlock_each{ss}{useN1N2});
      grind_audtlock_stageN1N2_save{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_audtlock_each{ss}{useN1N2});
      grind_nultlock_stageN1N2_save{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_nulttlock_each{ss}{useN1N2});
      grind_MStlock_stageN1N2_save{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_MStlock_each{ss}{useN1N2});
    end
    
    cfg=[];
    cfg.channel=chanuse;
    grave_tacPaud_stageAll{ss}=ft_timelockgrandaverage(cfg,tlock_tacPaud_each{ss}{useAll});
    grave_tacMSpN_stageAll{ss}=ft_timelockgrandaverage(cfg,tlock_tacMSpN_each{ss}{useAll});
    grave_tactlock_stageAll{ss}=ft_timelockgrandaverage(cfg,tlock_tactlock_each{ss}{useAll});
    grave_audtlock_stageAll{ss}=ft_timelockgrandaverage(cfg,tlock_audtlock_each{ss}{useAll});
    grave_nultlock_stageAll{ss}=ft_timelockgrandaverage(cfg,tlock_nulttlock_each{ss}{useAll});
    grave_MStlock_stageAll{ss}=ft_timelockgrandaverage(cfg,tlock_MStlock_each{ss}{useAll});
    
    grave_tacVSnul_stageAll{ss}=ft_timelockgrandaverage(cfg,tlock_tacVSnul_each{ss}{useAll});
    grave_audVSnul_stageAll{ss}=ft_timelockgrandaverage(cfg,tlock_audVSnul_each{ss}{useAll});
    grave_msVSnul_stageAll{ss}=ft_timelockgrandaverage(cfg,tlock_msVSnul_each{ss}{useAll});
    
    if ss>10
      grave_tacPaud_stageN1N2{ss}=ft_timelockgrandaverage(cfg,tlock_tacPaud_each{ss}{useN1N2});
      grave_tacMSpN_stageN1N2{ss}=ft_timelockgrandaverage(cfg,tlock_tacMSpN_each{ss}{useN1N2});
      grave_tactlock_stageN1N2{ss}=ft_timelockgrandaverage(cfg,tlock_tactlock_each{ss}{useN1N2});
      grave_audtlock_stageN1N2{ss}=ft_timelockgrandaverage(cfg,tlock_audtlock_each{ss}{useN1N2});
      grave_nultlock_stageN1N2{ss}=ft_timelockgrandaverage(cfg,tlock_nulttlock_each{ss}{useN1N2});
      grave_MStlock_stageN1N2{ss}=ft_timelockgrandaverage(cfg,tlock_MStlock_each{ss}{useN1N2});
      
      grave_tacVSnul_stageN1N2{ss}=ft_timelockgrandaverage(cfg,tlock_tacVSnul_each{ss}{useN1N2});
      grave_audVSnul_stageN1N2{ss}=ft_timelockgrandaverage(cfg,tlock_audVSnul_each{ss}{useN1N2});
      grave_msVSnul_stageN1N2{ss}=ft_timelockgrandaverage(cfg,tlock_msVSnul_each{ss}{useN1N2});
    end
    
    finduseAll=find(useAll);
    finduseN1N2=find(useN1N2);
    cfg=[];
    cfg.operation='subtract';
    cfg.parameter='avg';
    cfg.channel=chanuse;
    for ii=1:length(find(useAll))
      tlock_TPA_MSPN_stageAll{ss}{ii}=ft_math(cfg,tlock_tacPaud_each{ss}{finduseAll(ii)},tlock_tacMSpN_each{ss}{finduseAll(ii)});
    end
    if ss>10
      for ii=1:length(find(useN1N2))
        tlock_TPA_MSPN_stageN1N2{ss}{ii}=ft_math(cfg,tlock_tacPaud_each{ss}{finduseN1N2(ii)},tlock_tacMSpN_each{ss}{finduseN1N2(ii)});
      end
    end
    cfg=[];
    cfg.channel=chanuse;
    grave_TPA_MSPN_stageAll{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_TPA_MSPN_stageAll{ss}{:});
    if ss>10
      grave_TPA_MSPN_stageN1N2{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_TPA_MSPN_stageN1N2{ss}{:});
    end
    cfg.keepindividual='yes';
    grind_TPA_MSPN_stageAll{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_TPA_MSPN_stageAll{ss}{:});
    if ss>10
      grind_TPA_MSPN_stageN1N2{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_TPA_MSPN_stageN1N2{ss}{:});
    end
    
    %     cfg=[];
    %     cfg.latency=[-.5 1];
    %     grind_tacPaud_save{ll,tt,ss}=ft_selectdata(cfg,grind_tacPaud{ss});
    %     grind_tacMSpN_save{ll,tt,ss}=ft_selectdata(cfg,grind_tacMSpN{ss});
    %     grind_tacPaud_save{ll,tt,ss}=ft_selectdata(cfg,grind_tacPaud{ss});
    %     grind_tacMSpN_save{ll,tt,ss}=ft_selectdata(cfg,grind_tacMSpN{ss});
    
  end %ss
  
  
  if statsflag
    load eeg1010_neighb
    nsuball=length(find(sum(isnan(subuse(10:12,:)))==0));
    nsubN1N2=size(grind_TPA_MSPN_stageN1N2{ll,tt,11}.individual,1);
    
    
    cfg=[];
    if ll==1 || ll==3 || ll==4 || ll==5
      cfg.latency=[.1 .45];
    elseif ll==6
      cfg.latency=[.12 .47];
    elseif ll==7
      cfg.latency=[.17 .52];
    elseif ll==9
      cfg.latency=[.6 .95];
    end
    cfg.channel=chanuse;
    cfg.neighbours=neighbours;
    cfg.parameter='individual';
    cfg.method='montecarlo';
    cfg.numrandomization=2000;
    cfg.correctm='cluster';
    cfg.clusteralpha = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan = 2;
    cfg.statistic='depsamplesT';
    cfg.ivar=1;
    cfg.uvar=2;
    
    cfg.design=zeros(2,2*nsuball);
    cfg.design(1,:)=[ones(1,nsuball) 2*ones(1,nsuball)];
    cfg.design(2,:)=[1:nsuball 1:nsuball];
    
    statt_mc_WN1{ll}=ft_timelockstatistics(cfg,grind_TPA_MSPN_stageAll{ll,tt,11},grind_TPA_MSPN_stageAll{ll,tt,10});
    statt_ms_early_WN1{ll}=ft_timelockstatistics(cfg,grind_MStlock_stageAll_save{ll,tt,11},grind_MStlock_stageAll_save{ll,tt,10});
    
    cfg.design=zeros(2,2*nsubN1N2);
    cfg.design(1,:)=[ones(1,nsubN1N2) 2*ones(1,nsubN1N2)];
    cfg.design(2,:)=[1:nsubN1N2 1:nsubN1N2];
    
    statt_mc_N1N2{ll}=ft_timelockstatistics(cfg,grind_TPA_MSPN_stageN1N2{ll,tt,11},grind_TPA_MSPN_stageN1N2{ll,tt,12});
    statt_ms_early_N1N2{ll}=ft_timelockstatistics(cfg,grind_MStlock_stageN1N2_save{ll,tt,11},grind_MStlock_stageN1N2_save{ll,tt,12});
    
    if ll==1
      cfg.latency=[.1 .45]-.5;
    elseif ll==3
      cfg.latency=[.1 .45]-.07;
    elseif ll==4
      cfg.latency=[.1 .45]-.02;
    elseif ll==5
      cfg.latency=[.1 .45];
    elseif ll==6
      cfg.latency=[.12 .47];
    elseif ll==7
      cfg.latency=[.17 .52];
    elseif ll==9
      cfg.latency=[.6 .95];
    end
    
    statt_aud_early_N1N2{ll}=ft_timelockstatistics(cfg,grind_audtlock_stageN1N2_save{ll,tt,11},grind_audtlock_stageN1N2_save{ll,tt,12});
    
    cfg.design=zeros(2,2*nsuball);
    cfg.design(1,:)=[ones(1,nsuball) 2*ones(1,nsuball)];
    cfg.design(2,:)=[1:nsuball 1:nsuball];
    
    statt_aud_early_WN1{ll}=ft_timelockstatistics(cfg,grind_audtlock_stageAll_save{ll,tt,11},grind_audtlock_stageAll_save{ll,tt,10});
    
    cfg.latency=[.1 .45];
    
    statt_tac_early_WN1{ll}=ft_timelockstatistics(cfg,grind_tactlock_stageAll_save{ll,tt,11},grind_tactlock_stageAll_save{ll,tt,10});
    
    cfg.design=zeros(2,2*nsubN1N2);
    cfg.design(1,:)=[ones(1,nsubN1N2) 2*ones(1,nsubN1N2)];
    cfg.design(2,:)=[1:nsubN1N2 1:nsubN1N2];
    
    statt_tac_early_N1N2{ll}=ft_timelockstatistics(cfg,grind_tactlock_stageN1N2_save{ll,tt,11},grind_tactlock_stageN1N2_save{ll,tt,12});
    
    save([edir 'tlock_statmc_STAGECOMP_sleep' num2str(sleep) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'stat*');
    save([edir 'tlock_grind_STAGECOMP_sleep' num2str(sleep) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'grind*');
    
    cfg=[];
    cfg.neighbours=neighbours;
    cfg.parameter='individual';
    cfg.method='montecarlo';
    cfg.numrandomization=1000;
    cfg.correctm='cluster';
    cfg.clusteralpha = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan = 2;
    cfg.statistic='depsamplesFunivariate';
    cfg.ivar=1;
    cfg.uvar=2;
    cfg.tail=1;
    cfg.design=zeros(2,3*nsuball);
    cfg.design(1,:)=[ones(1,nsuball) 2*ones(1,nsuball) 3*ones(1,nsuball)];
    cfg.design(2,:)=[1:nsuball 1:nsuball 1:nsuball];
    
    if ll==1 || ll==3 || ll==4 || ll==5
      cfg.latency=[.1 .45];
    elseif ll==6
      cfg.latency=[.12 .47];
    elseif ll==7
      cfg.latency=[.17 .52];
    elseif ll==9
      cfg.latency=[.6 .95];
    end
    
    statt_mc_All{ll}=ft_timelockstatistics(cfg,grind_TPA_MSPN_stageAll{ll,tt,10},grind_TPA_MSPN_stageAll{ll,tt,11},grind_TPA_MSPN_stageAll{ll,tt,12});
    statt_ms_early_All{ll}=ft_timelockstatistics(cfg,grind_MStlock_stageAll_save{ll,tt,10},grind_MStlock_stageAll_save{ll,tt,11},grind_MStlock_stageAll_save{ll,tt,12});
    
    if ll==1
      cfg.latency=[.1 .45]-.5;
    elseif ll==3
      cfg.latency=[.1 .45]-.07;
    elseif ll==4
      cfg.latency=[.1 .45]-.02;
    elseif ll==5
      cfg.latency=[.1 .45];
    elseif ll==6
      cfg.latency=[.12 .47];
    elseif ll==7
      cfg.latency=[.17 .52];
    elseif ll==9
      cfg.latency=[.6 .95];
    end
    
    statt_aud_early_All{ll}=ft_timelockstatistics(cfg,grind_audtlock_stageAll_save{ll,tt,10},grind_audtlock_stageAll_save{ll,tt,11},grind_audtlock_stageAll_save{ll,tt,12});
    
    cfg.latency=[.1 .45];
    
    statt_tac_early_All{ll}=ft_timelockstatistics(cfg,grind_tactlock_stageAll_save{ll,tt,10},grind_tactlock_stageAll_save{ll,tt,11},grind_tactlock_stageAll_save{ll,tt,12});
    
    save([edir 'tlock_statmc_STAGECOMP_sleep' num2str(sleep) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'stat*');
  end
  
end  %ll




%% Assess reproducibility of stats output
if runagain==0
  stat1=load('tlock_statmc_sleep0');
  stat2=load('tlock_statmc_sleep0_rerun.mat');
  stat3=load('tlock_statmc_sleep0_rerun2.mat')
elseif runagain==1
  stat1=load('tlock_statmc_sleep0_runagain1');
  stat2=load('tlock_statmc_sleep0_runagain1_rerun.mat');
  stat3=load('tlock_statmc_sleep0_runagain1_rerun2.mat')
end

tt=3;ss=10;soalist=[1 3 4 5 6 7 9];
mcposp=nan(9,3);
mcnegp=nan(9,3);
snposp=nan(4,3);
snnegp=nan(4,3);
for ll=soalist
  try mcposp(ll,1)=stat1.statt_mc{ll,tt,ss}.posclusters(1).prob;end
  try mcnegp(ll,1)=stat1.statt_mc{ll,tt,ss}.negclusters(1).prob;end
  try mcposp(ll,2)=stat2.statt_mc{ll,tt,ss}.posclusters(1).prob;end
  try mcnegp(ll,2)=stat2.statt_mc{ll,tt,ss}.negclusters(1).prob;end
  try mcposp(ll,3)=stat3.statt_mc{ll,tt,ss}.posclusters(1).prob;end
  try mcnegp(ll,3)=stat3.statt_mc{ll,tt,ss}.negclusters(1).prob;end
  if ll<5
    try snposp(ll,1)=stat1.statt_synch{ll,tt,ss}.posclusters(1).prob;end
    try snnegp(ll,1)=stat1.statt_synch{ll,tt,ss}.negclusters(1).prob;end
    try snposp(ll,2)=stat2.statt_synch{ll,tt,ss}.posclusters(1).prob;end
    try snnegp(ll,2)=stat2.statt_synch{ll,tt,ss}.negclusters(1).prob;end
    try snposp(ll,3)=stat3.statt_synch{ll,tt,ss}.posclusters(1).prob;end
    try snnegp(ll,3)=stat3.statt_synch{ll,tt,ss}.negclusters(1).prob;end
  end
end

% regular contrast (TPA - MSPN),  pos and neg clusters
% for runagain==0
% mcposp
%     0.2879    0.2889    0.2969
%     0.0150    0.0150    0.0130    %% AT70          *
%     0.1879    0.2009    0.1879
%     0.0135    0.0205    0.0220    %% Simult
%        NaN       NaN       NaN
%     0.0610    0.0560    0.0550    % TA70 (trend)
%        NaN       NaN       NaN
%
% mcnegp
%        NaN       NaN       NaN
%     0.1304    0.1129    0.1219
%     0.0065    0.0070    0.0055    %% AT20          **
%     0.3088    0.2879    0.2929
%     0.2559    0.2559    0.2549
%     0.2014    0.2084    0.2079
%     0.1879    0.1784    0.1869
%
% for runagain==1
% mcposp
%     0.3068    0.2829    0.3038
%     0.0130    0.0105    0.0110    %% AT70          *
%     0.2949    0.2984    0.2964
%     0.2169    0.2054    0.2294
%     0.2604    0.2674    0.2704
%        NaN       NaN       NaN
%        NaN       NaN       NaN
%
% mcnegp
%        NaN       NaN       NaN
%     0.1174    0.1474    0.1409
%     0.0140    0.0150    0.0070   %% AT20          **
%     0.0600    0.0670    0.0755   % Simult trend
%     0.0055    0.0025    0.0040   %% AT20
%        NaN       NaN       NaN
%        NaN       NaN       NaN

% conclusion:
% 1) stats/pvalues sufficiently stable over reruns on given dataset
% 2) datasets differ sufficiently that it matters.



% MS-synch vs asynch contrast
% for runagain==0
% snposp
%        NaN       NaN       NaN
%     0.0115    0.0130    0.0080   %% AT70 and TA70
%        NaN       NaN       NaN
%
% snnegp
%     0.0120    0.0090    0.0115   %% AT500 and TA500
%     0.0005    0.0005    0.0005   %% AT70 and TA70
%     0.0030    0.0010    0.0020   %% AT20 and TA20

% for runagain==0
% snposp
%        NaN       NaN       NaN
%     0.0250    0.0245    0.0260   %% AT70 and TA70
%        NaN       NaN       NaN
%
% snnegp
%     0.0385    0.0500    0.0530   %% AT500 and TA500
%     0.0010    0.0005    0.0005   %% AT70 and TA70
%     0.0035    0.0055    0.0020   %% AT20 and TA20

mcposinclus(:,:,ll)=and(mcposinclus(:,:,ll),mcposclusmat(:,:,ll,iter));
%% Examining iterations of trial-comparisons
iterflag=1;
plotflag=0;
sleep=1;
tt=3;
ss=12;
[mcposp,mcnegp,snposp,snnegp]=deal(nan(9,20));
[mcposclusmat,mcnegclusmat,snposclusmat,snnegclusmat]= deal(zeros(62,351,9,20));
for iter=11
  %   while ~exist([edir 'tlock_statmc_sleep' num2str(sleep) '_iter' num2str(iter) '.mat'],'file')
  %     pause(600);
  %   end
  if sleep
    load([edir 'tlock_statmc_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '.mat']);
  else
    load([edir 'tlock_statmc_sleep' num2str(sleep) '_iter' num2str(iter) '.mat']);
  end
  for ll=[1 3 4 5 6 7 9]
    try mcposp(ll,iter)=statt_mc{ll,tt,ss}.posclusters(1).prob; catch
      mcposp(ll,iter)=nan; end
    try mcnegp(ll,iter)=statt_mc{ll,tt,ss}.negclusters(1).prob;  catch
      mcnegp(ll,iter)=nan; end
    try snposp(ll,iter)=statt_synch{ll,tt,ss}.posclusters(1).prob;  catch
      snposp(ll,iter)=nan; end
    try snnegp(ll,iter)=statt_synch{ll,tt,ss}.negclusters(1).prob;  catch
      snnegp(ll,iter)=nan; end
    if plotflag
      figure(ll);
      subplot(4,10,iter);try imagesc(statt_mc{ll,tt,ss}.posclusterslabelmat); caxis([0 5]);end
      subplot(4,10,iter+10);try imagesc(statt_mc{ll,tt,ss}.negclusterslabelmat); caxis([0 5]);end
      subplot(4,10,iter+20);try imagesc(statt_synch{ll,tt,ss}.posclusterslabelmat); caxis([0 5]);end
      subplot(4,10,iter+30);try imagesc(statt_synch{ll,tt,ss}.negclusterslabelmat); caxis([0 5]);end
    end
    try mcposclusmat(:,:,ll,iter)=statt_mc{ll,tt,ss}.posclusterslabelmat;end
    try mcnegclusmat(:,:,ll,iter)=statt_mc{ll,tt,ss}.negclusterslabelmat;end
    try snposclusmat(:,:,ll,iter)=statt_synch{ll,tt,ss}.posclusterslabelmat;end
    try snnegclusmat(:,:,ll,iter)=statt_synch{ll,tt,ss}.negclusterslabelmat;end
    
    for ii=1:length(statt_mc{ll,tt,ss}.cfg.previous{1}.previous)
      numtr(ll,tt,ss,ii,iter)=length(statt_mc{ll,tt,ss}.cfg.previous{1}.previous{ii}.previous{1}.previous.previous.trials);
    end
  end
end
save('pos_neg_stats.mat','*pos*','*neg*')

load('pos_neg_stats.mat')
% --> simple inclusive-intersection won't work: not always 'top' cluster
soalist=[1 3 4 5 6 7 9];
for ll=1:7
  figure(1);
  subplot(2,7,ll);  imagesc(mean(mcposclusmat(:,:,soalist(ll),:)==1,4));
  subplot(2,7,ll+7);imagesc(mean(mcnegclusmat(:,:,soalist(ll),:)==1,4));
  figure(2);
  subplot(2,7,ll);  imagesc(mean(mcposclusmat(:,:,soalist(ll),:)==2,4));
  subplot(2,7,ll+7);imagesc(mean(mcnegclusmat(:,:,soalist(ll),:)==2,4));
end


[mcposinclus,mcneginclus]=deal(zeros(62,351,9));
for ll=[1 3 4 5 6 7 9]
  for iter=1:10
    %     if all(all(isnan(mcposclusmat(:,:,ll,iter))))
    %       mcposclusmat(:,:,ll,iter)=zeros(62,351);
    %       mcnegclusmat(:,:,ll,iter)=zeros(62,351);
    %     end
    mcposinclus(:,:,ll)=and(mcposinclus(:,:,ll),mcposclusmat(:,:,ll,iter));
    mcneginclus(:,:,ll)=and(mcneginclus(:,:,ll),mcnegclusmat(:,:,ll,iter));
    
    
    mcnegclusmat(:,:,ll,iter)==1
  end
  mean(logical(mcnegclusmat(:,:,ll,:)),4)==1
end

cfg=[];
% for iter=1:10
%   load([edir 'tlock_grind_sleep' num2str(sleep) '_iter' num2str(iter) '.mat']);
% end

% ****
% % randperm(10)
% %      1     9     6    10     4     2     8     5     3     7
% % % USE iter1  % ************************
% this is most similar to iter27

% mcposp([1 3 4 5 6 7 9],:)
%        NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN
% *   0.0915    0.2769    0.0095    0.0595    0.0070    0.0440       NaN    0.0240    0.1259    0.0230
%     0.2679    0.1629    0.2364    0.0770       NaN    0.0710       NaN    0.3128       NaN       NaN
%     0.2779    0.2429    0.2429    0.0680    0.1669    0.2884    0.2254    0.1944    0.0585    0.0605
%        NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN
%     0.0870    0.2659    0.0530    0.0975    0.1289    0.1004    0.0855    0.0670    0.3358    0.1094
%        NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN

% mcnegp([1 3 4 5 6 7 9],:)
%        NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN
%     0.0885    0.0620    0.1949    0.0750    0.1054    0.1364    0.0550    0.2024    0.0885    0.0625
% *   0.0155    0.0095    0.0020    0.0220    0.0225    0.0355    0.0110    0.0275    0.0130    0.0020
% ?   0.0930    0.0135    0.0550    0.0595    0.0465    0.0560    0.0850    0.0445    0.0770    0.2754
% *   0.0220    0.0775    0.0150    0.0315    0.0880    0.2059    0.0240    0.0165    0.0410    0.1529
%     0.1874       NaN    0.0550    0.1594    0.2949       NaN    0.1569    0.1979       NaN    0.1094
%        NaN    0.3088    0.2669    0.1124       NaN       NaN    0.1579    0.1879    0.0645       NaN

% nanmean(mcposp([1 3 4 5 6 7 9],:),2)
%              p<0.05   p<0.10
%        NaN
%     0.0735   5/10 and 7/10;               median 0.044
%     0.1880
%     0.1826
%        NaN
%     0.1330
%        NaN
% nanmean(mcnegp([1 3 4 5 6 7 9],:),2)
%        NaN
%     0.1070
%     0.0160   10/10 and 10/10              median 0.0142
%     0.0805   3/10 and 9/10                median 0.0577
%     0.0674   6/10 and 8/10                median 0.0362            * inflated by not all same cluster
%     0.1658
%     0.1831


%
% snposp([1 3 4],:)
%        NaN       NaN       NaN       NaN    0.2569       NaN       NaN       NaN       NaN       NaN
%     0.0085    0.0130    0.0100    0.0395    0.0125    0.0100    0.0160    0.0120    0.0240    0.0120
%        NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN
%
% snnegp([1 3 4],:)
%     0.0145    0.0515    0.0205    0.2919    0.0075    0.0100    0.0115    0.0275    0.0140    0.0065
%     0.0005    0.0005    0.0015    0.0005    0.0005    0.0005    0.0005    0.0005    0.0005    0.0005
%     0.0045    0.0015    0.0020    0.0035    0.0015    0.0025    0.0025    0.0035    0.0020    0.0065

% nanmean(snposp([1 3 4],:),2)
%     0.2569
% *   0.0157   10/10                        median 0.0122
%        NaN
%
% nanmean(snnegp([1 3 4],:),2)
% *   0.0455    8/10 and 9/10               median 0.0142
% *   0.0006    10/10                       median 0.0005
% *   0.0030    10/10                       median 0.0025

%% Plotting above ERP sensor results


% First, plotting the final stats/ difference of conditions
% Below, plotting conditions on their own.

clearvars -except sub *dir
% chanplot{1}={'Fz' 'Cz' 'F1' 'F2' 'FC1' 'FC2' 'C1' 'C2'}; % frontocentral
% chanplot{2}={'CP5' 'P5' 'PO3' 'POz' 'PO4' 'P6' 'CP6' 'Pz' 'P3' 'P4'}; % occipital
% chanplot{3}={'FC6' 'FC4' 'F4' 'F6' 'F8' 'FT8' 'T8' 'C6' 'C4'}; % right frontotemporal
% chanlabel{1}='Frontocentral electrodes';
% chanlabel{2}='Occipital-parietal electrodes';
% chanlabel{3}='Right frontotemporal electrodes';
chanplot{1}={'Fz' 'Cz' 'F1' 'F2' 'FC1' 'FC2' 'C1' 'C2'}; % frontocentral
chanplot{2}={'CP5' 'POz' 'Pz' 'P3' 'P4' 'C4' 'O1' 'O2' 'P7' 'PO7'}; % occipital
chanlabel{1}='Frontocentral electrodes';
chanlabel{2}='Occipital-parietal electrodes';
iterflag=1;
sleep=0;
tt=3;
statwinorig=0;  % =1 means the latency range of .1 to .45 (first way);  =0 means 0 to .5 (final /better way)

if sleep
  iter=11;
  ss=12;
  trialkc=0;
else
  iter=27;
  ss=10;
  trialkc=-1;
end

if iterflag
  if sleep
    try
      load([edir 'tlock_statmc_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '_statwinorig' num2str(statwinorig) '.mat']);
      load([edir 'tlock_grind_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '_statwinorig' num2str(statwinorig) '.mat']);
    catch
      load([edir 'tlock_statmc_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat']);
      load([edir 'tlock_grind_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat']);
    end
  else
    try
      load([edir 'tlock_statmc_sleep' num2str(sleep) '_iter' num2str(iter) '_statwinorig' num2str(statwinorig) '.mat']);
      load([edir 'tlock_grind_sleep' num2str(sleep) '_iter' num2str(iter) '_statwinorig' num2str(statwinorig) '.mat']);
    catch
      load([edir 'tlock_statmc_sleep' num2str(sleep) '_iter' num2str(iter) '.mat']);
      load([edir 'tlock_grind_sleep' num2str(sleep) '_iter' num2str(iter) '.mat']);
    end
  end
else
  load([edir 'tlock_statmc_sleep' num2str(sleep) '.mat']);
  load([edir 'tlock_grind_sleep' num2str(sleep) '.mat']);
end

soalist=[1 3 4 5 6 7 9];
timwinstatflag=0; % =1 for using only stat.time, =0 for using [-0.5 1.0];
timwin=[-0.5 1];
soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];

scalediff=1;

coloruse=varycolor(10);
% optimsied to be maximally apart and distinct
% 1  TacPAud
% 2  Nul
% 3  MS_synch
% 4  Tac
% 5  diff_TPAMSPN
% 6  diff_synchAsynch
% 7  MS
% 8  MS_asynch
% 9  Aud
% 10 MSpN

close all
clear tmp*
for ll=soalist
  cfg=[];
  cfg.operation='subtract';
  cfg.parameter='individual';
  grind_TPA_MSPN{ll,tt,ss}=ft_math(cfg,grind_tacPaud_save{ll,tt,ss},grind_tacMSpN_save{ll,tt,ss});
  
  cfg=[];
  if timwinstatflag==1
    cfg.latency=[statt_mc{ll,tt,ss}.time(1) statt_mc{ll,tt,ss}.time(end)];
  elseif timwinstatflag==0
    cfg.latency=timwin;
    stattimwin=[statt_mc{ll,tt,ss}.time(1) statt_mc{ll,tt,ss}.time(end)];
  end
  cfg.channel=statt_mc{ll,tt,ss}.label;
  
  tmpn2=ft_selectdata(cfg,grind_nultlock_save{ll,tt,ss});
  tmpn2.dimord='chan_time';
  tmpn2.avg=squeeze(mean(tmpn2.individual,1));
  tmpn2=rmfield(tmpn2,'individual');
  
  tmpt4=ft_selectdata(cfg,grind_tactlock_save{ll,tt,ss});
  tmpt4.dimord='chan_time';
  tmpt4.avg=squeeze(mean(tmpt4.individual,1));
  tmpt4=rmfield(tmpt4,'individual');
  
  tmpa9=ft_selectdata(cfg,grind_audtlock_save{ll,tt,ss});
  tmpa9.dimord='chan_time';
  tmpa9.avg=squeeze(mean(tmpa9.individual,1));
  tmpa9=rmfield(tmpa9,'individual');
  
  tmpb7=ft_selectdata(cfg,grind_MStlock_save{ll,tt,ss});
  tmpb7.dimord='chan_time';
  tmpb7.avg=squeeze(mean(tmpb7.individual,1));
  tmpb7=rmfield(tmpb7,'individual');
  
  tmpu1=ft_selectdata(cfg,grind_tacPaud_save{ll,tt,ss});
  tmpu1.dimord='chan_time';
  tmpu1.avg=squeeze(mean(tmpu1.individual,1));
  tmpu1=rmfield(tmpu1,'individual');
  tmpu1.mask=statt_mc{ll,tt,ss}.mask;
  
  tmpm10=ft_selectdata(cfg,grind_tacMSpN_save{ll,tt,ss});
  tmpm10.dimord='chan_time';
  tmpm10.avg=squeeze(mean(tmpm10.individual,1));
  tmpm10=rmfield(tmpm10,'individual');
  tmpm10.mask=statt_mc{ll,tt,ss}.mask;
  
  tmpd5=ft_selectdata(cfg,grind_TPA_MSPN{ll,tt,ss});
  tmpd5.dimord='chan_time';
  tmpd5.avg=scalediff*squeeze(mean(tmpd5.individual,1));
  tmpd5=rmfield(tmpd5,'individual');
  tmpd5.mask=statt_mc{ll,tt,ss}.mask;
  
  if timwinstatflag==0
    %       tmpmask=mean(statt_mc{ll,tt,ss}.mask(match_str(tmpd5.label,cfg.channel),:),1);
    %       tmpd5.mask=zeros(1,length(tmpd5.time));
    %       tmpd5.mask(dsearchn(tmpd5.time',stattimwin(1)):dsearchn(tmpd5.time',stattimwin(end)))=tmpmask;
    tmpmask=statt_mc{ll,tt,ss}.mask;
    tmpu1.mask=zeros(size(tmpu1.avg,1),length(tmpu1.time));
    tmpu1.mask(:,dsearchn(tmpu1.time',stattimwin(1)):dsearchn(tmpu1.time',stattimwin(end)))=tmpmask;
    tmpm10.mask=zeros(size(tmpm10.avg,1),length(tmpm10.time));
    tmpm10.mask(:,dsearchn(tmpm10.time',stattimwin(1)):dsearchn(tmpm10.time',stattimwin(end)))=tmpmask;
    tmpd5.mask=zeros(size(tmpd5.avg,1),length(tmpd5.time));
    tmpd5.mask(:,dsearchn(tmpd5.time',stattimwin(1)):dsearchn(tmpd5.time',stattimwin(end)))=tmpmask;
  end
  
  %   cfg=[];
  %   if timwinstatflag==1
  %     cfg.latency=[statt_mc{ll,tt,ss}.time(1) statt_mc{ll,tt,ss}.time(end)];
  %   elseif timwinstatflag==0
  %     cfg.latency=[-0.5 1];
  %     stattimwin=[statt_mc{ll,tt,ss}.time(1) statt_mc{ll,tt,ss}.time(end)];
  %   end
  %   tmp=ft_selectdata(cfg,grave_TPA_MSPN{ll,tt,ss});
  %   tmp.mask=statt_mc{ll,tt,ss}.mask;
  %   tmpm=ft_selectdata(cfg,grind_tacMSpN_save{ll,tt,ss});
  %   tmpu=ft_selectdata(cfg,grind_tacPaud_save{ll,tt,ss});
  %   tmpm.mask=statt_mc{ll,tt,ss}.mask;
  %   tmpu.mask=statt_mc{ll,tt,ss}.mask;
  %
  %   tmpm.avg=squeeze(mean(tmpm.individual,1));
  %   tmpu.avg=squeeze(mean(tmpu.individual,1));
  %   tmpm.var=squeeze(var(tmpm.individual,1));
  %   tmpu.var=squeeze(var(tmpu.individual,1));
  %   tmpm.dimord='chan_time';
  %   tmpu.dimord='chan_time';
  %   tmpm=rmfield(tmpm,'individual');
  %   tmpu=rmfield(tmpu,'individual');
  
  %   figure(10+ll);
  %   for cg=1:length(chanplot)
  %     cfg=[];
  %     cfg.parameter='avg';
  %     cfg.layout='elec1010.lay';
  %     cfg.ylim=[-3 7];
  %     cfg.linewidth=3;
  %     cfg.xlim=timwin;
  %     cfg.channel=chanplot{cg};
  %     cfg.graphcolor=coloruse([4 9 1],:);
  %     subplot(3,1,cg)
  %     ft_singleplotER(cfg,tmp4,tmp9,tmp1);
  %   end
  %   legend('Tactile','Auditory','Sum Unisensory')
  %   print(10+ll,[fdir 'erp_TacAudTacPAud_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  %
  %   figure(20+ll);
  %   for cg=1:length(chanplot)
  %     cfg=[];
  %     cfg.parameter='avg';
  %     cfg.layout='elec1010.lay';
  %     cfg.ylim=[-3 8];
  %     cfg.linewidth=3;
  %     cfg.xlim=timwin;
  %     cfg.channel=chanplot{cg};
  %     cfg.graphcolor=coloruse([2 7 10],:);
  %     subplot(3,1,cg)
  %     ft_singleplotER(cfg,tmp2,tmp7,tmp10);
  %   end
  %   legend('Null','Multisensory','MultSens + Null')
  %   print(20+ll,[fdir 'erp_NulMSandMSPN_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  
  %   figure(ll);
  %   for cg=1:length(chanplot)+1
  for cg=1:length(chanplot)
    cfg=[];
    cfg.parameter='avg';
    cfg.layout='elec1010.lay';
    cfg.ylim=[-3 7];
    cfg.linewidth=3;
    cfg.xlim=timwin;
    if cg>length(chanplot)
      cfg.channel=tmpd5.label(any(tmpd5.mask,2));
    else
      cfg.channel=chanplot{cg};
    end
    cfg.graphcolor=coloruse([2 4 9 7],:);
    cfg.interactive='no';
    %     subplot(4,1,cg)
    figure(ll+10*(cg-1))
    ft_singleplotER(cfg,tmpn2,tmpt4,tmpa9,tmpb7);
    hold on;plot(tmpn2.time,0,'k');
    set(gca,'XTick',[-.5:.1:1])
    set(gca,'XTickLabel',{'-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'})
    set(gca,'FontSize',30)
    title([])
    plot([0 0],cfg.ylim,'Color',coloruse(4,:),'LineWidth',6)
    plot([soades(ll) soades(ll)],cfg.ylim,'Color',coloruse(9,:),'LineWidth',6)
    axis([-0.55 1 cfg.ylim(1) cfg.ylim(2)])
    if cg==3
      legend('Null','Tactile','Auditory','Multisensory')
    end
  end
  %   print(ll,[fdir 'erp_presum_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  print(ll,[fdir 'erp_presum_FC_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  print(ll+10,[fdir 'erp_presum_OP_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  
  %   figure(30+ll);
  %   for cg=1:length(chanplot)+1
  for cg=1:length(chanplot)
    cfg=[];
    cfg.parameter='avg';
    cfg.layout='elec1010.lay';
    cfg.ylim=[-5 8];
    cfg.linewidth=3;
    cfg.xlim=timwin;
    if cg>length(chanplot)
      cfg.channel=tmpd5.label(any(tmpd5.mask,2));
    else
      cfg.channel=chanplot{cg};
    end
    cfg.graphcolor=coloruse([1 10 5],:);
    cfg.interactive='no';
    cfg.maskparameter='mask';
    cfg.maskstyle='box'; % default
    %     subplot(4,1,cg)
    %     subplot(1,2,cg)
    figure(ll+10*(cg+1))
    ft_singleplotER(cfg,tmpu1,tmpm10,tmpd5);
    hold on;plot(tmpn2.time,0,'k');
    set(gca,'XTick',[-.5:.1:1])
    set(gca,'XTickLabel',{'-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'})
    set(gca,'FontSize',30)
    title([])
    plot([stattimwin(1) stattimwin(1)],cfg.ylim,'k--','LineWidth',6)
    plot([stattimwin(2) stattimwin(2)],cfg.ylim,'k--','LineWidth',6)
    plot([0 0],cfg.ylim,'Color',coloruse(4,:),'LineWidth',6)
    plot([soades(ll) soades(ll)],cfg.ylim,'Color',coloruse(9,:),'LineWidth',6)
    axis([-0.55 1 cfg.ylim(1) cfg.ylim(2)])
    if cg==3
      legend('Sum Unisensory','MultSens + Null','SumUnisens - MultsensNull')
    end
  end
  %   print(30+ll,[fdir 'erp_tacPaud_MSpN_diff_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  print(ll+20,[fdir 'erp_tacPaud_MSpN_diff_FC_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  print(ll+30,[fdir 'erp_tacPaud_MSpN_diff_OP_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  
  %   cfg=[];
  %   cfg.avgoverchan='yes';
  %   cfg.channel=chanplot{1};
  %   tmp1=ft_selectdata(cfg,tmp);
  %   if timwinstatflag==0
  %     tmp1mask=mean(tmp.mask(match_str(tmp.label,cfg.channel),:),1);
  %     tmp1.mask=zeros(1,length(tmp1.time));
  %     tmp1.mask(dsearchn(tmp1.time',stattimwin(1)):dsearchn(tmp1.time',stattimwin(end)))=tmp1mask;
  %   end
  %   tmp1.mask=logical(ceil(tmp1.mask));
  %   tmpm1=ft_selectdata(cfg,tmpm);
  %   tmpm1.mask=logical(ceil(tmpm1.mask));
  %   tmpu1=ft_selectdata(cfg,tmpu);
  %   tmpu1.mask=logical(ceil(tmpu1.mask));
  %
  %   figure(10*ll);
  %   cfg=[];
  %   cfg.parameter='avg';
  %   cfg.layout='elec1010.lay';
  %   cfg.maskparameter='mask';
  %   cfg.maskstyle='box'; % default
  %   cfg.ylim=[-3 7];
  %   cfg.linewidth=3;
  %   ft_singleplotER(cfg,tmpu1,tmpm1,tmp1);
  %   print(10*ll,[fdir 'erp_final_FC_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  %
  %   cfg=[];
  %   cfg.avgoverchan='yes';
  %   cfg.channel=chanplot{2};
  %   tmp1=ft_selectdata(cfg,tmp);
  %   if timwinstatflag==0
  %     tmp1mask=mean(tmp.mask(match_str(tmp.label,cfg.channel),:),1);
  %     tmp1.mask=zeros(1,length(tmp1.time));
  %     tmp1.mask(dsearchn(tmp1.time',stattimwin(1)):dsearchn(tmp1.time',stattimwin(end)))=tmp1mask;
  %   end
  %   tmp1.mask=logical(ceil(tmp1.mask));
  %   tmpm1=ft_selectdata(cfg,tmpm);
  %   tmpm1.mask=logical(ceil(tmpm1.mask));
  %   tmpu1=ft_selectdata(cfg,tmpu);
  %   tmpu1.mask=logical(ceil(tmpu1.mask));
  %
  %   figure(10*ll+1);
  %   cfg=[];
  %   cfg.parameter='avg';
  %   cfg.layout='elec1010.lay';
  %   cfg.maskparameter='mask';
  %   cfg.maskstyle='box'; % default
  %   cfg.ylim=[-3 3];
  %   cfg.linewidth=3;
  %   ft_singleplotER(cfg,tmpu1,tmpm1,tmp1);
  %   print(10*ll+1,[fdir 'erp_final_OP_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  
  %   if [tt==2 && any(ll==[3 4 5 6])] || [tt==3 && any(ll==[4 6 7])]
  if any(statt_mc{ll,tt,ss}.mask(:))
    masktime=find(any(tmpd5.mask,1));
    cfg=[];
    cfg.parameter='avg';
    cfg.layout='elec1010.lay';
    cfg.maskalpha=0.5;
    cfg.zlim=[-5 5];
    cfg.highlight='on';
    cfg.highlightsize=12;
    cfg.xlim=[tmpd5.time(masktime(1)) tmpd5.time(masktime(end))];
    cfg.comment='no';
    %     if ll==3
    %       cfg.xlim=[.1 .35];
    %     elseif ll==4
    %       cfg.xlim=[.06 .42];
    %     elseif ll==5
    %       cfg.xlim=[-.04 .36];
    %     elseif ll==6
    %       cfg.xlim=[.1 .44];
    %     end
    sigchannels=tmpd5.label(find(ceil(mean(tmpd5.mask(:,dsearchn(tmpd5.time',cfg.xlim(1)):dsearchn(tmpd5.time',cfg.xlim(2))),2))));
    cfg.highlightchannel=sigchannels;
    figure(100+ll);
    ft_topoplotER(cfg,tmpu1);
    print(100+ll,[fdir 'erp_topoU_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    figure(110+ll);
    ft_topoplotER(cfg,tmpm10);
    print(110+ll,[fdir 'erp_topoM_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    figure(120+ll);
    ft_topoplotER(cfg,tmpd5);
    print(120+ll,[fdir 'erp_topoDiff_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    
    %     figure(50+ll);
    %     %   for cg=1:length(chanplot)
    %     cfg=[];
    %     cfg.parameter='avg';
    %     cfg.layout='elec1010.lay';
    %     cfg.ylim=[-3 8];
    %     cfg.linewidth=3;
    %     cfg.xlim=timwin;
    %     cfg.channel=sigchannels;
    %     cfg.graphcolor=coloruse([1 10 5],:);
    %     cfg.interactive='no';
    %     cfg.maskparameter='mask';
    %     cfg.maskstyle='box'; % default
    %     subplot(3,1,3)
    %     if timwinstatflag==0
    %       %       tmpmask=mean(statt_mc{ll,tt,ss}.mask(match_str(tmpd5.label,cfg.channel),:),1);
    %       %       tmpd5.mask=zeros(1,length(tmpd5.time));
    %       %       tmpd5.mask(dsearchn(tmpd5.time',stattimwin(1)):dsearchn(tmpd5.time',stattimwin(end)))=tmpmask;
    %       tmpmask=statt_mc{ll,tt,ss}.mask;
    %       tmpd5.mask=zeros(size(tmpd5.avg,1),length(tmpd5.time));
    %       tmpd5.mask(:,dsearchn(tmpd5.time',stattimwin(1)):dsearchn(tmpd5.time',stattimwin(end)))=tmpmask;
    %     end
    %     ft_singleplotER(cfg,tmpu1,tmpm10,tmpd5);
    %     hold on;plot(tmpn2.time,0,'k');
    %     set(gca,'XTick',[cfg.xlim(1):.1:cfg.xlim(2)])
    %     plot([stattimwin(1) stattimwin(1)],cfg.ylim,'k--')
    %     plot([stattimwin(2) stattimwin(2)],cfg.ylim,'k--')
    %     legend('Sum Unisensory','MultSens + Null','SumUnisens - MultsensNull')
    %     print(50+ll,[fdir 'erp_tacPaud_MSpN_diff_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  end
end






close all;
% funny 'temporal' (in)congruency contrast
for ll=soalist
  if ll<5 
    cfg=[];
    cfg.operation='subtract';
    cfg.parameter='individual';
    grind_diffMSsynch{ll,tt,ss}=ft_math(cfg,grind_tMSsynch_save{ll,tt,ss},grind_tMSasynch_save{ll,tt,ss});
    cfg.channel=statt_synch{ll,tt,ss}.label;
    
    tmpd6=ft_selectdata(cfg,grind_diffMSsynch{ll,tt,ss});
    tmpd6.dimord='chan_time';
    tmpd6.avg=scalediff*squeeze(mean(tmpd6.individual,1));
    tmpd6=rmfield(tmpd6,'individual');
    tmpd6.mask=statt_synch{ll,tt,ss}.mask;
    
    tmps3=ft_selectdata(cfg,grind_tMSsynch_save{ll,tt,ss});
    tmps3.dimord='chan_time';
    tmps3.avg=squeeze(mean(tmps3.individual,1));
    tmps3=rmfield(tmps3,'individual');
    tmps3.mask=statt_synch{ll,tt,ss}.mask;
    
    tmpa8=ft_selectdata(cfg,grind_tMSasynch_save{ll,tt,ss});
    tmpa8.dimord='chan_time';
    tmpa8.avg=squeeze(mean(tmpa8.individual,1));
    tmpa8=rmfield(tmpa8,'individual');
    tmpa8.mask=statt_synch{ll,tt,ss}.mask;
    
    stattimwin=[statt_synch{ll,tt,ss}.time(1) statt_synch{ll,tt,ss}.time(end)];
    stattime=statt_synch{ll,tt,ss}.time;
    
    
    if timwinstatflag==0
      %       tmpmask=mean(statt_mc{ll,tt,ss}.mask(match_str(tmpd5.label,cfg.channel),:),1);
      %       tmpd5.mask=zeros(1,length(tmpd5.time));
      %       tmpd5.mask(dsearchn(tmpd5.time',stattimwin(1)):dsearchn(tmpd5.time',stattimwin(end)))=tmpmask;
      if stattimwin(end)>tmps3.time(end)
        stattimwin(end)=tmps3.time(end);
      end
      tmpmask=statt_synch{ll,tt,ss}.mask;
      tmps3.mask=zeros(size(tmps3.avg,1),length(tmps3.time));
      tmps3.mask(:,dsearchn(tmps3.time',stattimwin(1)):dsearchn(tmps3.time',stattimwin(end)))=tmpmask(:,1:dsearchn(stattime',stattimwin(end)));
      tmpa8.mask=zeros(size(tmpa8.avg,1),length(tmpa8.time));
      tmpa8.mask(:,dsearchn(tmpa8.time',stattimwin(1)):dsearchn(tmpa8.time',stattimwin(end)))=tmpmask(:,1:dsearchn(stattime',stattimwin(end)));
      tmpd6.mask=zeros(size(tmpd6.avg,1),length(tmpd6.time));
      tmpd6.mask(:,dsearchn(tmpd6.time',stattimwin(1)):dsearchn(tmpd6.time',stattimwin(end)))=tmpmask(:,1:dsearchn(stattime',stattimwin(end)));
    end
    
    %     cfg=[];
    %     if timwinstatflag==1
    %       cfg.latency=stattimwin;
    %     elseif timwinstatflag==0
    %       cfg.latency=[-0.5 1];
    %     end
    %     tmp=ft_selectdata(cfg,grave_TMSs_TMSa{ll,tt,ss});
    %     tmp.mask=statt_synch{ll,tt,ss}.mask;
    %     tmps=ft_selectdata(cfg,grind_tMSsynch_save{ll,tt,ss});
    %     tmpa=ft_selectdata(cfg,grind_tMSasynch_save{ll,tt,ss});
    %     tmps.mask=statt_synch{ll,tt,ss}.mask;
    %     tmpa.mask=statt_synch{ll,tt,ss}.mask;
    %
    %     tmps.avg=squeeze(mean(tmps.individual,1));
    %     tmpa.avg=squeeze(mean(tmpa.individual,1));
    %     tmps.var=squeeze(var(tmps.individual,1));
    %     tmpa.var=squeeze(var(tmpa.individual,1));
    %     tmps.dimord='chan_time';
    %     tmpa.dimord='chan_time';
    %     tmps=rmfield(tmps,'individual');
    %     tmpa=rmfield(tmpa,'individual');
    %
    %     cfg=[];
    %     cfg.avgoverchan='yes';
    %     cfg.channel=chanplot{1};
    %     tmp1=ft_selectdata(cfg,tmp);
    %     tmps1=ft_selectdata(cfg,tmps);
    %     tmpa1=ft_selectdata(cfg,tmpa);
    %     if timwinstatflag==0
    %       tmp1mask=mean(tmp.mask(match_str(tmp.label,cfg.channel),:),1);
    %       tmp1.mask=zeros(1,length(tmp1.time));
    %       tmp1.mask(dsearchn(tmp1.time',stattimwin(1)):dsearchn(tmp1.time',stattimwin(end)))=tmp1mask;
    %       tmps1mask=mean(tmps.mask(match_str(tmps.label,cfg.channel),:),1);
    %       tmps1.mask=zeros(1,length(tmps1.time));
    %       tmps1.mask(dsearchn(tmps1.time',stattimwin(1)):dsearchn(tmps1.time',stattimwin(end)))=tmps1mask;
    %       tmpa1mask=mean(tmpa.mask(match_str(tmpa.label,cfg.channel),:),1);
    %       tmpa1.mask=zeros(1,length(tmpa1.time));
    %       tmpa1.mask(dsearchn(tmpa1.time',stattimwin(1)):dsearchn(tmpa1.time',stattimwin(end)))=tmpa1mask;
    %     end
    %     tmp1.mask=logical(ceil(tmp1.mask));
    %     tmps1.mask=logical(ceil(tmps1.mask));
    %     tmpa1.mask=logical(ceil(tmpa1.mask));
    
    %     figure(40+ll);
    %     for cg=1:length(chanplot)+1
    for cg=1:length(chanplot)
      cfg=[];
      cfg.parameter='avg';
      cfg.layout='elec1010.lay';
      cfg.ylim=[-5 14];
      cfg.linewidth=3;
      cfg.xlim=timwin;
      if cg>length(chanplot)
        cfg.channel=tmpd6.label(any(tmpd6.mask,2));
      else
        cfg.channel=chanplot{cg};
      end
      cfg.graphcolor=coloruse([3 8 6],:);
      cfg.interactive='no';
      cfg.maskparameter='mask';
      cfg.maskstyle='box'; % default
      %       subplot(4,1,cg)
      %       subplot(1,2,cg)
      if cg==1
        figure(ll+40);
      elseif cg==2
        figure(ll+45);
      end
      ft_singleplotER(cfg,tmps3,tmpa8,tmpd6);
      hold on;plot(tmps3.time,0,'k');
      set(gca,'XTick',[cfg.xlim(1):.1:cfg.xlim(2)])
      plot([stattimwin(1) stattimwin(1)],cfg.ylim,'k--')
      plot([stattimwin(2) stattimwin(2)],cfg.ylim,'k--')
      plot([0 0],cfg.ylim,'Color',coloruse(4,:))
      plot([soades(10-ll) soades(10-ll)],cfg.ylim,'Color',coloruse(9,:))
      axis([-0.55 1 cfg.ylim(1) cfg.ylim(2)])
      if cg==1
        legend('MultSens Synchronous','MultSens Aynchronous','MSsynch - MSasynch')
      end
    end
    %     print(40+ll,[fdir 'erp_MSsynch_asynch_diff_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    print(ll+40,[fdir 'erp_MSsynch_asynch_diff_FC_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    print(ll+45,[fdir 'erp_MSsynch_asynch_diff_OP_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    
    %   figure(40+ll);
    %     for cg=1:length(chanplot)
    %       cfg=[];
    %       cfg.parameter='avg';
    %       cfg.layout='elec1010.lay';
    %       cfg.ylim=[-5 14];
    %       cfg.linewidth=3;
    %       cfg.xlim=timwin;
    %       cfg.channel=chanplot{cg};
    %       cfg.graphcolor=coloruse([3 8 6],:);
    %       subplot(3,1,cg)
    %       ft_singleplotER(cfg,tmp3,tmp8,tmp6);
    %     end
    %     legend('MultSens Synchronous','MultSens Aynchronous','MSsynch - MSasynch')
    %     print(40+ll,[fdir 'erp_MSsynch_asynch_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    
    
    %     figure(10*ll);
    %     cfg=[];
    %     cfg.parameter='avg';
    %     cfg.layout='elec1010.lay';
    %     cfg.maskparameter='mask';
    %     cfg.maskstyle='box'; % default
    %     cfg.ylim=[-5 14];
    %     cfg.linewidth=3;
    %     ft_singleplotER(cfg,tmps1,tmpa1,tmp1);
    %     print(10*ll,[fdir 'erp_synch_final_FC_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    
    %     cfg=[];
    %     cfg.avgoverchan='yes';
    %     cfg.channel=chanplot{2};
    %     tmp1=ft_selectdata(cfg,tmp);
    %     tmps1=ft_selectdata(cfg,tmps);
    %     tmpa1=ft_selectdata(cfg,tmpa);
    %     if timwinstatflag==0
    %       tmp1mask=mean(tmp.mask(match_str(tmp.label,cfg.channel),:),1);
    %       tmp1.mask=zeros(1,length(tmp1.time));
    %       tmp1.mask(dsearchn(tmp1.time',stattimwin(1)):dsearchn(tmp1.time',stattimwin(end)))=tmp1mask;
    %       tmps1mask=mean(tmps.mask(match_str(tmps.label,cfg.channel),:),1);
    %       tmps1.mask=zeros(1,length(tmps1.time));
    %       tmps1.mask(dsearchn(tmps1.time',stattimwin(1)):dsearchn(tmps1.time',stattimwin(end)))=tmps1mask;
    %       tmpa1mask=mean(tmpa.mask(match_str(tmpa.label,cfg.channel),:),1);
    %       tmpa1.mask=zeros(1,length(tmpa1.time));
    %       tmpa1.mask(dsearchn(tmpa1.time',stattimwin(1)):dsearchn(tmpa1.time',stattimwin(end)))=tmpa1mask;
    %     end
    %     tmp1.mask=logical(ceil(tmp1.mask));
    %     tmps1.mask=logical(ceil(tmps1.mask));
    %     tmpa1.mask=logical(ceil(tmpa1.mask));
    
    %     figure(10*ll+1);
    %     cfg=[];
    %     cfg.parameter='avg';
    %     cfg.layout='elec1010.lay';
    %     cfg.maskparameter='mask';
    %     cfg.maskstyle='box'; % default
    %     cfg.ylim=[-5 3];
    %     cfg.linewidth=3;
    %     ft_singleplotER(cfg,tmps1,tmpa1,tmp1);
    %     print(10*ll+1,[fdir 'erp_synch_final_OP_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    
    %     if [tt==2 && any(ll==[1 4])] || [tt==3 && any(ll==[1 4])]
    if [tt==2 && any(ll==[1 4])] || [tt==3 && any(ll==[1])]
      %       masktime=find(any(tmp.mask,1));
      masktime=find(any(tmpd6.mask,1));
      if find(diff(masktime)>2)
        error('stats gave different answer');
      end
      cfg=[];
      cfg.parameter='avg';
      cfg.layout='elec1010.lay';
      cfg.maskalpha=0.5;
      cfg.zlim=[-6 6];
      cfg.highlight='on';
      cfg.highlightsize=12;
      cfg.xlim=[tmpd6.time(masktime(1)) tmpd6.time(masktime(end))];
      %       cfg.highlightchannel=tmp.label(find(ceil(mean(tmp.mask(:,dsearchn(stattime',cfg.xlim(1)):dsearchn(stattime',cfg.xlim(2))),2))));
      %       cfg.highlightchannel=tmpd6.label(find(ceil(mean(tmpd6.mask(:,dsearchn(stattime',cfg.xlim(1)):dsearchn(stattime',cfg.xlim(2))),2))));
      cfg.highlightchannel=tmpd6.label(find(ceil(mean(tmpd6.mask(:,dsearchn(tmpd6.time',cfg.xlim(1)):dsearchn(tmpd6.time',cfg.xlim(2))),2))));
      figure(50+ll);
      ft_topoplotER(cfg,tmps3);
      print(50+ll,[fdir 'erp_synch_topoS_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
      figure(60+ll);
      ft_topoplotER(cfg,tmpa8);
      print(60+ll,[fdir 'erp_synch_topoA_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
      figure(70+ll);
      ft_topoplotER(cfg,tmpd6);
      print(70+ll,[fdir 'erp_synch_topoDiff_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    end
    if [tt==2 && any(ll==[3])] || [tt==3 && any(ll==[3]) && ss~=11]
      %       masktime=find(any(tmp.mask,1));
      masktime=find(any(tmpd6.mask,1));
      % two different effects going on in these 2 different time windows
      if length(find(diff(masktime)>2))~=1
        error('stats gave something different');
      end
      masktime1=masktime(1:find(diff(masktime)>2));
      masktime2=masktime(find(diff(masktime)>2)+1:end);
      cfg=[];
      cfg.parameter='avg';
      cfg.layout='elec1010.lay';
      cfg.maskalpha=0.5;
      cfg.zlim=[-6 6];
      cfg.highlight='on';
      cfg.highlightsize=12;
      %       cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
      cfg.xlim=[tmpd6.time(masktime1(1)) tmpd6.time(masktime1(end))];
      %       cfg.highlightchannel=tmp.label(find(ceil(mean(tmp.mask(:,dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2))),2))));
      %       cfg.highlightchannel=tmp.label(find(ceil(mean(tmp.mask(:,dsearchn(stattime',cfg.xlim(1)):dsearchn(stattime',cfg.xlim(2))),2))));
      cfg.highlightchannel=tmpd6.label(find(ceil(mean(tmpd6.mask(:,dsearchn(tmpd6.time',cfg.xlim(1)):dsearchn(tmpd6.time',cfg.xlim(2))),2))));
      figure(50+ll);
      ft_topoplotER(cfg,tmps3);
      print(50+ll,[fdir 'erp_synchm1_topoS_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
      figure(60+ll);
      ft_topoplotER(cfg,tmpa8);
      print(60+ll,[fdir 'erp_synchm1_topoA_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
      figure(70+ll);
      ft_topoplotER(cfg,tmpd6);
      print(70+ll,[fdir 'erp_synchm1_topoDiff_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
      cfg.xlim=[tmpd6.time(masktime2(1)) tmpd6.time(masktime2(end))];
      %       cfg.highlightchannel=tmp.label(find(ceil(mean(tmp.mask(:,dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2))),2))));
      cfg.highlightchannel=tmpd6.label(find(ceil(mean(tmpd6.mask(:,dsearchn(tmpd6.time',cfg.xlim(1)):dsearchn(tmpd6.time',cfg.xlim(2))),2))));
      figure(80+ll);
      ft_topoplotER(cfg,tmps3);
      print(80+ll,[fdir 'erp_synchm2_topoS_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
      figure(90+ll);
      ft_topoplotER(cfg,tmpa8);
      print(90+ll,[fdir 'erp_synchm2_topoA_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
      figure(100+ll);
      ft_topoplotER(cfg,tmpd6);
      print(100+ll,[fdir 'erp_synchm2_topoDiff_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    end
  end % ll
  
end



%%  Comparing awake to asleep

printflag=0;
plotflag=0;
dostats=0;

for tt=2
  for ll=[1 3 4 5 6 7 9]
    clearvars -except ll tt sub edir ddir sdir iiuse plotflag printflag dostats fstats*
    
    for sleep=[0 1]
      
      subuse=iiuse;
      
      submin=subuse(1)-1;
      subuseind=0;
      for ii=subuse
        %       for ii=setdiff(subuse,[8 9 10 12 14 15 16 17 18])
        cd([edir sub{ii} ])
        load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '.mat'])
        
        tk0=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(0) '.mat']);
        tk1=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '.mat']);
        
        % THis is preliminary...how best to included all stages later on?
        if sleep==0
          ss=10; % awake
        elseif sleep==1
          ss=23; % this is concatenation of N2 and N3
        end
        
        %         for ss=ssuse
        if sleep==0 %load for both at sleep0 then will still be in memory for sleep1
          numtrt(ll,tt,10,ii-submin)=numt_trials(ll,tt,10); % does this need to be per 'sleep01' as well?
          if audtacflag
            numtra(ll,tt,10,ii-submin)=numa_trials(ll,tt,10);
          end
          load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(1) '.mat'],'num*trials')
          numtrt(ll,tt,23,ii-submin)=numt_trials(ll,tt,23); % does this need to be per 'sleep01' as well?
          if audtacflag
            numtra(ll,tt,23,ii-submin)=numa_trials(ll,tt,23);
          end
        end
        %         end
        
        % discard from both if either sleep/wake doesn't have 'enough' (what is enough? 20?)
        if numtrt(ll,tt,10,ii-submin)<20 || numtrt(ll,tt,23,ii-submin)<20
          subuse=setdiff(subuse,ii);
        else
          subuseind=subuseind+1;
          tlock_tacPaud_each{subuseind}=tlock_tacPaud{ll,tt,ss};
          tlock_audPtac_each{subuseind}=tlock_audPtac{ll,tt,ss};
          tlock_tacMSpN_each{subuseind}=tlock_tacMSpN{ll,tt,ss};
          tlock_audMSpN_each{subuseind}=tlock_audMSpN{ll,tt,ss};
          if sleep==1
            fstats_tpa(:,ll,tt,:,subuseind)=featurestats_tacPaud(:,ll,tt,[10 23],ii);
            fstats_apt(:,ll,tt,:,subuseind)=featurestats_audPtac(:,ll,tt,[10 23],ii);
            fstats_tmspn(:,ll,tt,:,subuseind)=featurestats_tacMSpN(:,ll,tt,[10 23],ii);
            fstats_amspn(:,ll,tt,:,subuseind)=featurestats_audMSpN(:,ll,tt,[10 23],ii);
          end
        end
        %       clear *_s0
        clear tlock*N tlock*tac tlock*aud
      end
      subuseindfinal=subuseind;
      
      for ii=1:subuseindfinal
        cfg=[];
        cfg.operation='subtract';
        cfg.parameter='avg';
        tlock_TPA_MSPN{ii}=ft_math(cfg,tlock_tacPaud_each{ii},tlock_tacMSpN_each{ii});
        tlock_APT_MSPN{ii}=ft_math(cfg,tlock_audPtac_each{ii},tlock_audMSpN_each{ii});
      end
      cfg=[];
      cfg.keepindividual='yes';
      grind_TPA_MSPN{ss}=ft_timelockgrandaverage(cfg,tlock_TPA_MSPN{:});
      grind_APT_MSPN{ss}=ft_timelockgrandaverage(cfg,tlock_APT_MSPN{:});
      cfg=[];
      grave_TPA_MSPN{ss}=ft_timelockgrandaverage(cfg,tlock_TPA_MSPN{:});
      grave_APT_MSPN{ss}=ft_timelockgrandaverage(cfg,tlock_APT_MSPN{:});
      
      if plotflag
        figure(20);
        for ii=1:subuseindfinal
          if ss==23
            subplot(3,subuseindfinal,ii);               imagesc(tlock_TPA_MSPN{ii}.time,1:63,tlock_TPA_MSPN{ii}.avg);caxis([-6 6])
          elseif ss==10
            subplot(3,subuseindfinal,subuseindfinal+ii);imagesc(tlock_TPA_MSPN{ii}.time,1:63,tlock_TPA_MSPN{ii}.avg);caxis([-6 6])
          end
        end
        figure(21);
        for ii=1:subuseindfinal
          if ss==23
            subplot(3,subuseindfinal,ii);               imagesc(tlock_APT_MSPN{ii}.time,1:63,tlock_APT_MSPN{ii}.avg);caxis([-6 6])
          elseif ss==10
            subplot(3,subuseindfinal,subuseindfinal+ii);imagesc(tlock_APT_MSPN{ii}.time,1:63,tlock_APT_MSPN{ii}.avg);caxis([-6 6])
          end
        end
      end
      
      
    end % end ss
    
    cfg=[];
    cfg.operation='subtract';
    cfg.parameter='individual';
    grind_TPA_MSPN_sleepdiff=ft_math(cfg,grind_TPA_MSPN{23},grind_TPA_MSPN{10})
    grind_APT_MSPN_sleepdiff=ft_math(cfg,grind_APT_MSPN{23},grind_APT_MSPN{10})
    
    if plotflag
      figure(20);
      for ii=1:subuseindfinal
        subplot(3,subuseindfinal,2*subuseindfinal+ii);imagesc(grind_TPA_MSPN_sleepdiff.time,1:63,squeeze(grind_TPA_MSPN_sleepdiff.individual(ii,:,:)));caxis([-6 6])
      end
      figure(21);
      for ii=1:subuseindfinal
        subplot(3,subuseindfinal,2*subuseindfinal+ii);imagesc(grind_APT_MSPN_sleepdiff.time,1:63,squeeze(grind_APT_MSPN_sleepdiff.individual(ii,:,:)));caxis([-6 6])
      end
      if printflag
        print(20,['D:\audtac\figs\indiv_tpa_mspn_N23_W_diff_cond' num2str(ll) num2str(tt) '.png'],'-dpng')
        print(21,['D:\audtac\figs\indiv_apt_mspn_N23_W_diff_cond' num2str(ll) num2str(tt) '.png'],'-dpng')
      end
    end
    
    if plotflag
      topoplot_highlight(11,grave_TPA_MSPN{23},[-0.5 0.6],[]);
      topoplot_highlight(12,grave_TPA_MSPN{10},[-0.5 0.6],[]);
      topoplot_highlight(13,grave_APT_MSPN{23},[-0.1 1.1],[]);
      topoplot_highlight(14,grave_APT_MSPN{10},[-0.1 1.1],[]);
      
      if printflag
        print(11,['D:\audtac\figs\grave_tpamspnN23_topoOverTime_cond' num2str(ll) num2str(tt) '.png'],'-dpng')
        print(12,['D:\audtac\figs\grave_tpamspnW_topoOverTime_cond' num2str(ll) num2str(tt) '.png'],'-dpng')
        print(13,['D:\audtac\figs\grave_aptmspnN23_topoOverTime_cond' num2str(ll) num2str(tt) '.png'],'-dpng')
        print(14,['D:\audtac\figs\grave_aptmspnW_topoOverTime_cond' num2str(ll) num2str(tt) '.png'],'-dpng')
      end
    end
    
    if dostats
      load eeg1010_neighb
      
      nsub=subuseindfinal;
      
      cfg=[];
      cfg.latency=[-.1 .5];
      cfg.neighbours=neighbours;
      % cfg.parameter='avg';
      cfg.parameter='individual';
      cfg.method='montecarlo';
      % cfg.method='analytic';
      cfg.numrandomization=2000;
      % cfg.correctm='holm';
      cfg.correctm='cluster';
      cfg.clusteralpha = 0.05;
      cfg.clusterstatistic = 'maxsum';
      cfg.minnbchan = 2;
      cfg.statistic='depsamplesT';
      % cfg.statistic='indepsamplesregrT';
      % cfg.statistic='indepsamplesT';
      cfg.design=zeros(2,2*nsub);
      cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
      cfg.design(2,:)=[1:nsub 1:nsub];
      cfg.ivar=1;
      cfg.uvar=2;
      %     statt_mc=ft_timelockstatistics(cfg, grind_tacPaud, grind_tacMSpN);
      %     stata_mc=ft_timelockstatistics(cfg, grind_audPtac, grind_audMSpN);
      statt_mc=ft_timelockstatistics(cfg, grind_TPA_MSPN{23}, grind_TPA_MSPN{10});
      stata_mc=ft_timelockstatistics(cfg, grind_APT_MSPN{23}, grind_APT_MSPN{10});
      
      grave_TPA_MSPN_sleepdiff.avg=squeeze(mean(grind_TPA_MSPN_sleepdiff.individual,1));
      grave_TPA_MSPN_sleepdiff.time=grind_TPA_MSPN_sleepdiff.time;
      grave_TPA_MSPN_sleepdiff.label=grind_TPA_MSPN_sleepdiff.label;
      grave_TPA_MSPN_sleepdiff.dimord='chan_time';
      grave_APT_MSPN_sleepdiff.avg=squeeze(mean(grind_APT_MSPN_sleepdiff.individual,1));
      grave_APT_MSPN_sleepdiff.time=grind_APT_MSPN_sleepdiff.time;
      grave_APT_MSPN_sleepdiff.label=grind_APT_MSPN_sleepdiff.label;
      grave_APT_MSPN_sleepdiff.dimord='chan_time';
      
      if plotflag
        topoplot_highlight(22,grave_TPA_MSPN_sleepdiff,[statt_mc.time(1) statt_mc.time(end)],statt_mc);
        topoplot_highlight(23,grave_APT_MSPN_sleepdiff,[stata_mc.time(1) stata_mc.time(end)],stata_mc);
        if printflag
          print(22,['D:\audtac\figs\grSleepdiff_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) '.png'],'-dpng')
          print(23,['D:\audtac\figs\grSleepdiff_topoOverTime_mc_at_cond' num2str(ll) num2str(tt) '.png'],'-dpng')
        end
      end
    end
    
  end
  % do some processing on fstats* here, over all 'll' but within a 'tt'
  
  fstats_tpa(3,:,tt,:,:)
  
end
save([edir 'tlockSLEEP01_numtrlltt.mat'],'numtr*','grind*','fstats*');

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


%% Across subjects combining stats

% for ii=2:4
%   cd([edir sub{ii} ])
%   stat{ii}=load(['stat_erp_' sub{ii} '.mat']);
%   mmsu{ii}=load(['tlock_diffs_' sub{ii} '.mat']);
%   if ii==2  % find channels in common to all subjects
%     labelkeep=stat{ii}.stat1_s0{3}.label;
%   else
%     labelkeep=intersect(labelkeep,stat{ii}.stat1_s0{3}.label);
%   end
% end
%
% chanuse=match_str(labelkeep,'Fz');
%
%
% for ii=2:4
%   for ll=3:7
%     mask1_s0(:,:,ll,ii)=stat{ii}.stat1_s0{ll}.mask(match_str(stat{ii}.stat1_s0{ll}.label,labelkeep),:);
%     mask2_s0(:,:,ll,ii)=stat{ii}.stat2_s0{ll}.mask(match_str(stat{ii}.stat2_s0{ll}.label,labelkeep),:);
%     mask1_sall(:,:,ll,ii)=stat{ii}.stat1_sall{ll}.mask(match_str(stat{ii}.stat1_sall{ll}.label,labelkeep),:);
%     mask2_sall(:,:,ll,ii)=stat{ii}.stat2_sall{ll}.mask(match_str(stat{ii}.stat2_sall{ll}.label,labelkeep),:);
%     stat1_s0(:,:,ll,ii)=stat{ii}.stat1_s0{ll}.stat(match_str(stat{ii}.stat1_s0{ll}.label,labelkeep),:);
%     stat2_s0(:,:,ll,ii)=stat{ii}.stat2_s0{ll}.stat(match_str(stat{ii}.stat2_s0{ll}.label,labelkeep),:);
%     stat1_sall(:,:,ll,ii)=stat{ii}.stat1_sall{ll}.stat(match_str(stat{ii}.stat1_sall{ll}.label,labelkeep),:);
%     stat2_sall(:,:,ll,ii)=stat{ii}.stat2_sall{ll}.stat(match_str(stat{ii}.stat2_sall{ll}.label,labelkeep),:);
%
%     tpa_s0(:,:,ll,ii)=mmsu{ii}.tlock_tpa_mtamn_s0{ll}.avg(match_str(stat{ii}.stat1_s0{ll}.label,labelkeep),dsearchn(mmsu{ii}.tlock_tpa_mtamn_s0{ll}.time',-.2):dsearchn(mmsu{ii}.tlock_tpa_mtamn_s0{ll}.time',.5));
%     tpa_sall(:,:,ll,ii)=mmsu{ii}.tlock_tpa_mtamn_sall{ll}.avg(match_str(stat{ii}.stat1_s0{ll}.label,labelkeep),dsearchn(mmsu{ii}.tlock_tpa_mtamn_sall{ll}.time',-.2):dsearchn(mmsu{ii}.tlock_tpa_mtamn_sall{ll}.time',.5));
%     apt_s0(:,:,ll,ii)=mmsu{ii}.tlock_apt_matmn_s0{ll}.avg(match_str(stat{ii}.stat1_s0{ll}.label,labelkeep),dsearchn(mmsu{ii}.tlock_apt_matmn_s0{ll}.time',-.2):dsearchn(mmsu{ii}.tlock_apt_matmn_s0{ll}.time',.5));
%     apt_sall(:,:,ll,ii)=mmsu{ii}.tlock_apt_matmn_sall{ll}.avg(match_str(stat{ii}.stat1_s0{ll}.label,labelkeep),dsearchn(mmsu{ii}.tlock_apt_matmn_sall{ll}.time',-.2):dsearchn(mmsu{ii}.tlock_apt_matmn_sall{ll}.time',.5));
%
%   end
% end
%
% figure;
% for ll=3:7
%   subplot(1,5,ll-2);imagesc((nanmean(mask2_sall(:,:,ll,2:4),4)));caxis([0 1]);
% end
% figure;
% for ll=3:7
%   subplot(1,5,ll-2);imagesc((nanmean(mask1_sall(:,:,ll,2:4),4)));caxis([0 1]);
% end
% figure;
% for ll=3:7
%   subplot(1,5,ll-2);imagesc((nanmean(stat2_sall(:,:,ll,2:4),4)));caxis([-3 3]);
% end
% figure;
% for ll=3:7
%   subplot(1,5,ll-2);imagesc((nanmean(stat1_sall(:,:,ll,2:4),4)));caxis([-3 3]);
% end
%
% figure;for ll=3:7,subplot(1,5,ll-2);imagesc(-.2:.001:.5,1:length(labelkeep),mean(tpa_s0(:,:,ll,2:4),4));caxis([-4 4]);end
% figure;for ll=3:7,subplot(1,5,ll-2);imagesc(-.2:.001:.5,1:length(labelkeep),mean(tpa_sall(:,:,ll,2:4),4));caxis([-4 4]);end
% figure;for ll=3:7,subplot(1,5,ll-2);imagesc(-.2:.001:.5,1:length(labelkeep),mean(apt_s0(:,:,ll,2:4),4));caxis([-4 4]);end
% figure;for ll=3:7,subplot(1,5,ll-2);imagesc(-.2:.001:.5,1:length(labelkeep),mean(apt_sall(:,:,ll,2:4),4));caxis([-4 4]);end
%
%
% figure(111); % (Sum of Unisensory) minus (Multisensory plus Nul)
% for ll=3:7
%   subplot(2,5,ll-2);plot(-.2:.001:.5,mean(tpa_sall(chanuse,:,ll,2:4),4),'k');axis([-.2 .5 -7 7]);
%   legend('(T+A)-(TA-N), time0 is tactile')
%   subplot(2,5,ll-2+5);plot(-.2:.001:.5,mean(apt_sall(chanuse,:,ll,2:4),4),'k');axis([-.2 .5 -7 7])
%   legend('(A+T)-(AT-N), time0 is auditory')
% end


%%  Individual t-score

% ignore this for now.

time=-0.1:0.001:0.4;

ccc=nan(length(time),63,9,4,28);

for ii=5:18
  cd([edir sub{ii} ])
  load(['tlock_diffs_' sub{ii} '.mat']);
  if ii>7
    soalist=[1 3 4 5 6 7 9];
  else
    soalist=[3 4 5 6 7];
  end
  
  if ii==5
    labels=tlock_tac{5,1}.label;
  end
  
  %   chanuse=match_str(labels,tlock_tac{5,1}.label);
  %
  %   for ll=soalist
  %     for tt=1:4
  %
  %
  %       for teatime=1:length(time)
  %
  %         audtac=tlock_tac{ll,tt}.trial(:,:,dsearchn(tlock_tac{ll,tt}.time',time(teatime)));
  %         nul=tlock_nul_s0{ll+40,tt}.trial(:,:,dsearchn(tlock_nul_s0{ll+40,tt}.time',time(teatime)));
  %         tac=tlock_tac{10,tt}.trial(ttrialkept{ll,tt},:,dsearchn(tlock_tac{10,tt}.time',time(teatime)));
  %         aud=tlock_aud_s0{ll+40,tt}.trial(:,:,dsearchn(tlock_aud_s0{ll+40,tt}.time',time(teatime)));
  %
  %
  %         data=[tac; aud; audtac; nul];
  %         numtr=size(tac,1);
  %         design=zeros(size(data,1),5);
  %         design(:,1)=[zeros(0*numtr,1); ones(numtr,1); zeros(3*numtr,1)];
  %         design(:,2)=[zeros(1*numtr,1); ones(numtr,1); zeros(2*numtr,1)];
  %         design(:,3)=[zeros(2*numtr,1); ones(numtr,1); zeros(1*numtr,1)];
  %         design(:,4)=[zeros(3*numtr,1); ones(numtr,1); zeros(0*numtr,1)];
  %         design(:,5)=[zeros(0*numtr,1); ones(4*numtr,1); zeros(0*numtr,1)];
  %         cfg=[];
  %         cfg.glm.statistic='beta';
  %         cfg.glm.standardise=0;
  %         stat=ft_statfun_glm(cfg,data',design');
  %         beta=reshape(stat.stat,[size(design,2) size(data,2)]);
  %         con(teatime,chanuse,ll,tt,ii)=[1 1 -1 -1 0]*beta;
  %       end
  %     end
  %   end
  
  
  numtr=size(tlock_tac_s0{ll,tt}.trialinfo,1);
  mu=repmat(nanmean(tlock_tac_s0{ll,tt}.trial,1),[numtr 1 1]);
  sigma=repmat((nansum((tlock_tac_s0{ll,tt}.trial-mu).^2,1)/numtr).^0.5,[numtr 1 1]);
  zvalue=(tlock_tac_s0{ll,tt}.trial-mu)./sigma;
  zvalue=tlock_tac_s0{ll,tt}.trial./sigma;
  zvalue=squeeze(mu(1,:,:)./sigma(1,:,:));
  
  cfg.design=[];
  cfg.design(:,1)=[ ones(1,size(tlock_aud_s0{10,tt,ss}.trial,1)) ones(1,size(tlock_tac_s0{30+ll}.trial,1)) 2*ones(1,size(tlock_tac_s0{ll}.trial,1)) 2*ones(1,size(tlock_nul_s0{10,tt,ss}.trial,1))]';
  %   cfg.design(:,2)=[ ones(1,size(tlock_aud{10,tt,ss}.trial,1)) 2*ones(1,size(tlock_tac{30+ll}.trial,1)) 3*ones(1,size(tlock_tac{ll}.trial,1)) 4*ones(1,size(tlock_nul{10,tt,ss}.trial,1))]';
  stat1_s0{ll} = ft_timelockstatistics(cfg, tlock_aud_s0{10,tt,ss}, tlock_tac_s0{30+ll}, tlock_tac_s0{ll}, tlock_nul_s0{10,tt,ss})
  tlock_tpa_mtamn_s0{ll}.mask=stat1_s0{ll}.mask;
  
  cfg.design=[];
  cfg.design(:,1)=[ ones(1,size(tlock_aud_sall{10,tt,ss}.trial,1)) ones(1,size(tlock_tac_sall{30+ll}.trial,1)) 2*ones(1,size(tlock_tac_sall{ll}.trial,1)) 2*ones(1,size(tlock_nul_sall{10,tt,ss}.trial,1))]';
  %   cfg.design(:,2)=[ ones(1,size(tlock_aud{10,tt,ss}.trial,1)) 2*ones(1,size(tlock_tac{30+ll}.trial,1)) 3*ones(1,size(tlock_tac{ll}.trial,1)) 4*ones(1,size(tlock_nul{10,tt,ss}.trial,1))]';
  stat1_sall{ll} = ft_timelockstatistics(cfg, tlock_aud_sall{10,tt,ss}, tlock_tac_sall{30+ll}, tlock_tac_sall{ll}, tlock_nul_sall{10,tt,ss})
  tlock_tpa_mtamn_sall{ll}.mask=stat1_sall{ll}.mask;
end

% end
% thus, stat1_* is with Aud-alone at time zero, shifted tac, and AT with aud-first for 3, and tac-first for 7

for ll=soalist
  cfg.design=[];
  cfg.design(:,1)=[ ones(1,size(tlock_tac_s0{10,tt,ss}.trial,1)) ones(1,size(tlock_aud_s0{30+ll}.trial,1)) 2*ones(1,size(tlock_aud_s0{ll}.trial,1)) 2*ones(1,size(tlock_nul_s0{10,tt,ss}.trial,1))]';
  stat2_s0{ll} = ft_timelockstatistics(cfg, tlock_tac_s0{10,tt,ss}, tlock_aud_s0{30+ll}, tlock_aud_s0{ll}, tlock_nul_s0{10,tt,ss})
  
  cfg.design=[];
  cfg.design(:,1)=[ ones(1,size(tlock_tac_sall{10,tt,ss}.trial,1)) ones(1,size(tlock_aud_sall{30+ll}.trial,1)) 2*ones(1,size(tlock_aud_sall{ll}.trial,1)) 2*ones(1,size(tlock_nul_sall{10,tt,ss}.trial,1))]';
  stat2_sall{ll} = ft_timelockstatistics(cfg, tlock_tac_sall{10,tt,ss}, tlock_aud_sall{30+ll}, tlock_aud_sall{ll}, tlock_nul_sall{10,tt,ss})
end
% thus, stat2_* is with tac-alone at time zero, shifted aud, and AT with aud-first for 3, and tac-first for 7

save(['stat_erp_' sub{ii} '.mat'],'stat*')

figure;for ll=soalist,subplot(1,7,ll);imagesc(stat1_sall{ll}.time,1:62,stat1_sall{ll}.mask);end
figure;for ll=soalist,subplot(1,7,ll);imagesc(stat2_sall{ll}.time,1:62,stat2_sall{ll}.mask);end
figure;for ll=soalist,subplot(1,7,ll);imagesc(stat1_s0{ll}.time,1:62,stat1_s0{ll}.mask);end
figure;for ll=soalist,subplot(1,7,ll);imagesc(stat2_s0{ll}.time,1:62,stat2_s0{ll}.mask);end

cfg=[];
cfg.interactive='yes';
cfg.layout='EEG1010.lay';
ft_multiplotER(cfg,tlock_tpa_mtamn_sall{7});
ft_multiplotER(cfg,tlock_tpa_mtamn_sall{5});


clear *_s0
% end

%%  Stats on phase-sorting of unisensory ERP: attempt 1 (ignore now)

% related to phaset0=1 and phaset0use>0
plotflag=1;
printflag=0;

soalist=[1 3 4 5 6 7 9];
tt=3;
sleep=0;
tacaud=1;
if sleep
  ss=12;
  iter=11;
  iiuse=setdiff(iiBuse,1:7);
  trialkc=0;
else
  ss=10;
  iter=27;
  iiuse=setdiff(iiSuse,1:7);
  trialkc=-1;
end

subuseind=1;
for ii=iiuse
  % for ii=8:18
  cd([edir sub{ii}])
  load(['tlock_ERPphasebin_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat']);
  for bb=1:4
    tlockavg_tac5_phasebin{bb}{subuseind}=tlock_tactlock_phasebin{5,tt,ss,bb};
    tlockavg_aud5_phasebin{bb}{subuseind}=tlock_audtlock_phasebin{5,tt,ss,bb};
    tlockavg_nul5_phasebin{bb}{subuseind}=tlock_nultlock_phasebin{5,tt,ss,bb};
    tlockavg_ms15_phasebin{bb}{subuseind}=tlock_ms1tlock_phasebin{5,tt,ss,bb};
    
    tac5_dof(bb,subuseind)=tlockavg_tac5_phasebin{bb}{subuseind}.dof(18,1000);
    aud5_dof(bb,subuseind)=tlockavg_aud5_phasebin{bb}{subuseind}.dof(18,1000);
    nul5_dof(bb,subuseind)=tlockavg_nul5_phasebin{bb}{subuseind}.dof(18,1000);
    ms15_dof(bb,subuseind)=tlockavg_ms15_phasebin{bb}{subuseind}.dof(18,1000);
  end
  
  subuseind=subuseind+1;
  clear tlock*tlock_phasebin
end % ii

for bb=1:4,
  cfg=[];
  grave_tac_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_tac5_phasebin{bb}{:})
  grave_aud_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_aud5_phasebin{bb}{:})
  grave_nul_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_nul5_phasebin{bb}{:})
  grave_ms1_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_ms15_phasebin{bb}{:})
  
  cfg=[];
  cfg.keepindividual='yes';
  grind_tac_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_tac5_phasebin{bb}{:})
  grind_aud_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_aud5_phasebin{bb}{:})
  grind_nul_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_nul5_phasebin{bb}{:})
  grind_ms1_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_ms15_phasebin{bb}{:})
end

if plotflag
  if sleep
  else
    for bb=1:4
      topoplot_highlight(100+bb,grave_tac_pb{bb},[.0 .37],[]);
      if printflag
        %       print(111,['D:\audtac\figs\gravetacPaud_topoOverSTime_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
      end
    end
  end
  
  chanplot{1}={'Fz' 'FC1' 'FC2' 'F1' 'F2' 'C1' 'C2' 'Cz'};
  %   subplot(1,7,figind);
  cfg=[];
  cfg.channel=chanplot{1};
  cfg.xlim=[-0.5 1.1];
  cfg.ylim=[-7 7];
  figure;
  ft_singleplotER(cfg, grave_tac_pb{:})
  legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
  %   xlabel(['Tactile at time 0, ' sleepcond])
  %   ylabel(chanlabel{cc})
  title('Tactile Alone')
  figure;
  ft_singleplotER(cfg, grave_aud_pb{:})
  legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
  title('Auditory Alone')
  figure;
  ft_singleplotER(cfg, grave_nul_pb{:})
  legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
  title('Null Alone')
  figure;
  ft_singleplotER(cfg, grave_ms1_pb{:})
  legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
  title('Multisensory AT0')
  
  
end


load eeg1010_neighb
nsub=subuseind-1;
cfg=[];
cfg.latency=[.0 .4];
% cfg.channel=chanuse;
cfg.neighbours=neighbours;
cfg.parameter='individual';
cfg.method='montecarlo';
cfg.numrandomization=1000;
cfg.correctm='cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.statistic='depsamplesFunivariate';
cfg.design=zeros(2,4*nsub);
cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub) 3*ones(1,nsub) 4*ones(1,nsub)];
cfg.design(2,:)=[1:nsub 1:nsub 1:nsub 1:nsub];
cfg.ivar=1;
cfg.uvar=2;
cfg.tail=1;
stat_tacpb=ft_timelockstatistics(cfg, grind_tac_pb{:});
stat_audpb=ft_timelockstatistics(cfg, grind_aud_pb{:});
stat_nulpb=ft_timelockstatistics(cfg, grind_nul_pb{:});
stat_ms1pb=ft_timelockstatistics(cfg, grind_ms1_pb{:});

save([edir 'stat_pb_uni.mat'],'stat*','grave*','*dof');
save([edir 'grind_pb_uni.mat'],'grind*');

figure;imagesc(stat_nulpb.time,1:63,stat_nulpb.mask);title('Significance F-test mask: Null');xlabel('Time (s)');ylabel('Channel')
figure;imagesc(stat_nulpb.time,1:63,stat_tacpb.mask);title('Significance F-test mask: Tactile');xlabel('Time (s)');ylabel('Channel')
figure;imagesc(stat_nulpb.time,1:63,stat_audpb.mask);title('Significance F-test mask: Auditory');xlabel('Time (s)');ylabel('Channel')
figure;imagesc(stat_nulpb.time,1:63,stat_ms1pb.mask);title('Significance F-test mask: AT0');xlabel('Time (s)');ylabel('Channel')

% Nul is significant up to 120ms.   Do post-hoc t-test on anything showing
% diff after 120ms (namely Tac and AT0).
cfg=[];
cfg.latency=[.12 .4];
cfg.neighbours=neighbours;
cfg.parameter='individual';
cfg.method='montecarlo';
cfg.numrandomization=1000;
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
stat_tacpb_posthocT=ft_timelockstatistics(cfg, grind_tac_pb{1}, grind_tac_pb{4});
stat_ms1pb_posthocT=ft_timelockstatistics(cfg, grind_ms1_pb{1}, grind_ms1_pb{4});
% Tactile does show sig diff: (but not AT0)
figure;imagesc(stat_tacpb_posthocT.time,1:63,stat_tacpb_posthocT.mask);title('Significance T-test mask: Tactile');xlabel('Time (s)');ylabel('Channel')


% chanplot{2}=stat_tacpb.label(find(stat_tacpb.mask(:,145)));
%   cfg=[];
%   cfg.channel=chanplot{2};
%   cfg.xlim=[-0.5 1.1];
%   cfg.ylim=[-7 7];
%   figure;
%   ft_singleplotER(cfg, grave_tac_pb{:})
%   legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
% %   xlabel(['Tactile at time 0, ' sleepcond])
% %   ylabel(chanlabel{cc})
%   title('Tactile Alone')
%   figure;
%   ft_singleplotER(cfg, grave_aud_pb{:})
%   legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
%   title('Auditory Alone')
%   figure;
%   ft_singleplotER(cfg, grave_nul_pb{:})
%   legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
%   title('Null Alone')
%   figure;
%   ft_singleplotER(cfg, grave_ms1_pb{:})
%   legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
%   title('Multisensory AT0')
%
%   chanplot{3}=stat_tacpb.label(find(stat_ms1pb.mask(:,145)));
%   cfg=[];
%   cfg.channel=chanplot{3};
%   cfg.xlim=[-0.5 1.1];
%   cfg.ylim=[-7 7];
%   figure;
%   ft_singleplotER(cfg, grave_tac_pb{:})
%   legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
% %   xlabel(['Tactile at time 0, ' sleepcond])
% %   ylabel(chanlabel{cc})
%   title('Tactile Alone')
%   figure;
%   ft_singleplotER(cfg, grave_aud_pb{:})
%   legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
%   title('Auditory Alone')
%   figure;
%   ft_singleplotER(cfg, grave_nul_pb{:})
%   legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
%   title('Null Alone')
%   figure;
%   ft_singleplotER(cfg, grave_ms1_pb{:})
%   legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
%   title('Multisensory AT0')

%%  Stats on phase-sorting of unisensory and multisensry ERP versus Nul (use)

% related to phaset0=1 and phaset0use>0
plotflag=0;
statsflag=1;
soalist=[1 3 4 5 6 7 9];
tt=3;
sleep=1;
tacaud=1;
if sleep
  ss=12;
  iter=11;
  iiuse=setdiff(iiBuse,1:7);
  trialkc=-1;
else
  ss=10;
  iter=27;
  iiuse=setdiff(iiSuse,1:7);
  trialkc=-1;
end

subuseind=1;
for ii=iiuse
  % for ii=8:18
  cd([edir sub{ii}])
  load(['tlock_ERPphasebin_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat']);
  for bb=1:4
    
    tlockavg_tac5con1_phasebin{bb}{subuseind}=tlock_tac4mscon1_phasebin{5,3,ss,bb};
    tlockavg_aud5con1_phasebin{bb}{subuseind}=tlock_aud4mscon1_phasebin{5,3,ss,bb};
    tlockavg_nul5con1_phasebin{bb}{subuseind}=tlock_nul4mscon1_phasebin{5,3,ss,bb};
    
    if sleep
      tlockavg_tac5con1_absbin{bb}{subuseind}=tlock_tac4mscon1_absbin{5,3,ss,bb};
      tlockavg_aud5con1_absbin{bb}{subuseind}=tlock_aud4mscon1_absbin{5,3,ss,bb};
      tlockavg_nul5con1_absbin{bb}{subuseind}=tlock_nul4mscon1_absbin{5,3,ss,bb};
    end
    
    for ll=soalist
      tlockavg_msllcon1_phasebin{ll,bb}{subuseind}=tlock_ms4mscon1_phasebin{ll,3,ss,bb};
      tlockavg_nulllcon1_phasebin{ll,bb}{subuseind}=tlock_nul4mscon1_phasebin{ll,3,ss,bb};
      tlockavg_msllcon2_phasebin{ll,bb}{subuseind}=tlock_ms4mscon2_phasebin{ll,3,ss,bb};
      tlockavg_nulllcon2_phasebin{ll,bb}{subuseind}=tlock_nul4mscon2_phasebin{ll,3,ss,bb};
      if sleep
        tlockavg_msllcon1_absbin{ll,bb}{subuseind}=tlock_ms4mscon1_absbin{ll,3,ss,bb};
        tlockavg_nulllcon1_absbin{ll,bb}{subuseind}=tlock_nul4mscon1_absbin{ll,3,ss,bb};
        tlockavg_msllcon2_absbin{ll,bb}{subuseind}=tlock_ms4mscon2_absbin{ll,3,ss,bb};
        tlockavg_nulllcon2_absbin{ll,bb}{subuseind}=tlock_nul4mscon2_absbin{ll,3,ss,bb};
      end
    end % ll
    
  end
  
  subuseind=subuseind+1;
  clear tlock*tlock_phasebin
end % ii

for bb=1:4,
  cfg=[];
  grave_tac_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_tac5con1_phasebin{bb}{:});
  grave_aud_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_aud5con1_phasebin{bb}{:});
  grave_nul_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_nul5con1_phasebin{bb}{:});
  
  cfg=[];
  cfg.keepindividual='yes';
  grind_tac_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_tac5con1_phasebin{bb}{:});
  grind_aud_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_aud5con1_phasebin{bb}{:});
  grind_nul_pb{bb}=ft_timelockgrandaverage(cfg,tlockavg_nul5con1_phasebin{bb}{:});
  
  cfg=[];
  cfg.parameter='individual';
  cfg.operation='subtract';
  grind_tacMnul_pb{bb}=ft_math(cfg,grind_tac_pb{bb},grind_nul_pb{bb});
  grind_audMnul_pb{bb}=ft_math(cfg,grind_aud_pb{bb},grind_nul_pb{bb});
  
  grave_tacMnul_pb{bb}=grind_tacMnul_pb{bb};
  grave_tacMnul_pb{bb}.avg=squeeze(mean(grave_tacMnul_pb{bb}.individual,1));
  grave_tacMnul_pb{bb}=rmfield(grave_tacMnul_pb{bb},'individual');
  grave_tacMnul_pb{bb}.dimord='chan_time';
  
  grave_audMnul_pb{bb}=grind_audMnul_pb{bb};
  grave_audMnul_pb{bb}.avg=squeeze(mean(grave_audMnul_pb{bb}.individual,1));
  grave_audMnul_pb{bb}=rmfield(grave_audMnul_pb{bb},'individual');
  grave_audMnul_pb{bb}.dimord='chan_time';
  
  if sleep
    cfg=[];
    grave_tac_ab{bb}=ft_timelockgrandaverage(cfg,tlockavg_tac5con1_absbin{bb}{:});
    grave_aud_ab{bb}=ft_timelockgrandaverage(cfg,tlockavg_aud5con1_absbin{bb}{:});
    grave_nul_ab{bb}=ft_timelockgrandaverage(cfg,tlockavg_nul5con1_absbin{bb}{:});
    
    cfg=[];
    cfg.keepindividual='yes';
    grind_tac_ab{bb}=ft_timelockgrandaverage(cfg,tlockavg_tac5con1_absbin{bb}{:});
    grind_aud_ab{bb}=ft_timelockgrandaverage(cfg,tlockavg_aud5con1_absbin{bb}{:});
    grind_nul_ab{bb}=ft_timelockgrandaverage(cfg,tlockavg_nul5con1_absbin{bb}{:});
    
    cfg=[];
    cfg.parameter='individual';
    cfg.operation='subtract';
    grind_tacMnul_ab{bb}=ft_math(cfg,grind_tac_ab{bb},grind_nul_ab{bb});
    grind_audMnul_ab{bb}=ft_math(cfg,grind_aud_ab{bb},grind_nul_ab{bb});
    
    grave_tacMnul_ab{bb}=grind_tacMnul_ab{bb};
    grave_tacMnul_ab{bb}.avg=squeeze(mean(grave_tacMnul_ab{bb}.individual,1));
    grave_tacMnul_ab{bb}=rmfield(grave_tacMnul_ab{bb},'individual');
    grave_tacMnul_ab{bb}.dimord='chan_time';
    
    grave_audMnul_ab{bb}=grind_audMnul_ab{bb};
    grave_audMnul_ab{bb}.avg=squeeze(mean(grave_audMnul_ab{bb}.individual,1));
    grave_audMnul_ab{bb}=rmfield(grave_audMnul_ab{bb},'individual');
    grave_audMnul_ab{bb}.dimord='chan_time';
  end
  
  for ll=soalist
    cfg=[];
    grave_msllcon1_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_msllcon1_phasebin{ll,bb}{:});
    grave_nulllcon1_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_nulllcon1_phasebin{ll,bb}{:});
    grave_msllcon2_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_msllcon2_phasebin{ll,bb}{:});
    grave_nulllcon2_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_nulllcon2_phasebin{ll,bb}{:});
    cfg=[];
    cfg.keepindividual='yes';
    grind_msllcon1_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_msllcon1_phasebin{ll,bb}{:});
    grind_nulllcon1_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_nulllcon1_phasebin{ll,bb}{:});
    grind_msllcon2_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_msllcon2_phasebin{ll,bb}{:});
    grind_nulllcon2_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_nulllcon2_phasebin{ll,bb}{:});
    
    cfg=[];
    cfg.parameter='individual';
    cfg.operation='subtract';
    grind_msMnul1_pb{ll,bb}=ft_math(cfg,grind_msllcon1_pb{ll,bb},grind_nulllcon1_pb{ll,bb});
    grind_msMnul2_pb{ll,bb}=ft_math(cfg,grind_msllcon2_pb{ll,bb},grind_nulllcon2_pb{ll,bb});
    
    grave_msMnul1_pb{ll,bb}=grind_msMnul1_pb{ll,bb};
    grave_msMnul1_pb{ll,bb}.avg=squeeze(mean(grave_msMnul1_pb{ll,bb}.individual,1));
    grave_msMnul1_pb{ll,bb}=rmfield(grave_msMnul1_pb{ll,bb},'individual');
    grave_msMnul1_pb{ll,bb}.dimord='chan_time';
    
    grave_msMnul2_pb{ll,bb}=grind_msMnul2_pb{ll,bb};
    grave_msMnul2_pb{ll,bb}.avg=squeeze(mean(grave_msMnul2_pb{ll,bb}.individual,1));
    grave_msMnul2_pb{ll,bb}=rmfield(grave_msMnul2_pb{ll,bb},'individual');
    grave_msMnul2_pb{ll,bb}.dimord='chan_time';
    
    if sleep
      cfg=[];
      grave_msllcon1_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_msllcon1_absbin{ll,bb}{:});
      grave_nulllcon1_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_nulllcon1_absbin{ll,bb}{:});
      grave_msllcon2_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_msllcon2_absbin{ll,bb}{:});
      grave_nulllcon2_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_nulllcon2_absbin{ll,bb}{:});
      cfg=[];
      cfg.keepindividual='yes';
      grind_msllcon1_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_msllcon1_absbin{ll,bb}{:});
      grind_nulllcon1_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_nulllcon1_absbin{ll,bb}{:});
      grind_msllcon2_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_msllcon2_absbin{ll,bb}{:});
      grind_nulllcon2_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_nulllcon2_absbin{ll,bb}{:});
      
      cfg=[];
      cfg.parameter='individual';
      cfg.operation='subtract';
      grind_msMnul1_ab{ll,bb}=ft_math(cfg,grind_msllcon1_ab{ll,bb},grind_nulllcon1_ab{ll,bb});
      grind_msMnul2_ab{ll,bb}=ft_math(cfg,grind_msllcon2_ab{ll,bb},grind_nulllcon2_ab{ll,bb});
      
      grave_msMnul1_ab{ll,bb}=grind_msMnul1_ab{ll,bb};
      grave_msMnul1_ab{ll,bb}.avg=squeeze(mean(grave_msMnul1_ab{ll,bb}.individual,1));
      grave_msMnul1_ab{ll,bb}=rmfield(grave_msMnul1_ab{ll,bb},'individual');
      grave_msMnul1_ab{ll,bb}.dimord='chan_time';
      
      grave_msMnul2_ab{ll,bb}=grind_msMnul2_ab{ll,bb};
      grave_msMnul2_ab{ll,bb}.avg=squeeze(mean(grave_msMnul2_ab{ll,bb}.individual,1));
      grave_msMnul2_ab{ll,bb}=rmfield(grave_msMnul2_ab{ll,bb},'individual');
      grave_msMnul2_ab{ll,bb}.dimord='chan_time';
    end
  end % ll
end % bb
% save([edir 'grind_pb_uninul.mat'],'grind*')
save([edir 'grind_pb_uninul_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'grind*');

if statsflag
  load eeg1010_neighb
  nsub=subuseind-1;
  cfg=[];
  cfg.neighbours=neighbours;
  cfg.parameter='individual';
  cfg.method='montecarlo';
  cfg.numrandomization=1000;
  cfg.correctm='cluster';
  cfg.clusteralpha = 0.05;
  cfg.clusterstatistic = 'maxsum';
  cfg.minnbchan = 2;
  cfg.ivar=1;
  cfg.uvar=2;
  cfg.statistic='depsamplesFunivariate';
  cfg.design=zeros(2,4*nsub);
  cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub) 3*ones(1,nsub) 4*ones(1,nsub)];
  cfg.design(2,:)=[1:nsub 1:nsub 1:nsub 1:nsub];
  cfg.tail=1;
  
  if sleep
    cfg.latency=[.1 .75];
  else
    cfg.latency=[.1 .45];
  end
  
  stat_tacMnul_pb=ft_timelockstatistics(cfg, grind_tacMnul_pb{:});
  stat_audMnul_pb=ft_timelockstatistics(cfg, grind_audMnul_pb{:});
  
  if sleep
    stat_tacMnul_ab=ft_timelockstatistics(cfg, grind_tacMnul_ab{:});
    stat_audMnul_ab=ft_timelockstatistics(cfg, grind_audMnul_ab{:});
  end
  
  for ll=soalist
    if sleep
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[.1 .75];
      elseif ll==6
        cfg.latency=[.12 .77];
      elseif ll==7
        cfg.latency=[.17 .82];
      elseif ll==9
        cfg.latency=[.6 1.25];
      end
    else
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[.1 .45];
      elseif ll==6
        cfg.latency=[.12 .47];
      elseif ll==7
        cfg.latency=[.17 .52];
      elseif ll==9
        cfg.latency=[.6 .95];
      end
    end
    stat_msMnul1_pb{ll}=ft_timelockstatistics(cfg, grind_msMnul1_pb{ll,:});
    stat_msMnul2_pb{ll}=ft_timelockstatistics(cfg, grind_msMnul2_pb{ll,:});
    if sleep
      stat_msMnul1_ab{ll}=ft_timelockstatistics(cfg, grind_msMnul1_ab{ll,:});
      stat_msMnul2_ab{ll}=ft_timelockstatistics(cfg, grind_msMnul2_ab{ll,:});
    end
    save([edir 'stat_pb_uninul_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'stat*')
  end % ll
  
end % statsflag




if plotflag
  if ~exist('grave_tacMnul_pb','var')
    if ~exist('grind_tacMnul_pb','var')
      load([edir 'grind_pb_uninul_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat']);
    end
    
    for bb=1:4
      
      grave_tacMnul_pb{bb}=grind_tacMnul_pb{bb};
      grave_tacMnul_pb{bb}.avg=squeeze(mean(grave_tacMnul_pb{bb}.individual,1))
      grave_tacMnul_pb{bb}=rmfield(grave_tacMnul_pb{bb},'individual');
      grave_tacMnul_pb{bb}.dimord='chan_time';
      
      grave_audMnul_pb{bb}=grind_audMnul_pb{bb};
      grave_audMnul_pb{bb}.avg=squeeze(mean(grave_audMnul_pb{bb}.individual,1))
      grave_audMnul_pb{bb}=rmfield(grave_audMnul_pb{bb},'individual');
      grave_audMnul_pb{bb}.dimord='chan_time';
      
      if sleep
        grave_tacMnul_ab{bb}=grind_tacMnul_ab{bb};
        grave_tacMnul_ab{bb}.avg=squeeze(mean(grave_tacMnul_ab{bb}.individual,1))
        grave_tacMnul_ab{bb}=rmfield(grave_tacMnul_ab{bb},'individual');
        grave_tacMnul_ab{bb}.dimord='chan_time';
        
        grave_audMnul_ab{bb}=grind_audMnul_ab{bb};
        grave_audMnul_ab{bb}.avg=squeeze(mean(grave_audMnul_ab{bb}.individual,1))
        grave_audMnul_ab{bb}=rmfield(grave_audMnul_ab{bb},'individual');
        grave_audMnul_ab{bb}.dimord='chan_time';
        
      end
      
      for ll=soalist
        grave_msMnul1_pb{ll,bb}=grind_msMnul1_pb{ll,bb};
        grave_msMnul1_pb{ll,bb}.avg=squeeze(mean(grave_msMnul1_pb{ll,bb}.individual,1))
        grave_msMnul1_pb{ll,bb}=rmfield(grave_msMnul1_pb{ll,bb},'individual');
        grave_msMnul1_pb{ll,bb}.dimord='chan_time';
        
        grave_msMnul2_pb{ll,bb}=grind_msMnul2_pb{ll,bb};
        grave_msMnul2_pb{ll,bb}.avg=squeeze(mean(grave_msMnul2_pb{ll,bb}.individual,1))
        grave_msMnul2_pb{ll,bb}=rmfield(grave_msMnul2_pb{ll,bb},'individual');
        grave_msMnul2_pb{ll,bb}.dimord='chan_time';
        
        if sleep
          grave_msMnul1_ab{ll,bb}=grind_msMnul1_ab{ll,bb};
          grave_msMnul1_ab{ll,bb}.avg=squeeze(mean(grave_msMnul1_ab{ll,bb}.individual,1))
          grave_msMnul1_ab{ll,bb}=rmfield(grave_msMnul1_ab{ll,bb},'individual');
          grave_msMnul1_ab{ll,bb}.dimord='chan_time';
          
          grave_msMnul2_ab{ll,bb}=grind_msMnul2_ab{ll,bb};
          grave_msMnul2_ab{ll,bb}.avg=squeeze(mean(grave_msMnul2_ab{ll,bb}.individual,1))
          grave_msMnul2_ab{ll,bb}=rmfield(grave_msMnul2_ab{ll,bb},'individual');
          grave_msMnul2_ab{ll,bb}.dimord='chan_time';
        end
      end
      
    end  % bb
  end
  
  chanplot{1}={'Fz' 'Cz' 'F1' 'F2' 'FC1' 'FC2' 'C1' 'C2'}; % frontocentral
  chanplot{2}={'CP5' 'POz' 'Pz' 'P3' 'P4' 'C4' 'O1' 'O2' 'P7' 'PO7'}; % occipital
  
  
  cfg=[];
  cfg.channel=chanplot{1};
  cfg.xlim=[-0.5 1.1];
  cfg.ylim=[-7 7];
  figure(2);
  ft_singleplotER(cfg, grave_tacMnul_pb{:})
  legend({'Peak' 'P to T' 'Trough' 'T to P'})
  title('T-N')
  
  cfg=[];
  cfg.channel=chanplot{1};
  cfg.xlim=[-0.5 1.1];
  cfg.ylim=[-7 7];
  figure(8);
  ft_singleplotER(cfg, grave_audMnul_pb{:})
  legend({'Peak' 'P to T' 'Trough' 'T to P'})
  title('A-N')
  
  if sleep
    cfg=[];
    cfg.channel=chanplot{1};
    cfg.xlim=[-0.5 1.1];
    cfg.ylim=[-7 7];
    figure(12);
    ft_singleplotER(cfg, grave_tacMnul_ab{:})
    legend({'Low' 'Mid-low' 'Mid-high' 'High'})
    title('T-N')
    
    cfg=[];
    cfg.channel=chanplot{1};
    cfg.xlim=[-0.5 1.1];
    cfg.ylim=[-7 7];
    figure(18);
    ft_singleplotER(cfg, grave_audMnul_ab{:})
    legend({'Low' 'Mid-low' 'Mid-high' 'High'})
    title('A-N')
  end
  
  for ll=soalist
    
    cfg=[];
    cfg.channel=chanplot{1};
    cfg.xlim=[-0.5 1.1];
    cfg.ylim=[-7 7];
    figure(ll);
    ft_singleplotER(cfg, grave_msMnul1_pb{ll,:})
    legend({'Peak' 'P to T' 'Trough' 'T to P'})
    title('AT*-N')
    
    %     cfg=[];
    %     cfg.channel=chanplot{1};
    %     cfg.xlim=[-0.5 1.1];
    %     cfg.ylim=[-7 7];
    %     figure(100+ll);
    %     ft_singleplotER(cfg, grave_msMnul2_pb{ll,:})
    %     legend({'Peak' 'P to T' 'Trough' 'T to P'})
    %     title('AT*-N')
    
    if sleep
      cfg=[];
      cfg.channel=chanplot{1};
      cfg.xlim=[-0.5 1.1];
      cfg.ylim=[-7 7];
      figure(10+ll);
      ft_singleplotER(cfg, grave_msMnul1_ab{ll,:})
      legend({'Low' 'Mid-low' 'Mid-high' 'High'})
      title('AT*-N')
      
      %       cfg=[];
      %       cfg.channel=chanplot{1};
      %       cfg.xlim=[-0.5 1.1];
      %       cfg.ylim=[-7 7];
      %       figure(100+10+ll);
      %       ft_singleplotER(cfg, grave_msMnul2_ab{ll,:})
      %       legend({'Low' 'Mid-low' 'Mid-high' 'High'})
      %       title('AT*-N')
    end
  end
  
end

% followup stats
if statsflag
  load eeg1010_neighb
  nsub=subuseind-1;
  cfg=[];
  cfg.neighbours=neighbours;
  cfg.parameter='individual';
  cfg.method='montecarlo';
  cfg.numrandomization=1000;
  cfg.correctm='cluster';
  cfg.clusteralpha = 0.05;
  cfg.clusterstatistic = 'maxsum';
  cfg.minnbchan = 2;
  cfg.ivar=1;
  cfg.uvar=2;
  cfg.statistic='depsamplesFunivariate';
  cfg.design=zeros(2,4*nsub);
  cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub) 3*ones(1,nsub) 4*ones(1,nsub)];
  cfg.design(2,:)=[1:nsub 1:nsub 1:nsub 1:nsub];
  cfg.tail=1;
  
  if sleep
    cfg.latency=[.1 .75];
  else
    cfg.latency=[.1 .45];
  end
  
  stat_tacMnul_pb=ft_timelockstatistics(cfg, grind_tacMnul_pb{:});
  stat_audMnul_pb=ft_timelockstatistics(cfg, grind_audMnul_pb{:});
  
  if sleep
    stat_tacMnul_ab=ft_timelockstatistics(cfg, grind_tacMnul_ab{:});
    stat_audMnul_ab=ft_timelockstatistics(cfg, grind_audMnul_ab{:});
  end
  
  for ll=soalist
    if sleep
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[.1 .75];
      elseif ll==6
        cfg.latency=[.12 .77];
      elseif ll==7
        cfg.latency=[.17 .82];
      elseif ll==9
        cfg.latency=[.6 1.25];
      end
    else
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[.1 .45];
      elseif ll==6
        cfg.latency=[.12 .47];
      elseif ll==7
        cfg.latency=[.17 .52];
      elseif ll==9
        cfg.latency=[.6 .95];
      end
    end
    stat_msMnul1_pb{ll}=ft_timelockstatistics(cfg, grind_msMnul1_pb{ll,:});
    stat_msMnul2_pb{ll}=ft_timelockstatistics(cfg, grind_msMnul2_pb{ll,:});
    if sleep
      stat_msMnul1_ab{ll}=ft_timelockstatistics(cfg, grind_msMnul1_ab{ll,:});
      stat_msMnul2_ab{ll}=ft_timelockstatistics(cfg, grind_msMnul2_ab{ll,:});
    end
    save([edir 'stat_pb_uninul_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'stat*')
  end % ll
  
end % statsflag






%% Stats on phase-sorting MS ERP

plotflag=1;
printflag=1;
statsflag=0;
stats1flag=0;
sleepshortstatflag=1;

soalist=[1 3 4 5 6 7 9];
tt=3;
sleep=1;
tacaud=1;
if sleep
  ss=12;
  iter=11;
  iiuse=setdiff(iiBuse,1:7);
  trialkc=0;
  if trialkc==0
    iiuse=setdiff(iiuse,[18 24 ]);
  end
else
  ss=10;
  iter=27;
  iiuse=setdiff(iiSuse,1:7);
  trialkc=-1;
end

subuseind=1;
for ii=iiuse
  cd([edir sub{ii}])
  load(['tlock_ERPphasebin_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat']);
  for ll=soalist
    for bb=1:4
      tlockavg_tacPaud1_phasebin{ll,bb}{subuseind}=tlock_tacPaud_4mscon1_phasebin{ll,tt,ss,bb};
      tlockavg_tacMSpN1_phasebin{ll,bb}{subuseind}=tlock_tacMSpN_4mscon1_phasebin{ll,tt,ss,bb};
      if sleep
        tlockavg_tacPaud1_absbin{ll,bb}{subuseind}=tlock_tacPaud_4mscon1_absbin{ll,tt,ss,bb};
        tlockavg_tacMSpN1_absbin{ll,bb}{subuseind}=tlock_tacMSpN_4mscon1_absbin{ll,tt,ss,bb};
      else
        tlockavg_tacPaud1IAF_phasebin{ll,bb}{subuseind}=tlock_tacPaud_4mscon1IAF_phasebin{ll,tt,ss,bb};
        tlockavg_tacMSpN1IAF_phasebin{ll,bb}{subuseind}=tlock_tacMSpN_4mscon1IAF_phasebin{ll,tt,ss,bb};
      end
    end % bb
  end % ll
  freqsubfinal(subuseind)=freqsub(ii);
  try freqsleepfinal(subuseind)=freqsleep(ii); catch end
  
  subuseind=subuseind+1;
  clear tlock*tlock_phasebin
end  % ii
save([edir 'freqfinal_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'freq*final');

% subuseind=1;
% for ii=iiuse
%   cd([edir sub{ii}])
%   load(['tlock_ERPphasebin_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'freqsub');
%   freqsubfinal(subuseind)=freqsub(ii);
%
%   subuseind=subuseind+1;
%   clear tlock*tlock_phasebin
% end  % ii
figure(600);hist(freqsubfinal,7);

for ll=soalist
  for bb=1:4,
    %     cfg=[];
    %     cfg.parameter='avg';
    %     cfg.operation='subtract';
    %     tlockavg_TPAmMSPN1_pb{ll,bb}=ft_math(cfg,tlockavg_tacPaud1_phasebin{ll,bb}{:},tlockavg_tacMSpN1_phasebin{ll,bb}{:});
    
    cfg=[];
    grave_tacPaud1_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacPaud1_phasebin{ll,bb}{:});
    grave_tacMSpN1_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacMSpN1_phasebin{ll,bb}{:});
    if sleep
      grave_tacPaud1_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacPaud1_absbin{ll,bb}{:});
      grave_tacMSpN1_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacMSpN1_absbin{ll,bb}{:});
    else
      grave_tacPaud1IAF_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacPaud1IAF_phasebin{ll,bb}{:});
      grave_tacMSpN1IAF_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacMSpN1IAF_phasebin{ll,bb}{:});
    end
    %     grave_TPAmMSPN1_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_TPAmMSPN1_phasebin{ll,bb}{:});
    
    cfg=[];
    cfg.keepindividual='yes';
    grind_tacPaud1_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacPaud1_phasebin{ll,bb}{:});
    grind_tacMSpN1_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacMSpN1_phasebin{ll,bb}{:});
    if sleep
      grind_tacPaud1_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacPaud1_absbin{ll,bb}{:});
      grind_tacMSpN1_ab{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacMSpN1_absbin{ll,bb}{:});
    else
      grind_tacPaud1IAF_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacPaud1IAF_phasebin{ll,bb}{:});
      grind_tacMSpN1IAF_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_tacMSpN1IAF_phasebin{ll,bb}{:});
    end
    %     grind_TPAmMSPN1_pb{ll,bb}=ft_timelockgrandaverage(cfg,tlockavg_TPAmMSPN1_phasebin{ll,bb}{:});
    
    cfg=[];
    cfg.parameter='individual';
    cfg.operation='subtract';
    grind_TPAmMSPN1_pb{ll,bb}=ft_math(cfg,grind_tacPaud1_pb{ll,bb},grind_tacMSpN1_pb{ll,bb})
    if sleep
      grind_TPAmMSPN1_ab{ll,bb}=ft_math(cfg,grind_tacPaud1_ab{ll,bb},grind_tacMSpN1_ab{ll,bb})
    else
      grind_TPAmMSPN1IAF_pb{ll,bb}=ft_math(cfg,grind_tacPaud1IAF_pb{ll,bb},grind_tacMSpN1IAF_pb{ll,bb})
    end
    
    grave_TPAmMSPN1_pb{ll,bb}=grind_TPAmMSPN1_pb{ll,bb};
    grave_TPAmMSPN1_pb{ll,bb}.avg=squeeze(mean(grave_TPAmMSPN1_pb{ll,bb}.individual,1))
    grave_TPAmMSPN1_pb{ll,bb}=rmfield(grave_TPAmMSPN1_pb{ll,bb},'individual');
    grave_TPAmMSPN1_pb{ll,bb}.dimord='chan_time';
    
    if sleep
      grave_TPAmMSPN1_ab{ll,bb}=grind_TPAmMSPN1_ab{ll,bb};
      grave_TPAmMSPN1_ab{ll,bb}.avg=squeeze(mean(grave_TPAmMSPN1_ab{ll,bb}.individual,1))
      grave_TPAmMSPN1_ab{ll,bb}=rmfield(grave_TPAmMSPN1_ab{ll,bb},'individual');
      grave_TPAmMSPN1_ab{ll,bb}.dimord='chan_time';
    else
      grave_TPAmMSPN1IAF_pb{ll,bb}=grind_TPAmMSPN1IAF_pb{ll,bb};
      grave_TPAmMSPN1IAF_pb{ll,bb}.avg=squeeze(mean(grave_TPAmMSPN1IAF_pb{ll,bb}.individual,1))
      grave_TPAmMSPN1IAF_pb{ll,bb}=rmfield(grave_TPAmMSPN1IAF_pb{ll,bb},'individual');
      grave_TPAmMSPN1IAF_pb{ll,bb}.dimord='chan_time';
    end
  end
  cfg=[];
  cfg.parameter='individual';
  cfg.operation='subtract';
  grind_TPAmMSPN_4m2_pb{ll}=ft_math(cfg,grind_TPAmMSPN1_pb{ll,4},grind_TPAmMSPN1_pb{ll,2});
  grave_TPAmMSPN_4m2_pb{ll}=grind_TPAmMSPN_4m2_pb{ll};
  grave_TPAmMSPN_4m2_pb{ll}.avg=squeeze(mean(grave_TPAmMSPN_4m2_pb{ll}.individual,1));
  grave_TPAmMSPN_4m2_pb{ll}=rmfield(grave_TPAmMSPN_4m2_pb{ll},'individual');
  grave_TPAmMSPN_4m2_pb{ll}.dimord='chan_time';
  
  grind_TPAmMSPN_3m1_pb{ll}=ft_math(cfg,grind_TPAmMSPN1_pb{ll,3},grind_TPAmMSPN1_pb{ll,1});
  grave_TPAmMSPN_3m1_pb{ll}=grind_TPAmMSPN_3m1_pb{ll};
  grave_TPAmMSPN_3m1_pb{ll}.avg=squeeze(mean(grave_TPAmMSPN_3m1_pb{ll}.individual,1));
  grave_TPAmMSPN_3m1_pb{ll}=rmfield(grave_TPAmMSPN_3m1_pb{ll},'individual');
  grave_TPAmMSPN_3m1_pb{ll}.dimord='chan_time';
  
  if sleep
    grind_TPAmMSPN_4m1_ab{ll}=ft_math(cfg,grind_TPAmMSPN1_ab{ll,4},grind_TPAmMSPN1_ab{ll,1});
    grave_TPAmMSPN_4m1_ab{ll}=grind_TPAmMSPN_4m1_ab{ll};
    grave_TPAmMSPN_4m1_ab{ll}.avg=squeeze(mean(grave_TPAmMSPN_4m1_ab{ll}.individual,1));
    grave_TPAmMSPN_4m1_ab{ll}=rmfield(grave_TPAmMSPN_4m1_ab{ll},'individual');
    grave_TPAmMSPN_4m1_ab{ll}.dimord='chan_time';
  else
    grind_TPAmMSPNIAF_4m2_pb{ll}=ft_math(cfg,grind_TPAmMSPN1IAF_pb{ll,4},grind_TPAmMSPN1IAF_pb{ll,2});
    grave_TPAmMSPNIAF_4m2_pb{ll}=grind_TPAmMSPNIAF_4m2_pb{ll};
    grave_TPAmMSPNIAF_4m2_pb{ll}.avg=squeeze(mean(grave_TPAmMSPNIAF_4m2_pb{ll}.individual,1));
    grave_TPAmMSPNIAF_4m2_pb{ll}=rmfield(grave_TPAmMSPNIAF_4m2_pb{ll},'individual');
    grave_TPAmMSPNIAF_4m2_pb{ll}.dimord='chan_time';
    
    grind_TPAmMSPNIAF_3m1_pb{ll}=ft_math(cfg,grind_TPAmMSPN1IAF_pb{ll,3},grind_TPAmMSPN1IAF_pb{ll,1});
    grave_TPAmMSPNIAF_3m1_pb{ll}=grind_TPAmMSPNIAF_3m1_pb{ll};
    grave_TPAmMSPNIAF_3m1_pb{ll}.avg=squeeze(mean(grave_TPAmMSPNIAF_3m1_pb{ll}.individual,1));
    grave_TPAmMSPNIAF_3m1_pb{ll}=rmfield(grave_TPAmMSPNIAF_3m1_pb{ll},'individual');
    grave_TPAmMSPNIAF_3m1_pb{ll}.dimord='chan_time';
  end
end % ll


if plotflag
  for ll=soalist
    
    %     if sleep
    %     else
    %       for bb=1:4
    %         topoplot_highlight(100+bb,grave_tacMSpN1_pb{ll,bb},[.0 .37],[]);
    %         if printflag
    %           %       print(111,['D:\audtac\figs\gravetacPaud_topoOverSTime_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep)  '.png'],'-dpng')
    %         end
    %       end
    %     end
    
    chanplot{1}={'Fz' 'FC1' 'FC2' 'F1' 'F2' 'C1' 'C2' 'Cz'};
    %   subplot(1,7,figind);
    cfg=[];
    cfg.channel=chanplot{1};
    cfg.xlim=[-0.5 1.1];
    cfg.ylim=[-7 7];
    figure(ll);
    ft_singleplotER(cfg, grave_tacPaud1_pb{ll,:})
    legend({'Peak' 'P to T' 'Trough' 'T to P'})
    %     legend({'Trough' 'Low-mid phase' 'High-mid phase' 'Peak'})
    %   xlabel(['Tactile at time 0, ' sleepcond])
    %   ylabel(chanlabel{cc})
    title('T + A')
    
    figure(10+ll);
    ft_singleplotER(cfg, grave_tacMSpN1_pb{ll,:})
    legend({'Peak' 'P to T' 'Trough' 'T to P'})
    title('MS + N')
    
    figure(20+ll);
    ft_singleplotER(cfg, grave_TPAmMSPN1_pb{ll,:})
    legend({'Peak' 'P to T' 'Trough' 'T to P'})
    title('(T + A)- (MS + N)')
    
    figure(30+ll);
    ft_singleplotER(cfg, grave_TPAmMSPN_4m2_pb{ll}, grave_TPAmMSPN_3m1_pb{ll})
    title('(T + A)- (MS + N)')
    legend({'TtoP - PtoT' 'T - P'})
    
    if sleep
      figure(100+ll);
      ft_singleplotER(cfg, grave_tacPaud1_ab{ll,:})
      legend({'Low' 'Mid-low' 'Mid-high' 'High'})
      title('T + A')
      
      figure(100+10+ll);
      ft_singleplotER(cfg, grave_tacMSpN1_ab{ll,:})
      legend({'Low' 'Mid-low' 'Mid-high' 'High'})
      title('MS + N')
      
      figure(100+20+ll);
      ft_singleplotER(cfg, grave_TPAmMSPN1_ab{ll,:})
      legend({'Low' 'Mid-low' 'Mid-high' 'High'})
      title('(T + A)- (MS + N)')
      
      figure(100+30+ll);
      ft_singleplotER(cfg, grave_TPAmMSPN_4m1_ab{ll})
      title('(T + A)- (MS + N)')
      legend({'High vs Low'})
    else
      figure(100+ll);
      ft_singleplotER(cfg, grave_tacPaud1IAF_pb{ll,:})
      legend({'Peak' 'P to T' 'Trough' 'T to P'})
      title('T + A')
      
      figure(100+10+ll);
      ft_singleplotER(cfg, grave_tacMSpN1IAF_pb{ll,:})
      legend({'Peak' 'P to T' 'Trough' 'T to P'})
      title('MS + N')
      
      figure(100+20+ll);
      ft_singleplotER(cfg, grave_TPAmMSPN1IAF_pb{ll,:})
      legend({'Peak' 'P to T' 'Trough' 'T to P'})
      title('(T + A)- (MS + N)')
      
      figure(100+30+ll);
      ft_singleplotER(cfg, grave_TPAmMSPNIAF_4m2_pb{ll}, grave_TPAmMSPNIAF_3m1_pb{ll})
      title('(T + A)- (MS + N)')
      legend({'TtoP - PtoT' 'T - P'})
    end
    
    
  end
end


if statsflag
  load eeg1010_neighb
  nsub=subuseind-1;
  cfg=[];
  % cfg.channel=chanuse;
  cfg.neighbours=neighbours;
  cfg.parameter='individual';
  cfg.method='montecarlo';
  cfg.numrandomization=1000;
  cfg.correctm='cluster';
  cfg.clusteralpha = 0.05;
  cfg.clusterstatistic = 'maxsum';
  cfg.minnbchan = 2;
  cfg.ivar=1;
  cfg.uvar=2;
  cfg.statistic='depsamplesFunivariate';
  cfg.design=zeros(2,4*nsub);
  cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub) 3*ones(1,nsub) 4*ones(1,nsub)];
  cfg.design(2,:)=[1:nsub 1:nsub 1:nsub 1:nsub];
  cfg.tail=1;
  for ll=soalist
    if sleep && ~sleepshortstatflag
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[0 1];
      elseif ll==6
        cfg.latency=[0 1]+.02;
      elseif ll==7
        cfg.latency=[0 1]+.07;
      elseif ll==9
        cfg.latency=[0 1]+.5;
      end
    elseif sleep && sleepshortstatflag
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[0 .45];
      elseif ll==6
        cfg.latency=[0 .45]+.02;
      elseif ll==7
        cfg.latency=[0 .45]+.07;
      elseif ll==9
        cfg.latency=[0 .45]+.5;
      end
    else
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[.1 .45];
      elseif ll==6
        cfg.latency=[.12 .47];
      elseif ll==7
        cfg.latency=[.17 .52];
      elseif ll==9
        cfg.latency=[.6 .95];
      end
    end
    
    % these test F-test (1-way ANOVA) for any difference of the 4 bins against each other.
    stat_tacPaud1_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1_pb{ll,:});
    %   stat_tacPaud2_pb=ft_timelockstatistics(cfg, grind_tacPaud2_pb{ll,:});
    stat_tacMSpN1_pb{ll}=ft_timelockstatistics(cfg, grind_tacMSpN1_pb{ll,:});
    %   stat_tacMSpN2_pb=ft_timelockstatistics(cfg, grind_tacMSpN2_pb{ll,:});
    stat_TPAmMSPN1_pb{ll}=ft_timelockstatistics(cfg, grind_TPAmMSPN1_pb{ll,:});
    
    if sleep
      stat_tacPaud1_ab{ll}=ft_timelockstatistics(cfg, grind_tacPaud1_ab{ll,:});
      stat_tacMSpN1_ab{ll}=ft_timelockstatistics(cfg, grind_tacMSpN1_ab{ll,:});
      stat_TPAmMSPN1_ab{ll}=ft_timelockstatistics(cfg, grind_TPAmMSPN1_ab{ll,:});
    else
      stat_tacPaud1IAF_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1IAF_pb{ll,:});
      %   stat_tacPaud2_pb=ft_timelockstatistics(cfg, grind_tacPaud2_pb{ll,:});
      stat_tacMSpN1IAF_pb{ll}=ft_timelockstatistics(cfg, grind_tacMSpN1IAF_pb{ll,:});
      %   stat_tacMSpN2_pb=ft_timelockstatistics(cfg, grind_tacMSpN2_pb{ll,:});
      stat_TPAmMSPN1IAF_pb{ll}=ft_timelockstatistics(cfg, grind_TPAmMSPN1IAF_pb{ll,:});
    end
    %   save([edir 'stat_pb_mult.mat'],'stat*')
    %     save([edir 'stat_pb_mult_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'stat*')
    save([edir 'stat_pb_mult_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag) '.mat'],'stat*');
  end
  
  cfg.statistic='depsamplesT';
  cfg.design=zeros(2,2*nsub);
  cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
  cfg.design(2,:)=[1:nsub 1:nsub];
  cfg.tail=0;
  for ll=soalist
    if sleep && ~sleepshortstatflag
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[0 1];
      elseif ll==6
        cfg.latency=[0 1]+.02;
      elseif ll==7
        cfg.latency=[0 1]+.07;
      elseif ll==9
        cfg.latency=[0 1]+.5;
      end
    elseif sleep && sleepshortstatflag
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[0 .45];
      elseif ll==6
        cfg.latency=[0 .45]+.02;
      elseif ll==7
        cfg.latency=[0 .45]+.07;
      elseif ll==9
        cfg.latency=[0 .45]+.5;
      end
    else
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[.1 .45];
      elseif ll==6
        cfg.latency=[.12 .47];
      elseif ll==7
        cfg.latency=[.17 .52];
      elseif ll==9
        cfg.latency=[.6 .95];
      end
    end
    % these test t-test for any difference of the A+T vs MS+N for each bin separately
    stat_TPAmMSPN1_peak_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1_pb{ll,1}, grind_tacMSpN1_pb{ll,1});
    stat_TPAmMSPN1_ptot_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1_pb{ll,2}, grind_tacMSpN1_pb{ll,2});
    stat_TPAmMSPN1_trgh_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1_pb{ll,3}, grind_tacMSpN1_pb{ll,3});
    stat_TPAmMSPN1_ttop_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1_pb{ll,4}, grind_tacMSpN1_pb{ll,4});
    
    if sleep
      stat_TPAmMSPN1_low_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1_ab{ll,1}, grind_tacMSpN1_ab{ll,1});
      stat_TPAmMSPN1_midl_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1_ab{ll,2}, grind_tacMSpN1_ab{ll,2});
      stat_TPAmMSPN1_midh_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1_ab{ll,3}, grind_tacMSpN1_ab{ll,3});
      stat_TPAmMSPN1_high_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1_ab{ll,4}, grind_tacMSpN1_ab{ll,4});
    else
      stat_TPAmMSPN1IAF_peak_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1IAF_pb{ll,1}, grind_tacMSpN1IAF_pb{ll,1});
      stat_TPAmMSPN1IAF_ptot_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1IAF_pb{ll,2}, grind_tacMSpN1IAF_pb{ll,2});
      stat_TPAmMSPN1IAF_trgh_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1IAF_pb{ll,3}, grind_tacMSpN1IAF_pb{ll,3});
      stat_TPAmMSPN1IAF_ttop_pb{ll}=ft_timelockstatistics(cfg, grind_tacPaud1IAF_pb{ll,4}, grind_tacMSpN1IAF_pb{ll,4});
    end
  end
  
  %   save([edir 'stat_pb_mult.mat'],'stat*')
  %   save([edir 'stat_pb_mult_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'stat*')
  save([edir 'stat_pb_mult_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag) '.mat'],'stat*');
  clear grind*IAF* grind_*_*m*
  %   save([edir 'grind_pb_mult.mat'],'grind*')
  save([edir 'grind_pb_mult_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'grind*')
end

if stats1flag  % specific opposite phase or power contrasts (rather than full ANOVA)
  load eeg1010_neighb
  nsub=size(grind_TPAmMSPN1_pb{1,4}.individual,1);
  cfg=[];
  cfg.neighbours=neighbours;
  cfg.parameter='individual';
  cfg.method='montecarlo';
  cfg.numrandomization=1000;
  cfg.correctm='cluster';
  cfg.clusteralpha = 0.05;
  cfg.clusterstatistic = 'maxsum';
  cfg.minnbchan = 2;
  cfg.ivar=1;
  cfg.uvar=2;
  cfg.statistic='depsamplesT';
  cfg.design=zeros(2,2*nsub);
  cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
  cfg.design(2,:)=[1:nsub 1:nsub];
  cfg.tail=0;
  for ll=soalist
    if sleep==0
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[.1 .45];
      elseif ll==6
        cfg.latency=[.12 .47];
      elseif ll==7
        cfg.latency=[.17 .52];
      elseif ll==9
        cfg.latency=[.6 .95];
      end
    elseif sleep==1 && ~sleepshortstatflag
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[0 1];
      elseif ll==6
        cfg.latency=[0 1]+.02;
      elseif ll==7
        cfg.latency=[0 1]+.07;
      elseif ll==9
        cfg.latency=[0 1]+.5;
      end
    elseif sleep==1 && sleepshortstatflag
      if ll==1 || ll==3 || ll==4 || ll==5
        cfg.latency=[0 .45];
      elseif ll==6
        cfg.latency=[0 .45]+.02;
      elseif ll==7
        cfg.latency=[0 .45]+.07;
      elseif ll==9
        cfg.latency=[0 .45]+.5;
      end
    end
    stat_TPAmMSPN_ptotMttop_pb{ll}=ft_timelockstatistics(cfg, grind_TPAmMSPN1_pb{ll,2}, grind_TPAmMSPN1_pb{ll,4});
    stat_TPAmMSPN_peakMtrgh_pb{ll}=ft_timelockstatistics(cfg, grind_TPAmMSPN1_pb{ll,1}, grind_TPAmMSPN1_pb{ll,3});
    
    if sleep
      stat_TPAmMSPN_highMlow_ab{ll}=ft_timelockstatistics(cfg, grind_TPAmMSPN1_ab{ll,1}, grind_TPAmMSPN1_ab{ll,4});
    end
    
    if sleep==0 && (ll==3 || ll==5 || ll==7)  % conditions where previously there was significance
      aa=nan(2,2);
      try  aa(1,:)=[min(stat_TPAmMSPN1_peak_pb{ll}.time(find(mean(stat_TPAmMSPN1_peak_pb{ll}.mask,1)))) max(stat_TPAmMSPN1_peak_pb{ll}.time(find(mean(stat_TPAmMSPN1_peak_pb{ll}.mask,1))))]; catch end
      try  aa(2,:)=[min(stat_TPAmMSPN1_trgh_pb{ll}.time(find(mean(stat_TPAmMSPN1_trgh_pb{ll}.mask,1)))) max(stat_TPAmMSPN1_trgh_pb{ll}.time(find(mean(stat_TPAmMSPN1_trgh_pb{ll}.mask,1))))]; catch end
      cfg.latency=[min(aa(:,1)) max(aa(:,2))];
      stat_TPAmMSPN_peakMtrgh_pb_posthoctime{ll}=ft_timelockstatistics(cfg, grind_TPAmMSPN1_pb{ll,1}, grind_TPAmMSPN1_pb{ll,3});
      
      if ll==3 || ll==5
        aa=nan(2,2);
        try  aa(1,:)=[min(stat_TPAmMSPN1_ptot_pb{ll}.time(find(mean(stat_TPAmMSPN1_ptot_pb{ll}.mask,1)))) max(stat_TPAmMSPN1_ptot_pb{ll}.time(find(mean(stat_TPAmMSPN1_ptot_pb{ll}.mask,1))))]; catch end
        try  aa(2,:)=[min(stat_TPAmMSPN1_ttop_pb{ll}.time(find(mean(stat_TPAmMSPN1_ttop_pb{ll}.mask,1)))) max(stat_TPAmMSPN1_ttop_pb{ll}.time(find(mean(stat_TPAmMSPN1_ttop_pb{ll}.mask,1))))]; catch end
        cfg.latency=[min(aa(:,1)) max(aa(:,2))];
        stat_TPAmMSPN_ptotMttop_pb_posthoctime{ll}=ft_timelockstatistics(cfg, grind_TPAmMSPN1_pb{ll,2}, grind_TPAmMSPN1_pb{ll,4});
      end
    end
    
    
  end
  
  %   save([edir 'stat_pb_mult.mat'],'stat_TPAmMSPN_*M*','-append')
  %   save([edir 'stat_pb_mult_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat'],'stat_TPAmM*M*','-append')
  save([edir 'stat_pb_mult_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag) '.mat'],'stat_TPAmM*M*','-append');
end

% going through results
for ll=soalist
  anysig(ll,1)=length(find(stat_tacPaud1_pb{ll}.mask));
  anysig(ll,2)=length(find(stat_tacMSpN1_pb{ll}.mask));
  anysig(ll,3)=length(find(stat_TPAmMSPN1_pb{ll}.mask));
  anysig(ll,4)=length(find(stat_tacPaud1IAF_pb{ll}.mask));
  anysig(ll,5)=length(find(stat_tacMSpN1IAF_pb{ll}.mask));
  anysig(ll,6)=length(find(stat_TPAmMSPN1IAF_pb{ll}.mask));
  
  anysig1(ll,1)=length(find(stat_TPAmMSPN1_peak_pb{ll}.mask));
  anysig1(ll,2)=length(find(stat_TPAmMSPN1_ptot_pb{ll}.mask));
  anysig1(ll,3)=length(find(stat_TPAmMSPN1_trgh_pb{ll}.mask));
  anysig1(ll,4)=length(find(stat_TPAmMSPN1_ttop_pb{ll}.mask));
  anysig1(ll,5)=length(find(stat_TPAmMSPN1IAF_peak_pb{ll}.mask));
  anysig1(ll,6)=length(find(stat_TPAmMSPN1IAF_ptot_pb{ll}.mask));
  anysig1(ll,7)=length(find(stat_TPAmMSPN1IAF_trgh_pb{ll}.mask));
  anysig1(ll,8)=length(find(stat_TPAmMSPN1IAF_ttop_pb{ll}.mask));
  
end

%% Plotting phase-dependent results (mult)

% clearvars -except sub *dir
chanplot{1}={'Fz' 'Cz' 'F1' 'F2' 'FC1' 'FC2' 'C1' 'C2'}; % frontocentral
chanplot{2}={'CP5' 'POz' 'Pz' 'P3' 'P4' 'C4' 'O1' 'O2' 'P7' 'PO7'}; % occipital
chanlabel{1}='Frontocentral electrodes';
chanlabel{2}='Occipital-parietal electrodes';
iterflag=1;
tt=3;

sleep=1;
sleepshortstatflag=0;
if sleep
  iter=11;
  ss=12;
  trialkc=-1;
else
  iter=27;
  ss=10;
  trialkc=-1;
end

if iterflag
  if sleep
    try
      load([edir 'stat_pb_mult_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag) '.mat']);
    catch
      load([edir 'stat_pb_mult_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat']);
    end
    load([edir 'grind_pb_mult_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat']);
  else
    load([edir 'stat_pb_mult.mat']);
    load([edir 'grind_pb_mult.mat']);
  end
end

soalist=[1 3 4 5 6 7 9];
timwinstatflag=0; % =1 for using only stat.time, =0 for using [-0.5 1.0];
if sleepshortstatflag
  timwin=[-0.5 1];
else
  timwin=[-0.5 1.5];
end
soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];

scalediff=1;

bluered=varycolor(10);
coloruse=parula(8);

% For ANOVA results
for pb=1:2 % pb=1 means phase results, pb=2 means power results
  close all
  clear tmp*
  for ll=soalist
    
    cfg=[];
    if timwinstatflag==1
      cfg.latency=[stat_TPAmMSPN1_peak_pb{ll}.time(1) stat_TPAmMSPN1_peak_pb{ll}.time(end)];
    elseif timwinstatflag==0
      cfg.latency=timwin;
      stattimwin=[stat_TPAmMSPN1_peak_pb{ll}.time(1) stat_TPAmMSPN1_peak_pb{ll}.time(end)];
    end
    cfg.channel=stat_TPAmMSPN1_peak_pb{ll}.label;
    
    switch pb
      case 1
        tmp1=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,1});
        tmp1.dimord='chan_time';
        tmp1.avg=squeeze(mean(tmp1.individual,1));
        tmp1=rmfield(tmp1,'individual');
        
        tmp5=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,2});
        tmp5.dimord='chan_time';
        tmp5.avg=squeeze(mean(tmp5.individual,1));
        tmp5=rmfield(tmp5,'individual');
        
        tmp8=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,3});
        tmp8.dimord='chan_time';
        tmp8.avg=squeeze(mean(tmp8.individual,1));
        tmp8=rmfield(tmp8,'individual');
        
        tmp12=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,4});
        tmp12.dimord='chan_time';
        tmp12.avg=squeeze(mean(tmp12.individual,1));
        tmp12=rmfield(tmp12,'individual');
        
        tmpmask=stat_TPAmMSPN1_pb{ll}.mask;
      case 2
        tmp1=ft_selectdata(cfg,grind_TPAmMSPN1_ab{ll,1});
        tmp1.dimord='chan_time';
        tmp1.avg=squeeze(mean(tmp1.individual,1));
        tmp1=rmfield(tmp1,'individual');
        
        tmp5=ft_selectdata(cfg,grind_TPAmMSPN1_ab{ll,2});
        tmp5.dimord='chan_time';
        tmp5.avg=squeeze(mean(tmp5.individual,1));
        tmp5=rmfield(tmp5,'individual');
        
        tmp8=ft_selectdata(cfg,grind_TPAmMSPN1_ab{ll,3});
        tmp8.dimord='chan_time';
        tmp8.avg=squeeze(mean(tmp8.individual,1));
        tmp8=rmfield(tmp8,'individual');
        
        tmp12=ft_selectdata(cfg,grind_TPAmMSPN1_ab{ll,4});
        tmp12.dimord='chan_time';
        tmp12.avg=squeeze(mean(tmp12.individual,1));
        tmp12=rmfield(tmp12,'individual');
        
        tmpmask=stat_TPAmMSPN1_ab{ll}.mask;
    end
    
    
    if timwinstatflag==0
      tmp1.mask=zeros(size(tmp1.avg,1),length(tmp1.time));
      tmp1.mask(:,dsearchn(tmp1.time',stattimwin(1)):dsearchn(tmp1.time',stattimwin(end)))=tmpmask;
      tmp5.mask=zeros(size(tmp5.avg,1),length(tmp5.time));
      tmp5.mask(:,dsearchn(tmp5.time',stattimwin(1)):dsearchn(tmp5.time',stattimwin(end)))=tmpmask;
      tmp8.mask=zeros(size(tmp8.avg,1),length(tmp8.time));
      tmp8.mask(:,dsearchn(tmp8.time',stattimwin(1)):dsearchn(tmp8.time',stattimwin(end)))=tmpmask;
      tmp12.mask=zeros(size(tmp12.avg,1),length(tmp12.time));
      tmp12.mask(:,dsearchn(tmp12.time',stattimwin(1)):dsearchn(tmp12.time',stattimwin(end)))=tmpmask;
    else
      tmpu1{bb}.mask=tmpmask
      tmpm10{bb}.mask=tmpmask;
      tmpd5{bb}.mask=tmpmask;
    end
    
    for cg=1:length(chanplot)
      cfg=[];
      cfg.parameter='avg';
      cfg.layout='elec1010.lay';
      cfg.ylim=[-5 8];
      cfg.linewidth=3;
      cfg.xlim=timwin;
      if cg>length(chanplot)
        cfg.channel=tmp5.label(any(tmp5.mask,2));
      else
        cfg.channel=chanplot{cg};
      end
      cfg.graphcolor=coloruse([pb:2:8],:);
      cfg.interactive='no';
      cfg.maskparameter='mask';
      cfg.maskstyle='box'; % default
      figure(ll+10*(cg+1))
      ft_singleplotER(cfg,tmp1,tmp5,tmp8,tmp12);
      hold on;plot(tmp1.time,0,'k');
      set(gca,'XTick',[-.5:.1:1])
      set(gca,'XTickLabel',{'-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'})
      set(gca,'FontSize',30)
      title([])
      plot([0 0],cfg.ylim,'Color',bluered(4,:),'LineWidth',6)
      plot([soades(ll) soades(ll)],cfg.ylim,'Color',bluered(9,:),'LineWidth',6)
      plot([stattimwin(1) stattimwin(1)],cfg.ylim,'k--','LineWidth',6)
      plot([stattimwin(2) stattimwin(2)],cfg.ylim,'k--','LineWidth',6)
      axis([-0.55 1.5 cfg.ylim(1) cfg.ylim(2)])
      if cg==3
        legend('Sum Unisensory','MultSens + Null','SumUnisens - MultsensNull')
      end
    end
    switch pb
      case 1
        print(ll+20,[fdir 'erp_tacPaud_MSpN_diff_FC_pb_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag) '.png'],'-dpng')
        print(ll+30,[fdir 'erp_tacPaud_MSpN_diff_OP_pb_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag) '.png'],'-dpng')
      case 2
        print(ll+20,[fdir 'erp_tacPaud_MSpN_diff_FC_ab_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag) '.png'],'-dpng')
        print(ll+30,[fdir 'erp_tacPaud_MSpN_diff_OP_ab_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag) '.png'],'-dpng')
    end
    
    
    if any(tmpmask(:))
      masktime=find(any(tmp5.mask,1));
      cfg=[];
      cfg.parameter='avg';
      cfg.layout='elec1010.lay';
      cfg.maskalpha=0.5;
      cfg.zlim=[-5 5];
      cfg.highlight='on';
      cfg.highlightsize=12;
      cfg.xlim=[tmp5.time(masktime(1)) tmp5.time(masktime(end))];
      cfg.comment='no';
      sigchannels=tmp5.label(find(ceil(mean(tmp5.mask(:,dsearchn(tmp5.time',cfg.xlim(1)):dsearchn(tmp5.time',cfg.xlim(2))),2))));
      cfg.highlightchannel=sigchannels;
      figure(1000*bb+ll);
      ft_topoplotER(cfg,tmp1);
      switch pb
        case 1
          print(1000*bb+ll,[fdir 'erp_topoBin1_pb_' num2str(ll) num2str(tt) num2str(ss) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        case 2
          print(1000*bb+ll,[fdir 'erp_topoBin1_ab_' num2str(ll) num2str(tt) num2str(ss) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
      end
      figure(1000*bb+10+ll);
      ft_topoplotER(cfg,tmp5);
      switch pb
        case 1
          print(1000*bb+10+ll,[fdir 'erp_topoBin2_pb_' num2str(ll) num2str(tt) num2str(ss) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        case 2
          print(1000*bb+10+ll,[fdir 'erp_topoBin2_ab_' num2str(ll) num2str(tt) num2str(ss) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
      end
      figure(1000*bb+20+ll);
      ft_topoplotER(cfg,tmp8);
      switch pb
        case 1
          print(1000*bb+20+ll,[fdir 'erp_topoBin3_pb_' num2str(ll) num2str(tt) num2str(ss) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        case 2
          print(1000*bb+20+ll,[fdir 'erp_topoBin3_ab_' num2str(ll) num2str(tt) num2str(ss) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
      end
      figure(1000*bb+30+ll);
      ft_topoplotER(cfg,tmp12);
      switch pb
        case 1
          print(1000*bb+20+ll,[fdir 'erp_topoBin4_pb_' num2str(ll) num2str(tt) num2str(ss) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        case 2
          print(1000*bb+20+ll,[fdir 'erp_topoBin4_ab_' num2str(ll) num2str(tt) num2str(ss) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
      end
    end
    
    
  end % ll
end




% for individual bin results
for pb=1:2
  close all
  clear tmp*
  for ll=soalist
    for bb=1:4
      
      cfg=[];
      if timwinstatflag==1
        cfg.latency=[stat_TPAmMSPN1_peak_pb{ll}.time(1) stat_TPAmMSPN1_peak_pb{ll}.time(end)];
      elseif timwinstatflag==0
        cfg.latency=timwin;
        stattimwin=[stat_TPAmMSPN1_peak_pb{ll}.time(1) stat_TPAmMSPN1_peak_pb{ll}.time(end)];
      end
      cfg.channel=stat_TPAmMSPN1_peak_pb{ll}.label;
      
      switch pb
        case 1
          tmpu1{bb}=ft_selectdata(cfg,grind_tacPaud1_pb{ll,bb});
          tmpu1{bb}.dimord='chan_time';
          tmpu1{bb}.avg=squeeze(mean(tmpu1{bb}.individual,1));
          tmpu1{bb}=rmfield(tmpu1{bb},'individual');
          
          tmpm10{bb}=ft_selectdata(cfg,grind_tacMSpN1_pb{ll,bb});
          tmpm10{bb}.dimord='chan_time';
          tmpm10{bb}.avg=squeeze(mean(tmpm10{bb}.individual,1));
          tmpm10{bb}=rmfield(tmpm10{bb},'individual');
          
          tmpd5{bb}=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,bb});
          tmpd5{bb}.dimord='chan_time';
          tmpd5{bb}.avg=scalediff*squeeze(mean(tmpd5{bb}.individual,1));
          tmpd5{bb}=rmfield(tmpd5{bb},'individual');
          
          switch bb
            case 1
              tmpmask=stat_TPAmMSPN1_peak_pb{ll}.mask;
            case 2
              tmpmask=stat_TPAmMSPN1_ptot_pb{ll}.mask;
            case 3
              tmpmask=stat_TPAmMSPN1_trgh_pb{ll}.mask;
            case 4
              tmpmask=stat_TPAmMSPN1_ttop_pb{ll}.mask;
          end
        case 2
          tmpu1{bb}=ft_selectdata(cfg,grind_tacPaud1_ab{ll,bb});
          tmpu1{bb}.dimord='chan_time';
          tmpu1{bb}.avg=squeeze(mean(tmpu1{bb}.individual,1));
          tmpu1{bb}=rmfield(tmpu1{bb},'individual');
          
          tmpm10{bb}=ft_selectdata(cfg,grind_tacMSpN1_ab{ll,bb});
          tmpm10{bb}.dimord='chan_time';
          tmpm10{bb}.avg=squeeze(mean(tmpm10{bb}.individual,1));
          tmpm10{bb}=rmfield(tmpm10{bb},'individual');
          
          tmpd5{bb}=ft_selectdata(cfg,grind_TPAmMSPN1_ab{ll,bb});
          tmpd5{bb}.dimord='chan_time';
          tmpd5{bb}.avg=scalediff*squeeze(mean(tmpd5{bb}.individual,1));
          tmpd5{bb}=rmfield(tmpd5{bb},'individual');
          
          switch bb
            case 1
              tmpmask=stat_TPAmMSPN1_low_pb{ll}.mask;
            case 2
              tmpmask=stat_TPAmMSPN1_midl_pb{ll}.mask;
            case 3
              tmpmask=stat_TPAmMSPN1_midh_pb{ll}.mask;
            case 4
              tmpmask=stat_TPAmMSPN1_high_pb{ll}.mask;
          end
      end
      
      if timwinstatflag==0
        tmpu1{bb}.mask=zeros(size(tmpu1{bb}.avg,1),length(tmpu1{bb}.time));
        tmpu1{bb}.mask(:,dsearchn(tmpu1{bb}.time',stattimwin(1)):dsearchn(tmpu1{bb}.time',stattimwin(end)))=tmpmask;
        tmpm10{bb}.mask=zeros(size(tmpm10{bb}.avg,1),length(tmpm10{bb}.time));
        tmpm10{bb}.mask(:,dsearchn(tmpm10{bb}.time',stattimwin(1)):dsearchn(tmpm10{bb}.time',stattimwin(end)))=tmpmask;
        tmpd5{bb}.mask=zeros(size(tmpd5{bb}.avg,1),length(tmpd5{bb}.time));
        tmpd5{bb}.mask(:,dsearchn(tmpd5{bb}.time',stattimwin(1)):dsearchn(tmpd5{bb}.time',stattimwin(end)))=tmpmask;
      else
        tmpu1{bb}.mask=tmpmask
        tmpm10{bb}.mask=tmpmask;
        tmpd5{bb}.mask=tmpmask;
      end
      
      for cg=1:length(chanplot)
        cfg=[];
        cfg.parameter='avg';
        cfg.layout='elec1010.lay';
        cfg.ylim=[-5 8];
        cfg.linewidth=3;
        cfg.xlim=timwin;
        if cg>length(chanplot)
          cfg.channel=tmpd5{bb}.label(any(tmpd5{bb}.mask,2));
        else
          cfg.channel=chanplot{cg};
        end
        cfg.graphcolor=[bluered(1,:); bluered(10,:); coloruse(2*(bb-1)+pb,:)];
        cfg.interactive='no';
        cfg.maskparameter='mask';
        cfg.maskstyle='box'; % default
        figure(100*bb+ll+10*(cg+1))
        ft_singleplotER(cfg,tmpu1{bb},tmpm10{bb},tmpd5{bb});
        hold on;plot(tmpu1{bb}.time,0,'k');
        set(gca,'XTick',[-.5:.1:1])
        set(gca,'XTickLabel',{'-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'})
        set(gca,'FontSize',30)
        title([])
        plot([0 0],cfg.ylim,'Color',bluered(4,:),'LineWidth',6)
        plot([soades(ll) soades(ll)],cfg.ylim,'Color',bluered(9,:),'LineWidth',6)
        plot([stattimwin(1) stattimwin(1)],cfg.ylim,'k--','LineWidth',6)
        plot([stattimwin(2) stattimwin(2)],cfg.ylim,'k--','LineWidth',6)
        axis([-0.55 1.5 cfg.ylim(1) cfg.ylim(2)])
        if cg==3
          legend('Sum Unisensory','MultSens + Null','SumUnisens - MultsensNull')
        end
      end
      switch pb
        case 1
          print(100*bb+ll+20,[fdir 'erp_tacPaud_MSpN_diff_FC_pb_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_bin' num2str(bb)  '_shortstat' num2str(sleepshortstatflag) '.png'],'-dpng')
          print(100*bb+ll+30,[fdir 'erp_tacPaud_MSpN_diff_OP_pb_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_bin' num2str(bb)  '_shortstat' num2str(sleepshortstatflag) '.png'],'-dpng')
        case 2
          print(100*bb+ll+20,[fdir 'erp_tacPaud_MSpN_diff_FC_ab_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_bin' num2str(bb)  '_shortstat' num2str(sleepshortstatflag) '.png'],'-dpng')
          print(100*bb+ll+30,[fdir 'erp_tacPaud_MSpN_diff_OP_ab_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_bin' num2str(bb)  '_shortstat' num2str(sleepshortstatflag) '.png'],'-dpng')
      end
      
      if any(tmpmask(:))
        masktime=find(any(tmpd5{bb}.mask,1));
        cfg=[];
        cfg.parameter='avg';
        cfg.layout='elec1010.lay';
        cfg.maskalpha=0.5;
        cfg.zlim=[-5 5];
        cfg.highlight='on';
        cfg.highlightsize=12;
        cfg.xlim=[tmpd5{bb}.time(masktime(1)) tmpd5{bb}.time(masktime(end))];
        cfg.comment='no';
        sigchannels=tmpd5{bb}.label(find(ceil(mean(tmpd5{bb}.mask(:,dsearchn(tmpd5{bb}.time',cfg.xlim(1)):dsearchn(tmpd5{bb}.time',cfg.xlim(2))),2))));
        cfg.highlightchannel=sigchannels;
        figure(1000*bb+ll);
        ft_topoplotER(cfg,tmpu1{bb});
        switch pb
          case 1
            print(1000*bb+ll,[fdir 'erp_topoU_pb_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
          case 2
            print(1000*bb+ll,[fdir 'erp_topoU_ab_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        end
        figure(1000*bb+10+ll);
        ft_topoplotER(cfg,tmpm10{bb});
        switch pb
          case 1
            print(1000*bb+10+ll,[fdir 'erp_topoM_pb_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
          case 2
            print(1000*bb+10+ll,[fdir 'erp_topoM_ab_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        end
        figure(1000*bb+20+ll);
        ft_topoplotER(cfg,tmpd5{bb});
        switch pb
          case 1
            print(1000*bb+20+ll,[fdir 'erp_topoDiff_pb_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
          case 2
            print(1000*bb+20+ll,[fdir 'erp_topoDiff_ab_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        end
      end
      
    end % bb
  end % ll
end % pb


close all
clear tmp*
% Now, for specific ptotMttop and peakMtrgh and highMlow

for ll=soalist
  
  cfg=[];
  cfg.parameter='individual';
  cfg.operation='subtract';
  grind_TPAmMSPN_4m2_pb{ll}=ft_math(cfg,grind_TPAmMSPN1_pb{ll,4},grind_TPAmMSPN1_pb{ll,2});
  grave_TPAmMSPN_4m2_pb{ll}=grind_TPAmMSPN_4m2_pb{ll};
  grave_TPAmMSPN_4m2_pb{ll}.avg=squeeze(mean(grave_TPAmMSPN_4m2_pb{ll}.individual,1));
  grave_TPAmMSPN_4m2_pb{ll}=rmfield(grave_TPAmMSPN_4m2_pb{ll},'individual');
  grave_TPAmMSPN_4m2_pb{ll}.dimord='chan_time';
  
  grind_TPAmMSPN_3m1_pb{ll}=ft_math(cfg,grind_TPAmMSPN1_pb{ll,3},grind_TPAmMSPN1_pb{ll,1});
  grave_TPAmMSPN_3m1_pb{ll}=grind_TPAmMSPN_3m1_pb{ll};
  grave_TPAmMSPN_3m1_pb{ll}.avg=squeeze(mean(grave_TPAmMSPN_3m1_pb{ll}.individual,1));
  grave_TPAmMSPN_3m1_pb{ll}=rmfield(grave_TPAmMSPN_3m1_pb{ll},'individual');
  grave_TPAmMSPN_3m1_pb{ll}.dimord='chan_time';
  
  grind_TPAmMSPN_4m1_ab{ll}=ft_math(cfg,grind_TPAmMSPN1_ab{ll,4},grind_TPAmMSPN1_ab{ll,1});
  grave_TPAmMSPN_4m1_ab{ll}=grind_TPAmMSPN_4m1_ab{ll};
  grave_TPAmMSPN_4m1_ab{ll}.avg=squeeze(mean(grave_TPAmMSPN_4m1_ab{ll}.individual,1));
  grave_TPAmMSPN_4m1_ab{ll}=rmfield(grave_TPAmMSPN_4m1_ab{ll},'individual');
  grave_TPAmMSPN_4m1_ab{ll}.dimord='chan_time';
  
  for kk=1:3
    cfg=[];
    if timwinstatflag==1
      cfg.latency=[stat_TPAmMSPN1_peak_pb{ll}.time(1) stat_TPAmMSPN1_peak_pb{ll}.time(end)];
    elseif timwinstatflag==0
      cfg.latency=timwin;
      stattimwin=[stat_TPAmMSPN1_peak_pb{ll}.time(1) stat_TPAmMSPN1_peak_pb{ll}.time(end)];
    end
    cfg.channel=stat_TPAmMSPN1_peak_pb{ll}.label;
    
    switch kk
      case 1
        tmpu1=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,2});
        tmpu1.dimord='chan_time';
        tmpu1.avg=squeeze(mean(tmpu1.individual,1));
        tmpu1=rmfield(tmpu1,'individual');
        
        tmpm10=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,4});
        tmpm10.dimord='chan_time';
        tmpm10.avg=squeeze(mean(tmpm10.individual,1));
        tmpm10=rmfield(tmpm10,'individual');
        
        tmpd5=ft_selectdata(cfg,grind_TPAmMSPN_4m2_pb{ll});
        tmpd5.dimord='chan_time';
        tmpd5.avg=scalediff*squeeze(mean(tmpd5.individual,1));
        tmpd5=rmfield(tmpd5,'individual');
      case 2
        tmpu1=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,1});
        tmpu1.dimord='chan_time';
        tmpu1.avg=squeeze(mean(tmpu1.individual,1));
        tmpu1=rmfield(tmpu1,'individual');
        
        tmpm10=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,3});
        tmpm10.dimord='chan_time';
        tmpm10.avg=squeeze(mean(tmpm10.individual,1));
        tmpm10=rmfield(tmpm10,'individual');
        
        tmpd5=ft_selectdata(cfg,grind_TPAmMSPN_3m1_pb{ll});
        tmpd5.dimord='chan_time';
        tmpd5.avg=scalediff*squeeze(mean(tmpd5.individual,1));
        tmpd5=rmfield(tmpd5,'individual');
      case 3
        tmpu1=ft_selectdata(cfg,grind_TPAmMSPN1_ab{ll,1});
        tmpu1.dimord='chan_time';
        tmpu1.avg=squeeze(mean(tmpu1.individual,1));
        tmpu1=rmfield(tmpu1,'individual');
        
        tmpm10=ft_selectdata(cfg,grind_TPAmMSPN1_ab{ll,4});
        tmpm10.dimord='chan_time';
        tmpm10.avg=squeeze(mean(tmpm10.individual,1));
        tmpm10=rmfield(tmpm10,'individual');
        
        tmpd5=ft_selectdata(cfg,grind_TPAmMSPN_4m1_ab{ll});
        tmpd5.dimord='chan_time';
        tmpd5.avg=scalediff*squeeze(mean(tmpd5.individual,1));
        tmpd5=rmfield(tmpd5,'individual');
    end
    
    switch kk
      case 1
        tmpmask=stat_TPAmMSPN_ptotMttop_pb{ll}.mask;
      case 2
        tmpmask=stat_TPAmMSPN_peakMtrgh_pb{ll}.mask;
      case 3
        if trialkc==-1
          tmpmask=stat_TPAmMSPN_lowMhigh_ab{ll}.mask;
        elseif trialkc==0
          tmpmask=stat_TPAmMSPN_highMlow_ab{ll}.mask;
        end
    end
    
    if timwinstatflag==0
      tmpu1.mask=zeros(size(tmpu1.avg,1),length(tmpu1.time));
      tmpu1.mask(:,dsearchn(tmpu1.time',stattimwin(1)):dsearchn(tmpu1.time',stattimwin(end)))=tmpmask;
      tmpm10.mask=zeros(size(tmpm10.avg,1),length(tmpm10.time));
      tmpm10.mask(:,dsearchn(tmpm10.time',stattimwin(1)):dsearchn(tmpm10.time',stattimwin(end)))=tmpmask;
      tmpd5.mask=zeros(size(tmpd5.avg,1),length(tmpd5.time));
      tmpd5.mask(:,dsearchn(tmpd5.time',stattimwin(1)):dsearchn(tmpd5.time',stattimwin(end)))=tmpmask;
    else
      tmpu1.mask=tmpmask;
      tmpm10.mask=tmpmask;
      tmpd5.mask=tmpmask;
    end
    
    for cg=1:length(chanplot)
      cfg=[];
      cfg.parameter='avg';
      cfg.layout='elec1010.lay';
      cfg.ylim=[-5 8];
      cfg.linewidth=3;
      cfg.xlim=timwin;
      if cg>length(chanplot)
        cfg.channel=tmpd5.label(any(tmpd5.mask,2));
      else
        cfg.channel=chanplot{cg};
      end
      switch kk
        case 1
          cfg.graphcolor=[coloruse(3,:); coloruse(7,:);  bluered(10,:)]
        case 2
          cfg.graphcolor=[coloruse(1,:); coloruse(5,:);  bluered(10,:)]
        case 3
          cfg.graphcolor=[coloruse(2,:); coloruse(8,:);  bluered(10,:)]
      end
      cfg.interactive='no';
      cfg.maskparameter='mask';
      cfg.maskstyle='box'; % default
      figure(100*bb+ll+10*(cg+1))
      ft_singleplotER(cfg,tmpu1,tmpm10,tmpd5);
      hold on;plot(tmpu1.time,0,'k');
      set(gca,'XTick',[-.5:.1:1])
      set(gca,'XTickLabel',{'-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'})
      set(gca,'FontSize',30)
      title([])
      plot([0 0],cfg.ylim,'Color',bluered(4,:),'LineWidth',6)
      plot([soades(ll) soades(ll)],cfg.ylim,'Color',bluered(9,:),'LineWidth',6)
      plot([stattimwin(1) stattimwin(1)],cfg.ylim,'k--','LineWidth',6)
      plot([stattimwin(2) stattimwin(2)],cfg.ylim,'k--','LineWidth',6)
      axis([-0.55 1.5 cfg.ylim(1) cfg.ylim(2)])
      if cg==3
        legend('Sum Unisensory','MultSens + Null','SumUnisens - MultsensNull')
      end
    end
    %   print(30+ll,[fdir 'erp_tacPaud_MSpN_diff_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    switch kk
      case 1
        print(100*bb+ll+20,[fdir 'erp_tacPaud_MSpN_diff_FC_pb_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag)  '_bindiff42.png'],'-dpng')
        print(100*bb+ll+30,[fdir 'erp_tacPaud_MSpN_diff_OP_pb_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag)  '_bindiff42.png'],'-dpng')
      case 2
        print(100*bb+ll+20,[fdir 'erp_tacPaud_MSpN_diff_FC_pb_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag)  '_bindiff31.png'],'-dpng')
        print(100*bb+ll+30,[fdir 'erp_tacPaud_MSpN_diff_OP_pb_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag)  '_bindiff31.png'],'-dpng')
      case 3
        print(100*bb+ll+20,[fdir 'erp_tacPaud_MSpN_diff_FC_ab_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag)  '_bindiff41.png'],'-dpng')
        print(100*bb+ll+30,[fdir 'erp_tacPaud_MSpN_diff_OP_ab_' num2str(ll) num2str(tt) num2str(ss) '_trialkc' num2str(trialkc) '_shortstat' num2str(sleepshortstatflag)  '_bindiff41.png'],'-dpng')
    end
    
    if any(tmpmask(:))
      masktime=find(any(tmpm10.mask,1));
      cfg=[];
      cfg.parameter='avg';
      cfg.layout='elec1010.lay';
      cfg.maskalpha=0.5;
      cfg.zlim=[-5 5];
      cfg.highlight='on';
      cfg.highlightsize=12;
      cfg.xlim=[tmpm10.time(masktime(1)) tmpm10.time(masktime(end))];
      cfg.comment='no';
      sigchannels=tmpm10.label(find(ceil(mean(tmpm10.mask(:,dsearchn(tmpm10.time',cfg.xlim(1)):dsearchn(tmpm10.time',cfg.xlim(2))),2))));
      cfg.highlightchannel=sigchannels;
      figure(1000*bb+ll);
      ft_topoplotER(cfg,tmpu1);
      switch kk
        case 1
          print(1000*bb+ll,[fdir 'erp_topoU_pb_4m2_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        case 2
          print(1000*bb+ll,[fdir 'erp_topoU_pb_3m1_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        case 3
          print(1000*bb+ll,[fdir 'erp_topoU_ab_4m1_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
      end
      figure(1000*bb+10+ll);
      ft_topoplotER(cfg,tmpm10);
      switch kk
        case 1
          print(1000*bb+10+ll,[fdir 'erp_topoM_pb_4m2_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        case 2
          print(1000*bb+10+ll,[fdir 'erp_topoM_pb_3m1_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        case 3
          print(1000*bb+10+ll,[fdir 'erp_topoM_ab_4m1_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
      end
      figure(1000*bb+20+ll);
      ft_topoplotER(cfg,tmpd5);
      switch kk
        case 1
          print(1000*bb+20+ll,[fdir 'erp_topoDiff_pb_4m2_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        case 2
          print(1000*bb+20+ll,[fdir 'erp_topoDiff_pb_3m1_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
        case 3
          print(1000*bb+20+ll,[fdir 'erp_topoDiff_ab_4m1_' num2str(ll) num2str(tt) num2str(ss) '_bin' num2str(bb) '_shortstat' num2str(sleepshortstatflag)  '.png'],'-dpng')
      end
    end
    
  end % bb
end % ll



for ll=5 % 1vs3 and 2vs4 (ll=5 only significant finding)
  cfg=[];
  if timwinstatflag==1
    cfg.latency=[stat_TPAmMSPN_peakMtrgh_pb_posthoctime{ll}.time(1) stat_TPAmMSPN_peakMtrgh_pb_posthoctime{ll}.time(end)];
  elseif timwinstatflag==0
    cfg.latency=timwin;
    stattimwin=[stat_TPAmMSPN_peakMtrgh_pb_posthoctime{ll}.time(1) stat_TPAmMSPN_peakMtrgh_pb_posthoctime{ll}.time(end)];
  end
  cfg.channel=stat_TPAmMSPN_peakMtrgh_pb_posthoctime{ll}.label;
  
  % FINISH ME HERE!
  
  tmpu1=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,1});
  tmpu1.dimord='chan_time';
  tmpu1.avg=squeeze(mean(tmpu1.individual,1));
  tmpu1=rmfield(tmpu1,'individual');
  
  tmpm10=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,3});
  tmpm10.dimord='chan_time';
  tmpm10.avg=squeeze(mean(tmpm10.individual,1));
  tmpm10=rmfield(tmpm10,'individual');
  
  tmpd5=ft_selectdata(cfg,grind_TPAmMSPN_3m1_pb{ll});
  tmpd5.dimord='chan_time';
  tmpd5.avg=scalediff*squeeze(mean(tmpd5.individual,1));
  tmpd5=rmfield(tmpd5,'individual');
  
  
  tmpmask=stat_TPAmMSPN_peakMtrgh_pb_posthoctime{ll}.mask;
  if timwinstatflag==0
    tmpu1.mask=zeros(size(tmpu1.avg,1),length(tmpu1.time));
    tmpu1.mask(:,dsearchn(tmpu1.time',stattimwin(1)):dsearchn(tmpu1.time',stattimwin(end)))=tmpmask;
    tmpm10.mask=zeros(size(tmpm10.avg,1),length(tmpm10.time));
    tmpm10.mask(:,dsearchn(tmpm10.time',stattimwin(1)):dsearchn(tmpm10.time',stattimwin(end)))=tmpmask;
    tmpd5.mask=zeros(size(tmpd5.avg,1),length(tmpd5.time));
    tmpd5.mask(:,dsearchn(tmpd5.time',stattimwin(1)):dsearchn(tmpd5.time',stattimwin(end)))=tmpmask;
  else
    tmpu1{bb}.mask=tmpmask
    tmpm10{bb}.mask=tmpmask;
    tmpd5{bb}.mask=tmpmask;
  end
  
  for cg=1:length(chanplot)
    cfg=[];
    cfg.parameter='avg';
    cfg.layout='elec1010.lay';
    cfg.ylim=[-5 8];
    cfg.linewidth=3;
    cfg.xlim=timwin;
    if cg>length(chanplot)
      cfg.channel=tmpd5.label(any(tmpd5.mask,2));
    else
      cfg.channel=chanplot{cg};
    end
    cfg.graphcolor=coloruse([5 5 2],:);
    cfg.interactive='no';
    cfg.maskparameter='mask';
    cfg.maskstyle='box'; % default
    figure(ll+10*(cg+1))
    ft_singleplotER(cfg,tmpu1,tmpm10,tmpd5);
    hold on;plot(tmpu1.time,0,'k');
    set(gca,'XTick',[-.5:.1:1])
    set(gca,'XTickLabel',{'-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'})
    set(gca,'FontSize',30)
    title([])
    plot([stattimwin(1) stattimwin(1)],cfg.ylim,'k--','LineWidth',6)
    plot([stattimwin(2) stattimwin(2)],cfg.ylim,'k--','LineWidth',6)
    plot([0 0],cfg.ylim,'Color',coloruse(4,:),'LineWidth',6)
    plot([soades(ll) soades(ll)],cfg.ylim,'Color',coloruse(9,:),'LineWidth',6)
    axis([-0.55 1 cfg.ylim(1) cfg.ylim(2)])
    if cg==3
      legend('Sum Unisensory','MultSens + Null','SumUnisens - MultsensNull')
    end
  end
  %   print(30+ll,[fdir 'erp_tacPaud_MSpN_diff_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  print(ll+20,[fdir 'erp_tacPaud_MSpN_diff_PeakVsTrgh_FC_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  print(ll+30,[fdir 'erp_tacPaud_MSpN_diff_PeakVsTrgh_OP_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  
  if any(tmpmask(:))
    masktime=find(any(tmpd5.mask,1));
    cfg=[];
    cfg.parameter='avg';
    cfg.layout='elec1010.lay';
    cfg.maskalpha=0.5;
    cfg.zlim=[-5 5];
    cfg.highlight='on';
    cfg.highlightsize=12;
    cfg.xlim=[tmpd5.time(masktime(1)) tmpd5.time(masktime(end))];
    cfg.comment='no';
    sigchannels=tmpd5.label(find(ceil(mean(tmpd5.mask(:,dsearchn(tmpd5.time',cfg.xlim(1)):dsearchn(tmpd5.time',cfg.xlim(2))),2))));
    cfg.highlightchannel=sigchannels;
    figure(100+ll);
    ft_topoplotER(cfg,tmpu1);
    print(100+ll,[fdir 'erp_topoU_PeakVsTrgh_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    figure(100+10+ll);
    ft_topoplotER(cfg,tmpm10);
    print(100+10+ll,[fdir 'erp_topoM_PeakVsTrgh_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    figure(100+20+ll);
    ft_topoplotER(cfg,tmpd5);
    print(100+20+ll,[fdir 'erp_topoDiff_PeakVsTrgh_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  end
  
end

%% Plotting phase-dependent results (Uni vs Nul)

% clearvars -except sub *dir
chanplot{1}={'Fz' 'Cz' 'F1' 'F2' 'FC1' 'FC2' 'C1' 'C2'}; % frontocentral
chanplot{2}={'CP5' 'POz' 'Pz' 'P3' 'P4' 'C4' 'O1' 'O2' 'P7' 'PO7'}; % occipital
chanlabel{1}='Frontocentral electrodes';
chanlabel{2}='Occipital-parietal electrodes';
iterflag=1;
sleep=0;
tt=3;
if sleep
  iter=11;
  ss=12;
  trialkc=0;
else
  iter=27;
  ss=10;
  trialkc=-1;
end

if iterflag
  if sleep
  else
    load([edir 'stat_pb_uninul.mat']);
    load([edir 'grind_pb_uninul.mat']);
  end
end

soalist=[1 3 4 5 6 7 9];
timwinstatflag=0; % =1 for using only stat.time, =0 for using [-0.5 1.0];
timwin=[-0.5 1];
soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];

scalediff=2;

coloruse=varycolor(10);
% optimsied to be maximally apart and distinct
% 1  TacPAud
% 2  Nul
% 3  MS_synch
% 4  Tac
% 5  diff_TPAMSPN
% 6  diff_synchAsynch
% 7  MS
% 8  MS_asynch
% 9  Aud
% 10 MSpN

close all
clear tmp*

condname={'AT500' '' 'AT70' 'AT20' 'AT0' 'TA20' 'TA70' '' 'TA500' '' 'AUD' 'TAC'};
for ll=[soalist 11 12]
  
  cfg=[];
  cfg.channel=grind_tacMnul_pb{1}.label;
  
  switch ll
    case {1, 3, 4, 5, 6, 7, 9}
      if timwinstatflag==1
        cfg.latency=[stat_msMnul1_pb{ll}.time(1) stat_msMnul1_pb{ll}.time(end)];
      elseif timwinstatflag==0
        cfg.latency=timwin;
        stattimwin=[stat_msMnul1_pb{ll}.time(1) stat_msMnul1_pb{ll}.time(end)];
      end
      tmpu1=ft_selectdata(cfg,grind_msMnul1_pb{ll,1});
      tmpu1.dimord='chan_time';
      tmpu1.avg=squeeze(mean(tmpu1.individual,1));
      tmpu1=rmfield(tmpu1,'individual');
      
      tmpm10=ft_selectdata(cfg,grind_msMnul1_pb{ll,3});
      tmpm10.dimord='chan_time';
      tmpm10.avg=squeeze(mean(tmpm10.individual,1));
      tmpm10=rmfield(tmpm10,'individual');
      
      tmpd4=ft_selectdata(cfg,grind_msMnul1_pb{ll,2});
      tmpd4.dimord='chan_time';
      tmpd4.avg=squeeze(mean(tmpd4.individual,1));
      tmpd4=rmfield(tmpd4,'individual');
      
      tmpd7=ft_selectdata(cfg,grind_msMnul1_pb{ll,4});
      tmpd7.dimord='chan_time';
      tmpd7.avg=squeeze(mean(tmpd7.individual,1));
      tmpd7=rmfield(tmpd7,'individual');
      
      tmpmask=stat_msMnul1_pb{ll}.mask;
      
    case 11
      if timwinstatflag==1
        cfg.latency=[stat_audMnul_pb.time(1) stat_audMnul_pb.time(end)];
      elseif timwinstatflag==0
        cfg.latency=timwin;
        stattimwin=[stat_audMnul_pb.time(1) stat_audMnul_pb.time(end)];
      end
      tmpu1=ft_selectdata(cfg,grind_audMnul_pb{1});
      tmpu1.dimord='chan_time';
      tmpu1.avg=squeeze(mean(tmpu1.individual,1));
      tmpu1=rmfield(tmpu1,'individual');
      
      tmpm10=ft_selectdata(cfg,grind_audMnul_pb{3});
      tmpm10.dimord='chan_time';
      tmpm10.avg=squeeze(mean(tmpm10.individual,1));
      tmpm10=rmfield(tmpm10,'individual');
      
      tmpd4=ft_selectdata(cfg,grind_audMnul_pb{2});
      tmpd4.dimord='chan_time';
      tmpd4.avg=squeeze(mean(tmpd4.individual,1));
      tmpd4=rmfield(tmpd4,'individual');
      
      tmpd7=ft_selectdata(cfg,grind_audMnul_pb{4});
      tmpd7.dimord='chan_time';
      tmpd7.avg=squeeze(mean(tmpd7.individual,1));
      tmpd7=rmfield(tmpd7,'individual');
      
      tmpmask=stat_audMnul_pb.mask;
      
    case 12
      if timwinstatflag==1
        cfg.latency=[stat_tacMnul_pb.time(1) stat_tacMnul_pb.time(end)];
      elseif timwinstatflag==0
        cfg.latency=timwin;
        stattimwin=[stat_tacMnul_pb.time(1) stat_tacMnul_pb.time(end)];
      end
      tmpu1=ft_selectdata(cfg,grind_tacMnul_pb{1});
      tmpu1.dimord='chan_time';
      tmpu1.avg=squeeze(mean(tmpu1.individual,1));
      tmpu1=rmfield(tmpu1,'individual');
      
      tmpm10=ft_selectdata(cfg,grind_tacMnul_pb{3});
      tmpm10.dimord='chan_time';
      tmpm10.avg=squeeze(mean(tmpm10.individual,1));
      tmpm10=rmfield(tmpm10,'individual');
      
      tmpd4=ft_selectdata(cfg,grind_tacMnul_pb{2});
      tmpd4.dimord='chan_time';
      tmpd4.avg=squeeze(mean(tmpd4.individual,1));
      tmpd4=rmfield(tmpd4,'individual');
      
      tmpd7=ft_selectdata(cfg,grind_tacMnul_pb{4});
      tmpd7.dimord='chan_time';
      tmpd7.avg=squeeze(mean(tmpd7.individual,1));
      tmpd7=rmfield(tmpd7,'individual');
      
      tmpmask=stat_tacMnul_pb.mask;
      
  end
  
  
  
  
  if timwinstatflag==0
    tmpu1.mask=zeros(size(tmpu1.avg,1),length(tmpu1.time));
    tmpu1.mask(:,dsearchn(tmpu1.time',stattimwin(1)):dsearchn(tmpu1.time',stattimwin(end)))=tmpmask;
    tmpm10.mask=zeros(size(tmpm10.avg,1),length(tmpm10.time));
    tmpm10.mask(:,dsearchn(tmpm10.time',stattimwin(1)):dsearchn(tmpm10.time',stattimwin(end)))=tmpmask;
    tmpd4.mask=zeros(size(tmpd4.avg,1),length(tmpd4.time));
    tmpd4.mask(:,dsearchn(tmpd4.time',stattimwin(1)):dsearchn(tmpd4.time',stattimwin(end)))=tmpmask;
    tmpd7.mask=zeros(size(tmpd7.avg,1),length(tmpd7.time));
    tmpd7.mask(:,dsearchn(tmpd7.time',stattimwin(1)):dsearchn(tmpd7.time',stattimwin(end)))=tmpmask;
  else
    tmpu1.mask=tmpmask;
    tmpm10.mask=tmpmask;
    tmpd4.mask=tmpmask;
    tmpd7.mask=tmpmask;
  end
  
  
  for cg=1:length(chanplot)
    cfg=[];
    cfg.parameter='avg';
    cfg.layout='elec1010.lay';
    cfg.ylim=[-5 8];
    cfg.linewidth=3;
    cfg.xlim=timwin;
    if cg>length(chanplot)
      cfg.channel=tmpd4.label(any(tmpd4.mask,2));
    else
      cfg.channel=chanplot{cg};
    end
    cfg.graphcolor=coloruse([1 4 10 7],:);
    cfg.interactive='no';
    cfg.maskparameter='mask';
    cfg.maskstyle='box'; % default
    figure(ll*100+10*(cg+1))
    ft_singleplotER(cfg,tmpu1,tmpd7,tmpm10,tmpd4);
    hold on;plot(tmpu1.time,0,'k');
    set(gca,'XTick',[-.5:.1:1])
    set(gca,'XTickLabel',{'-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'})
    set(gca,'FontSize',30)
    title([])
    plot([stattimwin(1) stattimwin(1)],cfg.ylim,'k--','LineWidth',6)
    plot([stattimwin(2) stattimwin(2)],cfg.ylim,'k--','LineWidth',6)
    if ll<10 || ll==12
      plot([0 0],cfg.ylim,'Color',coloruse(4,:),'LineWidth',6)
    end
    if ll<10
      plot([soades(ll) soades(ll)],cfg.ylim,'Color',coloruse(9,:),'LineWidth',6)
    elseif ll==11
      plot([soades(5) soades(5)],cfg.ylim,'Color',coloruse(9,:),'LineWidth',6)
    end
    axis([-0.55 1 cfg.ylim(1) cfg.ylim(2)])
    if cg==3
      legend('Sum Unisensory','MultSens + Null','SumUnisens - MultsensNull')
    end
  end
  print(ll*100+20,[fdir 'erp_' condname{ll} 'vsNul_diff_PeakVsTrgh_FC_sleep' num2str(sleep) '_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  print(ll*100+30,[fdir 'erp_' condname{ll} 'vsNul_diff_PeakVsTrgh_OP_sleep' num2str(sleep) '_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  
  
end




for ll=5 % 1vs3 and 2vs4 (ll=5 only significant finding)
  cfg=[];
  if timwinstatflag==1
    cfg.latency=[stat_TPAmMSPN_peakMtrgh_pb_posthoctime{ll}.time(1) stat_TPAmMSPN_peakMtrgh_pb_posthoctime{ll}.time(end)];
  elseif timwinstatflag==0
    cfg.latency=timwin;
    stattimwin=[stat_TPAmMSPN_peakMtrgh_pb_posthoctime{ll}.time(1) stat_TPAmMSPN_peakMtrgh_pb_posthoctime{ll}.time(end)];
  end
  cfg.channel=stat_TPAmMSPN_peakMtrgh_pb_posthoctime{ll}.label;
  
  % FINISH ME HERE!
  
  tmpu1=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,1});
  tmpu1.dimord='chan_time';
  tmpu1.avg=squeeze(mean(tmpu1.individual,1));
  tmpu1=rmfield(tmpu1,'individual');
  
  tmpm10=ft_selectdata(cfg,grind_TPAmMSPN1_pb{ll,3});
  tmpm10.dimord='chan_time';
  tmpm10.avg=squeeze(mean(tmpm10.individual,1));
  tmpm10=rmfield(tmpm10,'individual');
  
  tmpd5=ft_selectdata(cfg,grind_TPAmMSPN_3m1_pb{ll});
  tmpd5.dimord='chan_time';
  tmpd5.avg=scalediff*squeeze(mean(tmpd5.individual,1));
  tmpd5=rmfield(tmpd5,'individual');
  
  
  tmpmask=stat_TPAmMSPN_peakMtrgh_pb_posthoctime{ll}.mask;
  if timwinstatflag==0
    tmpu1.mask=zeros(size(tmpu1.avg,1),length(tmpu1.time));
    tmpu1.mask(:,dsearchn(tmpu1.time',stattimwin(1)):dsearchn(tmpu1.time',stattimwin(end)))=tmpmask;
    tmpm10.mask=zeros(size(tmpm10.avg,1),length(tmpm10.time));
    tmpm10.mask(:,dsearchn(tmpm10.time',stattimwin(1)):dsearchn(tmpm10.time',stattimwin(end)))=tmpmask;
    tmpd5.mask=zeros(size(tmpd5.avg,1),length(tmpd5.time));
    tmpd5.mask(:,dsearchn(tmpd5.time',stattimwin(1)):dsearchn(tmpd5.time',stattimwin(end)))=tmpmask;
  else
    tmpu1{bb}.mask=tmpmask
    tmpm10{bb}.mask=tmpmask;
    tmpd5{bb}.mask=tmpmask;
  end
  
  for cg=1:length(chanplot)
    cfg=[];
    cfg.parameter='avg';
    cfg.layout='elec1010.lay';
    cfg.ylim=[-5 8];
    cfg.linewidth=3;
    cfg.xlim=timwin;
    if cg>length(chanplot)
      cfg.channel=tmpd5.label(any(tmpd5.mask,2));
    else
      cfg.channel=chanplot{cg};
    end
    cfg.graphcolor=coloruse([5 5 2],:);
    cfg.interactive='no';
    cfg.maskparameter='mask';
    cfg.maskstyle='box'; % default
    figure(ll+10*(cg+1))
    ft_singleplotER(cfg,tmpu1,tmpm10,tmpd5);
    hold on;plot(tmpu1.time,0,'k');
    set(gca,'XTick',[-.5:.1:1])
    set(gca,'XTickLabel',{'-0.5' '' ' ' '' ' ' '0' ' ' '' ' ' '' '0.5 ' '' ' ' ''  ' ' '1.0'})
    set(gca,'FontSize',30)
    title([])
    plot([stattimwin(1) stattimwin(1)],cfg.ylim,'k--','LineWidth',6)
    plot([stattimwin(2) stattimwin(2)],cfg.ylim,'k--','LineWidth',6)
    plot([0 0],cfg.ylim,'Color',coloruse(4,:),'LineWidth',6)
    plot([soades(ll) soades(ll)],cfg.ylim,'Color',coloruse(9,:),'LineWidth',6)
    axis([-0.55 1 cfg.ylim(1) cfg.ylim(2)])
    if cg==3
      legend('Sum Unisensory','MultSens + Null','SumUnisens - MultsensNull')
    end
  end
  %   print(30+ll,[fdir 'erp_tacPaud_MSpN_diff_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  print(ll+20,[fdir 'erp_tacPaud_MSpN_diff_PeakVsTrgh_FC_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  print(ll+30,[fdir 'erp_tacPaud_MSpN_diff_PeakVsTrgh_OP_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  
  if any(tmpmask(:))
    masktime=find(any(tmpd5.mask,1));
    cfg=[];
    cfg.parameter='avg';
    cfg.layout='elec1010.lay';
    cfg.maskalpha=0.5;
    cfg.zlim=[-5 5];
    cfg.highlight='on';
    cfg.highlightsize=12;
    cfg.xlim=[tmpd5.time(masktime(1)) tmpd5.time(masktime(end))];
    cfg.comment='no';
    sigchannels=tmpd5.label(find(ceil(mean(tmpd5.mask(:,dsearchn(tmpd5.time',cfg.xlim(1)):dsearchn(tmpd5.time',cfg.xlim(2))),2))));
    cfg.highlightchannel=sigchannels;
    figure(100+ll);
    ft_topoplotER(cfg,tmpu1);
    print(100+ll,[fdir 'erp_topoU_PeakVsTrgh_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    figure(100+10+ll);
    ft_topoplotER(cfg,tmpm10);
    print(100+10+ll,[fdir 'erp_topoM_PeakVsTrgh_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    figure(100+20+ll);
    ft_topoplotER(cfg,tmpd5);
    print(100+20+ll,[fdir 'erp_topoDiff_PeakVsTrgh_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  end
  
end

