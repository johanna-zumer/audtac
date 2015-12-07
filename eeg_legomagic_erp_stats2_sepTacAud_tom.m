% statistics of EEG awake data

clear all
close all
if ispc
  edir='F:\JohannaProject\jz_sleep_data\data\';%edir='D:\audtac\eeg_data\';
  ddir='G:\diaries\';%ddir='E:\JohannaProject\diaries\';%ddir='D:\audtac\legomagic\diaries\';
  bdir='D:\audtac\behav_data\';
  sdir='D:\audtac\spss_stuff\';
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
load('F:\JohannaProject\jz_sleep_data\data\iikeep.mat');
%load([edir 'iikeep.mat'])
soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];

%%

ttuse=[3];
% for ii=[iiuse]
for sleep=[1]
  %     for ii=[22 25 29 30 32]
  %   for ii=setdiff(union(iiSuse,iiBuse),[3:7])
  for ii= setdiff(iiBuse,[3:7 23 27])
    cd([edir sub{ii} ])
    clearvars -except ii sub edir ddir ii*use sleep featurestats* ttuse soades
    [raw_tac, raw_aud, raw_nul, artflag]=eeg_legomagic_epoching2(ii,sleep,1,0); % featfull=1, saveflag =0;
    %load(['raw3_each_rej_' sub{ii} '_sleep' num2str(sleep)],'raw_*'); % REMOVE THIS LINE
    
    %   try
    %     load(['tlock_trialSel_' sub{ii} '.mat']);
    %   catch
    
    %     [tlock_tac_s0,tlock_aud_s0,tlock_nul_s0,tr]=eeg_legomagic_trialSelection1(ii);
    %     save(['tlock_trialSel_' sub{ii} '.mat'],'tlock*s0','tr','-v7.3');
    %     [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection1_wakeSleep(ii,sleep);
    
    
    % For memory reasons, do TacPlusAud separately from AudPlusTac
    %% Do TacPlusAud first
    tacaud=1; % tacaud=1 means triggered on tactile; tacaud=0 means triggered on auditory
    for tt=ttuse % refers to lim=[no-limit .01 .005 .002];
      %   profile on
      %       [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection2_wakeSleep_sepTacAud(ii,sleep,tt,tacaud);
      [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection2_wakeSleep_sepTacAud(ii,sleep,tt,tacaud,raw_tac, raw_aud, raw_nul);
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
          fsample=1000;% 1/diff(tlock_tac{10,tt,ss}.time(1:2)); % sampling rate
        catch
          try
            fsample=1/diff(tlock_aud{10,tt,ss}.time(1:2)); % sampling rate
          catch
            fsample=1/diff(tlock_nul{10,tt,ss}.time(1:2)); % sampling rate
          end
        end
        %   for tt=1:4
        
        % Multisensory-shifted comparison
%         for ll=[1 3 4] % need to do this before other loops as the fields get deleted as it loops through ll
%           if ss==23
%             if length(tr.tmstrialkept1{ll,tt,12})>=2 && length(tr.tmstrialkept1{ll,tt,13})<2
%               cfg=[];
%               cfg.trials=tr.tmstrialkept1{ll,tt,12};
%               tlock_tacMSshift1{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{ll,tt,12});
%               cfg=[];
%               cfg.trials=tr.tmstrialkept2{ll,tt,12};
%               tlock_tacMSshift2{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{10-ll,tt,12});
%               cfg=[];
%               cfg.trials=tr.tmstrialkept3{ll,tt,12};
%               tlock_tacMSshift3{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{5,tt,12});
%               cfg=[];
%               cfg.trials=tr.tmstrialkept4{ll,tt,12};
%               tlock_tacMSshift4{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{5,tt,12});
%             elseif length(tr.tmstrialkept1{ll,tt,12})<2 && length(tr.tmstrialkept1{ll,tt,13})>=2
%               cfg=[];
%               cfg.trials=tr.tmstrialkept1{ll,tt,13};
%               tlock_tacMSshift1{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{ll,tt,13});
%               cfg=[];
%               cfg.trials=tr.tmstrialkept2{ll,tt,13};
%               tlock_tacMSshift2{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{10-ll,tt,13});
%               cfg=[];
%               cfg.trials=tr.tmstrialkept3{ll,tt,13};
%               tlock_tacMSshift3{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{5,tt,13});
%               cfg=[];
%               cfg.trials=tr.tmstrialkept4{ll,tt,13};
%               tlock_tacMSshift4{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{5,tt,13});
%             elseif length(tr.tmstrialkept1{ll,tt,12})>=2 && length(tr.tmstrialkept1{ll,tt,13})>=2
%               cfg=[];
%               cfg.trials=tr.tmstrialkept1{ll,tt,12};
%               tmp12=ft_selectdata(cfg,tlock_tac{ll,tt,12});
%               cfg=[];
%               cfg.trials=tr.tmstrialkept1{ll,tt,13};
%               tmp13=ft_selectdata(cfg,tlock_tac{ll,tt,13});
%               cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
%               cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
%               tlock_tacMSshift1{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
%               cfg=[];
%               cfg.trials=tr.tmstrialkept2{ll,tt,12};
%               tmp12=ft_selectdata(cfg,tlock_tac{10-ll,tt,12});
%               cfg=[];
%               cfg.trials=tr.tmstrialkept2{ll,tt,13};
%               tmp13=ft_selectdata(cfg,tlock_tac{10-ll,tt,13});
%               cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
%               cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
%               tlock_tacMSshift2{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
%               cfg=[];
%               cfg.trials=tr.tmstrialkept3{ll,tt,12};
%               tmp12=ft_selectdata(cfg,tlock_tac{5,tt,12});
%               cfg=[];
%               cfg.trials=tr.tmstrialkept3{ll,tt,13};
%               tmp13=ft_selectdata(cfg,tlock_tac{5,tt,13});
%               cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
%               cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
%               tlock_tacMSshift3{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
%               cfg=[];
%               cfg.trials=tr.tmstrialkept4{ll,tt,12};
%               tmp12=ft_selectdata(cfg,tlock_tac{5,tt,12});
%               cfg=[];
%               cfg.trials=tr.tmstrialkept4{ll,tt,13};
%               tmp13=ft_selectdata(cfg,tlock_tac{5,tt,13});
%               cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
%               cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
%               tlock_tacMSshift4{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
%             else
%               %               tlock_tac{ll,tt,ss}=[];
%               continue
%             end
%           else
%             if length(tr.tmstrialkept1{ll,tt,ss})>=2 && length(tr.tmstrialkept2{ll,tt,ss})>=2 && length(tr.tmstrialkept3{ll,tt,ss})>=2 && length(tr.tmstrialkept4{ll,tt,ss})>=2
%               cfg=[];
%               cfg.trials=tr.tmstrialkept1{ll,tt,ss};
%               tlock_tacMSshift1{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{ll,tt,ss});
%               cfg=[];
%               cfg.trials=tr.tmstrialkept2{ll,tt,ss};
%               tlock_tacMSshift2{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{10-ll,tt,ss});
%               cfg=[];
%               cfg.trials=tr.tmstrialkept3{ll,tt,ss};
%               tlock_tacMSshift3{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{5,tt,ss});
%               cfg=[];
%               cfg.trials=tr.tmstrialkept4{ll,tt,ss};
%               tlock_tacMSshift4{ll,tt,ss}=ft_selectdata(cfg,tlock_tac{5,tt,ss});
%             else
%               tlock_tacMSshift1{ll,tt,ss}=[];
%               tlock_tacMSshift2{ll,tt,ss}=[];
%               tlock_tacMSshift3{ll,tt,ss}=[];
%               tlock_tacMSshift4{ll,tt,ss}=[];
%             end
%           end
%           % do time-shifting here
%           if ~isempty(tlock_tacMSshift1{ll,tt,ss})
%             warning off
%             cfg=[];
%             cfg.offset=round(fsample*(-soades(ll)));
%             tlock_tacMSshift1{ll,tt,ss}=ft_redefinetrial(cfg,tlock_tacMSshift1{ll,tt,ss}); % Aud first shifted so that Aud at time 0
%             tlock_tacMSshift4{ll,tt,ss}=ft_redefinetrial(cfg,tlock_tacMSshift4{ll,tt,ss}); % Simult shifted so that both at time of second stim of 'll' condition
%             warning on
%           end
%         end %ll
%         
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
                continue
              end
%               if ll==max(soalist)
%                 tlock_tac{10,tt,12}=[];
%                 tlock_tac{10,tt,13}=[];
%               end
              
              %               cfg=[];
              %               cfg.trials=tr.a10trialkept{ll,tt,12};
              %               tmp12=ft_selectdata(cfg,tlock_aud{10,tt,12}); % create +20 here as 10 will be different for every ll,tt
              %               cfg.trials=tr.a10trialkept{ll,tt,13};
              %               tmp13=ft_selectdata(cfg,tlock_aud{10,tt,13}); % create +20 here as 10 will be different for every ll,tt
              %               cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
              %               cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
              %               tlock_aud{ll+20,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
%               if ll==max(soalist)
%                 tlock_aud{10,tt,12}=[];
%                 tlock_aud{10,tt,13}=[];
%               end
              
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
              
              
%               if ll==max(soalist) %% tw commented out - this was added to
%               save memory by JZ originally
%                 tlock_nul{10,tt,12}=[];
%                 tlock_nul{10,tt,13}=[];
%               end
%               
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
%                 if ll==max(soalist) && ss<12
%                   tlock_aud{10,tt,ss}=[];
%                 end
                
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
     %         tlock_tac{ll,tt,ss}=[];
              continue
            end
%             tlock_tac{ll,tt,12}=[];
%             tlock_tac{ll,tt,13}=[];
%             tlock_aud{ll,tt,12}=[];
%             tlock_aud{ll,tt,13}=[];
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
     %         tlock_aud{ll,tt,ss}=[];
              continue
            end
            
%             tlock_aud{ll,tt,12}=[];
%             tlock_aud{ll,tt,13}=[];
%             tlock_tac{ll,tt,12}=[];
%             tlock_tac{ll,tt,13}=[];
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
            else
              tlock_tactlock{ll,tt,ss}=[];
            end
          else % we for sure do not need ll>40 for tac here; it shouldn't even exist actually
 %           tlock_tac{ll,tt,ss}=[]; % clearing to help save memory as running
          end
%           if ss==12 || ss==13
%           else
%             tlock_tac{ll,tt,ss}=[]; % clearing to help save memory as running %% TW commented out!!
%           end
          
          if ll>40 % can we clear ll<40 earlier?
            if isfield(tlock_aud{ll,tt,ss},'trial') && size(tlock_aud{ll,tt,ss}.trial,1)
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
            else
              tlock_audtlock{ll,tt,ss}=[];
            end
            %           elseif ll<40 && ss~=12 && ss~=13
            %             tlock_aud{ll,tt,ss}=[];
          end
          if ss==12 || ss==13
          else
   %         tlock_aud{ll,tt,ss}=[];
          end
          %           catch
          %             disp('didnt work for this ll')
          %           end
        end % end ll
        

 
     
        
        
%         %         for ll=[soalist+50 soalist+60]
%         for ll=[soalist+50 ]
%           if length(tr.nllttrialkept{ll-50,tt,ss})>=2 && isfield(tlock_nul{ll,tt,ss},'trial') && size(tlock_nul{ll,tt,ss}.trial,1)
%             cfg=[];
%             cfg.lpfilter='yes';
%             cfg.lpfreq=40;
%             cfg.demean='yes';
%             %             cfg.baselinewindow=[-1.7 -0.6];
%             if ll==51
%               cfg.baselinewindow=[-.4 0]-.6;
%             elseif ll==53
%               cfg.baselinewindow=[-.4 0]-.17;
%             elseif ll==54
%               cfg.baselinewindow=[-.4 0]-.12;
%             elseif ll==55 || ll==56 || ll==57 || ll==59
%               cfg.baselinewindow=[-.4 0]-.1;
%             end
%             %           cfg.bpfilter='yes';
%             %           cfg.bpfreq=[1 40];
%             disp('ft_preprocessing nul')
%             tlock_nul{ll,tt,ss}=ft_preprocessing(cfg,tlock_nul{ll,tt,ss});
%             featurestats_nul(:,ll,tt,ss,ii)=mean(tlock_nul{ll,tt,ss}.trialinfo(:,12:19));
%             cfg=[];
%             cfg.covariance='yes';
%             cfg.covariancewindow=[-0.6 0.7]; % a full window valid for all ll
%             tlock_nultlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_nul{ll,tt,ss});
%           else
%    %         tlock_nultlock{ll,tt,ss}=[];
%           end
%     %      tlock_nul{ll,tt,ss}=[];
%         end
        
        
%         for ll=soalist
%           % create sum of unisensory conditions
%           cfg=[];
%           cfg.operation='add';
%           cfg.parameter='avg';
%           %           cfg.parameter={'avg' 'cov'};
%           % the call to ft_timelockanalysis is b/c .avg field not present after ft_redefinetrial
%           if numt_trials(ll,tt,ss)>=2
%             tlock_tacPaud{ll,tt,ss}=ft_math(cfg,tlock_tactlock{ll+20,tt,ss},tlock_audtlock{ll+40,tt,ss}); % keeping tac centred but jitter aud (to match tlock_tac{X,tt})
%             cfg.parameter='cov';
%             tmp=ft_math(cfg,tlock_tactlock{ll+20,tt,ss},tlock_audtlock{ll+40,tt,ss}); % keeping tac centred but jitter aud (to match tlock_tac{X,tt})
%             tlock_tacPaud{ll,tt,ss}.cov=tmp.cov;
%             featurestats_tacPaud(:,ll,tt,ss,ii)=mean([featurestats_tac(:,ll+20,tt,ss,ii), featurestats_aud(:,ll+40,tt,ss,ii)]');
%             tlock_tTacAlone{ll,tt,ss}=tlock_tactlock{ll+20,tt,ss}; % keeping tac centred but jitter aud (to match tlock_tac{X,tt})
%             tlock_tAudAlone{ll,tt,ss}=tlock_audtlock{ll+40,tt,ss}; % keeping tac centred but jitter aud (to match tlock_tac{X,tt})
%           else
%             tlock_tacPaud{ll,tt,ss}=[];
%             tlock_tTacAlone{ll,tt,ss}=[]; % keeping tac centred but jitter aud (to match tlock_tac{X,tt})
%             tlock_tAudAlone{ll,tt,ss}=[]; % keeping tac centred but jitter aud (to match tlock_tac{X,tt})
%           end
%           %           if numa_trials(ll,tt,ss)
%           %             tlock_audPtac{ll,tt,ss}=ft_math(cfg,tlock_audtlock{ll+20,tt,ss},tlock_tactlock{ll+40,tt,ss}); % keeping aud centred but jitter tac (to match tlock_aud{X,tt})
%           %             featurestats_audPtac(:,ll,tt,ss,ii)=mean([featurestats_aud(:,ll+20,tt,ss,ii), featurestats_tac(:,ll+40,tt,ss,ii)]');
%           %           else
%           %             tlock_audPtac{ll,tt,ss}=[];
%           %           end
%           
%           % create sum of multisensory plus nul conditions
%           cfg=[];
%           cfg.operation='add';
%           cfg.parameter='avg';
%           if numt_trials(ll,tt,ss)>=2
%             tlock_tacMSpN{ll,tt,ss}=ft_math(cfg,tlock_tactlock{ll,tt,ss},tlock_nultlock{ll+50,tt,ss}); % TA+N
%             cfg.parameter='cov';
%             tmp=ft_math(cfg,tlock_tactlock{ll,tt,ss},tlock_nultlock{ll+50,tt,ss}); % TA+N
%             tlock_tacMSpN{ll,tt,ss}.cov=tmp.cov;
%             featurestats_tacMSpN(:,ll,tt,ss,ii)=mean([featurestats_tac(:,ll,tt,ss,ii), featurestats_nul(:,ll+50,tt,ss,ii)]');
%             tlock_tMSAlone{ll,tt,ss}=tlock_tactlock{ll,tt,ss}; % TA+N
%             tlock_tNulAlone{ll,tt,ss}=tlock_nultlock{ll+50,tt,ss}; % TA+N
%           else
%             tlock_tacMSpN{ll,tt,ss}=[];
%             tlock_tMSAlone{ll,tt,ss}=[]; % TA+N
%             tlock_tNulAlone{ll,tt,ss}=[]; % TA+N
%           end
%           %           if numa_trials(ll,tt,ss)
%           %             tlock_audMSpN{ll,tt,ss}=ft_math(cfg,tlock_audtlock{ll,tt,ss},tlock_nultlock{ll+60,tt,ss}); % AT+N
%           %             featurestats_audMSpN(:,ll,tt,ss,ii)=mean([featurestats_aud(:,ll,tt,ss,ii), featurestats_nul(:,ll+60,tt,ss,ii)]');
%           %           else
%           %             tlock_audMSpN{ll,tt,ss}=[];
%           %           end
%         end % end ll
        
%% UTA - congruency analysis        
%         for ll=[1 3 4]
%           if ~isempty(tlock_tacMSshift1{ll,tt,ss}) && size(tlock_tacMSshift1{ll,tt,ss}.trial,1)
%             cfg=[];
%             cfg.lpfilter='yes';
%             cfg.lpfreq=40;
%             cfg.demean='yes';
%             cfg.baselinewindow=[-.4 0]-.1; % can use same baseline for all, as shifting has already occured
%             %               if ll==1
%             %                 cfg.baselinewindow=[-.4 0]-.6;
%             %               elseif ll==3
%             %                 cfg.baselinewindow=[-.4 0]-.17;
%             %               elseif ll==4
%             %                 cfg.baselinewindow=[-.4 0]-.12;
%             %               elseif ll==5
%             %                 cfg.baselinewindow=[-.4 0]-.1;
%             %               end
%             
%             disp('ft_preprocessing tac')
%             tlock_tacMSshift1{ll,tt,ss}=ft_preprocessing(cfg,tlock_tacMSshift1{ll,tt,ss});
%             tlock_tacMSshift2{ll,tt,ss}=ft_preprocessing(cfg,tlock_tacMSshift2{ll,tt,ss});
%             tlock_tacMSshift3{ll,tt,ss}=ft_preprocessing(cfg,tlock_tacMSshift3{ll,tt,ss});
%             tlock_tacMSshift4{ll,tt,ss}=ft_preprocessing(cfg,tlock_tacMSshift4{ll,tt,ss});
%             featurestats_tacMSshift1(:,ll,tt,ss,ii)=mean(tlock_tacMSshift1{ll,tt,ss}.trialinfo(:,12:19));
%             featurestats_tacMSshift2(:,ll,tt,ss,ii)=mean(tlock_tacMSshift2{ll,tt,ss}.trialinfo(:,12:19));
%             featurestats_tacMSshift3(:,ll,tt,ss,ii)=mean(tlock_tacMSshift3{ll,tt,ss}.trialinfo(:,12:19));
%             featurestats_tacMSshift4(:,ll,tt,ss,ii)=mean(tlock_tacMSshift4{ll,tt,ss}.trialinfo(:,12:19));
%             cfg=[];
%             cfg.covariance='yes';
%             cfg.covariancewindow=[-0.1 0.7]; % a full window valid for all ll
%             tlock_tacMSshift1tlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_tacMSshift1{ll,tt,ss});
%             tlock_tacMSshift2tlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_tacMSshift2{ll,tt,ss});
%             tlock_tacMSshift3tlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_tacMSshift3{ll,tt,ss});
%             tlock_tacMSshift4tlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_tacMSshift4{ll,tt,ss});
%           else
%             tlock_tacMSshift1tlock{ll,tt,ss}=[];
%             tlock_tacMSshift2tlock{ll,tt,ss}=[];
%             tlock_tacMSshift3tlock{ll,tt,ss}=[];
%             tlock_tacMSshift4tlock{ll,tt,ss}=[];
%           end
%           
%           if ~isempty(tlock_tacMSshift1tlock{ll,tt,ss})
%             % create sum of asynchronous multisensory conditions (e.g. AT70_shifted plus TA70)
%             cfg=[];
%             cfg.operation='add';
%             cfg.parameter='avg';
%             tlock_tMSasynch{ll,tt,ss}=ft_math(cfg,tlock_tacMSshift1tlock{ll,tt,ss},tlock_tacMSshift2tlock{ll,tt,ss});
%             cfg.parameter='cov';
%             tmp=ft_math(cfg,tlock_tacMSshift1tlock{ll,tt,ss},tlock_tacMSshift2tlock{ll,tt,ss});
%             tlock_tMSasynch{ll,tt,ss}.cov=tmp.cov;
%             featurestats_tMSasynch(:,ll,tt,ss,ii)=mean([featurestats_tacMSshift1(:,ll,tt,ss,ii), featurestats_tacMSshift2(:,ll,tt,ss,ii)]');
%             % create sum of synchronous multisensory conditions (e.g. TA0_shifted plus TA0)
%             cfg=[];
%             cfg.operation='add';
%             cfg.parameter='avg';
%             tlock_tMSsynch{ll,tt,ss}=ft_math(cfg,tlock_tacMSshift3tlock{ll,tt,ss},tlock_tacMSshift4tlock{ll,tt,ss});
%             cfg.parameter='cov';
%             tmp=ft_math(cfg,tlock_tacMSshift3tlock{ll,tt,ss},tlock_tacMSshift4tlock{ll,tt,ss});
%             tlock_tMSsynch{ll,tt,ss}.cov=tmp.cov;
%             featurestats_tMSsynch(:,ll,tt,ss,ii)=mean([featurestats_tacMSshift3(:,ll,tt,ss,ii), featurestats_tacMSshift4(:,ll,tt,ss,ii)]');
%           else
%             tlock_tMSasynch{ll,tt,ss}=[];
%             tlock_tMSsynch{ll,tt,ss}=[];
%           end
%         end % ll
       
        %% section added by tw
        for ll=1:101 % should it be 40?
     try
                if isfield(tlock_tac{ll,tt,ss},'trial') && size(tlock_tac{ll,tt,ss}.trial,1)
                cfg_alpha=[];
                cfg_alpha.bpfilter='yes';
                cfg_alpha.bpfreq=[8 12];
                cfg_alpha.demean='yes';
                tlock_tac{ll,tt,ss}.alpha=[];
                tlock_tac{ll,tt,ss}.alpha=ft_preprocessing(cfg_alpha,tlock_tac{ll,tt,ss});
                %cfg_alpha.baselinewindow=cfg.baselinewindow
                    for nTrials=1:length(tlock_tac{ll,tt,ss}.trial(:,1,1))
                        tmp.bp=[];
                        tmp.bp(1:length(tlock_tac{ll,tt,ss}.label),1:length(tlock_tac{ll,tt,ss}.trial))=tlock_tac{ll,tt,ss}.alpha.trial(nTrials,:,:);
                        tlock_tac{ll,tt,ss}.alpha.abs(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'abs');
                        tlock_tac{ll,tt,ss}.alpha.complex(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'complex');
                        tlock_tac{ll,tt,ss}.alpha.real(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'real');
                        tlock_tac{ll,tt,ss}.alpha.absreal(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'absreal');
                        tlock_tac{ll,tt,ss}.alpha.absimag(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'absimag');
                        tlock_tac{ll,tt,ss}.alpha.angle(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'angle');

                        tlock_tac{ll,tt,ss}.hilbert.alphaT0.abs(nTrials,:)=tlock_tac{ll,tt,ss}.alpha.abs(:,find(abs(tlock_tac{ll,tt,ss}.time)==min(abs(tlock_tac{ll,tt,ss}.time))),nTrials);
                        % calculate abs (amplitude envelope over 1000 ms pre-stimulus)
                        tlock_tac{ll,tt,ss}.hilbert.alphaT0.absMean(nTrials,:)= ...
                        mean(tlock_tac{ll,tt,ss}.alpha.abs(:,find(abs(tlock_tac{ll,tt,ss}.time)==min(abs(tlock_tac{ll,tt,ss}.time))):find(abs(tlock_tac{ll,tt,ss}.time)==min(abs(tlock_tac{ll,tt,ss}.time)))+1000,nTrials));
                        
                        tlock_tac{ll,tt,ss}.hilbert.alphaT0.complex(nTrials,:)=tlock_tac{ll,tt,ss}.alpha.complex(:,find(abs(tlock_tac{ll,tt,ss}.time)==min(abs(tlock_tac{ll,tt,ss}.time))),nTrials);
                        tlock_tac{ll,tt,ss}.hilbert.alphaT0.real(nTrials,:)=tlock_tac{ll,tt,ss}.alpha.real(:,find(abs(tlock_tac{ll,tt,ss}.time)==min(abs(tlock_tac{ll,tt,ss}.time))),nTrials);
                        tlock_tac{ll,tt,ss}.hilbert.alphaT0.absreal(nTrials,:)=tlock_tac{ll,tt,ss}.alpha.absreal(:,find(abs(tlock_tac{ll,tt,ss}.time)==min(abs(tlock_tac{ll,tt,ss}.time))),nTrials);
                        tlock_tac{ll,tt,ss}.hilbert.alphaT0.absimag(nTrials,:)=tlock_tac{ll,tt,ss}.alpha.absimag(:,find(abs(tlock_tac{ll,tt,ss}.time)==min(abs(tlock_tac{ll,tt,ss}.time))),nTrials);
                        tlock_tac{ll,tt,ss}.hilbert.alphaT0.angle(nTrials,:)=tlock_tac{ll,tt,ss}.alpha.angle(:,find(abs(tlock_tac{ll,tt,ss}.time)==min(abs(tlock_tac{ll,tt,ss}.time))),nTrials);
                    end
                    
                    % and for theta
                    cfg_theta=[];
                    cfg_theta.bpfilter='yes';
                    cfg_theta.bpfreq=[4 8];
                    cfg_theta.demean='yes';
                    tlock_tac{ll,tt,ss}.theta=[];
                    tlock_tac{ll,tt,ss}.theta=ft_preprocessing(cfg_theta,tlock_tac{ll,tt,ss});
                    %cfg_theta.baselinewindow=cfg.baselinewindow
                        for nTrials=1:length(tlock_tac{ll,tt,ss}.trial(:,1,1))
                            tmp.bp=[];
                            tmp.bp(1:length(tlock_tac{ll,tt,ss}.label),1:length(tlock_tac{ll,tt,ss}.trial))=tlock_tac{ll,tt,ss}.theta.trial(nTrials,:,:);
                            tlock_tac{ll,tt,ss}.theta.abs(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'abs');
                            tlock_tac{ll,tt,ss}.theta.complex(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'complex');
                            tlock_tac{ll,tt,ss}.theta.real(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'real');
                            tlock_tac{ll,tt,ss}.theta.absreal(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'absreal');
                            tlock_tac{ll,tt,ss}.theta.absimag(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'absimag');
                            tlock_tac{ll,tt,ss}.theta.angle(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'angle');

                            tlock_tac{ll,tt,ss}.hilbert.thetaT0.abs(nTrials,:)=tlock_tac{ll,tt,ss}.theta.abs(:,find(abs(tlock_tac{ll,tt,ss}.time)==min(abs(tlock_tac{ll,tt,ss}.time))),nTrials);
                            % calculate abs (amplitude envelope over 1000 ms pre-stimulus)
                            tlock_tac{ll,tt,ss}.hilbert.thetaT0.absMean(nTrials,:)= ...
                            mean(tlock_tac{ll,tt,ss}.theta.abs(:,find(abs(tlock_tac{ll,tt,ss}.time)==min(abs(tlock_tac{ll,tt,ss}.time))):find(abs(tlock_tac{ll,tt,ss}.time)==min(abs(tlock_tac{ll,tt,ss}.time)))+1000,nTrials));
                    
                            tlock_tac{ll,tt,ss}.hilbert.thetaT0.complex(nTrials,:)=tlock_tac{ll,tt,ss}.theta.complex(:,find(abs(tlock_tac{ll,tt,ss}.time)==min(abs(tlock_tac{ll,tt,ss}.time))),nTrials);
                            tlock_tac{ll,tt,ss}.hilbert.thetaT0.real(nTrials,:)=tlock_tac{ll,tt,ss}.theta.real(:,find(abs(tlock_tac{ll,tt,ss}.time)==min(abs(tlock_tac{ll,tt,ss}.time))),nTrials);
                            tlock_tac{ll,tt,ss}.hilbert.thetaT0.absreal(nTrials,:)=tlock_tac{ll,tt,ss}.theta.absreal(:,find(abs(tlock_tac{ll,tt,ss}.time)==min(abs(tlock_tac{ll,tt,ss}.time))),nTrials);
                            tlock_tac{ll,tt,ss}.hilbert.thetaT0.absimag(nTrials,:)=tlock_tac{ll,tt,ss}.theta.absimag(:,find(abs(tlock_tac{ll,tt,ss}.time)==min(abs(tlock_tac{ll,tt,ss}.time))),nTrials);
                            tlock_tac{ll,tt,ss}.hilbert.thetaT0.angle(nTrials,:)=tlock_tac{ll,tt,ss}.theta.angle(:,find(abs(tlock_tac{ll,tt,ss}.time)==min(abs(tlock_tac{ll,tt,ss}.time))),nTrials);
                        end
                end
     catch; end
     try
                if isfield(tlock_nul{ll,tt,ss},'trial') && size(tlock_nul{ll,tt,ss}.trial,1)
                cfg_alpha=[];
                cfg_alpha.bpfilter='yes';
                cfg_alpha.bpfreq=[8 12];
                cfg_alpha.demean='yes';
                tlock_nul{ll,tt,ss}.alpha=[];
                tlock_nul{ll,tt,ss}.alpha=ft_preprocessing(cfg_alpha,tlock_nul{ll,tt,ss});
                %cfg_alpha.baselinewindow=cfg.baselinewindow
                    for nTrials=1:length(tlock_nul{ll,tt,ss}.trial(:,1,1))
                        tmp.bp=[];
                        tmp.bp(1:length(tlock_nul{ll,tt,ss}.label),1:length(tlock_nul{ll,tt,ss}.trial))=tlock_nul{ll,tt,ss}.alpha.trial(nTrials,:,:);
                        tlock_nul{ll,tt,ss}.alpha.abs(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'abs');
                        tlock_nul{ll,tt,ss}.alpha.complex(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'complex');
                        tlock_nul{ll,tt,ss}.alpha.real(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'real');
                        tlock_nul{ll,tt,ss}.alpha.absreal(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'absreal');
                        tlock_nul{ll,tt,ss}.alpha.absimag(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'absimag');
                        tlock_nul{ll,tt,ss}.alpha.angle(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'angle');

                        tlock_nul{ll,tt,ss}.hilbert.alphaT0.abs(nTrials,:)=tlock_nul{ll,tt,ss}.alpha.abs(:,find(abs(tlock_nul{ll,tt,ss}.time)==min(abs(tlock_nul{ll,tt,ss}.time))),nTrials);
                        % calculate abs (amplitude envelope over 1000 ms pre-stimulus)
                        tlock_nul{ll,tt,ss}.hilbert.alphaT0.absMean(nTrials,:)= ...
                        mean(tlock_nul{ll,tt,ss}.alpha.abs(:,find(abs(tlock_nul{ll,tt,ss}.time)==min(abs(tlock_nul{ll,tt,ss}.time))):find(abs(tlock_nul{ll,tt,ss}.time)==min(abs(tlock_nul{ll,tt,ss}.time)))+1000,nTrials));
                    
                        tlock_nul{ll,tt,ss}.hilbert.alphaT0.complex(nTrials,:)=tlock_nul{ll,tt,ss}.alpha.complex(:,find(abs(tlock_nul{ll,tt,ss}.time)==min(abs(tlock_nul{ll,tt,ss}.time))),nTrials);
                        tlock_nul{ll,tt,ss}.hilbert.alphaT0.real(nTrials,:)=tlock_nul{ll,tt,ss}.alpha.real(:,find(abs(tlock_nul{ll,tt,ss}.time)==min(abs(tlock_nul{ll,tt,ss}.time))),nTrials);
                        tlock_nul{ll,tt,ss}.hilbert.alphaT0.absreal(nTrials,:)=tlock_nul{ll,tt,ss}.alpha.absreal(:,find(abs(tlock_nul{ll,tt,ss}.time)==min(abs(tlock_nul{ll,tt,ss}.time))),nTrials);
                        tlock_nul{ll,tt,ss}.hilbert.alphaT0.absimag(nTrials,:)=tlock_nul{ll,tt,ss}.alpha.absimag(:,find(abs(tlock_nul{ll,tt,ss}.time)==min(abs(tlock_nul{ll,tt,ss}.time))),nTrials);
                        tlock_nul{ll,tt,ss}.hilbert.alphaT0.angle(nTrials,:)=tlock_nul{ll,tt,ss}.alpha.angle(:,find(abs(tlock_nul{ll,tt,ss}.time)==min(abs(tlock_nul{ll,tt,ss}.time))),nTrials);
                    end

                    % and for theta
                    cfg_theta=[];
                    cfg_theta.bpfilter='yes';
                    cfg_theta.bpfreq=[4 8];
                    cfg_theta.demean='yes';
                    tlock_nul{ll,tt,ss}.theta=[];
                    tlock_nul{ll,tt,ss}.theta=ft_preprocessing(cfg_theta,tlock_nul{ll,tt,ss});
                    %cfg_theta.baselinewindow=cfg.baselinewindow
                        for nTrials=1:length(tlock_nul{ll,tt,ss}.trial(:,1,1))
                            tmp.bp=[];
                            tmp.bp(1:length(tlock_nul{ll,tt,ss}.label),1:length(tlock_nul{ll,tt,ss}.trial))=tlock_nul{ll,tt,ss}.theta.trial(nTrials,:,:);
                            tlock_nul{ll,tt,ss}.theta.abs(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'abs');
                            tlock_nul{ll,tt,ss}.theta.complex(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'complex');
                            tlock_nul{ll,tt,ss}.theta.real(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'real');
                            tlock_nul{ll,tt,ss}.theta.absreal(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'absreal');
                            tlock_nul{ll,tt,ss}.theta.absimag(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'absimag');
                            tlock_nul{ll,tt,ss}.theta.angle(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'angle');

                            tlock_nul{ll,tt,ss}.hilbert.thetaT0.abs(nTrials,:)=tlock_nul{ll,tt,ss}.theta.abs(:,find(abs(tlock_nul{ll,tt,ss}.time)==min(abs(tlock_nul{ll,tt,ss}.time))),nTrials);
                            % calculate abs (amplitude envelope over 1000 ms pre-stimulus)
                            tlock_nul{ll,tt,ss}.hilbert.thetaT0.absMean(nTrials,:)= ...
                            mean(tlock_nul{ll,tt,ss}.theta.abs(:,find(abs(tlock_nul{ll,tt,ss}.time)==min(abs(tlock_nul{ll,tt,ss}.time))):find(abs(tlock_nul{ll,tt,ss}.time)==min(abs(tlock_nul{ll,tt,ss}.time)))+1000,nTrials));
                    
                            tlock_nul{ll,tt,ss}.hilbert.thetaT0.complex(nTrials,:)=tlock_nul{ll,tt,ss}.theta.complex(:,find(abs(tlock_nul{ll,tt,ss}.time)==min(abs(tlock_nul{ll,tt,ss}.time))),nTrials);
                            tlock_nul{ll,tt,ss}.hilbert.thetaT0.real(nTrials,:)=tlock_nul{ll,tt,ss}.theta.real(:,find(abs(tlock_nul{ll,tt,ss}.time)==min(abs(tlock_nul{ll,tt,ss}.time))),nTrials);
                            tlock_nul{ll,tt,ss}.hilbert.thetaT0.absreal(nTrials,:)=tlock_nul{ll,tt,ss}.theta.absreal(:,find(abs(tlock_nul{ll,tt,ss}.time)==min(abs(tlock_nul{ll,tt,ss}.time))),nTrials);
                            tlock_nul{ll,tt,ss}.hilbert.thetaT0.absimag(nTrials,:)=tlock_nul{ll,tt,ss}.theta.absimag(:,find(abs(tlock_nul{ll,tt,ss}.time)==min(abs(tlock_nul{ll,tt,ss}.time))),nTrials);
                            tlock_nul{ll,tt,ss}.hilbert.thetaT0.angle(nTrials,:)=tlock_nul{ll,tt,ss}.theta.angle(:,find(abs(tlock_nul{ll,tt,ss}.time)==min(abs(tlock_nul{ll,tt,ss}.time))),nTrials);
                        end
                end
     catch; end
        end
        %%
        
      end
        eeg_legomagic_absKc
             
    
     % end ss
      
      
      %clear tlock_tac tlock_aud tlock_nul tr tlock*tlock
     % save(['tlock_hilbert_' sub{ii} '.mat'],'tlock*','tr','-v7.3');
%       save(['tlock_hilbert_tt3_to_save' sub{ii} '.mat'],'tlock_tac','tlock_nul','tr','-v7.3'); %
   %     save(['tlock_hilbert_tt3_to_save_3' sub{ii} '.mat'],'KC','Kc_Cdf', 'KC_ctr', '-v7.3'); %
         save(['tlock_hilbert_tt3_to_save_6_' sub{ii} '.mat'], 'abs*', '-v7.3'); %
      % save(['tlock_hilbert_tt3_to_save_4' sub{ii} '.mat'],'KC', '-v7.3'); %
      clear tlock_tac tlock_aud tlock_nul
      
    end % end tt
    
 

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
          fsample=1000;%/diff(tlock_tac{10,tt,ss}.time(1:2)); % sampling rate
        catch
          try
            fsample=1/diff(tlock_aud{10,tt,ss}.time(1:2)); % sampling rate
          catch
            fsample=1/diff(tlock_nul{10,tt,ss}.time(1:2)); % sampling rate
          end
        end
        
        %   for tt=1:4
        
%         % Multisensory-shifted comparison
%          for ll=[1 3 4] % need to do this before other loops as the fields get deleted as it loops through ll
%           if ss==23
%             if length(tr.amstrialkept1{ll,tt,12})>=2 && length(tr.amstrialkept1{ll,tt,13})<2
%               cfg=[];
%               cfg.trials=tr.amstrialkept1{ll,tt,12};
%               tlock_audMSshift1{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{ll,tt,12});
%               cfg=[];
%               cfg.trials=tr.amstrialkept2{ll,tt,12};
%               tlock_audMSshift2{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{10-ll,tt,12});
%               cfg=[];
%               cfg.trials=tr.amstrialkept3{ll,tt,12};
%               tlock_audMSshift3{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{5,tt,12});
%               cfg=[];
%               cfg.trials=tr.amstrialkept4{ll,tt,12};
%               tlock_audMSshift4{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{5,tt,12});
%             elseif length(tr.amstrialkept1{ll,tt,12})<2 && length(tr.amstrialkept1{ll,tt,13})>=2
%               cfg=[];
%               cfg.trials=tr.amstrialkept1{ll,tt,13};
%               tlock_audMSshift1{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{ll,tt,13});
%               cfg=[];
%               cfg.trials=tr.amstrialkept2{ll,tt,13};
%               tlock_audMSshift2{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{10-ll,tt,13});
%               cfg=[];
%               cfg.trials=tr.amstrialkept3{ll,tt,13};
%               tlock_audMSshift3{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{5,tt,13});
%               cfg=[];
%               cfg.trials=tr.amstrialkept4{ll,tt,13};
%               tlock_audMSshift4{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{5,tt,13});
%             elseif length(tr.amstrialkept1{ll,tt,12})>=2 && length(tr.amstrialkept1{ll,tt,13})>=2
%               cfg=[];
%               cfg.trials=tr.amstrialkept1{ll,tt,12};
%               tmp12=ft_selectdata(cfg,tlock_aud{ll,tt,12});
%               cfg=[];
%               cfg.trials=tr.amstrialkept1{ll,tt,13};
%               tmp13=ft_selectdata(cfg,tlock_aud{ll,tt,13});
%               cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
%               cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
%               tlock_audMSshift1{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
%               cfg=[];
%               cfg.trials=tr.amstrialkept2{ll,tt,12};
%               tmp12=ft_selectdata(cfg,tlock_aud{10-ll,tt,12});
%               cfg=[];
%               cfg.trials=tr.amstrialkept2{ll,tt,13};
%               tmp13=ft_selectdata(cfg,tlock_aud{10-ll,tt,13});
%               cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
%               cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
%               tlock_audMSshift2{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
%               cfg=[];
%               cfg.trials=tr.amstrialkept3{ll,tt,12};
%               tmp12=ft_selectdata(cfg,tlock_aud{5,tt,12});
%               cfg=[];
%               cfg.trials=tr.amstrialkept3{ll,tt,13};
%               tmp13=ft_selectdata(cfg,tlock_aud{5,tt,13});
%               cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
%               cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
%               tlock_audMSshift3{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
%               cfg=[];
%               cfg.trials=tr.amstrialkept4{ll,tt,12};
%               tmp12=ft_selectdata(cfg,tlock_aud{5,tt,12});
%               cfg=[];
%               cfg.trials=tr.amstrialkept4{ll,tt,13};
%               tmp13=ft_selectdata(cfg,tlock_aud{5,tt,13});
%               cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
%               cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
%               tlock_audMSshift4{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
%             else
%               %               tlock_aud{ll,tt,ss}=[];
%               continue
%             end
%           else
%             if length(tr.amstrialkept1{ll,tt,ss})>=2 && length(tr.amstrialkept2{ll,tt,ss})>=2 && length(tr.amstrialkept3{ll,tt,ss})>=2 && length(tr.amstrialkept4{ll,tt,ss})>=2
%               cfg=[];
%               cfg.trials=tr.amstrialkept1{ll,tt,ss};
%               tlock_audMSshift1{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{ll,tt,ss});
%               cfg=[];
%               cfg.trials=tr.amstrialkept2{ll,tt,ss};
%               tlock_audMSshift2{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{10-ll,tt,ss});
%               cfg=[];
%               cfg.trials=tr.amstrialkept3{ll,tt,ss};
%               tlock_audMSshift3{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{5,tt,ss});
%               cfg=[];
%               cfg.trials=tr.amstrialkept4{ll,tt,ss};
%               tlock_audMSshift4{ll,tt,ss}=ft_selectdata(cfg,tlock_aud{5,tt,ss});
%             else
%               tlock_audMSshift1{ll,tt,ss}=[];
%               tlock_audMSshift2{ll,tt,ss}=[];
%               tlock_audMSshift3{ll,tt,ss}=[];
%               tlock_audMSshift4{ll,tt,ss}=[];
%             end
%           end
%           % do time-shifting here
%           if ~isempty(tlock_audMSshift1{ll,tt,ss})
%             warning off
%             cfg=[];
%             cfg.offset=round(fsample*(-soades(ll)));
%             tlock_audMSshift2{ll,tt,ss}=ft_redefinetrial(cfg,tlock_audMSshift2{ll,tt,ss}); % Aud first shifted so that Aud at time 0
%             tlock_audMSshift4{ll,tt,ss}=ft_redefinetrial(cfg,tlock_audMSshift4{ll,tt,ss}); % Simult shifted so that both at time of second stim of 'll' condition
%             warning on
%           end
%         end %ll
        
        for ll=[soalist]% soalist+20 soalist+40] %%% CHANGED BY TOM %%%%
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
%                 tlock_tac{10,tt,12}=[];
%                 tlock_tac{10,tt,13}=[];
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
%                 tlock_aud{10,tt,12}=[];
%                 tlock_aud{10,tt,13}=[];
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
%                 tlock_nul{10,tt,12}=[];
%                 tlock_nul{10,tt,13}=[];
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
%                   tlock_tac{10,tt,ss}=[];
                end
                
                cfg=[];
                cfg.trials=tr.a10trialkept{ll,tt,ss};
                if length(cfg.trials)<2
                  tlock_aud{ll+20,tt,ss}=[]; % create +20 here as 10 will be different for every ll,tt
                else
                  tlock_aud{ll+20,tt,ss}=ft_selectdata(cfg,tlock_aud{10,tt,ss}); % create +20 here as 10 will be different for every ll,tt
                end
                if ll==max(soalist) && ss<12
%                   tlock_aud{10,tt,ss}=[];
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
%                   tlock_nul{10,tt,ss}=[];
                end
              end
              
            end
          end
          
          if ll<10
            %             numt_trials(ll,tt,ss)=size(tlock_tac{ll+20,tt,ss}.trial,1); % this should be the same for tacll20, audll40, tacll, nulll50
            if isempty(tlock_aud{ll+20,tt,ss})
%               numa_trials(ll,tt,ss)=0;
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
%             tlock_tac{ll,tt,12}=[];
%             tlock_tac{ll,tt,13}=[];
%             
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
%               tlock_aud{ll,tt,ss}=[];
              continue
            end
%             tlock_aud{ll,tt,12}=[];
%             tlock_aud{ll,tt,13}=[];
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
%               tlock_tac{ll,tt,ss}=[];
              continue
            end
%             tlock_tac{ll,tt,12}=[];
%             tlock_tac{ll,tt,13}=[];
            
            %             tmp12=tlock_aud{ll,tt,12};
            %             tmp13=tlock_aud{ll,tt,13};
            %             cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
            %             cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
            %             tlock_aud{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
%             tlock_aud{ll,tt,12}=[];
%             tlock_aud{ll,tt,13}=[];
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
%           % 40Hz lowpass filter (already 0.2Hz highpass on non-epoched data)
%           % also do baseline correction here
%           if ll>40
%             if isfield(tlock_tac{ll,tt,ss},'trial') && size(tlock_tac{ll,tt,ss}.trial,1)
%               cfg=[];
%               cfg.lpfilter='yes';
%               cfg.lpfreq=40;
%               cfg.demean='yes';
%               %               cfg.baselinewindow=[-1.7 -0.6];
%               if ll==49
%                 cfg.baselinewindow=[-.4 0]-.6;
%               elseif ll==47
%                 cfg.baselinewindow=[-.4 0]-.17;
%               elseif ll==46
%                 cfg.baselinewindow=[-.4 0]-.12;
%               elseif ll==45 || ll==44 || ll==43 || ll==41
%                 cfg.baselinewindow=[-.4 0]-.1;
%               end
%               disp('ft_preprocessing tac')
%               tlock_tac{ll,tt,ss}=ft_preprocessing(cfg,tlock_tac{ll,tt,ss});
%               featurestats_tac(:,ll,tt,ss,ii)=mean(tlock_tac{ll,tt,ss}.trialinfo(:,12:19));
%               cfg=[];
%               cfg.covariance='yes';
%               cfg.covariancewindow=[-0.6 0.7]; % a full window valid for all ll
%               tlock_tactlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
%             else
%               tlock_tactlock{ll,tt,ss}=[];
%             end
%           end
%           if ss==12 || ss==13
%           else
%          %   tlock_tac{ll,tt,ss}=[]; % clearing to help save memory as running
%           end
          
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
          %  tlock_aud{ll,tt,ss}=[];
          end
          %           catch
          %             disp('didnt work for this ll')
          %           end
        end % end ll
        
        
%         %         for ll=[soalist+50 soalist+60]
%         for ll=[soalist+60]
%           if length(tr.nllatrialkept{ll-60,tt,ss})>=2 && isfield(tlock_nul{ll,tt,ss},'trial') && size(tlock_nul{ll,tt,ss}.trial,1)
%             cfg=[];
%             cfg.lpfilter='yes';
%             cfg.lpfreq=40;
%             cfg.demean='yes';
%             %             cfg.baselinewindow=[-1.7 -0.6];
%             if ll==69
%               cfg.baselinewindow=[-.4 0]-.6;
%             elseif ll==67
%               cfg.baselinewindow=[-.4 0]-.17;
%             elseif ll==66
%               cfg.baselinewindow=[-.4 0]-.12;
%             elseif ll==65 || ll==64 || ll==63 || ll==61
%               cfg.baselinewindow=[-.4 0]-.1;
%             end
%             %           cfg.bpfilter='yes';
%             %           cfg.bpfreq=[1 40];
%             disp('ft_preprocessing nul')
%             tlock_nul{ll,tt,ss}=ft_preprocessing(cfg,tlock_nul{ll,tt,ss});
%             featurestats_nul(:,ll,tt,ss,ii)=mean(tlock_nul{ll,tt,ss}.trialinfo(:,12:19));
%             cfg=[];
%             cfg.covariance='yes';
%             cfg.covariancewindow=[-0.6 0.7]; % a full window valid for all ll
%             tlock_nultlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_nul{ll,tt,ss});
%           else
%             tlock_nultlock{ll,tt,ss}=[];
%           end
%           tlock_nul{ll,tt,ss}=[];
%         end
        
        
%         for ll=soalist
%           % create sum of unisensory conditions
%           cfg=[];
%           cfg.operation='add';
%           cfg.parameter='avg';
%           % the call to ft_timelockanalysis is b/c .avg field not present after ft_redefinetrial
%           %           if numt_trials(ll,tt,ss)
%           %             tlock_tacPaud{ll,tt,ss}=ft_math(cfg,tlock_tactlock{ll+20,tt,ss},tlock_audtlock{ll+40,tt,ss}); % keeping tac centred but jitter aud (to match tlock_tac{X,tt})
%           %             featurestats_tacPaud(:,ll,tt,ss,ii)=mean([featurestats_tac(:,ll+20,tt,ss,ii), featurestats_aud(:,ll+40,tt,ss,ii)]');
%           %           else
%           %             tlock_tacPaud{ll,tt,ss}=[];
%           %           end
%           if numa_trials(ll,tt,ss)>=2
%             tlock_audPtac{ll,tt,ss}=ft_math(cfg,tlock_audtlock{ll+20,tt,ss},tlock_tactlock{ll+40,tt,ss}); % keeping aud centred but jitter tac (to match tlock_aud{X,tt})
%             cfg.parameter='cov';
%             tmp=ft_math(cfg,tlock_audtlock{ll+20,tt,ss},tlock_tactlock{ll+40,tt,ss}); % keeping aud centred but jitter tac (to match tlock_aud{X,tt})
%             tlock_audPtac{ll,tt,ss}.cov=tmp.cov;
%             featurestats_audPtac(:,ll,tt,ss,ii)=mean([featurestats_aud(:,ll+20,tt,ss,ii), featurestats_tac(:,ll+40,tt,ss,ii)]');
%             tlock_aAudAlone{ll,tt,ss}=tlock_audtlock{ll+20,tt,ss}; % keeping aud centred but jitter tac (to match tlock_aud{X,tt})
%             tlock_aTacAlone{ll,tt,ss}=tlock_tactlock{ll+40,tt,ss}; % keeping aud centred but jitter tac (to match tlock_aud{X,tt})
%           else
%             tlock_audPtac{ll,tt,ss}=[];
%             tlock_aAudAlone{ll,tt,ss}=[]; % keeping aud centred but jitter tac (to match tlock_aud{X,tt})
%             tlock_aTacAlone{ll,tt,ss}=[]; % keeping aud centred but jitter tac (to match tlock_aud{X,tt})
%           end
%           
%           % create sum of multisensory plus nul conditions
%           cfg=[];
%           cfg.operation='add';
%           cfg.parameter='avg';
%           %           if numt_trials(ll,tt,ss)
%           %             tlock_tacMSpN{ll,tt,ss}=ft_math(cfg,tlock_tactlock{ll,tt,ss},tlock_nultlock{ll+50,tt,ss}); % TA+N
%           %             featurestats_tacMSpN(:,ll,tt,ss,ii)=mean([featurestats_tac(:,ll,tt,ss,ii), featurestats_nul(:,ll+50,tt,ss,ii)]');
%           %           else
%           %             tlock_tacMSpN{ll,tt,ss}=[];
%           %           end
%           if numa_trials(ll,tt,ss)>=2
%             tlock_audMSpN{ll,tt,ss}=ft_math(cfg,tlock_audtlock{ll,tt,ss},tlock_nultlock{ll+60,tt,ss}); % AT+N
%             cfg.parameter='cov';
%             tmp=ft_math(cfg,tlock_audtlock{ll,tt,ss},tlock_nultlock{ll+60,tt,ss}); % AT+N
%             tlock_audMSpN{ll,tt,ss}.cov=tmp.cov;
%             featurestats_audMSpN(:,ll,tt,ss,ii)=mean([featurestats_aud(:,ll,tt,ss,ii), featurestats_nul(:,ll+60,tt,ss,ii)]');
%             tlock_aMSAlone{ll,tt,ss}=tlock_audtlock{ll,tt,ss}; % AT+N
%             tlock_aNulAlone{ll,tt,ss}=tlock_nultlock{ll+60,tt,ss}; % AT+N
%           else
%             tlock_audMSpN{ll,tt,ss}=[];
%             tlock_aMSAlone{ll,tt,ss}=[]; % AT+N
%             tlock_aNulAlone{ll,tt,ss}=[]; % AT+N
%           end
%         end % end ll
%         
%         tlock_tacAll{tt,ss}=tlock_tac{100,tt,ss};
%         tlock_audAll{tt,ss}=tlock_aud{100,tt,ss};
%         tlock_tac19T{tt,ss}=tlock_tac{101,tt,ss};
%         tlock_aud19A{tt,ss}=tlock_aud{101,tt,ss};
%% UTA - congruency analysis        
%         for ll=[1 3 4]
%           if ~isempty(tlock_audMSshift1{ll,tt,ss})&& size(tlock_audMSshift1{ll,tt,ss}.trial,1)
%             cfg=[];
%             cfg.lpfilter='yes';
%             cfg.lpfreq=40;
%             cfg.demean='yes';
%             cfg.baselinewindow=[-.4 0]-.1; % can use same baseline for all, as shifting has already occured
%             %               if ll==1
%             %                 cfg.baselinewindow=[-.4 0]-.6;
%             %               elseif ll==3
%             %                 cfg.baselinewindow=[-.4 0]-.17;
%             %               elseif ll==4
%             %                 cfg.baselinewindow=[-.4 0]-.12;
%             %               elseif ll==5
%             %                 cfg.baselinewindow=[-.4 0]-.1;
%             %               end
%             
%             disp('ft_preprocessing aud')
%             tlock_audMSshift1{ll,tt,ss}=ft_preprocessing(cfg,tlock_audMSshift1{ll,tt,ss});
%             tlock_audMSshift2{ll,tt,ss}=ft_preprocessing(cfg,tlock_audMSshift2{ll,tt,ss});
%             tlock_audMSshift3{ll,tt,ss}=ft_preprocessing(cfg,tlock_audMSshift3{ll,tt,ss});
%             tlock_audMSshift4{ll,tt,ss}=ft_preprocessing(cfg,tlock_audMSshift4{ll,tt,ss});
%             featurestats_audMSshift1(:,ll,tt,ss,ii)=mean(tlock_audMSshift1{ll,tt,ss}.trialinfo(:,12:19));
%             featurestats_audMSshift2(:,ll,tt,ss,ii)=mean(tlock_audMSshift2{ll,tt,ss}.trialinfo(:,12:19));
%             featurestats_audMSshift3(:,ll,tt,ss,ii)=mean(tlock_audMSshift3{ll,tt,ss}.trialinfo(:,12:19));
%             featurestats_audMSshift4(:,ll,tt,ss,ii)=mean(tlock_audMSshift4{ll,tt,ss}.trialinfo(:,12:19));
%             cfg=[];
%             cfg.covariance='yes';
%             cfg.covariancewindow=[-0.1 0.7]; % a full window valid for all ll
%             tlock_audMSshift1tlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_audMSshift1{ll,tt,ss});
%             tlock_audMSshift2tlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_audMSshift2{ll,tt,ss});
%             tlock_audMSshift3tlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_audMSshift3{ll,tt,ss});
%             tlock_audMSshift4tlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_audMSshift4{ll,tt,ss});
%           else
%             tlock_audMSshift1tlock{ll,tt,ss}=[];
%             tlock_audMSshift2tlock{ll,tt,ss}=[];
%             tlock_audMSshift3tlock{ll,tt,ss}=[];
%             tlock_audMSshift4tlock{ll,tt,ss}=[];
%           end
%           
%           if ~isempty(tlock_audMSshift1tlock{ll,tt,ss})
%             % create sum of asynchronous multisensory conditions (e.g. AT70_shifted plus TA70)
%             cfg=[];
%             cfg.operation='add';
%             cfg.parameter='avg';
%             tlock_aMSasynch{ll,tt,ss}=ft_math(cfg,tlock_audMSshift1tlock{ll,tt,ss},tlock_audMSshift2tlock{ll,tt,ss});
%             cfg.parameter='cov';
%             tmp=ft_math(cfg,tlock_audMSshift1tlock{ll,tt,ss},tlock_audMSshift2tlock{ll,tt,ss});
%             tlock_aMSasynch{ll,tt,ss}.cov=tmp.cov;
%             featurestats_aMSasynch(:,ll,tt,ss,ii)=mean([featurestats_audMSshift1(:,ll,tt,ss,ii), featurestats_audMSshift2(:,ll,tt,ss,ii)]');
%             % create sum of synchronous multisensory conditions (e.g. TA0_shifted plus TA0)
%             cfg=[];
%             cfg.operation='add';
%             cfg.parameter='avg';
%             tlock_aMSsynch{ll,tt,ss}=ft_math(cfg,tlock_audMSshift3tlock{ll,tt,ss},tlock_audMSshift4tlock{ll,tt,ss});
%             cfg.parameter='cov';
%             tmp=ft_math(cfg,tlock_audMSshift3tlock{ll,tt,ss},tlock_audMSshift4tlock{ll,tt,ss});
%             tlock_aMSsynch{ll,tt,ss}.cov=tmp.cov;
%             featurestats_aMSsynch(:,ll,tt,ss,ii)=mean([featurestats_audMSshift3(:,ll,tt,ss,ii), featurestats_audMSshift4(:,ll,tt,ss,ii)]');
%           else
%             tlock_aMSasynch{ll,tt,ss}=[];
%             tlock_aMSsynch{ll,tt,ss}=[];
%           end
%         end % ll
        
 %% tw added    
    for ll=1:101
     try
                if isfield(tlock_aud{ll,tt,ss},'trial') && size(tlock_aud{ll,tt,ss}.trial,1)
                cfg_alpha=[];
                cfg_alpha.bpfilter='yes';
                cfg_alpha.bpfreq=[8 12];
                cfg_alpha.demean='yes';
                tlock_aud{ll,tt,ss}.alpha=[];
                tlock_aud{ll,tt,ss}.alpha=ft_preprocessing(cfg_alpha,tlock_aud{ll,tt,ss});
                %cfg_alpha.baselinewindow=cfg.baselinewindow
                    for nTrials=1:length(tlock_aud{ll,tt,ss}.trial(:,1,1))
                        tmp.bp=[];
                        tmp.bp(1:length(tlock_aud{ll,tt,ss}.label),1:length(tlock_aud{ll,tt,ss}.trial))=tlock_aud{ll,tt,ss}.alpha.trial(nTrials,:,:);
                        tlock_aud{ll,tt,ss}.alpha.abs(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'abs');
                        tlock_aud{ll,tt,ss}.alpha.complex(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'complex');
                        tlock_aud{ll,tt,ss}.alpha.real(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'real');
                        tlock_aud{ll,tt,ss}.alpha.absreal(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'absreal');
                        tlock_aud{ll,tt,ss}.alpha.absimag(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'absimag');
                        tlock_aud{ll,tt,ss}.alpha.angle(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'angle');

                        tlock_aud{ll,tt,ss}.hilbert.alphaT0.abs(nTrials,:)=tlock_aud{ll,tt,ss}.alpha.abs(:,find(abs(tlock_aud{ll,tt,ss}.time)==min(abs(tlock_aud{ll,tt,ss}.time))),nTrials);
                        % calculate abs (amplitude envelope over 1000 ms pre-stimulus)
                        tlock_aud{ll,tt,ss}.hilbert.alphaT0.absMean(nTrials,:)= ...
                        mean(tlock_aud{ll,tt,ss}.alpha.abs(:,find(abs(tlock_aud{ll,tt,ss}.time)==min(abs(tlock_aud{ll,tt,ss}.time))):find(abs(tlock_aud{ll,tt,ss}.time)==min(abs(tlock_aud{ll,tt,ss}.time)))+1000,nTrials));
                    
                        tlock_aud{ll,tt,ss}.hilbert.alphaT0.complex(nTrials,:)=tlock_aud{ll,tt,ss}.alpha.complex(:,find(abs(tlock_aud{ll,tt,ss}.time)==min(abs(tlock_aud{ll,tt,ss}.time))),nTrials);
                        tlock_aud{ll,tt,ss}.hilbert.alphaT0.real(nTrials,:)=tlock_aud{ll,tt,ss}.alpha.real(:,find(abs(tlock_aud{ll,tt,ss}.time)==min(abs(tlock_aud{ll,tt,ss}.time))),nTrials);
                        tlock_aud{ll,tt,ss}.hilbert.alphaT0.absreal(nTrials,:)=tlock_aud{ll,tt,ss}.alpha.absreal(:,find(abs(tlock_aud{ll,tt,ss}.time)==min(abs(tlock_aud{ll,tt,ss}.time))),nTrials);
                        tlock_aud{ll,tt,ss}.hilbert.alphaT0.absimag(nTrials,:)=tlock_aud{ll,tt,ss}.alpha.absimag(:,find(abs(tlock_aud{ll,tt,ss}.time)==min(abs(tlock_aud{ll,tt,ss}.time))),nTrials);
                        tlock_aud{ll,tt,ss}.hilbert.alphaT0.angle(nTrials,:)=tlock_aud{ll,tt,ss}.alpha.angle(:,find(abs(tlock_aud{ll,tt,ss}.time)==min(abs(tlock_aud{ll,tt,ss}.time))),nTrials);
                    end
% and for theta
                    cfg_theta=[];
                    cfg_theta.bpfilter='yes';
                    cfg_theta.bpfreq=[4 8];
                    cfg_theta.demean='yes';
                    tlock_aud{ll,tt,ss}.theta=[];
                    tlock_aud{ll,tt,ss}.theta=ft_preprocessing(cfg_theta,tlock_aud{ll,tt,ss});
                    %cfg_theta.baselinewindow=cfg.baselinewindow
                        for nTrials=1:length(tlock_aud{ll,tt,ss}.trial(:,1,1))
                            tmp.bp=[];
                            tmp.bp(1:length(tlock_aud{ll,tt,ss}.label),1:length(tlock_aud{ll,tt,ss}.trial))=tlock_aud{ll,tt,ss}.theta.trial(nTrials,:,:);
                            tlock_aud{ll,tt,ss}.theta.abs(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'abs');
                            tlock_aud{ll,tt,ss}.theta.complex(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'complex');
                            tlock_aud{ll,tt,ss}.theta.real(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'real');
                            tlock_aud{ll,tt,ss}.theta.absreal(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'absreal');
                            tlock_aud{ll,tt,ss}.theta.absimag(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'absimag');
                            tlock_aud{ll,tt,ss}.theta.angle(:,:,nTrials)=ft_preproc_hilbert(tmp.bp, 'angle');

                            tlock_aud{ll,tt,ss}.hilbert.thetaT0.abs(nTrials,:)=tlock_aud{ll,tt,ss}.theta.abs(:,find(abs(tlock_aud{ll,tt,ss}.time)==min(abs(tlock_aud{ll,tt,ss}.time))),nTrials);
                            % calculate abs (amplitude envelope over 1000 ms pre-stimulus)
                            tlock_aud{ll,tt,ss}.hilbert.thetaT0.absMean(nTrials,:)= ...
                            mean(tlock_aud{ll,tt,ss}.theta.abs(:,find(abs(tlock_aud{ll,tt,ss}.time)==min(abs(tlock_aud{ll,tt,ss}.time))):find(abs(tlock_aud{ll,tt,ss}.time)==min(abs(tlock_aud{ll,tt,ss}.time)))+1000,nTrials));
                    
                            tlock_aud{ll,tt,ss}.hilbert.thetaT0.complex(nTrials,:)=tlock_aud{ll,tt,ss}.theta.complex(:,find(abs(tlock_aud{ll,tt,ss}.time)==min(abs(tlock_aud{ll,tt,ss}.time))),nTrials);
                            tlock_aud{ll,tt,ss}.hilbert.thetaT0.real(nTrials,:)=tlock_aud{ll,tt,ss}.theta.real(:,find(abs(tlock_aud{ll,tt,ss}.time)==min(abs(tlock_aud{ll,tt,ss}.time))),nTrials);
                            tlock_aud{ll,tt,ss}.hilbert.thetaT0.absreal(nTrials,:)=tlock_aud{ll,tt,ss}.theta.absreal(:,find(abs(tlock_aud{ll,tt,ss}.time)==min(abs(tlock_aud{ll,tt,ss}.time))),nTrials);
                            tlock_aud{ll,tt,ss}.hilbert.thetaT0.absimag(nTrials,:)=tlock_aud{ll,tt,ss}.theta.absimag(:,find(abs(tlock_aud{ll,tt,ss}.time)==min(abs(tlock_aud{ll,tt,ss}.time))),nTrials);
                            tlock_aud{ll,tt,ss}.hilbert.thetaT0.angle(nTrials,:)=tlock_aud{ll,tt,ss}.theta.angle(:,find(abs(tlock_aud{ll,tt,ss}.time)==min(abs(tlock_aud{ll,tt,ss}.time))),nTrials);
                        end
                    
                end
     catch; end
    end       
        
    
        
        
      end  % end ss
      eeg_legomagic_phaseKc_distribution4   
        
      %clear tlock_tac tlock_aud tlock_nul tr tlock*tlock
      %save(['alock_hilbert_' sub{ii} '.mat'],'tlock*','tr','-v7.3');
   %   save(['alock_hilbert_tt3_to_save' sub{ii} '.mat'],'tlock_aud_to_save','tr','-v7.3');
        save(['alock_hilbert_tt3_to_save_5_' sub{ii} '.mat'], 'KC_ctr', '-v7.3'); %
      % save(['alock_hilbert_tt3_to_save_3' sub{ii} '.mat'],'KC','Kc_Cdf', 'KC_ctr', '-v7.3'); %
     % save(['alock_hilbert_tt3_to_save_4' sub{ii} '.mat'],'KC', '-v7.3'); %
       clear tlock_tac tlock_aud tlock_nul
    end % end tt
    

    %%
    
    
    
    %   save(['tlock_diffs_averef_' sub{ii} '.mat'],'*trialkept','tlock_*P*','tlock_*N*')
    % save(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '.mat'],'tlock_*P*','tlock_*N*','tlock_t*','tlock_a*','num*trials','featurestats_*')
    %save(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep)
    %'.mat'],'tlock_*P*','tlock_*N*','tlock_t*','tlock_a*','tlock*tlock','num*trials','featurestats_*')
    %%%%%%%%%TW commented last line%%%%%%
    
 
  end % ii
end % sleep

return

%% Group level: awake and asleep testing tactile alone response
% This is not for the main finding, but rather testing initial strange
% finding of no seeming Tactile response for N23 bed data.

plotflag=1;
printflag=0;
statsflag=0;

tt=3;
clearvars -except tt sub edir ddir iiuse plotflag printflag statsflag
% subuse=[5 8 9  10    12    14    15    16      18 21 23];
subuse=iiuse;
sleepss_tacAll=zeros(13,11,max(subuse));
sleepss_tac19T=zeros(13,11,max(subuse));
sleepss_tacAlone=zeros(13,11,max(subuse));
subind=1;
for ii=subuse
  cd([edir sub{ii} ])
  for sleep=[0 1]
    %   load(['tlock_diffs_' sub{ii} '.mat']);
    %       load(['tlock_diffs_averef_' sub{ii} '.mat']);
    load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '.mat'])
    %         load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '.mat'],'tr')
    tk0=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(0) '.mat']);
    tk1=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '.mat']);
    if sleep
      ssuse=tk1.tr.stageuse;
    else
      ssuse=tk0.tr.stageuse;
    end
    
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



%% Group level: awake and asleep separately, for each asynchrony separately

soadesc={'Aud first by 500ms' '' 'Aud first by 70ms' 'Aud first by 20ms' 'Simultaneous' 'Tac first by 20ms' 'Tac first by 70ms' '' 'Tac first by 500ms'};
plotflag=1;
printflag=1;
statsflag=1;
audtacflag=0;
soalist=[1 3 4 5 6 7 9];

chanuse_sleep0={'all' '-F4'};
chanuse_sleep1={'all' '-AF7' '-AF3' '-Fp1'};

% iiSuse=setdiff(iiSuse, 11);

% allcond_sameN=1; % means using ii>=8, but we doing that anyway with iiuse
% for tacaud=[1 0]
for sleep=[0]
  if sleep
    chanuse=chanuse_sleep1;
  else
    chanuse=chanuse_sleep0;
  end
  
  for tt=[2 3]
    close all
    figind=1;
    
    for ll=soalist
      % for ll=[3 4 5 6 7]
      %   for tt=1:4
      clearvars -except ll tt sub edir ddir ii* sleep *flag figind soadesc soalist chanuse* stat* grave*T* grind_t* grind_a* plv
      
      %     if ll==1 | ll==9
      %       subuse=8:32;
      %     else
      %       subuse=5:32;
      %     end
      %     if allcond_sameN
      %       subuse=8:32;
      %     end
      if sleep
        subuseall=setdiff(iiBuse,[3:7]);
      else
        subuseall=iiSuse;
      end
      %       subuseall=setdiff(iiuse,[ 27    28    30    31]);
      
      submin=subuseall(1)-1;
      subuseind=0;
      
      for ii=subuseall
        %       for ii=[8 9 10 12 14 15  16 17 18]
        cd([edir sub{ii} ])
        %   load(['tlock_diffs_' sub{ii} '.mat']);
        %       load(['tlock_diffs_averef_' sub{ii} '.mat']);
        load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '.mat'])
        %         load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '.mat'],'tr')
        tk0=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(0) '.mat']);
        tk1=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '.mat']);
        
        if 1
          % THis is preliminary...how best to included all stages later on?
          if sleep==0
            ss=10; % awake
            sleepcond='Awake';
          elseif sleep==1
            ssuse=12; % this is concatenation of N2 and N3
            sleepcond='Sleep N2';
            %             ssuse=23; % this is concatenation of N2 and N3
            %             sleepcond='Sleep (N2+N3)';
          end
        else
          if sleep
            ssuse=tk1.tr.stageuse;
          else
            ssuse=tk0.tr.stageuse;
          end
        end
        
        
        %                 for ss=ssuse
        subuse=subuseall; % reset to all for each sleep stage
        numtrt(ll,tt,ss,ii-submin)=numt_trials(ll,tt,ss); % does this need to be per 'sleep01' as well?
        if audtacflag
          numtra(ll,tt,ss,ii-submin)=numa_trials(ll,tt,ss);
        end
        
        if min(numtrt(ll,tt,ss,ii-submin),[],1)<20 % what is best number to use here?
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
          if ll<5
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
        %       clear *_s0
        clear tlock*N tlock*tac tlock*aud
      end % ii
      subuseindfinal=subuseind;
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
      grind_tactlock=ft_timelockgrandaverage(cfg,tlock_tactlock_each{:});
      grind_audtlock=ft_timelockgrandaverage(cfg,tlock_audtlock_each{:});
      grind_nultlock=ft_timelockgrandaverage(cfg,tlock_nulttlock_each{:});
      grind_MStlock=ft_timelockgrandaverage(cfg,tlock_MStlock_each{:});
      if ll<5
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
        if ll<5
      grave_tMSsynch=ft_timelockgrandaverage(cfg,tlock_tMSsynch_each{:});
      grave_tMSasynch=ft_timelockgrandaverage(cfg,tlock_tMSasynch_each{:});
        end
      if audtacflag
        grave_audPtac=ft_timelockgrandaverage(cfg,tlock_audPtac_each{:});
        grave_audMSpN=ft_timelockgrandaverage(cfg,tlock_audMSpN_each{:});
        if ll<5
        grave_aMSsynch=ft_timelockgrandaverage(cfg,tlock_aMSsynch_each{:});
        grave_aMSasynch=ft_timelockgrandaverage(cfg,tlock_aMSasynch_each{:});
        end
      end
      grave_tacVSnul=ft_timelockgrandaverage(cfg,tlock_tacVSnul_each{:});
      grave_audVSnul=ft_timelockgrandaverage(cfg,tlock_audVSnul_each{:});
      grave_msVSnul=ft_timelockgrandaverage(cfg,tlock_msVSnul_each{:});
      
      
      if plotflag
        if sleep
          topoplot_highlight(111,grave_tacPaud,[.07 .77],[]);
          topoplot_highlight(113,grave_tacMSpN,[.07 .77],[]);
          if printflag
            print(111,['D:\audtac\figs\gravetacPaud_topoOverSTime_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
            print(113,['D:\audtac\figs\gravetacMSpN_topoOverSTime_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          end
        else
          topoplot_highlight(111,grave_tacPaud,[.07 .37],[]);
          topoplot_highlight(113,grave_tacMSpN,[.07 .37],[]);
          if printflag
            print(111,['D:\audtac\figs\gravetacPaud_topoOverSTime_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
            print(113,['D:\audtac\figs\gravetacMSpN_topoOverSTime_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
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
      for cc=1:length(chanplot)
        if cc==2 % posterior, compute phase of alpha at time where signif findings
%           keyboard
          grindnames=whos('grind*');
          cfg=[];
          cfg.hilbert='complex';
          cfg.bpfilter='yes';
          cfg.bpfreq=[8 13];
          for gg=1:length(grindnames)
            hilgrind.(grindnames(gg).name)=ft_preprocessing(cfg,eval(grindnames(gg).name));
            tuse=dsearchn(hilgrind.(grindnames(gg).name).time',0):dsearchn(hilgrind.(grindnames(gg).name).time',.7);
            ccuse=match_str(hilgrind.(grindnames(gg).name).label,chanplot{cc});
            plv{ll,tt,ss}.(grindnames(gg).name)=squeeze(mean(hilgrind.(grindnames(gg).name).trial(:,ccuse,tuse)./abs(hilgrind.(grindnames(gg).name).trial(:,ccuse,tuse)),2));
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
        
        if plotflag
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

          if ll<5
            if sleep
              figure(140+cc*10); % 150, 160, 170
            else
              figure(170+cc*10); % 180, 190, 200
            end
            subplot(1,7,figind);
            cfg=[];
            cfg.channel=chanplot{cc}
            cfg.xlim=[-0.5 1.1];
            cfg.ylim=[-7 7];
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
            
            if ll<5
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
          
          
          
          
          
        end % plotflag
      end % cc
      if plotflag && printflag
        if sleep
          print(60,['D:\audtac\figs\gravetacMSpN_FC_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          print(70,['D:\audtac\figs\gravetacMSpN_OP_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          print(80,['D:\audtac\figs\gravetacMSpN_RT_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          print(90,['D:\audtac\figs\graveUniSens_FC_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          print(100,['D:\audtac\figs\graveUniSens_OP_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          print(110,['D:\audtac\figs\graveUniSens_RT_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          print(150,['D:\audtac\figs\gravetacMSsynch_FC_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          print(160,['D:\audtac\figs\gravetacMSsynch_OP_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          print(170,['D:\audtac\figs\gravetacMSsynch_RT_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          if audtacflag
            print(41,['D:\audtac\figs\graveaudMSpN_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          print(210,['D:\audtac\figs\graveaudMSsynch_FC_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          print(220,['D:\audtac\figs\graveaudMSsynch_OP_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          print(230,['D:\audtac\figs\graveaudMSsynch_RT_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          end
        else
          print(30,['D:\audtac\figs\gravetacMSpN_FC_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          print(40,['D:\audtac\figs\gravetacMSpN_OP_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          print(50,['D:\audtac\figs\gravetacMSpN_RT_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          print(120,['D:\audtac\figs\graveUniSens_FC_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          print(130,['D:\audtac\figs\graveUniSens_OP_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          print(140,['D:\audtac\figs\graveUniSens_RT_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          print(180,['D:\audtac\figs\gravetacMSsynch_FC_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          print(190,['D:\audtac\figs\gravetacMSsynch_OP_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          print(200,['D:\audtac\figs\gravetacMSsynch_RT_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          if audtacflag
            print(31,['D:\audtac\figs\graveaudMSpN_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          print(240,['D:\audtac\figs\graveaudMSsynch_FC_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          print(250,['D:\audtac\figs\graveaudMSsynch_OP_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
          print(260,['D:\audtac\figs\graveaudMSsynch_RT_erp_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
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
        if ll<5
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
        if ll<5
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
      
      if statsflag
        load eeg1010_neighb
        
        nsub=length(tlock_tacMSpN_each);
        
        cfg=[];
        if sleep
          %           cfg.latency=[.1 .8]; % longer to allow for Kc
          if ll==1 || ll==3 || ll==4 || ll==5
            cfg.latency=[.1 .8];
          elseif ll==6
            cfg.latency=[.12 .82];
          elseif ll==7
            cfg.latency=[.17 .87];
          elseif ll==9
            cfg.latency=[.6 .13];
          end
        else
          %           cfg.latency=[.1 .5]; % previously [-.1 .5]
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
        statt_mc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tacPaud, grind_tacMSpN);
        if audtacflag
          stata_mc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_audPtac, grind_audMSpN);
        end
        
%         % save out for brain-behaviour correlations
%         cfg=[];
%         cfg.latency=[statt_mc{ll,tt,ss}.time(1) statt_mc{ll,tt,ss}.time(end)];
%         grind_tacPaud_save{ll,tt,ss}=ft_selectdata(cfg,grind_tacPaud);
%         grind_tacMSpN_save{ll,tt,ss}=ft_selectdata(cfg,grind_tacMSpN);
        
        % late component
        cfg=[];
        if audtacflag==0
          if sleep
            keyboard % alter here
            %           cfg.latency=[.1 .8]; % longer to allow for Kc
            if ll==1 || ll==3 || ll==4 || ll==5
              cfg.latency=[.1 .8];
            elseif ll==6
              cfg.latency=[.12 .82];
            elseif ll==7
              cfg.latency=[.17 .87];
            elseif ll==9
              cfg.latency=[.6 .13];
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
        statt_latemc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tacPaud, grind_tacMSpN);
        if audtacflag
          stata_latemc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_audPtac, grind_audMSpN);
        end
        
        cfg.avgovertime='yes';
        statt_late{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tacPaud, grind_tacMSpN);
        if audtacflag
          stata_late{ll,tt,ss}=ft_timelockstatistics(cfg, grind_audPtac, grind_audMSpN);
        end
        
        % early and late all together
        cfg=[];
        if audtacflag==0
          if sleep
            keyboard % alter here
            %           cfg.latency=[.1 .8]; % longer to allow for Kc
            if ll==1 || ll==3 || ll==4 || ll==5
              cfg.latency=[.1 .8];
            elseif ll==6
              cfg.latency=[.12 .82];
            elseif ll==7
              cfg.latency=[.17 .87];
            elseif ll==9
              cfg.latency=[.6 .13];
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
        statt_allmc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tacPaud, grind_tacMSpN);
        statt_tacmc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tactlock, grind_nultlock);
        statt_audmc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_audtlock, grind_nultlock);
        statt_msmc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_MStlock,  grind_nultlock);
        if audtacflag
          stata_allmc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_audPtac, grind_audMSpN);
          stata_tacmc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tactlock, grind_nultlock);
          stata_audmc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_audtlock, grind_nultlock);
          stata_msmc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_MStlock,  grind_nultlock);
        end
        
        % Funny/Uta temporal (in)congruent contrast: (AT70 + TA70) vs (Simult + Simult_shift70)
        if ll<5
          cfg=[];
          if audtacflag==0
            if sleep
              if ll==4
                cfg.latency=[.12 .82];
              elseif ll==3
                cfg.latency=[.17 .87];
              elseif ll==1
                cfg.latency=[.6 .13];
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
          statt_synch{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tMSsynch, grind_tMSasynch);
          if audtacflag
            stata_synch{ll,tt,ss}=ft_timelockstatistics(cfg, grind_aMSsynch, grind_aMSasynch);
          end
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
            print(22,['D:\audtac\figs\grdiff_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(24,['D:\audtac\figs\grdiff_topoOverTime_late_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(26,['D:\audtac\figs\grdiff_topoOverTime_latemc_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(28,['D:\audtac\figs\grdiff_topoOverTime_allmc_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(10,['D:\audtac\figs\grdiff_topoOverTime_tacmc_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(12,['D:\audtac\figs\grdiff_topoOverTime_audmc_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(14,['D:\audtac\figs\grdiff_topoOverTime_msmc_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            if ll<5, print(32,['D:\audtac\figs\grdiff_topoOverTime_synch_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng');end
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
              print(23,['D:\audtac\figs\grdiff_topoOverTime_mc_at_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
              print(25,['D:\audtac\figs\grdiff_topoOverTime_late_at_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
              print(27,['D:\audtac\figs\grdiff_topoOverTime_latemc_at_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
              print(29,['D:\audtac\figs\grdiff_topoOverTime_allmc_at_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
              print(11,['D:\audtac\figs\grdiff_topoOverTime_tacmc_at_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
              print(13,['D:\audtac\figs\grdiff_topoOverTime_audmc_at_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
              print(15,['D:\audtac\figs\grdiff_topoOverTime_msmc_at_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
              if ll<5, print(33,['D:\audtac\figs\grdiff_topoOverTime_synch_at_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng');end
            end
          end
          
        end % if plotflag
        %       else
        %         print(23,['D:\audtac\figs\grdiff_topoOverTime_mc_at_cond' num2str(ll) num2str(tt) 'a.png'],'-dpng')
        %       end
      end % if statsflag
      
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
  end % tt
  %   save([edir 'tlock_statmc.mat'],'stat*','grave*'); % no point, as grave* isn't ll,tt,ss dependent
  %   save([edir 'tlock_statmc_sleep' num2str(sleep) '.mat'],'stat*');
  save([edir 'tlock_statmc_sleep' num2str(sleep) '.mat'],'stat*','grave*T*','grind_t*','plv');
  if audtacflag
    save([edir 'tlock_statmc_sleep' num2str(sleep) '.mat'],'stat*','grave*T*','grind_a*','plv');
  end
end % sleep
% end
% save([edir 'tlock_numtrlltt.mat'],'numtr*','grind*'); % no point as grind* isn't ll,tt,ss dependent
save([edir 'tlock_numtrlltt.mat'],'numtr*');

%% Plotting above ERP sensor results


% First, plotting the final stats/ difference of conditions
% Below, plotting conditions on their own.

clearvars -except sub *dir
chanplot{1}={'Fz' 'Cz' 'F1' 'F2' 'FC1' 'FC2' 'C1' 'C2'}; % frontocentral
chanplot{2}={'CP5' 'P5' 'PO3' 'POz' 'PO4' 'P6' 'CP6' 'Pz' 'P3' 'P4'}; % occipital
chanplot{3}={'FC6' 'FC4' 'F4' 'F6' 'F8' 'FT8' 'T8' 'C6' 'C4'}; % right frontotemporal
chanlabel{1}='Frontocentral electrodes';
chanlabel{2}='Occipital-parietal electrodes';
chanlabel{3}='Right frontotemporal electrodes';
sleep=0;
load([edir 'tlock_statmc_sleep' num2str(sleep) '.mat']);
soalist=[1 3 4 5 6 7 9];
tt=2;
ss=10;
timwinstatflag=0; % =1 for using only stat.time, =0 for using [-0.5 1.0];

for ll=soalist
  cfg=[];
  if timwinstatflag==1
    cfg.latency=[statt_mc{ll,tt,ss}.time(1) statt_mc{ll,tt,ss}.time(end)];
  elseif timwinstatflag==0
    cfg.latency=[-0.5 1];
    stattimwin=[statt_mc{ll,tt,ss}.time(1) statt_mc{ll,tt,ss}.time(end)];
  end
  tmp=ft_selectdata(cfg,grave_TPA_MSPN{ll,tt,ss});
  tmp.mask=statt_mc{ll,tt,ss}.mask;
  tmpm=ft_selectdata(cfg,grind_tacMSpN_save{ll,tt,ss});
  tmpu=ft_selectdata(cfg,grind_tacPaud_save{ll,tt,ss});
  tmpm.mask=statt_mc{ll,tt,ss}.mask;
  tmpu.mask=statt_mc{ll,tt,ss}.mask;
  
  tmpm.avg=squeeze(mean(tmpm.individual,1));
  tmpu.avg=squeeze(mean(tmpu.individual,1));
  tmpm.var=squeeze(var(tmpm.individual,1));
  tmpu.var=squeeze(var(tmpu.individual,1));
  tmpm.dimord='chan_time';
  tmpu.dimord='chan_time';
  tmpm=rmfield(tmpm,'individual');
  tmpu=rmfield(tmpu,'individual');
  
  cfg=[];
  cfg.avgoverchan='yes';
  cfg.channel=chanplot{1};
  tmp1=ft_selectdata(cfg,tmp);
  if timwinstatflag==0
    tmp1mask=mean(tmp.mask(match_str(tmp.label,cfg.channel),:),1);
    tmp1.mask=zeros(1,length(tmp1.time));
    tmp1.mask(dsearchn(tmp1.time',stattimwin(1)):dsearchn(tmp1.time',stattimwin(end)))=tmp1mask;
  end
  tmp1.mask=logical(ceil(tmp1.mask));
  tmpm1=ft_selectdata(cfg,tmpm);
  tmpm1.mask=logical(ceil(tmpm1.mask));
  tmpu1=ft_selectdata(cfg,tmpu);
  tmpu1.mask=logical(ceil(tmpu1.mask));
  
  figure(10*ll);
  cfg=[];
  cfg.parameter='avg';
  cfg.layout='elec1010.lay';
  cfg.maskparameter='mask';
  cfg.maskstyle='box'; % default
  cfg.ylim=[-3 7];
  cfg.linewidth=3;
  ft_singleplotER(cfg,tmpu1,tmpm1,tmp1);
  print(10*ll,[fdir 'erp_final_FC_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  
  cfg=[];
  cfg.avgoverchan='yes';
  cfg.channel=chanplot{2};
  tmp1=ft_selectdata(cfg,tmp);
  if timwinstatflag==0
    tmp1mask=mean(tmp.mask(match_str(tmp.label,cfg.channel),:),1);
    tmp1.mask=zeros(1,length(tmp1.time));
    tmp1.mask(dsearchn(tmp1.time',stattimwin(1)):dsearchn(tmp1.time',stattimwin(end)))=tmp1mask;
  end
  tmp1.mask=logical(ceil(tmp1.mask));
  tmpm1=ft_selectdata(cfg,tmpm);
  tmpm1.mask=logical(ceil(tmpm1.mask));
  tmpu1=ft_selectdata(cfg,tmpu);
  tmpu1.mask=logical(ceil(tmpu1.mask));
  
  figure(10*ll+1);
  cfg=[];
  cfg.parameter='avg';
  cfg.layout='elec1010.lay';
  cfg.maskparameter='mask';
  cfg.maskstyle='box'; % default
  cfg.ylim=[-3 3];
  cfg.linewidth=3;
  ft_singleplotER(cfg,tmpu1,tmpm1,tmp1);
  print(10*ll+1,[fdir 'erp_final_OP_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  
  if [tt==2 && any(ll==[3 4 5 6])] || [tt==3 && any(ll==[4 6 7])]
    masktime=find(any(tmp.mask,1));
    cfg=[];
    cfg.parameter='avg';
    cfg.layout='elec1010.lay';
    cfg.maskalpha=0.5;
    cfg.zlim=[-5 5];
    cfg.highlight='on';
    cfg.highlightsize=12;
    cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
    %     if ll==3
    %       cfg.xlim=[.1 .35];
    %     elseif ll==4
    %       cfg.xlim=[.06 .42];
    %     elseif ll==5
    %       cfg.xlim=[-.04 .36];
    %     elseif ll==6
    %       cfg.xlim=[.1 .44];
    %     end
    cfg.highlightchannel=tmp.label(find(ceil(mean(tmp.mask(:,dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2))),2))));
    figure(10*ll+2);
    ft_topoplotER(cfg,tmpu);
    print(10*ll+2,[fdir 'erp_topoU_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    figure(10*ll+3);
    ft_topoplotER(cfg,tmpm);
    print(10*ll+3,[fdir 'erp_topoM_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    figure(10*ll+4);
    ft_topoplotER(cfg,tmp);
    print(10*ll+4,[fdir 'erp_topoDiff_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  end
  
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
