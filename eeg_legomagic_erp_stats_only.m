
eeg_legomagic_preamble

%% Group level: awake and asleep separately, for each asynchrony separately

soadesc={'Aud first by 500ms' '' 'Aud first by 70ms' 'Aud first by 20ms' 'Simultaneous' 'Tac first by 20ms' 'Tac first by 70ms' '' 'Tac first by 500ms'};
plotflag=0;
printflag=0;
statsflag=1;
medsplitstatsflag=0;
loadprevstatsflag=0;
tacnulmsaudstatsflag=0;
% tacaudmsnul_earlyflag=1;
savegrindflag=1;
iterflag=1; % multiple iterations of random trial assignment
audtacflag=0;
runagainflag=0; % 2nd run of erp results (resampling of trials into ERP)
trialkcflag=1; % 0 for ignore trialkc (old results; all trials), or 1 for use trialkc value
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

mcseed=13;

chanuse_sleep0={'all' '-F4'};
chanuse_sleep1={'all' '-AF7' '-AF3' '-Fp1'};

% iiSuse=setdiff(iiSuse, 11);

% allcond_sameN=1; % means using ii>=8, but we doing that anyway with iiuse
% for tacaud=[1 0]
% for sleep=[1]
sleep=1;
if sleep
  chanuse=chanuse_sleep1;
  iteruse=11;
  trialkc=-1;  % change me: 0 (no Kc) or 1 (only Kc) or -1 (all trials)
  usetr=1;
else
  chanuse=chanuse_sleep0;
  %   warning('change me back to 27!!')
  %   iteruse=27;
  iteruse=32;
  trialkc=-1;
  usetr=3;
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
    %           for ii=[8 9]
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
        if usetr==2
          % load usetr=0 here; then later down load usetr=2
          tk1=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter)  '_usetr' num2str(0) '_trialkc' num2str(trialkc) '.mat']);
        else
          tk1=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter)  '_usetr' num2str(usetr) '_trialkc' num2str(trialkc) '.mat']);
        end
        if sleep
          load(['tlock_diffs_features_' sub{ii} '_sleep' num2str(sleep) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'],'featstruct*');
        end          
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
        %         ss=11; % N1
        %         ss=10; % W
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
        tlock_tacPaud_each{subuseind,1}=tlock_tacPaud{ll,tt,ss};
        tlock_tacMSpN_each{subuseind,1}=tlock_tacMSpN{ll,tt,ss};
        tlock_MStlock_each{subuseind,1}=tlock_tMSAlone{ll,tt,ss}; % ll is tac plus auditory multisensory, tac-locked
        tlock_tactlock_each{subuseind,1}=tlock_tTacAlone{ll,tt,ss}; % ll+20 is tac alone, tac-locked
        tlock_audtlock_each{subuseind,1}=tlock_tAudAlone{ll,tt,ss}; % ll+20 is aud alone, aud-locked
        tlock_nulttlock_each{subuseind,1}=tlock_tNulAlone{ll,tt,ss}; % ll+50 is nul, tac-locked
        if audtacflag
          tlock_nulatlock_each{subuseind}=tlock_aNulAlone{ll,tt,ss}; % ll+60 is nul, aud-locked
          tlock_audPtac_each{subuseind}=tlock_audPtac{ll,tt,ss};
          tlock_audMSpN_each{subuseind}=tlock_audMSpN{ll,tt,ss};
        end
        
        if sleep
          fn=fieldnames(featstruct_tacMSpN{ll,tt,ss});
          for ff=1:length(fn)
            tlock_tacPaud_each{subuseind}.(fn{ff})=featstruct_tacPaud{ll,tt,ss}.(fn{ff});
            tlock_tacMSpN_each{subuseind}.(fn{ff})=featstruct_tacMSpN{ll,tt,ss}.(fn{ff});
          end
          KcDuring_tacPaud{subuseind}=zeros(1,numt_trials(ll,tt,ss));
          KcPre_tacPaud{subuseind}=zeros(1,numt_trials(ll,tt,ss));
          KcEvoked_tacPaud{subuseind}=zeros(1,numt_trials(ll,tt,ss));
          for kk=1:numt_trials(ll,tt,ss)
            numKc=length(tlock_tacPaud_each{subuseind}.kc{kk});
            for nk=1:numKc
              if tlock_tacPaud_each{subuseind}.kc{kk}(nk).down<0+soades(ll) && tlock_tacPaud_each{subuseind}.kc{kk}(nk).upend>0+soades(ll)
                KcDuring_tacPaud{subuseind}(kk)=1;
              end
              if tlock_tacPaud_each{subuseind}.kc{kk}(nk).upend<0+soades(ll) && tlock_tacPaud_each{subuseind}.kc{kk}(nk).upend>-0.5+soades(ll)
                KcPre_tacPaud{subuseind}(kk)=1;
              end
              if tlock_tacPaud_each{subuseind}.kc{kk}(nk).negmax<0.8+soades(ll) && tlock_tacPaud_each{subuseind}.kc{kk}(nk).negmax>0.1+soades(ll)  %  Dang Vu 2011
                KcEvoked_tacPaud{subuseind}(kk)=1;
              end
            end % nk
            numKc=length(tlock_tacPaud_each{subuseind}.sw{kk});
            for nk=1:numKc
              if tlock_tacPaud_each{subuseind}.sw{kk}(nk).down<0+soades(ll) && tlock_tacPaud_each{subuseind}.sw{kk}(nk).upend>0+soades(ll)
                KcDuring_tacPaud{subuseind}(kk)=1;
              end
              if tlock_tacPaud_each{subuseind}.sw{kk}(nk).upend<0+soades(ll) && tlock_tacPaud_each{subuseind}.sw{kk}(nk).upend>-0.5+soades(ll)
                KcPre_tacPaud{subuseind}(kk)=1;
              end
              if tlock_tacPaud_each{subuseind}.sw{kk}(nk).negmax<0.8+soades(ll) && tlock_tacPaud_each{subuseind}.sw{kk}(nk).negmax>0.1+soades(ll)  %  Dang Vu 2011
                KcEvoked_tacPaud{subuseind}(kk)=1;
              end
            end % nk
          end  % kk
          KcDuring_tacMSpN{subuseind}=zeros(1,numt_trials(ll,tt,ss));
          KcPre_tacMSpN{subuseind}=zeros(1,numt_trials(ll,tt,ss));
          KcEvoked_tacMSpN{subuseind}=zeros(1,numt_trials(ll,tt,ss));
          for kk=1:numt_trials(ll,tt,ss)
            numKc=length(tlock_tacMSpN_each{subuseind}.kc{kk});
            for nk=1:numKc
              if tlock_tacMSpN_each{subuseind}.kc{kk}(nk).upend<0+soades(ll) && tlock_tacMSpN_each{subuseind}.kc{kk}(nk).upend>-0.5+soades(ll)
                KcPre_tacMSpN{subuseind}(kk)=1;
              end
              if tlock_tacMSpN_each{subuseind}.kc{kk}(nk).down<0+soades(ll) && tlock_tacMSpN_each{subuseind}.kc{kk}(nk).upend>0+soades(ll)
                KcDuring_tacMSpN{subuseind}(kk)=1;
              end
              if tlock_tacMSpN_each{subuseind}.kc{kk}(nk).negmax<0.8+soades(ll) && tlock_tacMSpN_each{subuseind}.kc{kk}(nk).negmax>0.1+soades(ll) %  Dang Vu 2011
                KcEvoked_tacMSpN{subuseind}(kk)=1;
              end
            end % nk
            numKc=length(tlock_tacMSpN_each{subuseind}.sw{kk});
            for nk=1:numKc
              if tlock_tacMSpN_each{subuseind}.sw{kk}(nk).upend<0+soades(ll) && tlock_tacMSpN_each{subuseind}.sw{kk}(nk).upend>-0.5+soades(ll)
                KcPre_tacMSpN{subuseind}(kk)=1;
              end
              if tlock_tacMSpN_each{subuseind}.sw{kk}(nk).down<0+soades(ll) && tlock_tacMSpN_each{subuseind}.sw{kk}(nk).upend>0+soades(ll)
                KcDuring_tacMSpN{subuseind}(kk)=1;
              end
              if tlock_tacMSpN_each{subuseind}.sw{kk}(nk).negmax<0.8+soades(ll) && tlock_tacMSpN_each{subuseind}.sw{kk}(nk).negmax>0.1+soades(ll)  %  Dang Vu 2011
                KcEvoked_tacMSpN{subuseind}(kk)=1;
              end
            end % nk
          end  % kk
          % Now spindles
          SpDuring_tacPaud{subuseind}=zeros(1,numt_trials(ll,tt,ss));
          SpPre_tacPaud{subuseind}=zeros(1,numt_trials(ll,tt,ss));
          SpEvoked_tacPaud{subuseind}=zeros(1,numt_trials(ll,tt,ss));
          for kk=1:numt_trials(ll,tt,ss)
            numSp=length(tlock_tacPaud_each{subuseind}.sp_fast{kk});
            for nk=1:numSp
              if tlock_tacPaud_each{subuseind}.sp_fast{kk}(nk).endtime<0+soades(ll) && tlock_tacPaud_each{subuseind}.sp_fast{kk}(nk).endtime>-0.5+soades(ll)
                SpPre_tacPaud{subuseind}(kk)=1;
              end
              if tlock_tacPaud_each{subuseind}.sp_fast{kk}(nk).starttime<0+soades(ll) && tlock_tacPaud_each{subuseind}.sp_fast{kk}(nk).endtime>0+soades(ll)
                SpDuring_tacPaud{subuseind}(kk)=1;
              end
              if tlock_tacPaud_each{subuseind}.sp_fast{kk}(nk).starttime<3+soades(ll) && tlock_tacPaud_each{subuseind}.sp_fast{kk}(nk).starttime>0.1+soades(ll)  %  Sato et al 2007
                SpEvoked_tacPaud{subuseind}(kk)=1;
              end
            end % nk
            numSp=length(tlock_tacPaud_each{subuseind}.sp_slow{kk});
            for nk=1:numSp
              if tlock_tacPaud_each{subuseind}.sp_slow{kk}(nk).endtime<0+soades(ll) && tlock_tacPaud_each{subuseind}.sp_slow{kk}(nk).endtime>-0.5+soades(ll)
                SpPre_tacPaud{subuseind}(kk)=1;
              end
              if tlock_tacPaud_each{subuseind}.sp_slow{kk}(nk).starttime<0+soades(ll) && tlock_tacPaud_each{subuseind}.sp_slow{kk}(nk).endtime>0+soades(ll)
                SpDuring_tacPaud{subuseind}(kk)=1;
              end
              if tlock_tacPaud_each{subuseind}.sp_slow{kk}(nk).starttime<3+soades(ll) && tlock_tacPaud_each{subuseind}.sp_slow{kk}(nk).starttime>0.1+soades(ll)  %  Sato et al 2007
                SpEvoked_tacPaud{subuseind}(kk)=1;
              end
            end % nk
          end  % kk
          SpDuring_tacMSpN{subuseind}=zeros(1,numt_trials(ll,tt,ss));
          SpPre_tacMSpN{subuseind}=zeros(1,numt_trials(ll,tt,ss));
          SpEvoked_tacMSpN{subuseind}=zeros(1,numt_trials(ll,tt,ss));
          for kk=1:numt_trials(ll,tt,ss)
            numSp=length(tlock_tacMSpN_each{subuseind}.sp_fast{kk});
            for nk=1:numSp
              if tlock_tacMSpN_each{subuseind}.sp_fast{kk}(nk).endtime<0+soades(ll) && tlock_tacMSpN_each{subuseind}.sp_fast{kk}(nk).endtime>-0.5+soades(ll)
                SpPre_tacMSpN{subuseind}(kk)=1;
              end
              if tlock_tacMSpN_each{subuseind}.sp_fast{kk}(nk).starttime<0+soades(ll) && tlock_tacMSpN_each{subuseind}.sp_fast{kk}(nk).endtime>0+soades(ll)
                SpDuring_tacMSpN{subuseind}(kk)=1;
              end
              if tlock_tacMSpN_each{subuseind}.sp_fast{kk}(nk).starttime<3+soades(ll) && tlock_tacMSpN_each{subuseind}.sp_fast{kk}(nk).starttime>0.1+soades(ll)  %  Sato et al 2007
                SpEvoked_tacMSpN{subuseind}(kk)=1;
              end
            end % nk
            numSp=length(tlock_tacMSpN_each{subuseind}.sp_slow{kk});
            for nk=1:numSp
              if tlock_tacMSpN_each{subuseind}.sp_slow{kk}(nk).endtime<0+soades(ll) && tlock_tacMSpN_each{subuseind}.sp_slow{kk}(nk).endtime>-0.5+soades(ll)
                SpPre_tacMSpN{subuseind}(kk)=1;
              end
              if tlock_tacMSpN_each{subuseind}.sp_slow{kk}(nk).starttime<0+soades(ll) && tlock_tacMSpN_each{subuseind}.sp_slow{kk}(nk).endtime>0+soades(ll)
                SpDuring_tacMSpN{subuseind}(kk)=1;
              end
              if tlock_tacMSpN_each{subuseind}.sp_slow{kk}(nk).starttime<3+soades(ll) && tlock_tacMSpN_each{subuseind}.sp_slow{kk}(nk).starttime>0.1+soades(ll)  %  Sato et al 2007
                SpEvoked_tacMSpN{subuseind}(kk)=1;
              end
            end % nk
          end  % kk

          
        end
        
        cfg=[];
        cfg.operation='subtract';
        cfg.parameter='avg';
        tlock_tacVSnul_each{subuseind,1}=ft_math(cfg,tlock_tTacAlone{ll,tt,ss},tlock_tNulAlone{ll,tt,ss});
        tlock_audVSnul_each{subuseind,1}=ft_math(cfg,tlock_tAudAlone{ll,tt,ss},tlock_tNulAlone{ll,tt,ss});
        tlock_msVSnul_each{subuseind,1}=ft_math(cfg,tlock_tMSAlone{ll,tt,ss},tlock_tNulAlone{ll,tt,ss});
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
    
    if usetr==2
      load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tt' num2str(tt) '_tacaud' num2str(tacaud) '_iter' num2str(iter+1) '_trialkc' num2str(trialkc) '.mat'])
      tk1=load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter+1)  '_usetr' num2str(usetr) '_trialkc' num2str(trialkc) '.mat']);
      tlock_tacPaud_each{subuseind,2}=tlock_tacPaud{ll,tt,ss};
      tlock_tacMSpN_each{subuseind,2}=tlock_tacMSpN{ll,tt,ss};
      tlock_MStlock_each{subuseind,2}=tlock_tMSAlone{ll,tt,ss}; % ll is tac plus auditory multisensory, tac-locked
      tlock_tactlock_each{subuseind,2}=tlock_tTacAlone{ll,tt,ss}; % ll+20 is tac alone, tac-locked
      tlock_audtlock_each{subuseind,2}=tlock_tAudAlone{ll,tt,ss}; % ll+20 is aud alone, aud-locked
      tlock_nulttlock_each{subuseind,2}=tlock_tNulAlone{ll,tt,ss}; % ll+50 is nul, tac-locked
      cfg=[];
      cfg.operation='subtract';
      cfg.parameter='avg';
      tlock_tacVSnul_each{subuseind,2}=ft_math(cfg,tlock_tTacAlone{ll,tt,ss},tlock_tNulAlone{ll,tt,ss});
      tlock_audVSnul_each{subuseind,2}=ft_math(cfg,tlock_tAudAlone{ll,tt,ss},tlock_tNulAlone{ll,tt,ss});
      tlock_msVSnul_each{subuseind,2}=ft_math(cfg,tlock_tMSAlone{ll,tt,ss},tlock_tNulAlone{ll,tt,ss});
    end
    clear tlock*N tlock*tac tlock*aud
    
  end % ii
  subuseindfinal=subuseind
  %     end
  
  keyboard
  
  for ii=1:subuseindfinal
    cfg=[];
    cfg.operation='subtract';
    cfg.parameter='avg';
    cfg.channel=chanuse;
    if usetr==0
      iterinduse=1;
    end
    for iterind=1:iterinduse
      tlock_TPA_MSPN{ii,iterind}=ft_math(cfg,tlock_tacPaud_each{ii,iterind},tlock_tacMSpN_each{ii,iterind});
    end
    if synchasynch && ll<5
      tlock_TMSs_TMSa{ii}=ft_math(cfg,tlock_tMSsynch_each{ii},tlock_tMSasynch_each{ii});
    end
    if audtacflag
      tlock_APT_MSPN{ii}=ft_math(cfg,tlock_audPtac_each{ii},tlock_audMSpN_each{ii});
      if ll<5
        tlock_AMSs_AMSa{ii}=ft_math(cfg,tlock_aMSsynch_each{ii},tlock_aMSasynch_each{ii});
      end
    end
    if usetr==2
      cfg.operation='(x1+x2)/2';
      tlock_TPA_MSPN_usetr2{ii}=ft_math(cfg,tlock_TPA_MSPN{ii,:});
      tlock_tacPaud_tmp{ii}=ft_math(cfg,tlock_tacPaud_each{ii,:});
      tlock_tacMSpN_tmp{ii}=ft_math(cfg,tlock_tacMSpN_each{ii,:});
      tlock_tactlock_tmp{ii}=ft_math(cfg,tlock_tactlock_each{ii,:});
      tlock_audtlock_tmp{ii}=ft_math(cfg,tlock_audtlock_each{ii,:});
      tlock_nulttlock_tmp{ii}=ft_math(cfg,tlock_nulttlock_each{ii,:});
      tlock_MStlock_tmp{ii}=ft_math(cfg,tlock_MStlock_each{ii,:});
    end
  end
  if usetr==2
    tlock_tacPaud_each=tlock_tacPaud_tmp;
    tlock_tacMSpN_each=tlock_tacMSpN_tmp;
    tlock_tactlock_each=tlock_tactlock_tmp;
    tlock_audtlock_each=tlock_audtlock_tmp;
    tlock_nulttlock_each=tlock_nulttlock_tmp;
    tlock_MStlock_each=tlock_MStlock_tmp;
  end
  cfg=[];
  cfg.channel=chanuse;
  if usetr<2
    grave_TPA_MSPN{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_TPA_MSPN{:});
  else
    grave_TPA_MSPN{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_TPA_MSPN_usetr2{:});
  end
  if synchasynch && ll<5
    grave_TMSs_TMSa{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_TMSs_TMSa{:});
  end
  if audtacflag
    grave_APT_MSPN{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_APT_MSPN{:});
    if ll<5
      grave_AMSs_AMSa{ll,tt,ss}=ft_timelockgrandaverage(cfg,tlock_AMSs_AMSa{:});
    end
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
  if usetr==2
    grind_TPA_MSPN=ft_timelockgrandaverage(cfg,tlock_TPA_MSPN_usetr2{:});
  end
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
  if usetr<2
    grave_tacVSnul=ft_timelockgrandaverage(cfg,tlock_tacVSnul_each{:});
    grave_audVSnul=ft_timelockgrandaverage(cfg,tlock_audVSnul_each{:});
    grave_msVSnul=ft_timelockgrandaverage(cfg,tlock_msVSnul_each{:});
  end
  
  
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
    if usetr==2
      grind_TPA_MSPN_zeros=grind_TPA_MSPN;
      grind_TPA_MSPN_zeros.individual=zeros(size(grind_TPA_MSPN.individual));
      statt_mc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_TPA_MSPN, grind_TPA_MSPN_zeros);
    else
      statt_mc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tacPaud, grind_tacMSpN);
    end
    if audtacflag
      stata_mc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_audPtac, grind_audMSpN);
    end
    %           if tacaudmsnul_earlyflag
    if usetr<2
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
    end % usetr
    
    if sleep
      for ii=1:subuseind
        KcPreMSpN_mean(ll,ii)=mean(KcPre_tacMSpN{ii});
        KcDuringMSpN_mean(ll,ii)=mean(KcDuring_tacMSpN{ii});
        KcEvokedMSpN_mean(ll,ii)=mean(KcEvoked_tacMSpN{ii});
        KcPreTpA_mean(ll,ii)=mean(KcPre_tacPaud{ii});
        KcDuringTpA_mean(ll,ii)=mean(KcDuring_tacPaud{ii});
        KcEvokedTpA_mean(ll,ii)=mean(KcEvoked_tacPaud{ii});
        SpPreMSpN_mean(ll,ii)=mean(SpPre_tacMSpN{ii});
        SpDuringMSpN_mean(ll,ii)=mean(SpDuring_tacMSpN{ii});
        SpEvokedMSpN_mean(ll,ii)=mean(SpEvoked_tacMSpN{ii});
        SpPreTpA_mean(ll,ii)=mean(SpPre_tacPaud{ii});
        SpDuringTpA_mean(ll,ii)=mean(SpDuring_tacPaud{ii});
        SpEvokedTpA_mean(ll,ii)=mean(SpEvoked_tacPaud{ii});
      end
      [P_KcPre(ll),H_KcPre(ll),STATS_KcPre{ll}] = signtest(KcPreMSpN_mean(ll,:),KcPreTpA_mean(ll,:));  % no assumptions on distribution using signtest
      [P_KcDuring(ll),H_KcDuring(ll),STATS_KcDuring{ll}] = signtest(KcDuringMSpN_mean(ll,:),KcDuringTpA_mean(ll,:));
      [P_KcEvoked(ll),H_KcEvoked(ll),STATS_KcEvoked{ll}] = signtest(KcEvokedMSpN_mean(ll,:),KcEvokedTpA_mean(ll,:));
      [P_SpPre(ll),H_SpPre(ll),STATS_SpPre{ll}] = signtest(SpPreMSpN_mean(ll,:),SpPreTpA_mean(ll,:));
      [P_SpDuring(ll),H_SpDuring(ll),STATS_SpDuring{ll}] = signtest(SpDuringMSpN_mean(ll,:),SpDuringTpA_mean(ll,:));
      [P_SpEvoked(ll),H_SpEvoked(ll),STATS_SpEvoked{ll}] = signtest(SpEvokedMSpN_mean(ll,:),SpEvokedTpA_mean(ll,:));
    end
    
    
    %   save([edir 'tlock_statmc.mat'],'stat*','grave*'); % no point, as grave* isn't ll,tt,ss dependent
    if iterflag
      if sleep
        %               if trialkcflag
        save([edir 'tlock_statmc_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '_statwinorig' num2str(statwinorig) '_ftver' num2str(ftver) 'mcseed' num2str(mcseed) '.mat'],'stat*','Kc*mean','Sp*mean','P_*','H_*','STATS_*');
        %               else
        %                 save([edir 'tlock_statmc_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '.mat'],'stat*');
        %               end
      else
        save([edir 'tlock_statmc_sleep' num2str(sleep) '_iter' num2str(iter) '_usetr' num2str(usetr) '_statwinorig' num2str(statwinorig) '_ftver' num2str(ftver) 'mcseed' num2str(mcseed) '.mat'],'stat*');
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
    keyboard
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
  
  if medsplitstatsflag && nsub>1

    
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


    
    
    load eeg1010_neighb
    cfg=[];
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
    cfg.channel=chanuse;
    cfg.neighbours=neighbours;
    cfg.parameter='individual';
    cfg.method='montecarlo';
    cfg.numrandomization=2000;
    cfg.correctm='cluster';
    cfg.clusteralpha = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan = 2;
    cfg.ivar=1;
    cfg.randomseed=mcseed;

    cfg.statistic='depsamplesT';
    % cfg.statistic='indepsamplesT';
    cfg.design=zeros(2,2*nsub);
    cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)];
    cfg.design(2,:)=[1:nsub 1:nsub];
    cfg.uvar=2;

    statt_mc{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tacPaud, grind_tacMSpN);
    
  end
  
  
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
    %     break
  else
    if sleep
      save([edir 'tlock_grind_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '_statwinorig' num2str(statwinorig) '_ftver' num2str(ftver) '.mat'],'grave*T*','grind_*save','plv');
    else
      save([edir 'tlock_grind_sleep' num2str(sleep) '_iter' num2str(iter) '_usetr' num2str(usetr) '_statwinorig' num2str(statwinorig) '_ftver' num2str(ftver) '.mat'],'grave*T*','grind_*save','plv');
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

