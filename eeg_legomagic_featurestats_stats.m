
eeg_legomagic_preamble

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


