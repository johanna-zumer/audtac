eeg_legomagic_preamble

%%
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
      cfg.parameter={'powspctrm' 'plvspctrm' 'plvabs','stdpow','stdplv'};
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
      cfg.parameter={'powspctrm' 'plvspctrm' 'plvabs','stdpow','stdplv'};
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
        cfg.parameter={'powspctrm' 'plvspctrm' 'plvabs'};
        freqloall_TPA_MSPN_comb2{ii}=ft_math(cfg,freqloall_tacPaud_comb2{ii},freqloall_tacMSpN_comb2{ii});
        freqhiall_TPA_MSPN_comb2{ii}=ft_math(cfg,freqhiall_tacPaud_comb2{ii},freqhiall_tacMSpN_comb2{ii});
        
        if audtacflag
          freqloall_APT_MSPN_comb2{ii}=ft_math(cfg,freqloall_audPtac_comb2{ii},freqloall_audMSpN_comb2{ii});
          freqhiall_APT_MSPN_comb2{ii}=ft_math(cfg,freqhiall_audPtac_comb2{ii},freqhiall_audMSpN_comb2{ii});
        end
      end % end ii
      
      cfg=[];
      cfg.parameter={'powspctrm' 'plvspctrm' 'plvabs'};
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


