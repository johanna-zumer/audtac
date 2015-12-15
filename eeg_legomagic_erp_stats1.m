% statistics of EEG awake data

clear all
close all
if ispc
  edir='D:\audtac\eeg_data\';
  ddir='D:\audtac\legomagic\diaries\';
else
  edir='/mnt/hgfs/D/audtac/eeg_data/';
  ddir='/mnt/hgfs/D/audtac/legomagic/diaries/';
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
for ii=[5 iiuse]
  % for ii=5
  for sleep=[0 1]
    cd([edir sub{ii} ])
    clearvars -except ii sub edir ddir iiuse sleep featurestats*
    
    %   try
    %     load(['tlock_trialSel_' sub{ii} '.mat']);
    %   catch
    
    %     [tlock_tac_s0,tlock_aud_s0,tlock_nul_s0,tr]=eeg_legomagic_trialSelection1(ii);
    %     save(['tlock_trialSel_' sub{ii} '.mat'],'tlock*s0','tr','-v7.3');
    %     [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection1_wakeSleep(ii,sleep);
    
    for tt=1:4 % refers to lim=[no-limit .01 .005 .002];
      %   profile on
      [tlock_tac,tlock_aud,tlock_nul,tr]=eeg_legomagic_trialSelection2_wakeSleep(ii,sleep,tt);
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
              cfg=[];
              cfg.trials=tr.t10trialkept{ll,tt,12};
              tmp12=ft_selectdata(cfg,tlock_tac{10,tt,12}); % create +20 here as 10 will be different for every ll,tt
              cfg.trials=tr.t10trialkept{ll,tt,13};
              tmp13=ft_selectdata(cfg,tlock_tac{10,tt,13}); % create +20 here as 10 will be different for every ll,tt
              cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
              cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
              tlock_tac{ll+20,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
              if ll==max(soalist)
                tlock_tac{10,tt,12}=[];
                tlock_tac{10,tt,13}=[];
              end
              
              cfg=[];
              cfg.trials=tr.a10trialkept{ll,tt,12};
              tmp12=ft_selectdata(cfg,tlock_aud{10,tt,12}); % create +20 here as 10 will be different for every ll,tt
              cfg.trials=tr.a10trialkept{ll,tt,13};
              tmp13=ft_selectdata(cfg,tlock_aud{10,tt,13}); % create +20 here as 10 will be different for every ll,tt
              cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
              cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
              tlock_aud{ll+20,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
              if ll==max(soalist)
                tlock_aud{10,tt,12}=[];
                tlock_aud{10,tt,13}=[];
              end
              
              cfg=[];
              cfg.trials=tr.nllttrialkept{ll,tt,12};
              tmp12=ft_selectdata(cfg,tlock_nul{10,tt,12});
              cfg.trials=tr.nllttrialkept{ll,tt,13};
              tmp13=ft_selectdata(cfg,tlock_nul{10,tt,13});
              cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
              cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
              tlock_nul{ll+50,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
              
              cfg=[];
              cfg.trials=tr.nllatrialkept{ll,tt,12};
              tmp12=ft_selectdata(cfg,tlock_nul{10,tt,12});
              cfg.trials=tr.nllatrialkept{ll,tt,13};
              tmp13=ft_selectdata(cfg,tlock_nul{10,tt,13});
              cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
              cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
              tlock_nul{ll+60,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
              if ll==max(soalist)
                tlock_nul{10,tt,12}=[];
                tlock_nul{10,tt,13}=[];
              end
              
            else
              cfg=[];
              cfg.trials=tr.t10trialkept{ll,tt,ss};
              tlock_tac{ll+20,tt,ss}=ft_selectdata(cfg,tlock_tac{10,tt,ss}); % create +20 here as 10 will be different for every ll,tt
              if ll==max(soalist) && ss<12
                tlock_tac{10,tt,ss}=[];
              end
              
              cfg=[];
              cfg.trials=tr.a10trialkept{ll,tt,ss};
              tlock_aud{ll+20,tt,ss}=ft_selectdata(cfg,tlock_aud{10,tt,ss}); % create +20 here as 10 will be different for every ll,tt
              if ll==max(soalist) && ss<12
                tlock_aud{10,tt,ss}=[];
              end
              
              cfg=[];
              cfg.trials=tr.nllttrialkept{ll,tt,ss};
              tlock_nul{ll+50,tt,ss}=ft_selectdata(cfg,tlock_nul{10,tt,ss});
              cfg=[];
              cfg.trials=tr.nllatrialkept{ll,tt,ss};
              tlock_nul{ll+60,tt,ss}=ft_selectdata(cfg,tlock_nul{10,tt,ss});
              if ll==max(soalist) && ss<12
                tlock_nul{10,tt,ss}=[];
              end
              
            end
          end
          
          if ll<10
            numt_trials(ll,tt,ss)=size(tlock_tac{ll+20,tt,ss}.trial,1); % this should be the same for tacll20, audll40, tacll, nulll50
            numa_trials(ll,tt,ss)=size(tlock_aud{ll+20,tt,ss}.trial,1); % this should be the same for audll20, tacll40, audll, nulll60
          end
          
          if ss==23 && [ll<10 || ll>40]
            tmp12=tlock_tac{ll,tt,12};
            tmp13=tlock_tac{ll,tt,13};
            cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
            cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
            tlock_tac{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
            tlock_tac{ll,tt,12}=[];
            tlock_tac{ll,tt,13}=[];
            
            tmp12=tlock_aud{ll,tt,12};
            tmp13=tlock_aud{ll,tt,13};
            cfg=[];cfg.latency=[tmp12.time(1) tmp12.time(end)];tmp13=ft_selectdata(cfg,tmp13);
            cfg=[];cfg.latency=[tmp13.time(1) tmp13.time(end)];tmp12=ft_selectdata(cfg,tmp12);
            tlock_aud{ll,tt,ss}=ft_appendtimelock([],tmp12,tmp13);
            tlock_aud{ll,tt,12}=[];
            tlock_aud{ll,tt,13}=[];
          end
          
          if ll>40
            if ~isempty(tlock_tac{ll,tt,ss})
              tlock_tac{ll,tt,ss}=trim_nans(tlock_tac{ll,tt,ss});
            end
            if ~isempty(tlock_aud{ll,tt,ss})
              tlock_aud{ll,tt,ss}=trim_nans(tlock_aud{ll,tt,ss});
            end
          end
          
          % 40Hz lowpass filter (already 0.2Hz highpass on non-epoched data)
          % also do baseline correction here
          if isfield(tlock_tac{ll,tt,ss},'trial') && size(tlock_tac{ll,tt,ss}.trial,1)
            cfg=[];
            cfg.lpfilter='yes';
            cfg.lpfreq=40;
            cfg.demean='yes';
            cfg.baselinewindow=[-1.7 -0.6];
            disp('ft_preprocessing tac')
            tlock_tac{ll,tt,ss}=ft_preprocessing(cfg,tlock_tac{ll,tt,ss});
            featurestats_tac(:,ll,tt,ss,ii)=mean(tlock_tac{ll,tt,ss}.trialinfo(:,12:19));
            cfg=[];
            tlock_tactlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_tac{ll,tt,ss});
          else
            tlock_tactlock{ll,tt,ss}=[];
          end
          if ss==12 || ss==13
          else
            tlock_tac{ll,tt,ss}=[]; % clearing to help save memory as running
          end
          
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
            featurestats_aud(:,ll,tt,ss,ii)=mean(tlock_aud{ll,tt,ss}.trialinfo(:,12:19));
            cfg=[];
            tlock_audtlock{ll,tt,ss}=ft_timelockanalysis(cfg,tlock_aud{ll,tt,ss});
          else
            tlock_audtlock{ll,tt,ss}=[];
          end
          if ss==12 || ss==13
          else
            tlock_aud{ll,tt,ss}=[];
          end
          %           catch
          %             disp('didnt work for this ll')
          %           end
        end % end ll
        
        
        for ll=[soalist+50 soalist+60]
          if isfield(tlock_nul{ll,tt,ss},'trial') && size(tlock_nul{ll,tt,ss}.trial,1)
            cfg=[];
            cfg.lpfilter='yes';
            cfg.lpfreq=40;
            cfg.demean='yes';
            cfg.baselinewindow=[-1.7 -0.6];
            %           cfg.bpfilter='yes';
            %           cfg.bpfreq=[1 40];
            disp('ft_preprocessing nul')
            tlock_nul{ll,tt,ss}=ft_preprocessing(cfg,tlock_nul{ll,tt,ss});
            featurestats_nul(:,ll,tt,ss,ii)=mean(tlock_nul{ll,tt,ss}.trialinfo(:,12:19));
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
          if numt_trials(ll,tt,ss)
            tlock_tacPaud{ll,tt,ss}=ft_math(cfg,tlock_tactlock{ll+20,tt,ss},tlock_audtlock{ll+40,tt,ss}); % keeping tac centred but jitter aud (to match tlock_tac{X,tt})
            featurestats_tacPaud(:,ll,tt,ss,ii)=mean([featurestats_tac(:,ll+20,tt,ss,ii), featurestats_aud(:,ll+40,tt,ss,ii)]');
          else
            tlock_tacPaud{ll,tt,ss}=[];
          end
          if numa_trials(ll,tt,ss)
            tlock_audPtac{ll,tt,ss}=ft_math(cfg,tlock_audtlock{ll+20,tt,ss},tlock_tactlock{ll+40,tt,ss}); % keeping aud centred but jitter tac (to match tlock_aud{X,tt})
            featurestats_audPtac(:,ll,tt,ss,ii)=mean([featurestats_aud(:,ll+20,tt,ss,ii), featurestats_tac(:,ll+40,tt,ss,ii)]');
          else
            tlock_audPtac{ll,tt,ss}=[];
          end
          
          % create sum of multisensory plus nul conditions
          cfg=[];
          cfg.operation='add';
          cfg.parameter='avg';
          if numt_trials(ll,tt,ss)
            tlock_tacMSpN{ll,tt,ss}=ft_math(cfg,tlock_tactlock{ll,tt,ss},tlock_nultlock{ll+50,tt,ss}); % TA+N
            featurestats_tacMSpN(:,ll,tt,ss,ii)=mean([featurestats_tac(:,ll,tt,ss,ii), featurestats_nul(:,ll+50,tt,ss,ii)]');
          else
            tlock_tacMSpN{ll,tt,ss}=[];
          end
          if numa_trials(ll,tt,ss)
            tlock_audMSpN{ll,tt,ss}=ft_math(cfg,tlock_audtlock{ll,tt,ss},tlock_nultlock{ll+60,tt,ss}); % AT+N
            featurestats_audMSpN(:,ll,tt,ss,ii)=mean([featurestats_aud(:,ll,tt,ss,ii), featurestats_nul(:,ll+60,tt,ss,ii)]');
          else
            tlock_audMSpN{ll,tt,ss}=[];
          end
        end % end ll
      
             
        
      end  % end ss
      
      clear tlock_tac tlock_aud tlock_nul tr
    end
    %   save(['tlock_diffs_averef_' sub{ii} '.mat'],'*trialkept','tlock_*P*','tlock_*N*')
    save(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '.mat'],'tlock_*P*','tlock_*N*','num*trials','featurestats_*')
    
    
  end
end

return


%% Group level

allcond_sameN=1; % means using ii>=8, but we doing that anyway with iiuse
% for ll=[1 3 4 5 6 7 9]
for ll=[3 4 5 6 7]
  for tt=1:4
    for sleep=[0 1]
      clearvars -except allcond_sameN ll tt sub edir ddir iiuse sleep
      
      %     if ll==1 | ll==9
      %       subuse=8:32;
      %     else
      %       subuse=5:32;
      %     end
      %     if allcond_sameN
      %       subuse=8:32;
      %     end
      subuse=iiuse;
      
      submin=subuse(1)-1;
      subuseind=0;
      for ii=subuse
        cd([edir sub{ii} ])
        %   load(['tlock_diffs_' sub{ii} '.mat']);
        %       load(['tlock_diffs_averef_' sub{ii} '.mat']);
        load(['tlock_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '.mat'])
        load(['trialkept_tt' num2str(tt) '_sleep' num2str(sleep) '.mat'],'tr')
        
        % THis is preliminary...how best to included all stages later on?
        if sleep==0
          ss=10; % awake
        elseif sleep==1
          ss=23; % this is concatenation of N2 and N3
        end
        
        %         for ss=ssuse
        numtrt(ll,tt,ss,ii-submin)=numt_trials(ll,tt,ss); % does this need to be per 'sleep01' as well?
        numtra(ll,tt,ss,ii-submin)=numa_trials(ll,tt,ss);
        %         end
        
        
        if sum(numtrt(ll,tt,ss,ii-submin),3)<20 % what is best number to use here?
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
          tlock_tacPaud_each{subuseind}=tlock_tacPaud{ll,tt,ss};
          tlock_audPtac_each{subuseind}=tlock_audPtac{ll,tt,ss};
          tlock_tacMSpN_each{subuseind}=tlock_tacMSpN{ll,tt,ss};
          tlock_audMSpN_each{subuseind}=tlock_audMSpN{ll,tt,ss};
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
          
        end
        %       clear *_s0
        clear tlock*N tlock*tac tlock*aud
      end
      subuseindfinal=subuseind;
      
      figure(20);
      %     plotind=1;
      %     for ii=1:length(subuse)
      for ii=1:subuseindfinal
        cfg=[];
        cfg.parameter='avg';
        cfg.operation='subtract';
        diff{ii}=ft_math(cfg,tlock_tacPaud_each{ii},tlock_tacMSpN_each{ii});
        subplot(3,subuseindfinal,ii);imagesc(tlock_tacPaud_each{ii}.time,1:63,tlock_tacPaud_each{ii}.avg);caxis([-6 6])
        subplot(3,subuseindfinal,subuseindfinal+ii);imagesc(tlock_tacMSpN_each{ii}.time,1:63,tlock_tacMSpN_each{ii}.avg);caxis([-6 6])
        subplot(3,subuseindfinal,2*subuseindfinal+ii);imagesc(diff{ii}.time,1:63,diff{ii}.avg);caxis([-6 6])
        %       plotind=plotind+1;
      end
      %       if allcond_sameN
      print(20,['D:\audtac\figs\sumuni_mspn_diff_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
      %       else
      %         print(20,['D:\audtac\figs\sumuni_mspn_diff_ta_cond' num2str(ll) num2str(tt) 'a.png'],'-dpng')
      %       end
      figure(21);
      %     plotind=1;
      for ii=1:subuseindfinal
        cfg=[];
        cfg.parameter='avg';
        cfg.operation='subtract';
        diff{ii}=ft_math(cfg,tlock_audPtac_each{ii},tlock_audMSpN_each{ii});
        subplot(3,subuseindfinal,ii);imagesc(tlock_audPtac_each{ii}.time,1:63,tlock_audPtac_each{ii}.avg);caxis([-6 6])
        subplot(3,subuseindfinal,subuseindfinal+ii);imagesc(tlock_audMSpN_each{ii}.time,1:63,tlock_audMSpN_each{ii}.avg);caxis([-6 6])
        subplot(3,subuseindfinal,2*subuseindfinal+ii);imagesc(diff{ii}.time,1:63,diff{ii}.avg);caxis([-6 6])
        %       plotind=plotind+1;
      end
      %       if allcond_sameN
      print(21,['D:\audtac\figs\sumuni_mspn_diff_at_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
      %       else
      %         print(21,['D:\audtac\figs\sumuni_mspn_diff_at_cond' num2str(ll) num2str(tt) 'a.png'],'-dpng')
      %       end
      
      cfg=[];
      cfg.keepindividual='yes';
      grind_tacPaud=ft_timelockgrandaverage(cfg,tlock_tacPaud_each{:});
      grind_audPtac=ft_timelockgrandaverage(cfg,tlock_audPtac_each{:});
      grind_tacMSpN=ft_timelockgrandaverage(cfg,tlock_tacMSpN_each{:});
      grind_audMSpN=ft_timelockgrandaverage(cfg,tlock_audMSpN_each{:});
      
      
      cfg=[];
      grave_tacPaud=ft_timelockgrandaverage(cfg,tlock_tacPaud_each{:});
      grave_audPtac=ft_timelockgrandaverage(cfg,tlock_audPtac_each{:});
      grave_tacMSpN=ft_timelockgrandaverage(cfg,tlock_tacMSpN_each{:});
      grave_audMSpN=ft_timelockgrandaverage(cfg,tlock_audMSpN_each{:});
      
      topoplot_highlight(11,grave_tacPaud,[-0.5 0.6],[]);
      topoplot_highlight(12,grave_audPtac,[-0.1 1.1],[]);
      topoplot_highlight(13,grave_tacMSpN,[-0.5 0.6],[]);
      topoplot_highlight(14,grave_audMSpN,[-0.1 1.1],[]);
      
      %       if allcond_sameN
      print(11,['D:\audtac\figs\gravetacPaud_topoOverTime_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
      %       else
      %         print(11,['D:\audtac\figs\gravetacPaud_topoOverTime_cond' num2str(ll) num2str(tt)  '_28.png'],'-dpng')
      %       end
      %       if allcond_sameN
      print(12,['D:\audtac\figs\graveaudPtac_topoOverTime_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
      %       else
      %         print(12,['D:\audtac\figs\graveaudPtac_topoOverTime_cond' num2str(ll) num2str(tt)  '_28.png'],'-dpng')
      %       end
      %       if allcond_sameN
      print(13,['D:\audtac\figs\gravetacMSpN_topoOverTime_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
      %       else
      %         print(13,['D:\audtac\figs\gravetacMSpN_topoOverTime_cond' num2str(ll) num2str(tt)  '_28.png'],'-dpng')
      %       end
      %       if allcond_sameN
      print(14,['D:\audtac\figs\graveaudMSpN_topoOverTime_cond' num2str(ll) num2str(tt) num2str(sleep)  '.png'],'-dpng')
      %       else
      %         print(14,['D:\audtac\figs\graveaudMSpN_topoOverTime_cond' num2str(ll) num2str(tt)  '_28.png'],'-dpng')
      %       end
      
      
      
      for ii=1:subuseindfinal
        cfg=[];
        cfg.operation='subtract';
        cfg.parameter='avg';
        tlock_TPA_MSPN{ii}=ft_math(cfg,tlock_tacPaud_each{ii},tlock_tacMSpN_each{ii});
        tlock_APT_MSPN{ii}=ft_math(cfg,tlock_audPtac_each{ii},tlock_audMSpN_each{ii});
      end
      cfg=[];
      grave_TPA_MSPN=ft_timelockgrandaverage(cfg,tlock_TPA_MSPN{:});
      grave_APT_MSPN=ft_timelockgrandaverage(cfg,tlock_APT_MSPN{:});
      
      
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
      
      load eeg1010_neighb
      
      nsub=length(tlock_audMSpN_each);
      
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
      statt_mc=ft_timelockstatistics(cfg, grind_tacPaud, grind_tacMSpN);
      stata_mc=ft_timelockstatistics(cfg, grind_audPtac, grind_audMSpN);
      
      
      topoplot_highlight(22,grave_TPA_MSPN,[statt_mc.time(1) statt_mc.time(end)],statt_mc);
      
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
      %       %      keyboard
      %       cfg.comment = 'xlim';
      %       cfg.commentpos = 'title';
      %       cfg.layout = 'elec1010.lay';
      %       ft_topoplotER(cfg, grave_TPA_MSPN);
      %     end
      %       if allcond_sameN
      print(22,['D:\audtac\figs\grdiff_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
      %       else
      %         print(22,['D:\audtac\figs\grdiff_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) 'a.png'],'-dpng')
      %       end
      
      topoplot_highlight(23,grave_APT_MSPN,[stata_mc.time(1) stata_mc.time(end)],stata_mc);
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
      %       %      keyboard
      %       cfg.comment = 'xlim';
      %       cfg.commentpos = 'title';
      %       cfg.layout = 'elec1010.lay';
      %       ft_topoplotER(cfg, grave_APT_MSPN);
      %     end
      %       if allcond_sameN
      print(23,['D:\audtac\figs\grdiff_topoOverTime_mc_at_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
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
      %       %      keyboard
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
      %       %      keyboard
      %       cfg.comment = 'xlim';
      %       cfg.commentpos = 'title';
      %       cfg.layout = 'elec1010.lay';
      %       ft_topoplotER(cfg, grave_APT_MSPN);
      %     end
      %     print(25,['D:\audtac\figs\grdiff_topoOverTime_an_at_cond' num2str(ll) num2str(tt) 'a.png'],'-dpng')
    end
  end
end
save([edir 'tlock_numtrlltt.mat'],'numtr*','grind*');


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
