eeg_legomagic_preamble

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

