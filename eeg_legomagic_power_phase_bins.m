function [angbin,absbin,freqtr]=eeg_legomagic_power_phase_bins(ii,sub,sleep,ss,trialkc,phaset0use)
%
% INPUTS
%  ii  = subject index
%  sub  = subject structure index (e.g.  sub{8}='e08')
%  sleep =  1 for lying/bed dataset, 0 for sitting/awake dataset
%  ss   =  sleep stage (10 = W, 11=N1, 12=n2, 13=N3)
%  trialkc =  -1 for use all trials, 0 for non-Kc trials, and 1 for Kc-only trials
%  phaset0use = method for computing phase.   1,2,3 all use FFT, 4,5,6 use filter+Hilbert
%              1 and 4 use Cz sensor only
%              2 and 5 use a cluster of frontocentral channels
%              3 and 6 use a PCA component most correlated to Cz
%              I've settled on using "2" always now, but left all the other code in just in case.
%
% OUTPUTS
%  angbin = structure of binary index (size: trials X bins) indicating if that trial is in that phase bin or not
%  absbin = same as above but for power bins
%  freqtr = containing trial indices
%
% Note 1: this is called by Johanna in eeg_legomagic_erp_stats2_sepTacAud.m 
%       and called by Tom in .....

% Note 2:  See bottom after 'return' for additional info/tips


tt=3; % we are for sure only use tt=3 now
soalist=[1 3 4 5 6 7 9];

switch phaset0use
  case {1 2 3}
    if trialkc==0
      load(['freqtime0_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_tt' num2str(tt) '.mat'],'freq*');
      freqtr=load(['freqtime0_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_tt' num2str(tt) '.mat'],'tr');
    else
      try
        load(['freqtime0_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_tt' num2str(tt) '_trialkc' num2str(trialkc) '.mat'],'freq*');
        freqtr=load(['freqtime0_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_tt' num2str(tt) '_trialkc' num2str(trialkc) '.mat'],'tr');
      catch
        if sleep==0 && trialkc==-1
          load(['freqtime0_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_tt' num2str(tt) '.mat'],'freq*');
          freqtr=load(['freqtime0_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_tt' num2str(tt) '.mat'],'tr');
        end
      end
    end
  case {4 5 6}
    load(['tlock_unimultsensHilbert_' sub{ii} '_sleep0_tt' num2str(tt) '_tacaud1_iter' num2str(iter) '_trialkc-1.mat']);
  otherwise
    error('not valid phaset0use')
end


for ll=soalist
  switch phaset0use
    case 1
      audang{ll,tt,ss}=angle(freqlo_tAudAlone_time0{ll,tt,ss}.fourierspctrm(:,18,10,9));
    case 2
      
      % choose peak freq per subject
      if sleep==0
        freqtest=7:13;
        [mx,avemind]=max(mean([squeeze(mean(abs(freqlo_tac4mscon1_chanFC_time0{1,3,10}.fourierspctrm(:,:,7:13)),1)) squeeze(mean(abs(freqlo_aud4mscon1_chanFC_time0{1,3,10}.fourierspctrm(:,:,7:13)),1)) squeeze(mean(abs(freqlo_ms4mscon1_first_chanFC_time0{1,3,10}.fourierspctrm(:,:,7:13)),1)) squeeze(mean(abs(freqlo_nul4mscon1_chanFC_time0{1,3,10}.fourierspctrm(:,:,7:13)),1))],2));
        freqsub(ii)=freqtest(avemind);
        freqsleep(ii)=10;
        timeind=9;
        timesleep=1;
        timesub=1;
      elseif sleep==1
        freqtest=1:4;
        [mx,avemind1]=max(mean([squeeze(mean(abs(freqlo_tac4mscon1_chanFC_time0{1,3,12}.fourierspctrm(:,1,1:4,1)),1)) squeeze(mean(abs(freqlo_aud4mscon1_chanFC_time0{1,3,12}.fourierspctrm(:,1,1:4,1)),1)) squeeze(mean(abs(freqlo_ms4mscon1_first_chanFC_time0{1,3,12}.fourierspctrm(:,1,1:4,1)),1)) squeeze(mean(abs(freqlo_nul4mscon1_chanFC_time0{1,3,12}.fourierspctrm(:,1,1:4,1)),1))],2));
        
        %                     % based on simulations, we found this to be indicative of the predominant frequency
        %                     [mx,avemind2]=max(mean([squeeze(abs(mean(freqlo_tac4mscon1_chanFC_time0{1,3,12}.fourierspctrm(:,1,1:4,1),1))) squeeze(abs(mean(freqlo_aud4mscon1_chanFC_time0{1,3,12}.fourierspctrm(:,1,1:4,1),1))) squeeze(abs(mean(freqlo_ms4mscon1_first_chanFC_time0{1,3,12}.fourierspctrm(:,1,1:4,1),1))) squeeze(abs(mean(freqlo_nul4mscon1_chanFC_time0{1,3,12}.fourierspctrm(:,1,1:4,1),1)))],2));
        %                     freqsleep(ii)=freqtest(avemind2);
        freqsleep(ii)=2;
        freqsub(ii)=freqtest(avemind1);
        timeind=1;
        timesleep=2; % -250ms for 2Hz for phase
        timesleepabs=1; % -500ms for 2Hz for abs
        timesub=avemind1;
        timesubabs=max(1,avemind1-1);
      end
      
      nulang{ll,tt,ss}=angle(freqlo_tNulAlone_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timeind));
      tacang{ll,tt,ss}=angle(freqlo_tTacAlone_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timeind));
      audang{ll,tt,ss}=angle(freqlo_tAudAlone_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timeind));
      ms1ang{ll,tt,ss}=angle(freqlo_tMSAlone_first_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timeind));
      ms2ang{ll,tt,ss}=angle(freqlo_tMSAlone_second_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timeind));
      tac4mscon1ang{ll,tt,ss}=angle(freqlo_tac4mscon1_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timesleep));
      tac4mscon2ang{ll,tt,ss}=angle(freqlo_tac4mscon2_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timesleep));
      aud4mscon1ang{ll,tt,ss}=angle(freqlo_aud4mscon1_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timesleep));
      aud4mscon2ang{ll,tt,ss}=angle(freqlo_aud4mscon2_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timesleep));
      ms4mscon1ang{ll,tt,ss}=angle(freqlo_ms4mscon1_first_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timesleep));
      ms4mscon2ang{ll,tt,ss}=angle(freqlo_ms4mscon2_first_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timesleep));
      nul4mscon1ang{ll,tt,ss}=angle(freqlo_nul4mscon1_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timesleep));
      nul4mscon2ang{ll,tt,ss}=angle(freqlo_nul4mscon2_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timesleep));
      
% % % Do not delete, but not using anymore
%       nulangIAF{ll,tt,ss}=angle(freqlo_tNulAlone_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timeind));
%       tacangIAF{ll,tt,ss}=angle(freqlo_tTacAlone_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timeind));
%       audangIAF{ll,tt,ss}=angle(freqlo_tAudAlone_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timeind));
%       ms1angIAF{ll,tt,ss}=angle(freqlo_tMSAlone_first_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timeind));
%       ms2angIAF{ll,tt,ss}=angle(freqlo_tMSAlone_second_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timeind));
%       tac4mscon1angIAF{ll,tt,ss}=angle(freqlo_tac4mscon1_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timesub));
%       tac4mscon2angIAF{ll,tt,ss}=angle(freqlo_tac4mscon2_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timesub));
%       aud4mscon1angIAF{ll,tt,ss}=angle(freqlo_aud4mscon1_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timesub));
%       aud4mscon2angIAF{ll,tt,ss}=angle(freqlo_aud4mscon2_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timesub));
%       ms4mscon1angIAF{ll,tt,ss}=angle(freqlo_ms4mscon1_first_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timesub));
%       ms4mscon2angIAF{ll,tt,ss}=angle(freqlo_ms4mscon2_first_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timesub));
%       nul4mscon1angIAF{ll,tt,ss}=angle(freqlo_nul4mscon1_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timesub));
%       nul4mscon2angIAF{ll,tt,ss}=angle(freqlo_nul4mscon2_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timesub));
      
      % power
      nulabs{ll,tt,ss}=abs(freqlo_tNulAlone_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timeind));
      tacabs{ll,tt,ss}=abs(freqlo_tTacAlone_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timeind));
      audabs{ll,tt,ss}=abs(freqlo_tAudAlone_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timeind));
      ms1abs{ll,tt,ss}=abs(freqlo_tMSAlone_first_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timeind));
      ms2abs{ll,tt,ss}=abs(freqlo_tMSAlone_second_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timeind));
      tac4mscon1abs{ll,tt,ss}=abs(freqlo_tac4mscon1_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timesleepabs));
      tac4mscon2abs{ll,tt,ss}=abs(freqlo_tac4mscon2_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timesleepabs));
      aud4mscon1abs{ll,tt,ss}=abs(freqlo_aud4mscon1_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timesleepabs));
      aud4mscon2abs{ll,tt,ss}=abs(freqlo_aud4mscon2_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timesleepabs));
      ms4mscon1abs{ll,tt,ss}=abs(freqlo_ms4mscon1_first_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timesleepabs));
      ms4mscon2abs{ll,tt,ss}=abs(freqlo_ms4mscon2_first_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timesleepabs));
      nul4mscon1abs{ll,tt,ss}=abs(freqlo_nul4mscon1_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timesleepabs));
      nul4mscon2abs{ll,tt,ss}=abs(freqlo_nul4mscon2_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsleep(ii),timesleepabs));
      
% % % Do not delete, but not using anymore
%       nulabsIAF{ll,tt,ss}=abs(freqlo_tNulAlone_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timeind));
%       tacabsIAF{ll,tt,ss}=abs(freqlo_tTacAlone_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timeind));
%       audabsIAF{ll,tt,ss}=abs(freqlo_tAudAlone_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timeind));
%       ms1absIAF{ll,tt,ss}=abs(freqlo_tMSAlone_first_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timeind));
%       ms2absIAF{ll,tt,ss}=abs(freqlo_tMSAlone_second_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timeind));
%       tac4mscon1absIAF{ll,tt,ss}=abs(freqlo_tac4mscon1_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timesubabs));
%       tac4mscon2absIAF{ll,tt,ss}=abs(freqlo_tac4mscon2_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timesubabs));
%       aud4mscon1absIAF{ll,tt,ss}=abs(freqlo_aud4mscon1_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timesubabs));
%       aud4mscon2absIAF{ll,tt,ss}=abs(freqlo_aud4mscon2_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timesubabs));
%       ms4mscon1absIAF{ll,tt,ss}=abs(freqlo_ms4mscon1_first_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timesubabs));
%       ms4mscon2absIAF{ll,tt,ss}=abs(freqlo_ms4mscon2_first_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timesubabs));
%       nul4mscon1absIAF{ll,tt,ss}=abs(freqlo_nul4mscon1_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timesubabs));
%       nul4mscon2absIAF{ll,tt,ss}=abs(freqlo_nul4mscon2_chanFC_time0{ll,tt,ss}.fourierspctrm(:,1,freqsub(ii),timesubabs));
    case 3
      audang{ll,tt,ss}=angle(freqlo_tAudAlone_chanPCA_time0{ll,tt,ss}.fourierspctrm(:,1,10,9));
    case [4 5 6]
    otherwise
      error('not valid phaset0use')
  end
  
  % Bin by phase:
  angbin.tacbin{ll,tt,ss}=bin_my_angle(tacang{ll,tt,ss},4);
  angbin.audbin{ll,tt,ss}=bin_my_angle(audang{ll,tt,ss},4);
  angbin.nulbin{ll,tt,ss}=bin_my_angle(nulang{ll,tt,ss},4);
  angbin.ms1bin{ll,tt,ss}=bin_my_angle(ms1ang{ll,tt,ss},4);
  angbin.ms2bin{ll,tt,ss}=bin_my_angle(ms2ang{ll,tt,ss},4);
  angbin.tac4mscon1bin{ll,tt,ss}=bin_my_angle(tac4mscon1ang{ll,tt,ss},4);
  angbin.tac4mscon2bin{ll,tt,ss}=bin_my_angle(tac4mscon2ang{ll,tt,ss},4);
  angbin.aud4mscon1bin{ll,tt,ss}=bin_my_angle(aud4mscon1ang{ll,tt,ss},4);
  angbin.aud4mscon2bin{ll,tt,ss}=bin_my_angle(aud4mscon2ang{ll,tt,ss},4);
  angbin.nul4mscon1bin{ll,tt,ss}=bin_my_angle(nul4mscon1ang{ll,tt,ss},4);
  angbin.nul4mscon2bin{ll,tt,ss}=bin_my_angle(nul4mscon2ang{ll,tt,ss},4);
  angbin.ms4mscon1bin{ll,tt,ss}=bin_my_angle(ms4mscon1ang{ll,tt,ss},4);
  angbin.ms4mscon2bin{ll,tt,ss}=bin_my_angle(ms4mscon2ang{ll,tt,ss},4);
  
% % Do not delete, but not using anymore
%   tacIAFbin{ll,tt,ss}=bin_my_angle(tacangIAF{ll,tt,ss},4);
%   audIAFbin{ll,tt,ss}=bin_my_angle(audangIAF{ll,tt,ss},4);
%   nulIAFbin{ll,tt,ss}=bin_my_angle(nulangIAF{ll,tt,ss},4);
%   ms1IAFbin{ll,tt,ss}=bin_my_angle(ms1angIAF{ll,tt,ss},4);
%   ms2IAFbin{ll,tt,ss}=bin_my_angle(ms2angIAF{ll,tt,ss},4);
%   tac4mscon1IAFbin{ll,tt,ss}=bin_my_angle(tac4mscon1angIAF{ll,tt,ss},4);
%   tac4mscon2IAFbin{ll,tt,ss}=bin_my_angle(tac4mscon2angIAF{ll,tt,ss},4);
%   aud4mscon1IAFbin{ll,tt,ss}=bin_my_angle(aud4mscon1angIAF{ll,tt,ss},4);
%   aud4mscon2IAFbin{ll,tt,ss}=bin_my_angle(aud4mscon2angIAF{ll,tt,ss},4);
%   nul4mscon1IAFbin{ll,tt,ss}=bin_my_angle(nul4mscon1angIAF{ll,tt,ss},4);
%   nul4mscon2IAFbin{ll,tt,ss}=bin_my_angle(nul4mscon2angIAF{ll,tt,ss},4);
%   ms4mscon1IAFbin{ll,tt,ss}=bin_my_angle(ms4mscon1angIAF{ll,tt,ss},4);
%   ms4mscon2IAFbin{ll,tt,ss}=bin_my_angle(ms4mscon2angIAF{ll,tt,ss},4);
  

  % Bin by power:
  % Set same power limit within participant for all data going in to same contrast
  
  [sortabs,sortind]=sort([tacabs{ll,tt,ss}' audabs{ll,tt,ss}' nulabs{ll,tt,ss}' ms1abs{ll,tt,ss}' ms2abs{ll,tt,ss}']);
  numperbin=floor(length(sortabs)/4);
  for nn=2:4
    xbin(nn)=sortabs(numperbin*(nn-1)+1);
  end
  absbin.tacabsbin{ll,tt,ss}=bin_my_abs(tacabs{ll,tt,ss},xbin);
  absbin.audabsbin{ll,tt,ss}=bin_my_abs(audabs{ll,tt,ss},xbin);
  absbin.nulabsbin{ll,tt,ss}=bin_my_abs(nulabs{ll,tt,ss},xbin);
  absbin.ms1absbin{ll,tt,ss}=bin_my_abs(ms1abs{ll,tt,ss},xbin);
  absbin.ms2absbin{ll,tt,ss}=bin_my_abs(ms2abs{ll,tt,ss},xbin);
  
  [sortabs,sortind]=sort([tac4mscon1abs{ll,tt,ss}' aud4mscon1abs{ll,tt,ss}' nul4mscon1abs{ll,tt,ss}' ms4mscon1abs{ll,tt,ss}']);
  numperbin=floor(length(sortabs)/4);
  for nn=2:4
    xbin(nn)=sortabs(numperbin*(nn-1)+1);
  end
  absbin.tac4mscon1absbin{ll,tt,ss}=bin_my_abs(tac4mscon1abs{ll,tt,ss},xbin);
  absbin.aud4mscon1absbin{ll,tt,ss}=bin_my_abs(aud4mscon1abs{ll,tt,ss},xbin);
  absbin.nul4mscon1absbin{ll,tt,ss}=bin_my_abs(nul4mscon1abs{ll,tt,ss},xbin);
  absbin.ms4mscon1absbin{ll,tt,ss}=bin_my_abs(ms4mscon1abs{ll,tt,ss},xbin);  
  %               powdist{ii,ll}=sortabs;
  
  [sortabs,sortind]=sort([tac4mscon2abs{ll,tt,ss}' aud4mscon2abs{ll,tt,ss}' nul4mscon2abs{ll,tt,ss}' ms4mscon2abs{ll,tt,ss}']);
  numperbin=floor(length(sortabs)/4);
  for nn=2:4
    xbin(nn)=sortabs(numperbin*(nn-1)+1);
  end
  absbin.tac4mscon2absbin{ll,tt,ss}=bin_my_abs(tac4mscon2abs{ll,tt,ss},xbin);
  absbin.aud4mscon2absbin{ll,tt,ss}=bin_my_abs(aud4mscon2abs{ll,tt,ss},xbin);
  absbin.nul4mscon2absbin{ll,tt,ss}=bin_my_abs(nul4mscon2abs{ll,tt,ss},xbin);
  absbin.ms4mscon2absbin{ll,tt,ss}=bin_my_abs(ms4mscon2abs{ll,tt,ss},xbin);
  

%   % Interesting idea, but not used
%   [sortabs,sortind]=sort([tac4mscon1abs{5,tt,ss}' aud4mscon1abs{5,tt,ss}' nul4mscon1abs{5,tt,ss}' ms4mscon1abs{1,tt,ss}' ms4mscon1abs{3,tt,ss}' ms4mscon1abs{4,tt,ss}' ms4mscon1abs{5,tt,ss}' ms4mscon1abs{6,tt,ss}' ms4mscon1abs{7,tt,ss}' ms4mscon1abs{9,tt,ss}']);
% numperbin=floor(length(sortabs)/4);
% for nn=2:4
%   xbin(nn)=sortabs(numperbin*(nn-1)+1);
% end
% for ll=soalist
%   tac4mscon1absallbin{ll,tt,ss}=bin_my_abs(tac4mscon1abs{ll,tt,ss},xbin);
%   aud4mscon1absallbin{ll,tt,ss}=bin_my_abs(aud4mscon1abs{ll,tt,ss},xbin);
%   nul4mscon1absallbin{ll,tt,ss}=bin_my_abs(nul4mscon1abs{ll,tt,ss},xbin);
%   ms4mscon1absallbin{ll,tt,ss} =bin_my_abs(ms4mscon1abs{ll,tt,ss},xbin);
% end


% % % Do not delete, but not using anymore
%   [sortabs,sortind]=sort([tacabsIAF{ll,tt,ss}' audabsIAF{ll,tt,ss}' nulabsIAF{ll,tt,ss}' ms1absIAF{ll,tt,ss}' ms2absIAF{ll,tt,ss}']);
%   numperbin=floor(length(sortabs)/4);
%   for nn=2:4
%     xbin(nn)=sortabs(numperbin*(nn-1)+1);
%   end
%   tacabsIAFbin{ll,tt,ss}=bin_my_abs(tacabsIAF{ll,tt,ss},xbin);
%   audabsIAFbin{ll,tt,ss}=bin_my_abs(audabsIAF{ll,tt,ss},xbin);
%   nulabsIAFbin{ll,tt,ss}=bin_my_abs(nulabsIAF{ll,tt,ss},xbin);
%   ms1absIAFbin{ll,tt,ss}=bin_my_abs(ms1absIAF{ll,tt,ss},xbin);
%   ms2absIAFbin{ll,tt,ss}=bin_my_abs(ms2absIAF{ll,tt,ss},xbin);
%   
% % % Do not delete, but not using anymore
%   [sortabs,sortind]=sort([tac4mscon1absIAF{ll,tt,ss}' aud4mscon1absIAF{ll,tt,ss}' nul4mscon1absIAF{ll,tt,ss}' ms4mscon1absIAF{ll,tt,ss}']);
%   numperbin=floor(length(sortabs)/4);
%   for nn=2:4
%     xbin(nn)=sortabs(numperbin*(nn-1)+1);
%   end
%   tac4mscon1absIAFbin{ll,tt,ss}=bin_my_abs(tac4mscon1absIAF{ll,tt,ss},xbin);
%   aud4mscon1absIAFbin{ll,tt,ss}=bin_my_abs(aud4mscon1absIAF{ll,tt,ss},xbin);
%   nul4mscon1absIAFbin{ll,tt,ss}=bin_my_abs(nul4mscon1absIAF{ll,tt,ss},xbin);
%   ms4mscon1absIAFbin{ll,tt,ss}=bin_my_abs(ms4mscon1absIAF{ll,tt,ss},xbin);  
%   %               powIAFdist{ii,ll}=sortabs;
%   
% % % Do not delete, but not using anymore
%   [sortabs,sortind]=sort([tac4mscon2absIAF{ll,tt,ss}' aud4mscon2absIAF{ll,tt,ss}' nul4mscon2absIAF{ll,tt,ss}' ms4mscon2absIAF{ll,tt,ss}']);
%   numperbin=floor(length(sortabs)/4);
%   for nn=2:4
%     xbin(nn)=sortabs(numperbin*(nn-1)+1);
%   end
%   tac4mscon2absIAFbin{ll,tt,ss}=bin_my_abs(tac4mscon2absIAF{ll,tt,ss},xbin);
%   aud4mscon2absIAFbin{ll,tt,ss}=bin_my_abs(aud4mscon2absIAF{ll,tt,ss},xbin);
%   nul4mscon2absIAFbin{ll,tt,ss}=bin_my_abs(nul4mscon2absIAF{ll,tt,ss},xbin);
%   ms4mscon2absIAFbin{ll,tt,ss}=bin_my_abs(ms4mscon2absIAF{ll,tt,ss},xbin);

end % ll

return

%%   Use this code to un-struct the variables in code external to this one, if desired, after calling this one

% tacbin = angbin.tacbin;
% audbin = angbin.audbin;
% nulbin = angbin.nulbin;
% ms1bin = angbin.ms1bin;
% ms2bin = angbin.ms2bin;
% tac4mscon1bin = angbin.tac4mscon1bin;
% tac4mscon2bin = angbin.tac4mscon2bin;
% aud4mscon1bin = angbin.aud4mscon1bin;
% aud4mscon2bin = angbin.aud4mscon2bin;
% nul4mscon1bin = angbin.nul4mscon1bin;
% nul4mscon2bin = angbin.nul4mscon2bin;
% ms4mscon1bin = angbin.ms4mscon1bin;
% ms4mscon2bin = angbin.ms4mscon2bin;
% 
% tacabsbin = absbin.tacabsbin;
% audabsbin = absbin.audabsbin;
% nulabsbin = absbin.nulabsbin;
% ms1absbin = absbin.ms1absbin;
% ms2absbin = absbin.ms2absbin;
% tac4mscon1absbin = absbin.tac4mscon1absbin;
% aud4mscon1absbin = absbin.aud4mscon1absbin;
% nul4mscon1absbin = absbin.nul4mscon1absbin;
% ms4mscon1absbin = absbin.ms4mscon1absbin;
% tac4mscon2absbin = absbin.tac4mscon2absbin;
% aud4mscon2absbin = absbin.aud4mscon2absbin;
% nul4mscon2absbin = absbin.nul4mscon2absbin;
% ms4mscon2absbin = absbin.ms4mscon2absbin;


%% Which variables to use:

% All above were computed with tactile onset at time 0.

% These are with trial numbers equated and the location of 'time 0' equated across conditions.  
% Thus, e.g. tac4mscon1bin{1,tt,ss} (compare to AT500) has the phase computed prior to -500ms as it wil be compared to ms4mscon1bin{1,tt,ss} which has auditory at -500.
% But tac4mscon1bin{5,tt,ss} (compare to AT0) has the phase computed prior to 0ms.

% *con1* means the phase taken at -250ms (for 2 Hz presumably??) prior to onset of first stimulus (in MS case) or equivalent time for unisensory
% *con2* means the phase taken at -250ms (for 2 Hz presumably??) prior to onset of second stimulus

% power is taken -500ms prior to either first or second stimulus (con1 or con2)

% The variables without *con* mean that the time is relative to first stimulus onset, irrespective of later comparison (thus always -250ms prior to tactile or auditory alone)

% Thus, I am focussing only on *con1* variables, but include the others in case they are useful

% To know which trials they correspond to, please see freqtr.tr
% freqtr.tr.t10trialkept    indicates which trials to use for tac
% freqtr.tr.nllttrialkept   indicates which trials to use for nul
% freqtr.tr.all40trialkept  indicates which trials to use for aud
% freqtr.tr.tlltrialkept    indicates which trials to use for multisensory

% Be sure to use unique(freqtr.tr.t10trialkept{ll,ss,tt}), rather than just freqtr.tr.t10trialkept{ll,ss,tt},
% since in a few rare cases, there are repeated numbers.



