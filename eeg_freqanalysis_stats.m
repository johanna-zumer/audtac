eeg_legomagic_preamble

%%  Main stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load elec1010_neighb.mat

% allcond_sameN=1;
soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];
tacbasemax=min(-.15,soades-.15);
audbasemax=min(-.15,fliplr(soades)-.15);

ylimlo=[4 7; 8 12; 14 30];
ylimhi=[30 45; 55 80];

soadesc={'Aud first by 500ms' '' 'Aud first by 70ms' 'Aud first by 20ms' 'Simultaneous' 'Tac first by 20ms' 'Tac first by 70ms' '' 'Tac first by 500ms'};
statsflag=1;
itcdepsampflag=1;
plotflag=0;
printflag=0;
audtacflag=0;
comb2flag=1;
fftaddflag=0;
synchasynch=0;
mcseed=13;  % montecarlo cfg.randomseed
savegrindflag=1;   %=1 saves out 'grind' 

sleep=1;
if sleep
  trialkc=-1;
  usetr=1;
  subuseall=iiBuse;
  iter=11;  
else
  trialkc=-1;
  usetr=3; % 1 with 27, or 2 or 3 with 31 or 32
  subuseall=setdiff(iiSuse,[])
  %     iter=27;
  iter=31;
end
resetusetr=0;

soalist=[1 3 4 5 6 7 9];
% soalist=[3 4 5 6 7 9];

chanuse_sleep0={'all' '-F4'};
chanuse_sleep1={'all' '-AF7' '-AF3' '-Fp1'};

tt=3;
% for sleep=[1]
if sleep
  chanuse=chanuse_sleep1;
else
  chanuse=chanuse_sleep0;
end
%   for tt=[3]
figind=1;


for ll=soalist
  % for ll=[9];
  clearvars -except ll tt sub* *dir ii*use sleep *flag figind soadesc soalist chanuse* ylim* neigh* *basemax synch* mcseed *tr trialkc iter
  
  if resetusetr
    usetr=2;
  end
  
  submin=subuseall(1)-1;
  subuseind=0;
  subuseind_nKD=0;
  subuseind_nSD=0;
  subuseind_nKD_nSD=0;
  
  % Baseline correct each participant prior to entering to stats???? NO
  for ii=subuseall
    %   for ii=8:9
    cd([edir sub{ii} ])
    %       load(['freq_diffs_averef_' sub{ii} '.mat']);
    try
      if fftaddflag
        load(['freq_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_tt' num2str(tt) '.mat'])
      else
        if iter==27 && sleep==0
          load(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_tt' num2str(tt) '.mat'])
        elseif sleep==0
          load(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'])
        elseif sleep==1
          try
            load(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat'])
          catch
            load(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_tt' num2str(tt) '_trialkc' num2str(trialkc) '.mat'])
          end
        end
      end
      if audtacflag
        load(['freq_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(0) '_tt' num2str(tt) '.mat'])
      end
    catch
      if tt==2
        load(['freq_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '.mat'])
        if audtacflag
          load(['freq_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(0) '.mat'])
        end
      end
    end
    if audtacflag
      tka=load(['trialkeptTFR_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(0) '_iter' num2str(iter) '_usetr' num2str(usetr) '.mat']);
    end
    try
      if usetr==2
        % load usetr=0 here; then later down load usetr=2
        tkt=load(['trialkeptTFR_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter) '_usetr' num2str(3) '_trialkc' num2str(trialkc) '.mat']);
      else
        tkt=load(['trialkeptTFR_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter) '_usetr' num2str(usetr) '.mat']);
      end
    catch
      tkt=load(['trialkeptTFR_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter) '_usetr' num2str(usetr) '_trialkc' num2str(trialkc) '.mat']);
    end
    
    if 1
      % THis is preliminary...how best to included all stages later on?
      if sleep==0
        ssuse=10; % awake
        sleepcond='Awake W';
      elseif sleep==1
        ssuse=12; % N2
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
    numtrt(ll,tt,ss,ii-submin)=tkt.numcondtfinal(ll,tt,ss);
    if audtacflag
      numtra(ll,tt,ss,ii-submin)=tka.numcondafinal(ll,tt,ss);
    end
    
    
    if min(numtrt(ll,tt,ss,ii-submin),[],1)<20 % what is best number to use here?
      subuse=setdiff(subuse,ii);
    else
      subuseind=subuseind+1;
      freqloall_tacPaud_comb1{subuseind,1}=freqlo_tacPaud_comb{ll,tt,ss,1};
      freqloall_tacMSpN_comb1{subuseind,1}=freqlo_tacMSpN_comb{ll,tt,ss,1};
      freqhiall_tacPaud_comb1{subuseind,1}=freqhi_tacPaud_comb{ll,tt,ss,1};
      freqhiall_tacMSpN_comb1{subuseind,1}=freqhi_tacMSpN_comb{ll,tt,ss,1};
      freqhiall_tNulAlone_comb{subuseind,1}=freqhi_tNulAlone_comb{ll,tt,ss};
      freqloall_tNulAlone_comb{subuseind,1}=freqlo_tNulAlone_comb{ll,tt,ss};
      freqhiall_tTacAlone_comb{subuseind,1}=freqhi_tTacAlone_comb{ll,tt,ss};
      freqloall_tTacAlone_comb{subuseind,1}=freqlo_tTacAlone_comb{ll,tt,ss};
      freqhiall_tAudAlone_comb{subuseind,1}=freqhi_tAudAlone_comb{ll,tt,ss};
      freqloall_tAudAlone_comb{subuseind,1}=freqlo_tAudAlone_comb{ll,tt,ss};
      freqhiall_tMSAlone_comb{subuseind,1}=freqhi_tMSAlone_comb{ll,tt,ss};
      freqloall_tMSAlone_comb{subuseind,1}=freqlo_tMSAlone_comb{ll,tt,ss};
      freqloall_tacPaud_comb1{subuseind,1}.dimord='chan_freq_time';
      freqloall_tacMSpN_comb1{subuseind,1}.dimord='chan_freq_time';
      freqhiall_tacPaud_comb1{subuseind,1}.dimord='chan_freq_time';
      freqhiall_tacMSpN_comb1{subuseind,1}.dimord='chan_freq_time';
      freqhiall_tNulAlone_comb{subuseind,1}.dimord='chan_freq_time';
      freqloall_tNulAlone_comb{subuseind,1}.dimord='chan_freq_time';
      freqhiall_tTacAlone_comb{subuseind,1}.dimord='chan_freq_time';
      freqloall_tTacAlone_comb{subuseind,1}.dimord='chan_freq_time';
      freqhiall_tAudAlone_comb{subuseind,1}.dimord='chan_freq_time';
      freqloall_tAudAlone_comb{subuseind,1}.dimord='chan_freq_time';
      freqhiall_tMSAlone_comb{subuseind,1}.dimord='chan_freq_time';
      freqloall_tMSAlone_comb{subuseind,1}.dimord='chan_freq_time';
      if sleep
        if size(freqlo_tacPaud_nKD{ll,tt,ss}.cumtapcnt,1)>20
            subuseind_nKD=subuseind_nKD+1;
            freqloall_tacPaud_nKD{subuseind_nKD}=freqlo_tacPaud_nKD{ll,tt,ss};
            freqloall_tacMSpN_nKD{subuseind_nKD}=freqlo_tacMSpN_nKD{ll,tt,ss};
            freqhiall_tacPaud_nKD{subuseind_nKD}=freqhi_tacPaud_nKD{ll,tt,ss};
            freqhiall_tacMSpN_nKD{subuseind_nKD}=freqhi_tacMSpN_nKD{ll,tt,ss};
            freqloall_tacPaud_nKD{subuseind_nKD}.dimord='chan_freq_time';
            freqloall_tacMSpN_nKD{subuseind_nKD}.dimord='chan_freq_time';
            freqhiall_tacPaud_nKD{subuseind_nKD}.dimord='chan_freq_time';
            freqhiall_tacMSpN_nKD{subuseind_nKD}.dimord='chan_freq_time';
            freqhiall_tNulAlone_nKD{subuseind_nKD}=freqhi_tNulAlone_nKD{ll,tt,ss};
            freqloall_tNulAlone_nKD{subuseind_nKD}=freqlo_tNulAlone_nKD{ll,tt,ss};
            freqhiall_tTacAlone_nKD{subuseind_nKD}=freqhi_tTacAlone_nKD{ll,tt,ss};
            freqloall_tTacAlone_nKD{subuseind_nKD}=freqlo_tTacAlone_nKD{ll,tt,ss};
            freqhiall_tAudAlone_nKD{subuseind_nKD}=freqhi_tAudAlone_nKD{ll,tt,ss};
            freqloall_tAudAlone_nKD{subuseind_nKD}=freqlo_tAudAlone_nKD{ll,tt,ss};
            freqhiall_tMSAlone_nKD{subuseind_nKD}=freqhi_tMSAlone_nKD{ll,tt,ss};
            freqloall_tMSAlone_nKD{subuseind_nKD}=freqlo_tMSAlone_nKD{ll,tt,ss};
            freqhiall_tNulAlone_nKD{subuseind_nKD}.dimord='chan_freq_time';
            freqloall_tNulAlone_nKD{subuseind_nKD}.dimord='chan_freq_time';
            freqhiall_tTacAlone_nKD{subuseind_nKD}.dimord='chan_freq_time';
            freqloall_tTacAlone_nKD{subuseind_nKD}.dimord='chan_freq_time';
            freqhiall_tAudAlone_nKD{subuseind_nKD}.dimord='chan_freq_time';
            freqloall_tAudAlone_nKD{subuseind_nKD}.dimord='chan_freq_time';
            freqhiall_tMSAlone_nKD{subuseind_nKD}.dimord='chan_freq_time';
            freqloall_tMSAlone_nKD{subuseind_nKD}.dimord='chan_freq_time';
        end
        if size(freqlo_tacPaud_nSD{ll,tt,ss}.cumtapcnt,1)>20
            subuseind_nSD=subuseind_nSD+1;
            freqloall_tacPaud_nSD{subuseind_nSD}=freqlo_tacPaud_nSD{ll,tt,ss};
            freqloall_tacMSpN_nSD{subuseind_nSD}=freqlo_tacMSpN_nSD{ll,tt,ss};
            freqhiall_tacPaud_nSD{subuseind_nSD}=freqhi_tacPaud_nSD{ll,tt,ss};
            freqhiall_tacMSpN_nSD{subuseind_nSD}=freqhi_tacMSpN_nSD{ll,tt,ss};
            freqloall_tacPaud_nSD{subuseind_nSD}.dimord='chan_freq_time';
            freqloall_tacMSpN_nSD{subuseind_nSD}.dimord='chan_freq_time';
            freqhiall_tacPaud_nSD{subuseind_nSD}.dimord='chan_freq_time';
            freqhiall_tacMSpN_nSD{subuseind_nSD}.dimord='chan_freq_time';
            freqhiall_tNulAlone_nSD{subuseind_nSD}=freqhi_tNulAlone_nSD{ll,tt,ss};
            freqloall_tNulAlone_nSD{subuseind_nSD}=freqlo_tNulAlone_nSD{ll,tt,ss};
            freqhiall_tTacAlone_nSD{subuseind_nSD}=freqhi_tTacAlone_nSD{ll,tt,ss};
            freqloall_tTacAlone_nSD{subuseind_nSD}=freqlo_tTacAlone_nSD{ll,tt,ss};
            freqhiall_tAudAlone_nSD{subuseind_nSD}=freqhi_tAudAlone_nSD{ll,tt,ss};
            freqloall_tAudAlone_nSD{subuseind_nSD}=freqlo_tAudAlone_nSD{ll,tt,ss};
            freqhiall_tMSAlone_nSD{subuseind_nSD}=freqhi_tMSAlone_nSD{ll,tt,ss};
            freqloall_tMSAlone_nSD{subuseind_nSD}=freqlo_tMSAlone_nSD{ll,tt,ss};
            freqhiall_tNulAlone_nSD{subuseind_nSD}.dimord='chan_freq_time';
            freqloall_tNulAlone_nSD{subuseind_nSD}.dimord='chan_freq_time';
            freqhiall_tTacAlone_nSD{subuseind_nSD}.dimord='chan_freq_time';
            freqloall_tTacAlone_nSD{subuseind_nSD}.dimord='chan_freq_time';
            freqhiall_tAudAlone_nSD{subuseind_nSD}.dimord='chan_freq_time';
            freqloall_tAudAlone_nSD{subuseind_nSD}.dimord='chan_freq_time';
            freqhiall_tMSAlone_nSD{subuseind_nSD}.dimord='chan_freq_time';
            freqloall_tMSAlone_nSD{subuseind_nSD}.dimord='chan_freq_time';
        end
        if size(freqlo_tacPaud_nSD{ll,tt,ss}.cumtapcnt,1)>20
            subuseind_nKD_nSD=subuseind_nKD_nSD+1;
            freqloall_tacPaud_nKD_nSD{subuseind_nKD_nSD}=freqlo_tacPaud_nKD_nSD{ll,tt,ss};
            freqloall_tacMSpN_nKD_nSD{subuseind_nKD_nSD}=freqlo_tacMSpN_nKD_nSD{ll,tt,ss};
            freqhiall_tacPaud_nKD_nSD{subuseind_nKD_nSD}=freqhi_tacPaud_nKD_nSD{ll,tt,ss};
            freqhiall_tacMSpN_nKD_nSD{subuseind_nKD_nSD}=freqhi_tacMSpN_nKD_nSD{ll,tt,ss};
            freqloall_tacPaud_nKD_nSD{subuseind_nKD_nSD}.dimord='chan_freq_time';
            freqloall_tacMSpN_nKD_nSD{subuseind_nKD_nSD}.dimord='chan_freq_time';
            freqhiall_tacPaud_nKD_nSD{subuseind_nKD_nSD}.dimord='chan_freq_time';
            freqhiall_tacMSpN_nKD_nSD{subuseind_nKD_nSD}.dimord='chan_freq_time';
            freqhiall_tNulAlone_nKD_nSD{subuseind_nKD_nSD}=freqhi_tNulAlone_nKD_nSD{ll,tt,ss};
            freqloall_tNulAlone_nKD_nSD{subuseind_nKD_nSD}=freqlo_tNulAlone_nKD_nSD{ll,tt,ss};
            freqhiall_tTacAlone_nKD_nSD{subuseind_nKD_nSD}=freqhi_tTacAlone_nKD_nSD{ll,tt,ss};
            freqloall_tTacAlone_nKD_nSD{subuseind_nKD_nSD}=freqlo_tTacAlone_nKD_nSD{ll,tt,ss};
            freqhiall_tAudAlone_nKD_nSD{subuseind_nKD_nSD}=freqhi_tAudAlone_nKD_nSD{ll,tt,ss};
            freqloall_tAudAlone_nKD_nSD{subuseind_nKD_nSD}=freqlo_tAudAlone_nKD_nSD{ll,tt,ss};
            freqhiall_tMSAlone_nKD_nSD{subuseind_nKD_nSD}=freqhi_tMSAlone_nKD_nSD{ll,tt,ss};
            freqloall_tMSAlone_nKD_nSD{subuseind_nKD_nSD}=freqlo_tMSAlone_nKD_nSD{ll,tt,ss};
            freqhiall_tNulAlone_nKD_nSD{subuseind_nKD_nSD}.dimord='chan_freq_time';
            freqloall_tNulAlone_nKD_nSD{subuseind_nKD_nSD}.dimord='chan_freq_time';
            freqhiall_tTacAlone_nKD_nSD{subuseind_nKD_nSD}.dimord='chan_freq_time';
            freqloall_tTacAlone_nKD_nSD{subuseind_nKD_nSD}.dimord='chan_freq_time';
            freqhiall_tAudAlone_nKD_nSD{subuseind_nKD_nSD}.dimord='chan_freq_time';
            freqloall_tAudAlone_nKD_nSD{subuseind_nKD_nSD}.dimord='chan_freq_time';
            freqhiall_tMSAlone_nKD_nSD{subuseind_nKD_nSD}.dimord='chan_freq_time';
            freqloall_tMSAlone_nKD_nSD{subuseind_nKD_nSD}.dimord='chan_freq_time';
        end
      end
      if synchasynch && ll<5
        freqloall_tMSasynch_comb1{subuseind}=freqlo_tMSasynch_comb{ll,tt,ss};
        freqloall_tMSsynch_comb1{subuseind}=freqlo_tMSsynch_comb{ll,tt,ss};
        freqhiall_tMSasynch_comb1{subuseind}=freqhi_tMSasynch_comb{ll,tt,ss};
        freqhiall_tMSsynch_comb1{subuseind}=freqhi_tMSsynch_comb{ll,tt,ss};
        freqloall_tMSasynch_comb1{subuseind}.dimord='chan_freq_time';
        freqloall_tMSsynch_comb1{subuseind}.dimord='chan_freq_time';
        freqhiall_tMSasynch_comb1{subuseind}.dimord='chan_freq_time';
        freqhiall_tMSsynch_comb1{subuseind}.dimord='chan_freq_time';
      end
      if comb2flag
        freqloall_tacPaud_comb2{subuseind,1}=freqlo_tacPaud_comb{ll,tt,ss,2};
        freqloall_tacMSpN_comb2{subuseind,1}=freqlo_tacMSpN_comb{ll,tt,ss,2};
        freqhiall_tacPaud_comb2{subuseind,1}=freqhi_tacPaud_comb{ll,tt,ss,2};
        freqhiall_tacMSpN_comb2{subuseind,1}=freqhi_tacMSpN_comb{ll,tt,ss,2};
        freqloall_tacPaud_comb2{subuseind,1}.dimord='chan_freq_time';
        freqloall_tacMSpN_comb2{subuseind,1}.dimord='chan_freq_time';
        freqhiall_tacPaud_comb2{subuseind,1}.dimord='chan_freq_time';
        freqhiall_tacMSpN_comb2{subuseind,1}.dimord='chan_freq_time';
        if synchasynch ll<5
          freqloall_tMSasynch_comb2{subuseind}=freqlo_tMSasynch_comb{ll,tt,ss};
          freqloall_tMSsynch_comb2{subuseind}=freqlo_tMSsynch_comb{ll,tt,ss};
          freqhiall_tMSasynch_comb2{subuseind}=freqhi_tMSasynch_comb{ll,tt,ss};
          freqhiall_tMSsynch_comb2{subuseind}=freqhi_tMSsynch_comb{ll,tt,ss};
          freqloall_tMSasynch_comb2{subuseind}.dimord='chan_freq_time';
          freqloall_tMSsynch_comb2{subuseind}.dimord='chan_freq_time';
          freqhiall_tMSasynch_comb2{subuseind}.dimord='chan_freq_time';
          freqhiall_tMSsynch_comb2{subuseind}.dimord='chan_freq_time';
        end
      end
      if fftaddflag
        freqloall_tacPaud_fftadd{subuseind}=freqlo_tacPaud_fftadd{ll,tt,ss};
        freqloall_tacMSpN_fftadd{subuseind}=freqlo_tacMSpN_fftadd{ll,tt,ss};
        freqhiall_tacPaud_fftadd{subuseind}=freqhi_tacPaud_fftadd{ll,tt,ss};
        freqhiall_tacMSpN_fftadd{subuseind}=freqhi_tacMSpN_fftadd{ll,tt,ss};
        freqhiall_tNulAlone_fftadd{subuseind}=freqhi_tNulAlone_fftadd{ll,tt,ss};
        freqloall_tNulAlone_fftadd{subuseind}=freqlo_tNulAlone_fftadd{ll,tt,ss};
        freqhiall_tTacAlone_fftadd{subuseind}=freqhi_tTacAlone_fftadd{ll,tt,ss};
        freqloall_tTacAlone_fftadd{subuseind}=freqlo_tTacAlone_fftadd{ll,tt,ss};
        freqhiall_tAudAlone_fftadd{subuseind}=freqhi_tAudAlone_fftadd{ll,tt,ss};
        freqloall_tAudAlone_fftadd{subuseind}=freqlo_tAudAlone_fftadd{ll,tt,ss};
        freqhiall_tMSAlone_fftadd{subuseind}=freqhi_tMSAlone_fftadd{ll,tt,ss};
        freqloall_tMSAlone_fftadd{subuseind}=freqlo_tMSAlone_fftadd{ll,tt,ss};
        freqloall_tacPaud_fftadd{subuseind}.dimord='chan_freq_time';
        freqloall_tacMSpN_fftadd{subuseind}.dimord='chan_freq_time';
        freqhiall_tacPaud_fftadd{subuseind}.dimord='chan_freq_time';
        freqhiall_tacMSpN_fftadd{subuseind}.dimord='chan_freq_time';
        freqhiall_tNulAlone_fftadd{subuseind}.dimord='chan_freq_time';
        freqloall_tNulAlone_fftadd{subuseind}.dimord='chan_freq_time';
        freqhiall_tTacAlone_fftadd{subuseind}.dimord='chan_freq_time';
        freqloall_tTacAlone_fftadd{subuseind}.dimord='chan_freq_time';
        freqhiall_tAudAlone_fftadd{subuseind}.dimord='chan_freq_time';
        freqloall_tAudAlone_fftadd{subuseind}.dimord='chan_freq_time';
        freqhiall_tMSAlone_fftadd{subuseind}.dimord='chan_freq_time';
        freqloall_tMSAlone_fftadd{subuseind}.dimord='chan_freq_time';
        freqloall_tacPaud_fftadd{subuseind}.crsspctrm=[];
        freqloall_tacMSpN_fftadd{subuseind}.crsspctrm=[];
        freqhiall_tacPaud_fftadd{subuseind}.crsspctrm=[];
        freqhiall_tacMSpN_fftadd{subuseind}.crsspctrm=[];
        freqhiall_tNulAlone_fftadd{subuseind}.crsspctrm=[];
        freqloall_tNulAlone_fftadd{subuseind}.crsspctrm=[];
        freqhiall_tTacAlone_fftadd{subuseind}.crsspctrm=[];
        freqloall_tTacAlone_fftadd{subuseind}.crsspctrm=[];
        freqhiall_tAudAlone_fftadd{subuseind}.crsspctrm=[];
        freqloall_tAudAlone_fftadd{subuseind}.crsspctrm=[];
        freqhiall_tMSAlone_fftadd{subuseind}.crsspctrm=[];
        freqloall_tMSAlone_fftadd{subuseind}.crsspctrm=[];
      end
      
      if audtacflag
        freqloall_audPtac_comb1{subuseind}=freqlo_audPtac_comb{ll,tt,ss,1};
        freqloall_audMSpN_comb1{subuseind}=freqlo_audMSpN_comb{ll,tt,ss,1};
        freqhiall_audPtac_comb1{subuseind}=freqhi_audPtac_comb{ll,tt,ss,1};
        freqhiall_audMSpN_comb1{subuseind}=freqhi_audMSpN_comb{ll,tt,ss,1};
        freqhiall_aNulAlone_comb{subuseind}=freqhi_aNulAlone_comb{ll,tt,ss};
        freqloall_aNulAlone_comb{subuseind}=freqlo_aNulAlone_comb{ll,tt,ss};
        freqhiall_aTacAlone_comb{subuseind}=freqhi_aTacAlone_comb{ll,tt,ss};
        freqloall_aTacAlone_comb{subuseind}=freqlo_aTacAlone_comb{ll,tt,ss};
        freqhiall_aAudAlone_comb{subuseind}=freqhi_aAudAlone_comb{ll,tt,ss};
        freqloall_aAudAlone_comb{subuseind}=freqlo_aAudAlone_comb{ll,tt,ss};
        freqhiall_aMSAlone_comb{subuseind}=freqhi_aMSAlone_comb{ll,tt,ss};
        freqloall_aMSAlone_comb{subuseind}=freqlo_aMSAlone_comb{ll,tt,ss};
        freqloall_audPtac_comb1{subuseind}.dimord='chan_freq_time';
        freqloall_audMSpN_comb1{subuseind}.dimord='chan_freq_time';
        freqhiall_audPtac_comb1{subuseind}.dimord='chan_freq_time';
        freqhiall_audMSpN_comb1{subuseind}.dimord='chan_freq_time';
        if comb2flag
          freqloall_audPtac_comb2{subuseind}=freqlo_audPtac_comb{ll,tt,ss,2};
          freqloall_audMSpN_comb2{subuseind}=freqlo_audMSpN_comb{ll,tt,ss,2};
          freqhiall_audPtac_comb2{subuseind}=freqhi_audPtac_comb{ll,tt,ss,2};
          freqhiall_audMSpN_comb2{subuseind}=freqhi_audMSpN_comb{ll,tt,ss,2};
          freqloall_audPtac_comb2{subuseind}.dimord='chan_freq_time';
          freqloall_audMSpN_comb2{subuseind}.dimord='chan_freq_time';
          freqhiall_audPtac_comb2{subuseind}.dimord='chan_freq_time';
          freqhiall_audMSpN_comb2{subuseind}.dimord='chan_freq_time';
        end
        freqhiall_aNulAlone_comb{subuseind}.dimord='chan_freq_time';
        freqloall_aNulAlone_comb{subuseind}.dimord='chan_freq_time';
        freqhiall_aTacAlone_comb{subuseind}.dimord='chan_freq_time';
        freqloall_aTacAlone_comb{subuseind}.dimord='chan_freq_time';
        freqhiall_aAudAlone_comb{subuseind}.dimord='chan_freq_time';
        freqloall_aAudAlone_comb{subuseind}.dimord='chan_freq_time';
        freqhiall_aMSAlone_comb{subuseind}.dimord='chan_freq_time';
        freqloall_aMSAlone_comb{subuseind}.dimord='chan_freq_time';
        try
          freqloall_audPtac_fftadd{subuseind}=freqlo_audPtac_fftadd{ll,tt,ss};
          freqloall_audMSpN_fftadd{subuseind}=freqlo_audMSpN_fftadd{ll,tt,ss};
          freqhiall_audPtac_fftadd{subuseind}=freqhi_audPtac_fftadd{ll,tt,ss};
          freqhiall_audMSpN_fftadd{subuseind}=freqhi_audMSpN_fftadd{ll,tt,ss};
          freqhiall_aNulAlone_fftadd{subuseind}=freqhi_aNulAlone_fftadd{ll,tt,ss};
          freqloall_aNulAlone_fftadd{subuseind}=freqlo_aNulAlone_fftadd{ll,tt,ss};
          freqhiall_aTacAlone_fftadd{subuseind}=freqhi_aTacAlone_fftadd{ll,tt,ss};
          freqloall_aTacAlone_fftadd{subuseind}=freqlo_aTacAlone_fftadd{ll,tt,ss};
          freqhiall_aAudAlone_fftadd{subuseind}=freqhi_aAudAlone_fftadd{ll,tt,ss};
          freqloall_aAudAlone_fftadd{subuseind}=freqlo_aAudAlone_fftadd{ll,tt,ss};
          freqhiall_aMSAlone_fftadd{subuseind}=freqhi_aMSAlone_fftadd{ll,tt,ss};
          freqloall_aMSAlone_fftadd{subuseind}=freqlo_aMSAlone_fftadd{ll,tt,ss};
          freqloall_audPtac_fftadd{subuseind}.dimord='chan_freq_time';
          freqloall_audMSpN_fftadd{subuseind}.dimord='chan_freq_time';
          freqhiall_audPtac_fftadd{subuseind}.dimord='chan_freq_time';
          freqhiall_audMSpN_fftadd{subuseind}.dimord='chan_freq_time';
          freqhiall_aNulAlone_fftadd{subuseind}.dimord='chan_freq_time';
          freqloall_aNulAlone_fftadd{subuseind}.dimord='chan_freq_time';
          freqhiall_aTacAlone_fftadd{subuseind}.dimord='chan_freq_time';
          freqloall_aTacAlone_fftadd{subuseind}.dimord='chan_freq_time';
          freqhiall_aAudAlone_fftadd{subuseind}.dimord='chan_freq_time';
          freqloall_aAudAlone_fftadd{subuseind}.dimord='chan_freq_time';
          freqhiall_aMSAlone_fftadd{subuseind}.dimord='chan_freq_time';
          freqloall_aMSAlone_fftadd{subuseind}.dimord='chan_freq_time';
          freqloall_audPtac_fftadd{subuseind}.crsspctrm=[];
          freqloall_audMSpN_fftadd{subuseind}.crsspctrm=[];
          freqhiall_audPtac_fftadd{subuseind}.crsspctrm=[];
          freqhiall_audMSpN_fftadd{subuseind}.crsspctrm=[];
          freqhiall_aNulAlone_fftadd{subuseind}.crsspctrm=[];
          freqloall_aNulAlone_fftadd{subuseind}.crsspctrm=[];
          freqhiall_aTacAlone_fftadd{subuseind}.crsspctrm=[];
          freqloall_aTacAlone_fftadd{subuseind}.crsspctrm=[];
          freqhiall_aAudAlone_fftadd{subuseind}.crsspctrm=[];
          freqloall_aAudAlone_fftadd{subuseind}.crsspctrm=[];
          freqhiall_aMSAlone_fftadd{subuseind}.crsspctrm=[];
          freqloall_aMSAlone_fftadd{subuseind}.crsspctrm=[];
        end
      end
    end
    %         end % ss
    clear freqlo_* freqhi_*
    
    if usetr==2
      try
        load(['freqcomb_diffs_averef_' sub{ii} '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter+1) '_trialkc' num2str(trialkc) '.mat'])
        %       tkt=load(['trialkeptTFR_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter+1) '_usetr' num2str(2) '_trialkc' num2str(trialkc) '.mat']);
        
        freqloall_tacPaud_comb1{subuseind,2}=freqlo_tacPaud_comb{ll,tt,ss,1};
        freqloall_tacMSpN_comb1{subuseind,2}=freqlo_tacMSpN_comb{ll,tt,ss,1};
        freqloall_tNulAlone_comb{subuseind,2}=freqlo_tNulAlone_comb{ll,tt,ss};
        freqloall_tTacAlone_comb{subuseind,2}=freqlo_tTacAlone_comb{ll,tt,ss};
        freqloall_tAudAlone_comb{subuseind,2}=freqlo_tAudAlone_comb{ll,tt,ss};
        freqloall_tMSAlone_comb{subuseind,2}=freqlo_tMSAlone_comb{ll,tt,ss};
        freqloall_tacPaud_comb1{subuseind,2}.dimord='chan_freq_time';
        freqloall_tacMSpN_comb1{subuseind,2}.dimord='chan_freq_time';
        freqloall_tNulAlone_comb{subuseind,2}.dimord='chan_freq_time';
        freqloall_tTacAlone_comb{subuseind,2}.dimord='chan_freq_time';
        freqloall_tAudAlone_comb{subuseind,2}.dimord='chan_freq_time';
        freqloall_tMSAlone_comb{subuseind,2}.dimord='chan_freq_time';
        if comb2flag
          freqloall_tacPaud_comb2{subuseind,2}=freqlo_tacPaud_comb{ll,tt,ss,2};
          freqloall_tacMSpN_comb2{subuseind,2}=freqlo_tacMSpN_comb{ll,tt,ss,2};
          freqloall_tacPaud_comb2{subuseind,2}.dimord='chan_freq_time';
          freqloall_tacMSpN_comb2{subuseind,2}.dimord='chan_freq_time';
        end
        freqhiall_tacPaud_comb1{subuseind,2}=freqhi_tacPaud_comb{ll,tt,ss,1};
        freqhiall_tacMSpN_comb1{subuseind,2}=freqhi_tacMSpN_comb{ll,tt,ss,1};
        freqhiall_tNulAlone_comb{subuseind,2}=freqhi_tNulAlone_comb{ll,tt,ss};
        freqhiall_tTacAlone_comb{subuseind,2}=freqhi_tTacAlone_comb{ll,tt,ss};
        freqhiall_tAudAlone_comb{subuseind,2}=freqhi_tAudAlone_comb{ll,tt,ss};
        freqhiall_tMSAlone_comb{subuseind,2}=freqhi_tMSAlone_comb{ll,tt,ss};
        freqhiall_tacPaud_comb1{subuseind,2}.dimord='chan_freq_time';
        freqhiall_tacMSpN_comb1{subuseind,2}.dimord='chan_freq_time';
        freqhiall_tNulAlone_comb{subuseind,2}.dimord='chan_freq_time';
        freqhiall_tTacAlone_comb{subuseind,2}.dimord='chan_freq_time';
        freqhiall_tAudAlone_comb{subuseind,2}.dimord='chan_freq_time';
        freqhiall_tMSAlone_comb{subuseind,2}.dimord='chan_freq_time';
        if comb2flag
          freqhiall_tacPaud_comb2{subuseind,2}=freqhi_tacPaud_comb{ll,tt,ss,2};
          freqhiall_tacMSpN_comb2{subuseind,2}=freqhi_tacMSpN_comb{ll,tt,ss,2};
          freqhiall_tacPaud_comb2{subuseind,2}.dimord='chan_freq_time';
          freqhiall_tacMSpN_comb2{subuseind,2}.dimord='chan_freq_time';
        end
      catch
        if ~exist(['trialkeptTFR_tt' num2str(tt) '_sleep' num2str(sleep) '_tacaud' num2str(1) '_iter' num2str(iter) '_usetr' num2str(usetr) '_trialkc' num2str(trialkc) '.mat'],'file')
          error('something not right with this iter/usetr combination')
        end
      end
    end
  end % ii
  
  
  resetusetr=0;
  if usetr~=2
    iterinduse=1;
  elseif usetr==2 && size(freqloall_tacPaud_comb1,2)==1
    iterinduse=1;
    usetr=3;
    resetusetr=1;
  else
    iterinduse=2;
  end
  for ii=1:subuseind
    cfg=[];
    cfg.operation='subtract';
    cfg.parameter='powspctrm';
    if fftaddflag
      freqloall_TPA_MSPN_fftadd{ii}=ft_math(cfg,freqloall_tacPaud_fftadd{ii},freqloall_tacMSpN_fftadd{ii});
      freqhiall_TPA_MSPN_fftadd{ii}=ft_math(cfg,freqhiall_tacPaud_fftadd{ii},freqhiall_tacMSpN_fftadd{ii});
      if 0
        freqloall_TPA_MSPN_fftaddbn{ii}=ft_math(cfg,freqloall_tacPaud_fftaddbn{ii},freqloall_tacMSpN_fftaddbn{ii});
        freqhiall_TPA_MSPN_fftaddbn{ii}=ft_math(cfg,freqhiall_tacPaud_fftaddbn{ii},freqhiall_tacMSpN_fftaddbn{ii});
        freqloall_TPA_MSPN_comb1bn{ii}=ft_math(cfg,freqloall_tacPaud_comb1bn{ii},freqloall_tacMSpN_comb1bn{ii});
        freqhiall_TPA_MSPN_comb1bn{ii}=ft_math(cfg,freqhiall_tacPaud_comb1bn{ii},freqhiall_tacMSpN_comb1bn{ii});
      end
      
      if audtacflag
        freqloall_APT_MSPN_fftadd{ii}=ft_math(cfg,freqloall_audPtac_fftadd{ii},freqloall_audMSpN_fftadd{ii});
        freqhiall_APT_MSPN_fftadd{ii}=ft_math(cfg,freqhiall_audPtac_fftadd{ii},freqhiall_audMSpN_fftadd{ii});
        if 0
          freqloall_APT_MSPN_fftaddbn{ii}=ft_math(cfg,freqloall_audPtac_fftaddbn{ii},freqloall_audMSpN_fftaddbn{ii});
          freqhiall_APT_MSPN_fftaddbn{ii}=ft_math(cfg,freqhiall_audPtac_fftaddbn{ii},freqhiall_audMSpN_fftaddbn{ii});
          freqloall_APT_MSPN_comb1bn{ii}=ft_math(cfg,freqloall_audPtac_comb1bn{ii},freqloall_audMSpN_comb1bn{ii});
          freqhiall_APT_MSPN_comb1bn{ii}=ft_math(cfg,freqhiall_audPtac_comb1bn{ii},freqhiall_audMSpN_comb1bn{ii});
        end
      end
    end
    cfg.parameter={'powspctrm' 'plvspctrm' 'plvabs'};
    for iterind=1:iterinduse
      freqloall_tacPaud_comb1{ii,iterind}.plvabs=abs(freqloall_tacPaud_comb1{ii,iterind}.plvspctrm);
      freqloall_tacMSpN_comb1{ii,iterind}.plvabs=abs(freqloall_tacMSpN_comb1{ii,iterind}.plvspctrm);
      freqhiall_tacPaud_comb1{ii,iterind}.plvabs=abs(freqhiall_tacPaud_comb1{ii,iterind}.plvspctrm);
      freqhiall_tacMSpN_comb1{ii,iterind}.plvabs=abs(freqhiall_tacMSpN_comb1{ii,iterind}.plvspctrm);
      freqloall_TPA_MSPN_comb1{ii,iterind}=ft_math(cfg,freqloall_tacPaud_comb1{ii,iterind},freqloall_tacMSpN_comb1{ii,iterind});
      freqhiall_TPA_MSPN_comb1{ii,iterind}=ft_math(cfg,freqhiall_tacPaud_comb1{ii,iterind},freqhiall_tacMSpN_comb1{ii,iterind});
    end
    if sleep && ss==12
      if ii<=subuseind_nKD
        freqloall_tacPaud_nKD{ii}.plvabs=abs(freqloall_tacPaud_nKD{ii}.plvspctrm);
        freqloall_tacMSpN_nKD{ii}.plvabs=abs(freqloall_tacMSpN_nKD{ii}.plvspctrm);
        freqhiall_tacPaud_nKD{ii}.plvabs=abs(freqhiall_tacPaud_nKD{ii}.plvspctrm);
        freqhiall_tacMSpN_nKD{ii}.plvabs=abs(freqhiall_tacMSpN_nKD{ii}.plvspctrm);
        freqloall_TPA_MSPN_nKD{ii}=ft_math(cfg,freqloall_tacPaud_nKD{ii},freqloall_tacMSpN_nKD{ii});
        freqhiall_TPA_MSPN_nKD{ii}=ft_math(cfg,freqhiall_tacPaud_nKD{ii},freqhiall_tacMSpN_nKD{ii});
      end
      if ii<=subuseind_nSD
        freqloall_tacPaud_nSD{ii}.plvabs=abs(freqloall_tacPaud_nSD{ii}.plvspctrm);
        freqloall_tacMSpN_nSD{ii}.plvabs=abs(freqloall_tacMSpN_nSD{ii}.plvspctrm);
        freqhiall_tacPaud_nSD{ii}.plvabs=abs(freqhiall_tacPaud_nSD{ii}.plvspctrm);
        freqhiall_tacMSpN_nSD{ii}.plvabs=abs(freqhiall_tacMSpN_nSD{ii}.plvspctrm);
        freqloall_TPA_MSPN_nSD{ii}=ft_math(cfg,freqloall_tacPaud_nSD{ii},freqloall_tacMSpN_nSD{ii});
        freqhiall_TPA_MSPN_nSD{ii}=ft_math(cfg,freqhiall_tacPaud_nSD{ii},freqhiall_tacMSpN_nSD{ii});
      end
      if ii<=subuseind_nKD_nSD
        freqloall_tacPaud_nKD_nSD{ii}.plvabs=abs(freqloall_tacPaud_nKD_nSD{ii}.plvspctrm);
        freqloall_tacMSpN_nKD_nSD{ii}.plvabs=abs(freqloall_tacMSpN_nKD_nSD{ii}.plvspctrm);
        freqhiall_tacPaud_nKD_nSD{ii}.plvabs=abs(freqhiall_tacPaud_nKD_nSD{ii}.plvspctrm);
        freqhiall_tacMSpN_nKD_nSD{ii}.plvabs=abs(freqhiall_tacMSpN_nKD_nSD{ii}.plvspctrm);
        freqloall_TPA_MSPN_nKD_nSD{ii}=ft_math(cfg,freqloall_tacPaud_nKD_nSD{ii},freqloall_tacMSpN_nKD_nSD{ii});
        freqhiall_TPA_MSPN_nKD_nSD{ii}=ft_math(cfg,freqhiall_tacPaud_nKD_nSD{ii},freqhiall_tacMSpN_nKD_nSD{ii});
      end
    end
    if synchasynch && ll<5
      freqloall_TMSs_TMSa_comb1{ii}=ft_math(cfg,freqloall_tMSsynch_comb1{ii},freqloall_tMSasynch_comb1{ii});
      freqhiall_TMSs_TMSa_comb1{ii}=ft_math(cfg,freqhiall_tMSsynch_comb1{ii},freqhiall_tMSasynch_comb1{ii});
    end
    if comb2flag
      for iterind=1:iterinduse
        freqloall_tacPaud_comb2{ii,iterind}.plvabs=abs(freqloall_tacPaud_comb2{ii,iterind}.plvspctrm);
        freqloall_tacMSpN_comb2{ii,iterind}.plvabs=abs(freqloall_tacMSpN_comb2{ii,iterind}.plvspctrm);
        freqhiall_tacPaud_comb2{ii,iterind}.plvabs=abs(freqhiall_tacPaud_comb2{ii,iterind}.plvspctrm);
        freqhiall_tacMSpN_comb2{ii,iterind}.plvabs=abs(freqhiall_tacMSpN_comb2{ii,iterind}.plvspctrm);
        freqloall_TPA_MSPN_comb2{ii,iterind}=ft_math(cfg,freqloall_tacPaud_comb2{ii,iterind},freqloall_tacMSpN_comb2{ii,iterind});
        freqhiall_TPA_MSPN_comb2{ii,iterind}=ft_math(cfg,freqhiall_tacPaud_comb2{ii,iterind},freqhiall_tacMSpN_comb2{ii,iterind});
      end
      if synchasynch && ll<5
        freqloall_TMSs_TMSa_comb2{ii}=ft_math(cfg,freqloall_tMSsynch_comb2{ii},freqloall_tMSasynch_comb2{ii});
        freqhiall_TMSs_TMSa_comb2{ii}=ft_math(cfg,freqhiall_tMSsynch_comb2{ii},freqhiall_tMSasynch_comb2{ii});
      end
    end
    
    if audtacflag
      freqloall_APT_MSPN_comb1{ii}=ft_math(cfg,freqloall_audPtac_comb1{ii},freqloall_audMSpN_comb1{ii});
      freqhiall_APT_MSPN_comb1{ii}=ft_math(cfg,freqhiall_audPtac_comb1{ii},freqhiall_audMSpN_comb1{ii});
      if comb2flag
        freqloall_APT_MSPN_comb2{ii}=ft_math(cfg,freqloall_audPtac_comb2{ii},freqloall_audMSpN_comb2{ii});
        freqhiall_APT_MSPN_comb2{ii}=ft_math(cfg,freqhiall_audPtac_comb2{ii},freqhiall_audMSpN_comb2{ii});
      end
    end
    if usetr==2
      cfg.operation='(x1+x2)/2';
      cfg.parameter={'powspctrm'};
      
      freqloall_TPA_MSPN_comb1_usetr2{ii}=ft_math(cfg,freqloall_TPA_MSPN_comb1{ii,:});
      freqhiall_TPA_MSPN_comb1_usetr2{ii}=ft_math(cfg,freqhiall_TPA_MSPN_comb1{ii,:});
      freqloall_TPA_MSPN_comb2_usetr2{ii}=ft_math(cfg,freqloall_TPA_MSPN_comb2{ii,:});
      freqhiall_TPA_MSPN_comb2_usetr2{ii}=ft_math(cfg,freqhiall_TPA_MSPN_comb2{ii,:});
      freqloall_tacPaud_comb1tmp{ii}=ft_math(cfg,freqloall_tacPaud_comb1{ii,:});
      freqloall_tacPaud_comb2tmp{ii}=ft_math(cfg,freqloall_tacPaud_comb2{ii,:});
      freqhiall_tacPaud_comb1tmp{ii}=ft_math(cfg,freqhiall_tacPaud_comb1{ii,:});
      freqhiall_tacPaud_comb2tmp{ii}=ft_math(cfg,freqhiall_tacPaud_comb2{ii,:});
      freqloall_tacMSpN_comb1tmp{ii}=ft_math(cfg,freqloall_tacMSpN_comb1{ii,:});
      freqloall_tacMSpN_comb2tmp{ii}=ft_math(cfg,freqloall_tacMSpN_comb2{ii,:});
      freqhiall_tacMSpN_comb1tmp{ii}=ft_math(cfg,freqhiall_tacMSpN_comb1{ii,:});
      freqhiall_tacMSpN_comb2tmp{ii}=ft_math(cfg,freqhiall_tacMSpN_comb2{ii,:});
      freqloall_tTacAlone_combtmp{ii}=ft_math(cfg,freqloall_tTacAlone_comb{ii,:});
      freqloall_tAudAlone_combtmp{ii}=ft_math(cfg,freqloall_tAudAlone_comb{ii,:});
      freqloall_tNulAlone_combtmp{ii}=ft_math(cfg,freqloall_tNulAlone_comb{ii,:});
      freqloall_tMSAlone_combtmp{ii}=ft_math(cfg,freqloall_tMSAlone_comb{ii,:});
      freqhiall_tTacAlone_combtmp{ii}=ft_math(cfg,freqhiall_tTacAlone_comb{ii,:});
      freqhiall_tAudAlone_combtmp{ii}=ft_math(cfg,freqhiall_tAudAlone_comb{ii,:});
      freqhiall_tNulAlone_combtmp{ii}=ft_math(cfg,freqhiall_tNulAlone_comb{ii,:});
      freqhiall_tMSAlone_combtmp{ii}=ft_math(cfg,freqhiall_tMSAlone_comb{ii,:});
      cfg.parameter={'plvspctrm'};
      cfg.operation='(x1./abs(x1) + x2./abs(x2))/2'; % This is correct for getting plvavgang (later in ft_freqgrandaverage)
      tmp=ft_math(cfg,freqloall_TPA_MSPN_comb1{ii,:});
      freqloall_TPA_MSPN_comb1_usetr2{ii}.plvspctrm=tmp.plvspctrm;
      tmp=ft_math(cfg,freqhiall_TPA_MSPN_comb1{ii,:});
      freqhiall_TPA_MSPN_comb1_usetr2{ii}.plvspctrm=tmp.plvspctrm;
      tmp=ft_math(cfg,freqloall_TPA_MSPN_comb2{ii,:});
      freqloall_TPA_MSPN_comb2_usetr2{ii}.plvspctrm=tmp.plvspctrm;
      tmp=ft_math(cfg,freqhiall_TPA_MSPN_comb2{ii,:});
      freqhiall_TPA_MSPN_comb2_usetr2{ii}.plvspctrm=tmp.plvspctrm;
      tmp=ft_math(cfg,freqloall_tacPaud_comb1{ii,:});
      freqloall_tacPaud_comb1tmp{ii}.plvspctrm=tmp.plvspctrm;
      tmp=ft_math(cfg,freqloall_tacPaud_comb2{ii,:});
      freqloall_tacPaud_comb2tmp{ii}.plvspctrm=tmp.plvspctrm;
      tmp=ft_math(cfg,freqhiall_tacPaud_comb1{ii,:});
      freqhiall_tacPaud_comb1tmp{ii}.plvspctrm=tmp.plvspctrm;
      tmp=ft_math(cfg,freqhiall_tacPaud_comb2{ii,:});
      freqhiall_tacPaud_comb2tmp{ii}.plvspctrm=tmp.plvspctrm;
      tmp=ft_math(cfg,freqloall_tacMSpN_comb1{ii,:});
      freqloall_tacMSpN_comb1tmp{ii}.plvspctrm=tmp.plvspctrm;
      tmp=ft_math(cfg,freqloall_tacMSpN_comb2{ii,:});
      freqloall_tacMSpN_comb2tmp{ii}.plvspctrm=tmp.plvspctrm;
      tmp=ft_math(cfg,freqhiall_tacMSpN_comb1{ii,:});
      freqhiall_tacMSpN_comb1tmp{ii}.plvspctrm=tmp.plvspctrm;
      tmp=ft_math(cfg,freqhiall_tacMSpN_comb2{ii,:});
      freqhiall_tacMSpN_comb2tmp{ii}.plvspctrm=tmp.plvspctrm;
      tmp=ft_math(cfg,freqloall_tTacAlone_comb{ii,:});
      freqloall_tTacAlone_combtmp{ii}.plvspctrm=tmp.plvspctrm;
      tmp=ft_math(cfg,freqloall_tAudAlone_comb{ii,:});
      freqloall_tAudAlone_combtmp{ii}.plvspctrm=tmp.plvspctrm;
      tmp=ft_math(cfg,freqloall_tNulAlone_comb{ii,:});
      freqloall_tNulAlone_combtmp{ii}.plvspctrm=tmp.plvspctrm;
      tmp=ft_math(cfg,freqloall_tMSAlone_comb{ii,:});
      freqloall_tMSAlone_combtmp{ii}.plvspctrm=tmp.plvspctrm;
      tmp=ft_math(cfg,freqhiall_tTacAlone_comb{ii,:});
      freqhiall_tTacAlone_combtmp{ii}.plvspctrm=tmp.plvspctrm;
      tmp=ft_math(cfg,freqhiall_tAudAlone_comb{ii,:});
      freqhiall_tAudAlone_combtmp{ii}.plvspctrm=tmp.plvspctrm;
      tmp=ft_math(cfg,freqhiall_tNulAlone_comb{ii,:});
      freqhiall_tNulAlone_combtmp{ii}.plvspctrm=tmp.plvspctrm;
      tmp=ft_math(cfg,freqhiall_tMSAlone_comb{ii,:});
      freqhiall_tMSAlone_combtmp{ii}.plvspctrm=tmp.plvspctrm;
      
      cfg.parameter={'plvspctrm'};
      cfg.operation='(abs(x1) + abs(x2))/2'; % This is correct for getting plvavgabs (must alter ft_freqgrandaverage)
      tmp=ft_math(cfg,freqloall_TPA_MSPN_comb1{ii,:});
      freqloall_TPA_MSPN_comb1_usetr2{ii}.plvavgabs=tmp.plvspctrm;
      tmp=ft_math(cfg,freqhiall_TPA_MSPN_comb1{ii,:});
      freqhiall_TPA_MSPN_comb1_usetr2{ii}.plvavgabs=tmp.plvspctrm;
      tmp=ft_math(cfg,freqloall_TPA_MSPN_comb2{ii,:});
      freqloall_TPA_MSPN_comb2_usetr2{ii}.plvavgabs=tmp.plvspctrm;
      tmp=ft_math(cfg,freqhiall_TPA_MSPN_comb2{ii,:});
      freqhiall_TPA_MSPN_comb2_usetr2{ii}.plvavgabs=tmp.plvspctrm;
      tmp=ft_math(cfg,freqloall_tacPaud_comb1{ii,:});
      freqloall_tacPaud_comb1tmp{ii}.plvavgabs=tmp.plvspctrm;
      tmp=ft_math(cfg,freqloall_tacPaud_comb2{ii,:});
      freqloall_tacPaud_comb2tmp{ii}.plvavgabs=tmp.plvspctrm;
      tmp=ft_math(cfg,freqhiall_tacPaud_comb1{ii,:});
      freqhiall_tacPaud_comb1tmp{ii}.plvavgabs=tmp.plvspctrm;
      tmp=ft_math(cfg,freqhiall_tacPaud_comb2{ii,:});
      freqhiall_tacPaud_comb2tmp{ii}.plvavgabs=tmp.plvspctrm;
      tmp=ft_math(cfg,freqloall_tacMSpN_comb1{ii,:});
      freqloall_tacMSpN_comb1tmp{ii}.plvavgabs=tmp.plvspctrm;
      tmp=ft_math(cfg,freqloall_tacMSpN_comb2{ii,:});
      freqloall_tacMSpN_comb2tmp{ii}.plvavgabs=tmp.plvspctrm;
      tmp=ft_math(cfg,freqhiall_tacMSpN_comb1{ii,:});
      freqhiall_tacMSpN_comb1tmp{ii}.plvavgabs=tmp.plvspctrm;
      tmp=ft_math(cfg,freqhiall_tacMSpN_comb2{ii,:});
      freqhiall_tacMSpN_comb2tmp{ii}.plvavgabs=tmp.plvspctrm;
      tmp=ft_math(cfg,freqloall_tTacAlone_comb{ii,:});
      freqloall_tTacAlone_combtmp{ii}.plvavgabs=tmp.plvspctrm;
      tmp=ft_math(cfg,freqloall_tAudAlone_comb{ii,:});
      freqloall_tAudAlone_combtmp{ii}.plvavgabs=tmp.plvspctrm;
      tmp=ft_math(cfg,freqloall_tNulAlone_comb{ii,:});
      freqloall_tNulAlone_combtmp{ii}.plvavgabs=tmp.plvspctrm;
      tmp=ft_math(cfg,freqloall_tMSAlone_comb{ii,:});
      freqloall_tMSAlone_combtmp{ii}.plvavgabs=tmp.plvspctrm;
      tmp=ft_math(cfg,freqhiall_tTacAlone_comb{ii,:});
      freqhiall_tTacAlone_combtmp{ii}.plvavgabs=tmp.plvspctrm;
      tmp=ft_math(cfg,freqhiall_tAudAlone_comb{ii,:});
      freqhiall_tAudAlone_combtmp{ii}.plvavgabs=tmp.plvspctrm;
      tmp=ft_math(cfg,freqhiall_tNulAlone_comb{ii,:});
      freqhiall_tNulAlone_combtmp{ii}.plvavgabs=tmp.plvspctrm;
      tmp=ft_math(cfg,freqhiall_tMSAlone_comb{ii,:});
      freqhiall_tMSAlone_combtmp{ii}.plvavgabs=tmp.plvspctrm;
      
    end
  end % end ii
  if usetr==2
    freqloall_tacPaud_comb1=freqloall_tacPaud_comb1tmp;
    freqloall_tacPaud_comb2=freqloall_tacPaud_comb2tmp;
    freqhiall_tacPaud_comb1=freqhiall_tacPaud_comb1tmp;
    freqhiall_tacPaud_comb2=freqhiall_tacPaud_comb2tmp;
    freqloall_tacMSpN_comb1=freqloall_tacMSpN_comb1tmp;
    freqloall_tacMSpN_comb2=freqloall_tacMSpN_comb2tmp;
    freqhiall_tacMSpN_comb1=freqhiall_tacMSpN_comb1tmp;
    freqhiall_tacMSpN_comb2=freqhiall_tacMSpN_comb2tmp;
    freqloall_tTacAlone_comb=freqloall_tTacAlone_combtmp;
    freqloall_tAudAlone_comb=freqloall_tAudAlone_combtmp;
    freqloall_tNulAlone_comb=freqloall_tNulAlone_combtmp;
    freqloall_tMSAlone_comb=freqloall_tMSAlone_combtmp;
    freqhiall_tTacAlone_comb=freqhiall_tTacAlone_combtmp;
    freqhiall_tAudAlone_comb=freqhiall_tAudAlone_combtmp;
    freqhiall_tNulAlone_comb=freqhiall_tNulAlone_combtmp;
    freqhiall_tMSAlone_comb=freqhiall_tMSAlone_combtmp;
    clear freq*tmp
  end
  
  cfg=[];
  %   if fftaddflag
  %     gravelo_TPA_MSPN_fftadd=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_fftadd{:});
  %     gravehi_TPA_MSPN_fftadd=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_fftadd{:});
  %     if 0
  %       gravelo_TPA_MSPN_fftaddbn=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_fftaddbn{:});
  %       gravehi_TPA_MSPN_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_fftaddbn{:});
  %       gravelo_TPA_MSPN_comb1bn=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb1bn{:});
  %       gravehi_TPA_MSPN_comb1bn=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb1bn{:});
  %     end
  %
  %     if audtacflag
  %       gravelo_APT_MSPN_fftadd=ft_freqgrandaverage(cfg,freqloall_APT_MSPN_fftadd{:});
  %       gravehi_APT_MSPN_fftadd=ft_freqgrandaverage(cfg,freqhiall_APT_MSPN_fftadd{:});
  %       if 0
  %         gravelo_APT_MSPN_fftaddbn=ft_freqgrandaverage(cfg,freqloall_APT_MSPN_fftaddbn{:});
  %         gravehi_APT_MSPN_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_APT_MSPN_fftaddbn{:});
  %         gravelo_APT_MSPN_comb1bn=ft_freqgrandaverage(cfg,freqloall_APT_MSPN_comb1bn{:});
  %         gravehi_APT_MSPN_comb1bn=ft_freqgrandaverage(cfg,freqhiall_APT_MSPN_comb1bn{:});
  %       end
  %     end
  %   end
  cfg.parameter={'powspctrm' 'plvspctrm' 'plvabs'};
  if usetr==2
    gravelo_TPA_MSPN_comb1=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb1_usetr2{:});
    gravehi_TPA_MSPN_comb1=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb1_usetr2{:});
    cfg.parameter={'plvavgabs'};
    tmp=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb1_usetr2{:});
    gravelo_TPA_MSPN_comb1.plvavgabs=tmp.plvavgabs;
    tmp=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb1_usetr2{:});
    gravehi_TPA_MSPN_comb1.plvavgabs=tmp.plvavgabs;
    cfg.parameter={'powspctrm' 'plvspctrm' 'plvabs'};
  else
    cfg.parameter={'powspctrm' 'plvabs' 'plvspctrm'};
    gravelo_TPA_MSPN_comb1=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb1{:});
    gravehi_TPA_MSPN_comb1=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb1{:});
  end
  if sleep && ss==12
    cfg.parameter={'powspctrm' 'plvabs' 'plvspctrm'};
    gravelo_TPA_MSPN_nKD=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_nKD{:});
    gravehi_TPA_MSPN_nKD=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_nKD{:});
    gravelo_TPA_MSPN_nSD=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_nSD{:});
    gravehi_TPA_MSPN_nSD=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_nSD{:});
    gravelo_TPA_MSPN_nKD_nSD=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_nKD_nSD{:});
    gravehi_TPA_MSPN_nKD_nSD=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_nKD_nSD{:});
  end
  if synchasynch && ll<5
    gravelo_TMSs_TMSa_comb1=ft_freqgrandaverage(cfg,freqloall_TMSs_TMSa_comb1{:});
    gravehi_TMSs_TMSa_comb1=ft_freqgrandaverage(cfg,freqhiall_TMSs_TMSa_comb1{:});
  end
  if comb2flag
    if usetr==2
      gravelo_TPA_MSPN_comb2=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb2_usetr2{:});
      gravehi_TPA_MSPN_comb2=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb2_usetr2{:});
      cfg.parameter={'plvavgabs'};
      tmp=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb2_usetr2{:});
      gravelo_TPA_MSPN_comb2.plvavgabs=tmp.plvavgabs;
      tmp=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb2_usetr2{:});
      gravehi_TPA_MSPN_comb2.plvavgabs=tmp.plvavgabs;
      cfg.parameter={'powspctrm' 'plvspctrm' 'plvabs'};
    else
      gravelo_TPA_MSPN_comb2=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb2{:});
      gravehi_TPA_MSPN_comb2=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb2{:});
    end
    if synchasynch && ll<5
      gravelo_TMSs_TMSa_comb2=ft_freqgrandaverage(cfg,freqloall_TMSs_TMSa_comb2{:});
      gravehi_TMSs_TMSa_comb2=ft_freqgrandaverage(cfg,freqhiall_TMSs_TMSa_comb2{:});
    end
  end
  if audtacflag
    gravelo_APT_MSPN_comb1=ft_freqgrandaverage(cfg,freqloall_APT_MSPN_comb1{:});
    gravehi_APT_MSPN_comb1=ft_freqgrandaverage(cfg,freqhiall_APT_MSPN_comb1{:});
    if comb2flag
      gravelo_APT_MSPN_comb2=ft_freqgrandaverage(cfg,freqloall_APT_MSPN_comb2{:});
      gravehi_APT_MSPN_comb2=ft_freqgrandaverage(cfg,freqhiall_APT_MSPN_comb2{:});
      %       cfg.parameter={'plvavgabs'};
      %       tmp=ft_freqgrandaverage(cfg,freqloall_APT_MSPN_comb2{:});
      %       gravelo_APT_MSPN_comb2.plvavgabs=tmp.plvavgabs;
      %       tmp=ft_freqgrandaverage(cfg,freqhiall_APT_MSPN_comb2{:});
      %       gravehi_APT_MSPN_comb2.plvavgabs=tmp.plvavgabs;
      %       cfg.parameter={'powspctrm' 'plvspctrm' 'plvabs'};
    end
  end
  
  if 0
    figure;
    cfg=[];
    cfg.layout='elec1010.lay';
    cfg.baselinetype='relchange';
    cfg.baseline=[-1.2 tacbasemax(ll)];
    cfg.zlim='maxabs';
    cfg.parameter='powspctrm';
    ft_multiplotTFR(cfg,gravehi_TPA_MSPN_comb1);
    ft_multiplotTFR(cfg,gravelo_TPA_MSPN_comb1);
    cfg.parameter='plvavgabs';
    ft_multiplotTFR(cfg,gravehi_TPA_MSPN_comb1);
    ft_multiplotTFR(cfg,gravelo_TPA_MSPN_comb1);
  end
  
  
  
  % assessing unisensory
  if 0
    if plotflag
      
      for ii=1:length(subuse)
        figure;
        cfg=[];
        cfg.layout='elec1010.lay';
        cfg.baselinetype='relchange';
        cfg.baseline=[-1.2 tacbasemax(ll)];
        cfg.zlim=[-1.5 1.5];
        cfg.ylim=[4 24];
        subplot(2,2,1);ft_multiplotTFR(cfg,freqloall_tNulAlone_comb{ii});
        subplot(2,2,2);ft_multiplotTFR(cfg,freqloall_tTacAlone_comb{ii});
        subplot(2,2,3);ft_multiplotTFR(cfg,freqloall_tAudAlone_comb{ii});
        subplot(2,2,4);ft_multiplotTFR(cfg,freqloall_tMSAlone_comb{ii});
      end
    end
  end
  
  
  
  % assessing combination of conditions
  if 0
    if plotflag
      figure(20);
      for ii=1:length(subuse)
        cfg=[];
        cfg.parameter='powspctrm';
        cfg.operation='subtract';
        diff{ii}=ft_math(cfg,freqloall_tacPaud_comb1{ii},freqloall_tacMSpN_comb1{ii});
        subplot(3,length(subuse),ii);imagesc(freqloall_tacPaud_comb1{ii}.time,freqloall_tacPaud_comb1{ii}.freq,squeeze(mean(freqloall_tacPaud_comb1{ii}.powspctrm,1)));caxis([-6 6])
        subplot(3,length(subuse),length(subuse)+ii);imagesc(freqloall_tacMSpN_comb1{ii}.time,freqloall_tacPaud_comb1{ii}.freq,squeeze(mean(freqloall_tacMSpN_comb1{ii}.powspctrm,1)));caxis([-6 6])
        subplot(3,length(subuse),2*length(subuse)+ii);imagesc(diff{ii}.time,freqloall_tacPaud_comb1{ii}.freq,squeeze(mean(diff{ii}.powspctrm,1)));caxis([-6 6])
      end
      if printflag
        print(20,[fdir 'sumuni_mspn_fl_diff_ta_cond' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
      end
      figure(22);
      for ii=1:length(subuse)
        cfg=[];
        cfg.parameter='powspctrm';
        cfg.operation='subtract';
        diff{ii}=ft_math(cfg,freqhiall_tacPaud_comb1{ii},freqhiall_tacMSpN_comb1{ii});
        subplot(3,length(subuse),ii);imagesc(freqhiall_tacPaud_comb1{ii}.time,freqhiall_tacPaud_comb1{ii}.freq,squeeze(mean(freqhiall_tacPaud_comb1{ii}.powspctrm,1)));caxis([-1 1])
        subplot(3,length(subuse),length(subuse)+ii);imagesc(freqhiall_tacMSpN_comb1{ii}.time,freqhiall_tacPaud_comb1{ii}.freq,squeeze(mean(freqhiall_tacMSpN_comb1{ii}.powspctrm,1)));caxis([-1 1])
        subplot(3,length(subuse),2*length(subuse)+ii);imagesc(diff{ii}.time,freqhiall_tacPaud_comb1{ii}.freq,squeeze(mean(diff{ii}.powspctrm,1)));caxis([-.2 .2])
      end
      if printflag
        print(22,[fdir 'sumuni_mspn_fh_diff_ta_cond' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
      end
      if audtacflag
        figure(21);
        for ii=1:length(subuse)
          cfg=[];
          cfg.parameter='powspctrm';
          cfg.operation='subtract';
          diff{ii}=ft_math(cfg,freqloall_audPtac_comb1{ii},freqloall_audMSpN_comb1{ii});
          subplot(3,length(subuse),ii);imagesc(freqloall_audPtac_comb1{ii}.time,freqloall_tacPaud_comb1{ii}.freq,squeeze(mean(freqloall_audPtac_comb1{ii}.powspctrm,1)));caxis([-6 6])
          subplot(3,length(subuse),length(subuse)+ii);imagesc(freqloall_audMSpN_comb1{ii}.time,freqloall_tacPaud_comb1{ii}.freq,squeeze(mean(freqloall_audMSpN_comb1{ii}.powspctrm,1)));caxis([-6 6])
          subplot(3,length(subuse),2*length(subuse)+ii);imagesc(diff{ii}.time,freqloall_tacPaud_comb1{ii}.freq,squeeze(mean(diff{ii}.powspctrm,1)));caxis([-6 6])
        end
        if printflag
          print(21,[fdir 'sumuni_mspn_fl_diff_at_cond' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
        end
        figure(23);
        for ii=1:length(subuse)
          cfg=[];
          cfg.parameter='powspctrm';
          cfg.operation='subtract';
          diff{ii}=ft_math(cfg,freqhiall_audPtac_comb1{ii},freqhiall_audMSpN_comb1{ii});
          subplot(3,length(subuse),ii);imagesc(freqhiall_audPtac_comb1{ii}.time,freqhiall_tacPaud_comb1{ii}.freq,squeeze(mean(freqhiall_audPtac_comb1{ii}.powspctrm,1)));caxis([-1 1])
          subplot(3,length(subuse),length(subuse)+ii);imagesc(freqhiall_audMSpN_comb1{ii}.time,freqhiall_tacPaud_comb1{ii}.freq,squeeze(mean(freqhiall_audMSpN_comb1{ii}.powspctrm,1)));caxis([-1 1])
          subplot(3,length(subuse),2*length(subuse)+ii);imagesc(diff{ii}.time,freqhiall_tacPaud_comb1{ii}.freq,squeeze(mean(diff{ii}.powspctrm,1)));caxis([-.2 .2])
        end
        if printflag
          print(23,[fdir 'sumuni_mspn_fh_diff_at_cond' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
        end
      end
    end
  end
  
  
  if fftaddflag
    for ii=1:length(subuse)
      
      cfg=[];
      cfg.baselinetype='relchange';
      cfg.baseline=[-1.3 tacbasemax(ll)];
      freqloall_tacPaud_fftaddbn{ii}=ft_freqbaseline(cfg,freqloall_tacPaud_fftadd{ii});
      freqloall_tacMSpN_fftaddbn{ii}=ft_freqbaseline(cfg,freqloall_tacMSpN_fftadd{ii});
      freqhiall_tacPaud_fftaddbn{ii}=ft_freqbaseline(cfg,freqhiall_tacPaud_fftadd{ii});
      freqhiall_tacMSpN_fftaddbn{ii}=ft_freqbaseline(cfg,freqhiall_tacMSpN_fftadd{ii});
      freqloall_tacPaud_comb1bn{ii}=ft_freqbaseline(cfg,freqloall_tacPaud_comb1{ii});
      freqloall_tacMSpN_comb1bn{ii}=ft_freqbaseline(cfg,freqloall_tacMSpN_comb1{ii});
      freqhiall_tacPaud_comb1bn{ii}=ft_freqbaseline(cfg,freqhiall_tacPaud_comb1{ii});
      freqhiall_tacMSpN_comb1bn{ii}=ft_freqbaseline(cfg,freqhiall_tacMSpN_comb1{ii});
      
      if audtacflag
        cfg.baseline=[-1.3 audbasemax(ll)];
        freqloall_audPtac_fftaddbn{ii}=ft_freqbaseline(cfg,freqloall_audPtac_fftadd{ii});
        freqloall_audMSpN_fftaddbn{ii}=ft_freqbaseline(cfg,freqloall_audMSpN_fftadd{ii});
        freqhiall_audPtac_fftaddbn{ii}=ft_freqbaseline(cfg,freqhiall_audPtac_fftadd{ii});
        freqhiall_audMSpN_fftaddbn{ii}=ft_freqbaseline(cfg,freqhiall_audMSpN_fftadd{ii});
        freqloall_audPtac_comb1bn{ii}=ft_freqbaseline(cfg,freqloall_audPtac_comb1{ii});
        freqloall_audMSpN_comb1bn{ii}=ft_freqbaseline(cfg,freqloall_audMSpN_comb1{ii});
        freqhiall_audPtac_comb1bn{ii}=ft_freqbaseline(cfg,freqhiall_audPtac_comb1{ii});
        freqhiall_audMSpN_comb1bn{ii}=ft_freqbaseline(cfg,freqhiall_audMSpN_comb1{ii});
      end
      
      %       baselinetime=[dsearchn(freqloall_tacPaud_fftadd{ii}.time',-1.2) dsearchn(freqloall_tacPaud_fftadd{ii}.time',-0.7)];
      %       freqloall_tacPaud_fabn{ii}=freqloall_tacPaud_fftadd{ii};
      %       freqloall_tacPaud_fabn{ii}.powspctrm=freqloall_tacPaud_fabn{ii}.powspctrm./repmat(nanmean(freqloall_tacPaud_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_fabn{ii}.time)]);
      %       freqloall_tacMSpN_fabn{ii}=freqloall_tacMSpN_fftadd{ii};
      %       freqloall_tacMSpN_fabn{ii}.powspctrm=freqloall_tacMSpN_fabn{ii}.powspctrm./repmat(nanmean(freqloall_tacMSpN_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_fabn{ii}.time)]);
      %       freqloall_audPtac_fabn{ii}=freqloall_audPtac_fftadd{ii};
      %       freqloall_audPtac_fabn{ii}.powspctrm=freqloall_audPtac_fabn{ii}.powspctrm./repmat(nanmean(freqloall_audPtac_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_fabn{ii}.time)]);
      %       freqloall_audMSpN_fabn{ii}=freqloall_audMSpN_fftadd{ii};
      %       freqloall_audMSpN_fabn{ii}.powspctrm=freqloall_audMSpN_fabn{ii}.powspctrm./repmat(nanmean(freqloall_audMSpN_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_fabn{ii}.time)]);
      %       freqhiall_tacPaud_fabn{ii}=freqhiall_tacPaud_fftadd{ii};
      %       freqhiall_tacPaud_fabn{ii}.powspctrm=freqhiall_tacPaud_fabn{ii}.powspctrm./repmat(nanmean(freqhiall_tacPaud_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_fabn{ii}.time)]);
      %       freqhiall_tacMSpN_fabn{ii}=freqhiall_tacMSpN_fftadd{ii};
      %       freqhiall_tacMSpN_fabn{ii}.powspctrm=freqhiall_tacMSpN_fabn{ii}.powspctrm./repmat(nanmean(freqhiall_tacMSpN_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_fabn{ii}.time)]);
      %       freqhiall_audPtac_fabn{ii}=freqhiall_audPtac_fftadd{ii};
      %       freqhiall_audPtac_fabn{ii}.powspctrm=freqhiall_audPtac_fabn{ii}.powspctrm./repmat(nanmean(freqhiall_audPtac_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_fabn{ii}.time)]);
      %       freqhiall_audMSpN_fabn{ii}=freqhiall_audMSpN_fftadd{ii};
      %       freqhiall_audMSpN_fabn{ii}.powspctrm=freqhiall_audMSpN_fabn{ii}.powspctrm./repmat(nanmean(freqhiall_audMSpN_fabn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_fabn{ii}.time)]);
      %       baselinetime=[dsearchn(freqloall_tacPaud_comb1{ii}.time',-1.2) dsearchn(freqloall_tacPaud_comb1{ii}.time',-0.7)];
      %       freqloall_tacPaud_cbbn{ii}=freqloall_tacPaud_comb1{ii};
      %       freqloall_tacPaud_cbbn{ii}.powspctrm=freqloall_tacPaud_cbbn{ii}.powspctrm./repmat(nanmean(freqloall_tacPaud_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_cbbn{ii}.time)]);
      %       freqloall_tacMSpN_cbbn{ii}=freqloall_tacMSpN_comb1{ii};
      %       freqloall_tacMSpN_cbbn{ii}.powspctrm=freqloall_tacMSpN_cbbn{ii}.powspctrm./repmat(nanmean(freqloall_tacMSpN_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_cbbn{ii}.time)]);
      %       freqloall_audPtac_cbbn{ii}=freqloall_audPtac_comb1{ii};
      %       freqloall_audPtac_cbbn{ii}.powspctrm=freqloall_audPtac_cbbn{ii}.powspctrm./repmat(nanmean(freqloall_audPtac_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_cbbn{ii}.time)]);
      %       freqloall_audMSpN_cbbn{ii}=freqloall_audMSpN_comb1{ii};
      %       freqloall_audMSpN_cbbn{ii}.powspctrm=freqloall_audMSpN_cbbn{ii}.powspctrm./repmat(nanmean(freqloall_audMSpN_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqloall_tacPaud_cbbn{ii}.time)]);
      %       freqhiall_tacPaud_cbbn{ii}=freqhiall_tacPaud_comb1{ii};
      %       freqhiall_tacPaud_cbbn{ii}.powspctrm=freqhiall_tacPaud_cbbn{ii}.powspctrm./repmat(nanmean(freqhiall_tacPaud_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_cbbn{ii}.time)]);
      %       freqhiall_tacMSpN_cbbn{ii}=freqhiall_tacMSpN_comb1{ii};
      %       freqhiall_tacMSpN_cbbn{ii}.powspctrm=freqhiall_tacMSpN_cbbn{ii}.powspctrm./repmat(nanmean(freqhiall_tacMSpN_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_cbbn{ii}.time)]);
      %       freqhiall_audPtac_cbbn{ii}=freqhiall_audPtac_comb1{ii};
      %       freqhiall_audPtac_cbbn{ii}.powspctrm=freqhiall_audPtac_cbbn{ii}.powspctrm./repmat(nanmean(freqhiall_audPtac_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_cbbn{ii}.time)]);
      %       freqhiall_audMSpN_cbbn{ii}=freqhiall_audMSpN_comb1{ii};
      %       freqhiall_audMSpN_cbbn{ii}.powspctrm=freqhiall_audMSpN_cbbn{ii}.powspctrm./repmat(nanmean(freqhiall_audMSpN_cbbn{ii}.powspctrm(:,:,baselinetime(1):baselinetime(2)),3),[1 1 length(freqhiall_tacPaud_cbbn{ii}.time)]);
    end
  end
  
  cfg=[];
  cfg.keepindividual='yes';
  if fftaddflag
    grindlo_tacPaud_fftadd=ft_freqgrandaverage(cfg,freqloall_tacPaud_fftadd{:});
    grindlo_tacMSpN_fftadd=ft_freqgrandaverage(cfg,freqloall_tacMSpN_fftadd{:});
    grindhi_tacPaud_fftadd=ft_freqgrandaverage(cfg,freqhiall_tacPaud_fftadd{:});
    grindhi_tacMSpN_fftadd=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_fftadd{:});
    grindlo_tNulAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_tNulAlone_fftadd{:});
    grindhi_tNulAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_tNulAlone_fftadd{:});
    grindlo_tTacAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_tTacAlone_fftadd{:});
    grindhi_tTacAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_tTacAlone_fftadd{:});
    grindlo_tAudAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_tAudAlone_fftadd{:});
    grindhi_tAudAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_tAudAlone_fftadd{:});
    grindlo_tMSAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_tMSAlone_fftadd{:});
    grindhi_tMSAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_tMSAlone_fftadd{:});
    grindlo_tacPaud_fftaddbn=ft_freqgrandaverage(cfg,freqloall_tacPaud_fftaddbn{:});
    grindlo_tacMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqloall_tacMSpN_fftaddbn{:});
    grindhi_tacPaud_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_tacPaud_fftaddbn{:});
    grindhi_tacMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_fftaddbn{:});
    grindlo_tacPaud_comb1bn=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb1bn{:});
    grindlo_tacMSpN_comb1bn=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb1bn{:});
    grindhi_tacPaud_comb1bn=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb1bn{:});
    grindhi_tacMSpN_comb1bn=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb1bn{:});
    
    if audtacflag
      grindlo_audPtac_fftadd=ft_freqgrandaverage(cfg,freqloall_audPtac_fftadd{:});
      grindlo_audMSpN_fftadd=ft_freqgrandaverage(cfg,freqloall_audMSpN_fftadd{:});
      grindhi_audPtac_fftadd=ft_freqgrandaverage(cfg,freqhiall_audPtac_fftadd{:});
      grindhi_audMSpN_fftadd=ft_freqgrandaverage(cfg,freqhiall_audMSpN_fftadd{:});
      grindlo_audPtac_fftaddbn=ft_freqgrandaverage(cfg,freqloall_audPtac_fftaddbn{:});
      grindlo_audMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqloall_audMSpN_fftaddbn{:});
      grindhi_audPtac_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_audPtac_fftaddbn{:});
      grindhi_audMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_audMSpN_fftaddbn{:});
      grindlo_audPtac_comb1bn=ft_freqgrandaverage(cfg,freqloall_audPtac_comb1bn{:});
      grindlo_audMSpN_comb1bn=ft_freqgrandaverage(cfg,freqloall_audMSpN_comb1bn{:});
      grindhi_audPtac_comb1bn=ft_freqgrandaverage(cfg,freqhiall_audPtac_comb1bn{:});
      grindhi_audMSpN_comb1bn=ft_freqgrandaverage(cfg,freqhiall_audMSpN_comb1bn{:});
      grindlo_aNulAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_aNulAlone_fftadd{:});
      grindhi_aNulAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_aNulAlone_fftadd{:});
      grindlo_aTacAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_aTacAlone_fftadd{:});
      grindhi_aTacAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_aTacAlone_fftadd{:});
      grindlo_aAudAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_aAudAlone_fftadd{:});
      grindhi_aAudAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_aAudAlone_fftadd{:});
      grindlo_aMSAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_aMSAlone_fftadd{:});
      grindhi_aMSAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_aMSAlone_fftadd{:});
    end
  end
  
  cfg.parameter={'powspctrm' 'plvspctrm' 'plvabs'};
  grindlo_tacPaud_comb1=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb1{:});
  grindlo_tacMSpN_comb1=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb1{:});
  grindhi_tacPaud_comb1=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb1{:});
  grindhi_tacMSpN_comb1=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb1{:});
  if sleep && ss==12
    grindlo_tacPaud_nKD=ft_freqgrandaverage(cfg,freqloall_tacPaud_nKD{:});
    grindlo_tacMSpN_nKD=ft_freqgrandaverage(cfg,freqloall_tacMSpN_nKD{:});
    grindhi_tacPaud_nKD=ft_freqgrandaverage(cfg,freqhiall_tacPaud_nKD{:});
    grindhi_tacMSpN_nKD=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_nKD{:});
    grindlo_tacPaud_nSD=ft_freqgrandaverage(cfg,freqloall_tacPaud_nSD{:});
    grindlo_tacMSpN_nSD=ft_freqgrandaverage(cfg,freqloall_tacMSpN_nSD{:});
    grindhi_tacPaud_nSD=ft_freqgrandaverage(cfg,freqhiall_tacPaud_nSD{:});
    grindhi_tacMSpN_nSD=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_nSD{:});
    grindlo_tacPaud_nKD_nSD=ft_freqgrandaverage(cfg,freqloall_tacPaud_nKD_nSD{:});
    grindlo_tacMSpN_nKD_nSD=ft_freqgrandaverage(cfg,freqloall_tacMSpN_nKD_nSD{:});
    grindhi_tacPaud_nKD_nSD=ft_freqgrandaverage(cfg,freqhiall_tacPaud_nKD_nSD{:});
    grindhi_tacMSpN_nKD_nSD=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_nKD_nSD{:});
  end
  %   cfg.parameter={'plvavgabs'};
  %   tmp=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb1{:});
  %   grindlo_tacPaud_comb1.plvavgabs=tmp.plvavgabs;
  %   tmp=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb1{:});
  %   grindlo_tacMSpN_comb1.plvavgabs=tmp.plvavgabs;
  %   tmp=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb1{:});
  %   grindhi_tacPaud_comb1.plvavgabs=tmp.plvavgabs;
  %   tmp=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb1{:});
  %   grindhi_tacMSpN_comb1.plvavgabs=tmp.plvavgabs;
  %   cfg.parameter={'powspctrm' 'plvspctrm' 'plvabs'};
  if synchasynch && ll<5
    grindlo_tMSsynch_comb1=ft_freqgrandaverage(cfg,freqloall_tMSsynch_comb1{:});
    grindlo_tMSasynch_comb1=ft_freqgrandaverage(cfg,freqloall_tMSasynch_comb1{:});
    grindhi_tMSsynch_comb1=ft_freqgrandaverage(cfg,freqhiall_tMSsynch_comb1{:});
    grindhi_tMSasynch_comb1=ft_freqgrandaverage(cfg,freqhiall_tMSasynch_comb1{:});
  end
  if comb2flag
    grindlo_tacPaud_comb2=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb2{:});
    grindlo_tacMSpN_comb2=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb2{:});
    grindhi_tacPaud_comb2=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb2{:});
    grindhi_tacMSpN_comb2=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb2{:});
    %     cfg.parameter={'plvavgabs'};
    %     tmp=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb2{:});
    %     grindlo_tacPaud_comb2.plvavgabs=tmp.plvavgabs;
    %     tmp=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb2{:});
    %     grindlo_tacMSpN_comb2.plvavgabs=tmp.plvavgabs;
    %     tmp=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb2{:});
    %     grindhi_tacPaud_comb2.plvavgabs=tmp.plvavgabs;
    %     tmp=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb2{:});
    %     grindhi_tacMSpN_comb2.plvavgabs=tmp.plvavgabs;
    %     cfg.parameter={'powspctrm' 'plvspctrm' 'plvabs'};
    if synchasynch && ll<5
      grindlo_tMSsynch_comb2=ft_freqgrandaverage(cfg,freqloall_tMSsynch_comb2{:});
      grindlo_tMSasynch_comb2=ft_freqgrandaverage(cfg,freqloall_tMSasynch_comb2{:});
      grindhi_tMSsynch_comb2=ft_freqgrandaverage(cfg,freqhiall_tMSsynch_comb2{:});
      grindhi_tMSasynch_comb2=ft_freqgrandaverage(cfg,freqhiall_tMSasynch_comb2{:});
    end
  end
  if 0 % not needed really
    grindlo_tNulAlone_comb=ft_freqgrandaverage(cfg,freqloall_tNulAlone_comb{:});
    grindhi_tNulAlone_comb=ft_freqgrandaverage(cfg,freqhiall_tNulAlone_comb{:});
    grindlo_tTacAlone_comb=ft_freqgrandaverage(cfg,freqloall_tTacAlone_comb{:});
    grindhi_tTacAlone_comb=ft_freqgrandaverage(cfg,freqhiall_tTacAlone_comb{:});
    grindlo_tAudAlone_comb=ft_freqgrandaverage(cfg,freqloall_tAudAlone_comb{:});
    grindhi_tAudAlone_comb=ft_freqgrandaverage(cfg,freqhiall_tAudAlone_comb{:});
    grindlo_tMSAlone_comb=ft_freqgrandaverage(cfg,freqloall_tMSAlone_comb{:});
    grindhi_tMSAlone_comb=ft_freqgrandaverage(cfg,freqhiall_tMSAlone_comb{:});
  end
  %   cfg.parameter={'plvavgabs'};
  %   tmp=ft_freqgrandaverage(cfg,freqloall_tNulAlone_comb{:});
  %   grindlo_tNulAlone_comb.plvavgabs=tmp.plvavgabs;
  %   tmp=ft_freqgrandaverage(cfg,freqhiall_tNulAlone_comb{:});
  %   grindhi_tNulAlone_comb.plvavgabs=tmp.plvavgabs;
  %   tmp=ft_freqgrandaverage(cfg,freqloall_tTacAlone_comb{:});
  %   grindlo_tTacAlone_comb.plvavgabs=tmp.plvavgabs;
  %   tmp=ft_freqgrandaverage(cfg,freqhiall_tTacAlone_comb{:});
  %   grindhi_tTacAlone_comb.plvavgabs=tmp.plvavgabs;
  %   tmp=ft_freqgrandaverage(cfg,freqloall_tAudAlone_comb{:});
  %   grindlo_tAudAlone_comb.plvavgabs=tmp.plvavgabs;
  %   tmp=ft_freqgrandaverage(cfg,freqhiall_tAudAlone_comb{:});
  %   grindhi_tAudAlone_comb.plvavgabs=tmp.plvavgabs;
  %   tmp=ft_freqgrandaverage(cfg,freqloall_tMSAlone_comb{:});
  %   grindlo_tMSAlone_comb.plvavgabs=tmp.plvavgabs;
  %   tmp=ft_freqgrandaverage(cfg,freqhiall_tMSAlone_comb{:});
  %   grindhi_tMSAlone_comb.plvavgabs=tmp.plvavgabs;
  %   cfg.parameter={'powspctrm' 'plvspctrm' 'plvabs'};
  if usetr==2
    grindlo_TPA_MSPN_comb1=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb1_usetr2{:});
    grindhi_TPA_MSPN_comb1=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb1_usetr2{:});
    cfg.parameter={'plvavgabs'};
    tmp=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb1_usetr2{:});
    grindlo_TPA_MSPN_comb1.plvavgabs=tmp.plvavgabs;
    tmp=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb1_usetr2{:});
    grindhi_TPA_MSPN_comb1.plvavgabs=tmp.plvavgabs;
    cfg.parameter={'powspctrm' 'plvspctrm' 'plvabs'};
    if comb2flag
      grindlo_TPA_MSPN_comb2=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb2_usetr2{:});
      grindhi_TPA_MSPN_comb2=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb2_usetr2{:});
      cfg.parameter={'plvavgabs'};
      tmp=ft_freqgrandaverage(cfg,freqloall_TPA_MSPN_comb2_usetr2{:});
      grindlo_TPA_MSPN_comb2.plvavgabs=tmp.plvavgabs;
      tmp=ft_freqgrandaverage(cfg,freqhiall_TPA_MSPN_comb2_usetr2{:});
      grindhi_TPA_MSPN_comb2.plvavgabs=tmp.plvavgabs;
      cfg.parameter={'powspctrm' 'plvspctrm' 'plvabs'};
    end
  end
  
  if audtacflag
    grindlo_audPtac_comb1=ft_freqgrandaverage(cfg,freqloall_audPtac_comb1{:});
    grindlo_audMSpN_comb1=ft_freqgrandaverage(cfg,freqloall_audMSpN_comb1{:});
    grindhi_audPtac_comb1=ft_freqgrandaverage(cfg,freqhiall_audPtac_comb1{:});
    grindhi_audMSpN_comb1=ft_freqgrandaverage(cfg,freqhiall_audMSpN_comb1{:});
    if comb2flag
      grindlo_audPtac_comb2=ft_freqgrandaverage(cfg,freqloall_audPtac_comb2{:});
      grindlo_audMSpN_comb2=ft_freqgrandaverage(cfg,freqloall_audMSpN_comb2{:});
      grindhi_audPtac_comb2=ft_freqgrandaverage(cfg,freqhiall_audPtac_comb2{:});
      grindhi_audMSpN_comb2=ft_freqgrandaverage(cfg,freqhiall_audMSpN_comb2{:});
    end
    grindlo_aNulAlone_comb=ft_freqgrandaverage(cfg,freqloall_aNulAlone_comb{:});
    grindhi_aNulAlone_comb=ft_freqgrandaverage(cfg,freqhiall_aNulAlone_comb{:});
    grindlo_aTacAlone_comb=ft_freqgrandaverage(cfg,freqloall_aTacAlone_comb{:});
    grindhi_aTacAlone_comb=ft_freqgrandaverage(cfg,freqhiall_aTacAlone_comb{:});
    grindlo_aAudAlone_comb=ft_freqgrandaverage(cfg,freqloall_aAudAlone_comb{:});
    grindhi_aAudAlone_comb=ft_freqgrandaverage(cfg,freqhiall_aAudAlone_comb{:});
    grindlo_aMSAlone_comb=ft_freqgrandaverage(cfg,freqloall_aMSAlone_comb{:});
    grindhi_aMSAlone_comb=ft_freqgrandaverage(cfg,freqhiall_aMSAlone_comb{:});
  end
  
  if savegrindflag
    save([edir 'grindTFR_cond' num2str(ll) '_sleep' num2str(sleep) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '_usetr' num2str(usetr) '_mcseed' num2str(mcseed) '_itc' '.mat'],'grindlo*');
%     continue   % if want to do stats, then statsflag=1
  end
  
  
  
  
  cfg=[];
  if fftaddflag
    gravelo_tacPaud_fftadd=ft_freqgrandaverage(cfg,freqloall_tacPaud_fftadd{:});
    gravelo_tacMSpN_fftadd=ft_freqgrandaverage(cfg,freqloall_tacMSpN_fftadd{:});
    gravehi_tacPaud_fftadd=ft_freqgrandaverage(cfg,freqhiall_tacPaud_fftadd{:});
    gravehi_tacMSpN_fftadd=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_fftadd{:});
    gravelo_tNulAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_tNulAlone_fftadd{:});
    gravehi_tNulAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_tNulAlone_fftadd{:});
    gravelo_tTacAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_tTacAlone_fftadd{:});
    gravehi_tTacAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_tTacAlone_fftadd{:});
    gravelo_tAudAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_tAudAlone_fftadd{:});
    gravehi_tAudAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_tAudAlone_fftadd{:});
    gravelo_tMSAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_tMSAlone_fftadd{:});
    gravehi_tMSAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_tMSAlone_fftadd{:});
    gravelo_tacPaud_fftaddbn=ft_freqgrandaverage(cfg,freqloall_tacPaud_fftaddbn{:});
    gravelo_tacMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqloall_tacMSpN_fftaddbn{:});
    gravehi_tacPaud_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_tacPaud_fftaddbn{:});
    gravehi_tacMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_fftaddbn{:});
    gravelo_tacPaud_comb1bn=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb1bn{:});
    gravelo_tacMSpN_comb1bn=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb1bn{:});
    gravehi_tacPaud_comb1bn=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb1bn{:});
    gravehi_tacMSpN_comb1bn=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb1bn{:});
    
    if audtacflag
      gravelo_audPtac_fftadd=ft_freqgrandaverage(cfg,freqloall_audPtac_fftadd{:});
      gravelo_audMSpN_fftadd=ft_freqgrandaverage(cfg,freqloall_audMSpN_fftadd{:});
      gravehi_audPtac_fftadd=ft_freqgrandaverage(cfg,freqhiall_audPtac_fftadd{:});
      gravehi_audMSpN_fftadd=ft_freqgrandaverage(cfg,freqhiall_audMSpN_fftadd{:});
      gravelo_audPtac_fftaddbn=ft_freqgrandaverage(cfg,freqloall_audPtac_fftaddbn{:});
      gravelo_audMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqloall_audMSpN_fftaddbn{:});
      gravehi_audPtac_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_audPtac_fftaddbn{:});
      gravehi_audMSpN_fftaddbn=ft_freqgrandaverage(cfg,freqhiall_audMSpN_fftaddbn{:});
      gravelo_audPtac_comb1bn=ft_freqgrandaverage(cfg,freqloall_audPtac_comb1bn{:});
      gravelo_audMSpN_comb1bn=ft_freqgrandaverage(cfg,freqloall_audMSpN_comb1bn{:});
      gravehi_audPtac_comb1bn=ft_freqgrandaverage(cfg,freqhiall_audPtac_comb1bn{:});
      gravehi_audMSpN_comb1bn=ft_freqgrandaverage(cfg,freqhiall_audMSpN_comb1bn{:});
      gravelo_aNulAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_aNulAlone_fftadd{:});
      gravehi_aNulAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_aNulAlone_fftadd{:});
      gravelo_aTacAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_aTacAlone_fftadd{:});
      gravehi_aTacAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_aTacAlone_fftadd{:});
      gravelo_aAudAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_aAudAlone_fftadd{:});
      gravehi_aAudAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_aAudAlone_fftadd{:});
      gravelo_aMSAlone_fftadd=ft_freqgrandaverage(cfg,freqloall_aMSAlone_fftadd{:});
      gravehi_aMSAlone_fftadd=ft_freqgrandaverage(cfg,freqhiall_aMSAlone_fftadd{:});
    end
  end
  
  
  cfg.parameter={'powspctrm' 'plvspctrm' 'plvabs'};
  gravelo_tacPaud_comb1=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb1{:});
  gravelo_tacMSpN_comb1=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb1{:});
  gravehi_tacPaud_comb1=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb1{:});
  gravehi_tacMSpN_comb1=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb1{:});
  if sleep && ss==12
    gravelo_tacPaud_nKD=ft_freqgrandaverage(cfg,freqloall_tacPaud_nKD{:});
    gravelo_tacMSpN_nKD=ft_freqgrandaverage(cfg,freqloall_tacMSpN_nKD{:});
    gravehi_tacPaud_nKD=ft_freqgrandaverage(cfg,freqhiall_tacPaud_nKD{:});
    gravehi_tacMSpN_nKD=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_nKD{:});
    gravelo_tacPaud_nSD=ft_freqgrandaverage(cfg,freqloall_tacPaud_nSD{:});
    gravelo_tacMSpN_nSD=ft_freqgrandaverage(cfg,freqloall_tacMSpN_nSD{:});
    gravehi_tacPaud_nSD=ft_freqgrandaverage(cfg,freqhiall_tacPaud_nSD{:});
    gravehi_tacMSpN_nSD=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_nSD{:});
    gravelo_tacPaud_nKD_nSD=ft_freqgrandaverage(cfg,freqloall_tacPaud_nKD_nSD{:});
    gravelo_tacMSpN_nKD_nSD=ft_freqgrandaverage(cfg,freqloall_tacMSpN_nKD_nSD{:});
    gravehi_tacPaud_nKD_nSD=ft_freqgrandaverage(cfg,freqhiall_tacPaud_nKD_nSD{:});
    gravehi_tacMSpN_nKD_nSD=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_nKD_nSD{:});
  end
  %     cfg.parameter={'plvavgabs'};
  %     tmp=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb1{:});
  %     gravelo_tacPaud_comb1.plvavgabs=tmp.plvavgabs;
  %     tmp=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb1{:});
  %     gravelo_tacMSpN_comb1.plvavgabs=tmp.plvavgabs;
  %     tmp=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb1{:});
  %     gravehi_tacPaud_comb1.plvavgabs=tmp.plvavgabs;
  %     tmp=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb1{:});
  %     gravehi_tacMSpN_comb1.plvavgabs=tmp.plvavgabs;
  %     cfg.parameter={'powspctrm' 'plvspctrm' 'plvabs'};
  if synchasynch && ll<5
    gravelo_tMSsynch_comb1=ft_freqgrandaverage(cfg,freqloall_tMSsynch_comb1{:});
    gravelo_tMSasynch_comb1=ft_freqgrandaverage(cfg,freqloall_tMSasynch_comb1{:});
    gravehi_tMSsynch_comb1=ft_freqgrandaverage(cfg,freqhiall_tMSsynch_comb1{:});
    gravehi_tMSasynch_comb1=ft_freqgrandaverage(cfg,freqhiall_tMSasynch_comb1{:});
  end
  if comb2flag
    gravelo_tacPaud_comb2=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb2{:});
    gravelo_tacMSpN_comb2=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb2{:});
    gravehi_tacPaud_comb2=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb2{:});
    gravehi_tacMSpN_comb2=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb2{:});
    %     cfg.parameter={'plvavgabs'};
    %     tmp=ft_freqgrandaverage(cfg,freqloall_tacPaud_comb2{:});
    %     gravelo_tacPaud_comb2.plvavgabs=tmp.plvavgabs;
    %     tmp=ft_freqgrandaverage(cfg,freqloall_tacMSpN_comb2{:});
    %     gravelo_tacMSpN_comb2.plvavgabs=tmp.plvavgabs;
    %     tmp=ft_freqgrandaverage(cfg,freqhiall_tacPaud_comb2{:});
    %     gravehi_tacPaud_comb2.plvavgabs=tmp.plvavgabs;
    %     tmp=ft_freqgrandaverage(cfg,freqhiall_tacMSpN_comb2{:});
    %     gravehi_tacMSpN_comb2.plvavgabs=tmp.plvavgabs;
    %     cfg.parameter={'powspctrm' 'plvspctrm' 'plvabs'};
    if synchasynch && ll<5
      gravelo_tMSsynch_comb2=ft_freqgrandaverage(cfg,freqloall_tMSsynch_comb2{:});
      gravelo_tMSasynch_comb2=ft_freqgrandaverage(cfg,freqloall_tMSasynch_comb2{:});
      gravehi_tMSsynch_comb2=ft_freqgrandaverage(cfg,freqhiall_tMSsynch_comb2{:});
      gravehi_tMSasynch_comb2=ft_freqgrandaverage(cfg,freqhiall_tMSasynch_comb2{:});
    end
  end
  if 0
    gravelo_tNulAlone_comb=ft_freqgrandaverage(cfg,freqloall_tNulAlone_comb{:});
    gravehi_tNulAlone_comb=ft_freqgrandaverage(cfg,freqhiall_tNulAlone_comb{:});
    gravelo_tTacAlone_comb=ft_freqgrandaverage(cfg,freqloall_tTacAlone_comb{:});
    gravehi_tTacAlone_comb=ft_freqgrandaverage(cfg,freqhiall_tTacAlone_comb{:});
    gravelo_tAudAlone_comb=ft_freqgrandaverage(cfg,freqloall_tAudAlone_comb{:});
    gravehi_tAudAlone_comb=ft_freqgrandaverage(cfg,freqhiall_tAudAlone_comb{:});
    gravelo_tMSAlone_comb=ft_freqgrandaverage(cfg,freqloall_tMSAlone_comb{:});
    gravehi_tMSAlone_comb=ft_freqgrandaverage(cfg,freqhiall_tMSAlone_comb{:});
  end
  %     cfg.parameter={'plvavgabs'};
  %     tmp=ft_freqgrandaverage(cfg,freqloall_tNulAlone_comb{:});
  %     gravelo_tNulAlone_comb.plvavgabs=tmp.plvavgabs;
  %     tmp=ft_freqgrandaverage(cfg,freqhiall_tNulAlone_comb{:});
  %     gravehi_tNulAlone_comb.plvavgabs=tmp.plvavgabs;
  %     tmp=ft_freqgrandaverage(cfg,freqloall_tTacAlone_comb{:});
  %     gravelo_tTacAlone_comb.plvavgabs=tmp.plvavgabs;
  %     tmp=ft_freqgrandaverage(cfg,freqhiall_tTacAlone_comb{:});
  %     gravehi_tTacAlone_comb.plvavgabs=tmp.plvavgabs;
  %     tmp=ft_freqgrandaverage(cfg,freqloall_tAudAlone_comb{:});
  %     gravelo_tAudAlone_comb.plvavgabs=tmp.plvavgabs;
  %     tmp=ft_freqgrandaverage(cfg,freqhiall_tAudAlone_comb{:});
  %     gravehi_tAudAlone_comb.plvavgabs=tmp.plvavgabs;
  %     tmp=ft_freqgrandaverage(cfg,freqloall_tMSAlone_comb{:});
  %     gravelo_tMSAlone_comb.plvavgabs=tmp.plvavgabs;
  %     tmp=ft_freqgrandaverage(cfg,freqhiall_tMSAlone_comb{:});
  %     gravehi_tMSAlone_comb.plvavgabs=tmp.plvavgabs;
  %     cfg.parameter={'powspctrm' 'plvspctrm' 'plvabs'};
  
  if audtacflag
    gravelo_audPtac_comb1=ft_freqgrandaverage(cfg,freqloall_audPtac_comb1{:});
    gravelo_audMSpN_comb1=ft_freqgrandaverage(cfg,freqloall_audMSpN_comb1{:});
    gravehi_audPtac_comb1=ft_freqgrandaverage(cfg,freqhiall_audPtac_comb1{:});
    gravehi_audMSpN_comb1=ft_freqgrandaverage(cfg,freqhiall_audMSpN_comb1{:});
    if comb2flag
      gravelo_audPtac_comb2=ft_freqgrandaverage(cfg,freqloall_audPtac_comb2{:});
      gravelo_audMSpN_comb2=ft_freqgrandaverage(cfg,freqloall_audMSpN_comb2{:});
      gravehi_audPtac_comb2=ft_freqgrandaverage(cfg,freqhiall_audPtac_comb2{:});
      gravehi_audMSpN_comb2=ft_freqgrandaverage(cfg,freqhiall_audMSpN_comb2{:});
    end
    gravelo_aNulAlone_comb=ft_freqgrandaverage(cfg,freqloall_aNulAlone_comb{:});
    gravehi_aNulAlone_comb=ft_freqgrandaverage(cfg,freqhiall_aNulAlone_comb{:});
    gravelo_aTacAlone_comb=ft_freqgrandaverage(cfg,freqloall_aTacAlone_comb{:});
    gravehi_aTacAlone_comb=ft_freqgrandaverage(cfg,freqhiall_aTacAlone_comb{:});
    gravelo_aAudAlone_comb=ft_freqgrandaverage(cfg,freqloall_aAudAlone_comb{:});
    gravehi_aAudAlone_comb=ft_freqgrandaverage(cfg,freqhiall_aAudAlone_comb{:});
    gravelo_aMSAlone_comb=ft_freqgrandaverage(cfg,freqloall_aMSAlone_comb{:});
    gravehi_aMSAlone_comb=ft_freqgrandaverage(cfg,freqhiall_aMSAlone_comb{:});
  end
  
  if plotflag
    powalonemaxz=[0.5 0.4 0.3 0.2 0.2]; % one per frequency band
    plvabsmaxz=[.6 .4 .2 .1 .1];
    % individual conditions first
    
    % fftadd first
    if 0 % this is redundant with comb for non-added-together conditions
      topoplotTFR_highlight(61,gravelo_tNulAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(62,gravelo_tTacAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(63,gravelo_tAudAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(64,gravelo_tMSAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(71,gravelo_tNulAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(72,gravelo_tTacAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(73,gravelo_tAudAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(74,gravelo_tMSAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(81,gravelo_tNulAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(82,gravelo_tTacAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(83,gravelo_tAudAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(84,gravelo_tMSAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(66,gravehi_tNulAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(67,gravehi_tTacAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(68,gravehi_tAudAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(69,gravehi_tMSAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(76,gravehi_tNulAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(77,gravehi_tTacAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(78,gravehi_tAudAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05);
      topoplotTFR_highlight(79,gravehi_tMSAlone_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05);
      if printflag
        ylim=ylimlo(1,:);
        print(61,[fdir 'gravelo_tNulAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(62,[fdir 'gravelo_tTacAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(63,[fdir 'gravelo_tAudAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(64,[fdir 'gravelo_tMSAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        ylim=ylimlo(2,:);
        print(71,[fdir 'gravelo_tNulAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(72,[fdir 'gravelo_tTacAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(73,[fdir 'gravelo_tAudAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(74,[fdir 'gravelo_tMSAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        ylim=ylimlo(3,:);
        print(81,[fdir 'gravelo_tNulAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(82,[fdir 'gravelo_tTacAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(83,[fdir 'gravelo_tAudAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(84,[fdir 'gravelo_tMSAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        ylim=ylimhi(1,:);
        print(66,[fdir 'gravehi_tNulAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(67,[fdir 'gravehi_tTacAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(68,[fdir 'gravehi_tAudAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(69,[fdir 'gravehi_tMSAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        ylim=ylimhi(2,:);
        print(76,[fdir 'gravehi_tNulAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(77,[fdir 'gravehi_tTacAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(78,[fdir 'gravehi_tAudAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(79,[fdir 'gravehi_tMSAlone_fftadd_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      end
      close all
    end
    
    % comb second
    topoplotTFR_highlight(61,gravelo_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(1) powalonemaxz(1)]);
    topoplotTFR_highlight(62,gravelo_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(1) powalonemaxz(1)]);
    topoplotTFR_highlight(63,gravelo_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(1) powalonemaxz(1)]);
    topoplotTFR_highlight(64,gravelo_tMSAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(1) powalonemaxz(1)]);
    topoplotTFR_highlight(71,gravelo_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(2) powalonemaxz(2)]);
    topoplotTFR_highlight(72,gravelo_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(2) powalonemaxz(2)]);
    topoplotTFR_highlight(73,gravelo_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(2) powalonemaxz(2)]);
    topoplotTFR_highlight(74,gravelo_tMSAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(2) powalonemaxz(2)]);
    topoplotTFR_highlight(81,gravelo_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(3) powalonemaxz(3)]);
    topoplotTFR_highlight(82,gravelo_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(3) powalonemaxz(3)]);
    topoplotTFR_highlight(83,gravelo_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(3) powalonemaxz(3)]);
    topoplotTFR_highlight(84,gravelo_tMSAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(3) powalonemaxz(3)]);
    topoplotTFR_highlight(66,gravehi_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(4) powalonemaxz(4)]);
    topoplotTFR_highlight(67,gravehi_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(4) powalonemaxz(4)]);
    topoplotTFR_highlight(68,gravehi_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(4) powalonemaxz(4)]);
    topoplotTFR_highlight(69,gravehi_tMSAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(4) powalonemaxz(4)]);
    topoplotTFR_highlight(76,gravehi_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(5) powalonemaxz(5)]);
    topoplotTFR_highlight(77,gravehi_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(5) powalonemaxz(5)]);
    topoplotTFR_highlight(78,gravehi_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(5) powalonemaxz(5)]);
    topoplotTFR_highlight(79,gravehi_tMSAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(5) powalonemaxz(5)]);
    if printflag
      ylim=ylimlo(1,:);
      print(61,[fdir 'gravelo_tNulAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(62,[fdir 'gravelo_tTacAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(63,[fdir 'gravelo_tAudAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(64,[fdir 'gravelo_tMSAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(2,:);
      print(71,[fdir 'gravelo_tNulAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(72,[fdir 'gravelo_tTacAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(73,[fdir 'gravelo_tAudAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(74,[fdir 'gravelo_tMSAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(3,:);
      print(81,[fdir 'gravelo_tNulAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(82,[fdir 'gravelo_tTacAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(83,[fdir 'gravelo_tAudAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(84,[fdir 'gravelo_tMSAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(1,:);
      print(66,[fdir 'gravehi_tNulAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(67,[fdir 'gravehi_tTacAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(68,[fdir 'gravehi_tAudAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(69,[fdir 'gravehi_tMSAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(2,:);
      print(76,[fdir 'gravehi_tNulAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(77,[fdir 'gravehi_tTacAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(78,[fdir 'gravehi_tAudAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(79,[fdir 'gravehi_tMSAlone_comb_topoOverTime_mc_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
    end
    close all
    
    % comb with PLV
    topoplotTFR_highlight(141,gravelo_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(142,gravelo_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(143,gravelo_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(144,gravelo_tMSAlone_comb ,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(151,gravelo_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(152,gravelo_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(153,gravelo_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(154,gravelo_tMSAlone_comb ,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(161,gravelo_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(162,gravelo_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(163,gravelo_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(164,gravelo_tMSAlone_comb ,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(146,gravehi_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(147,gravehi_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(148,gravehi_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(149,gravehi_tMSAlone_comb ,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(156,gravehi_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(5)]);
    topoplotTFR_highlight(157,gravehi_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(5)]);
    topoplotTFR_highlight(158,gravehi_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(5)]);
    topoplotTFR_highlight(159,gravehi_tMSAlone_comb ,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(5)]);
    
    topoplotTFR_highlight(1141,gravelo_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(1142,gravelo_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(1143,gravelo_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(1144,gravelo_tMSAlone_comb ,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(1151,gravelo_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(1152,gravelo_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(1153,gravelo_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(1154,gravelo_tMSAlone_comb ,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(1161,gravelo_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(1162,gravelo_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(1163,gravelo_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(1164,gravelo_tMSAlone_comb ,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(1146,gravehi_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(1147,gravehi_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(1148,gravehi_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(1149,gravehi_tMSAlone_comb ,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(1156,gravehi_tNulAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(5)]);
    topoplotTFR_highlight(1157,gravehi_tTacAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(5)]);
    topoplotTFR_highlight(1158,gravehi_tAudAlone_comb,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(5)]);
    topoplotTFR_highlight(1159,gravehi_tMSAlone_comb ,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(5)]);
    if printflag
      ylim=ylimlo(1,:);
      print(141,[fdir 'gravelo_tNulAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(142,[fdir 'gravelo_tTacAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(143,[fdir 'gravelo_tAudAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(144,[fdir 'gravelo_tMSAlone_comb_plvabsavg_topoOverTime_ta_cond'  num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(2,:);
      print(151,[fdir 'gravelo_tNulAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(152,[fdir 'gravelo_tTacAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(153,[fdir 'gravelo_tAudAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(154,[fdir 'gravelo_tMSAlone_comb_plvabsavg_topoOverTime_ta_cond'  num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(3,:);
      print(161,[fdir 'gravelo_tNulAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(162,[fdir 'gravelo_tTacAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(163,[fdir 'gravelo_tAudAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(164,[fdir 'gravelo_tMSAlone_comb_plvabsavg_topoOverTime_ta_cond'  num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(1,:);
      print(146,[fdir 'gravehi_tNulAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(147,[fdir 'gravehi_tTacAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(148,[fdir 'gravehi_tAudAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(149,[fdir 'gravehi_tMSAlone_comb_plvabsavg_topoOverTime_ta_cond'  num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(2,:);
      print(156,[fdir 'gravehi_tNulAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(157,[fdir 'gravehi_tTacAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(158,[fdir 'gravehi_tAudAlone_comb_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(159,[fdir 'gravehi_tMSAlone_comb_plvabsavg_topoOverTime_ta_cond'  num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      
      ylim=ylimlo(1,:);
      print(1141,[fdir 'gravelo_tNulAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1142,[fdir 'gravelo_tTacAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1143,[fdir 'gravelo_tAudAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1144,[fdir 'gravelo_tMSAlone_comb_plvavgabs_topoOverTime_ta_cond'  num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(2,:);
      print(1151,[fdir 'gravelo_tNulAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1152,[fdir 'gravelo_tTacAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1153,[fdir 'gravelo_tAudAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1154,[fdir 'gravelo_tMSAlone_comb_plvavgabs_topoOverTime_ta_cond'  num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(3,:);
      print(1161,[fdir 'gravelo_tNulAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1162,[fdir 'gravelo_tTacAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1163,[fdir 'gravelo_tAudAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1164,[fdir 'gravelo_tMSAlone_comb_plvavgabs_topoOverTime_ta_cond'  num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(1,:);
      print(1146,[fdir 'gravehi_tNulAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1147,[fdir 'gravehi_tTacAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1148,[fdir 'gravehi_tAudAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1149,[fdir 'gravehi_tMSAlone_comb_plvavgabs_topoOverTime_ta_cond'  num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(2,:);
      print(1156,[fdir 'gravehi_tNulAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1157,[fdir 'gravehi_tTacAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1158,[fdir 'gravehi_tAudAlone_comb_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1159,[fdir 'gravehi_tMSAlone_comb_plvavgabs_topoOverTime_ta_cond'  num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
    end
    close all
    %          % single condition alone plvavgang and plvangavg ??
    
    
    % contrast of conditions second
    % power first
    if fftaddflag
      topoplotTFR_highlight(110,gravelo_tacPaud_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(1) powalonemaxz(1)]);
      topoplotTFR_highlight(111,gravelo_tacMSpN_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(1) powalonemaxz(1)]);
      topoplotTFR_highlight(112,gravelo_tacPaud_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(2) powalonemaxz(2)]);
      topoplotTFR_highlight(113,gravelo_tacMSpN_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(2) powalonemaxz(2)]);
      topoplotTFR_highlight(114,gravelo_tacPaud_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(3) powalonemaxz(3)]);
      topoplotTFR_highlight(115,gravelo_tacMSpN_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(3) powalonemaxz(3)]);
      topoplotTFR_highlight(116,gravehi_tacPaud_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(4) powalonemaxz(4)]);
      topoplotTFR_highlight(117,gravehi_tacMSpN_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(4) powalonemaxz(4)]);
      topoplotTFR_highlight(118,gravehi_tacPaud_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(5) powalonemaxz(5)]);
      topoplotTFR_highlight(119,gravehi_tacMSpN_fftadd,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(5) powalonemaxz(5)]);
    end
    
    topoplotTFR_highlight(120,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(1) powalonemaxz(1)]);
    topoplotTFR_highlight(121,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(1) powalonemaxz(1)]);
    topoplotTFR_highlight(122,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(2) powalonemaxz(2)]);
    topoplotTFR_highlight(123,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(2) powalonemaxz(2)]);
    topoplotTFR_highlight(124,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(3) powalonemaxz(3)]);
    topoplotTFR_highlight(125,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(3) powalonemaxz(3)]);
    topoplotTFR_highlight(126,gravehi_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(4) powalonemaxz(4)]);
    topoplotTFR_highlight(127,gravehi_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(4) powalonemaxz(4)]);
    topoplotTFR_highlight(128,gravehi_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(5) powalonemaxz(5)]);
    topoplotTFR_highlight(129,gravehi_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05,'powspctrm',[-powalonemaxz(5) powalonemaxz(5)]);
    
    if comb2flag
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
    end
    
    if printflag
      ylim=ylimlo(1,:);
      if fftaddflag
        print(110,[fdir 'gravelo_tacPaud_fftadd_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(111,[fdir 'gravelo_tacMSpN_fftadd_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      end
      print(120,[fdir 'gravelo_tacPaud_comb1_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(121,[fdir 'gravelo_tacMSpN_comb1_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(130,[fdir 'gravelo_tacPaud_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(131,[fdir 'gravelo_tacMSpN_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(2,:);
      if fftaddflag
        print(112,[fdir 'gravelo_tacPaud_fftadd_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(113,[fdir 'gravelo_tacMSpN_fftadd_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      end
      print(122,[fdir 'gravelo_tacPaud_comb1_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(123,[fdir 'gravelo_tacMSpN_comb1_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(132,[fdir 'gravelo_tacPaud_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(133,[fdir 'gravelo_tacMSpN_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(3,:);
      if fftaddflag
        print(114,[fdir 'gravelo_tacPaud_fftadd_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(115,[fdir 'gravelo_tacMSpN_fftadd_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      end
      print(124,[fdir 'gravelo_tacPaud_comb1_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(125,[fdir 'gravelo_tacMSpN_comb1_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(134,[fdir 'gravelo_tacPaud_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(135,[fdir 'gravelo_tacMSpN_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(1,:);
      if fftaddflag
        print(116,[fdir 'gravehi_tacPaud_fftadd_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(117,[fdir 'gravehi_tacMSpN_fftadd_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      end
      print(126,[fdir 'gravehi_tacPaud_comb1_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(127,[fdir 'gravehi_tacMSpN_comb1_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(136,[fdir 'gravehi_tacPaud_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(137,[fdir 'gravehi_tacMSpN_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(2,:);
      if fftaddflag
        print(118,[fdir 'gravehi_tacPaud_fftadd_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        print(119,[fdir 'gravehi_tacMSpN_fftadd_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      end
      print(128,[fdir 'gravehi_tacPaud_comb1_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(129,[fdir 'gravehi_tacMSpN_comb1_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(138,[fdir 'gravehi_tacPaud_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(139,[fdir 'gravehi_tacMSpN_comb2_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
    end
    close all
    
    % plv second
    topoplotTFR_highlight(220,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(221,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(222,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(223,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(224,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(225,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(226,gravehi_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(227,gravehi_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(228,gravehi_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(5)]);
    topoplotTFR_highlight(229,gravehi_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvabsavg',[0 plvabsmaxz(5)]);
    topoplotTFR_highlight(1220,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(1221,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(1)]);
    topoplotTFR_highlight(1222,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(1223,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(2)]);
    topoplotTFR_highlight(1224,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(1225,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(3)]);
    topoplotTFR_highlight(1226,gravehi_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(1227,gravehi_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(4)]);
    topoplotTFR_highlight(1228,gravehi_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(5)]);
    topoplotTFR_highlight(1229,gravehi_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgabs',[0 plvabsmaxz(5)]);
    
    topoplotTFR_highlight(250,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvangavg',[-pi pi]);
    topoplotTFR_highlight(250,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvangavg',[-pi pi]);
    topoplotTFR_highlight(251,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvangavg',[-pi pi]);
    topoplotTFR_highlight(252,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvangavg',[-pi pi]);
    topoplotTFR_highlight(253,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvangavg',[-pi pi]);
    topoplotTFR_highlight(254,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvangavg',[-pi pi]);
    topoplotTFR_highlight(255,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvangavg',[-pi pi]);
    topoplotTFR_highlight(256,gravehi_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvangavg',[-pi pi]);
    topoplotTFR_highlight(257,gravehi_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvangavg',[-pi pi]);
    topoplotTFR_highlight(258,gravehi_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvangavg',[-pi pi]);
    topoplotTFR_highlight(259,gravehi_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvangavg',[-pi pi]);
    topoplotTFR_highlight(1250,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgang',[-pi pi]);
    topoplotTFR_highlight(1251,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(1,:),[],[],0.05,'plvavgang',[-pi pi]);
    topoplotTFR_highlight(1252,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgang',[-pi pi]);
    topoplotTFR_highlight(1253,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[],[],0.05,'plvavgang',[-pi pi]);
    topoplotTFR_highlight(1254,gravelo_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgang',[-pi pi]);
    topoplotTFR_highlight(1255,gravelo_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(3,:),[],[],0.05,'plvavgang',[-pi pi]);
    topoplotTFR_highlight(1256,gravehi_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgang',[-pi pi]);
    topoplotTFR_highlight(1257,gravehi_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(1,:),[],[],0.05,'plvavgang',[-pi pi]);
    topoplotTFR_highlight(1258,gravehi_tacPaud_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgang',[-pi pi]);
    topoplotTFR_highlight(1259,gravehi_tacMSpN_comb1,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[],[],0.05,'plvavgang',[-pi pi]);
    
    if comb2flag
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
    end
    
    if printflag
      ylim=ylimlo(1,:);
      print(220,[fdir 'gravelo_tacPaud_comb1_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(221,[fdir 'gravelo_tacMSpN_comb1_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(250,[fdir 'gravelo_tacPaud_comb1_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(251,[fdir 'gravelo_tacMSpN_comb1_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(2,:);
      print(222,[fdir 'gravelo_tacPaud_comb1_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(223,[fdir 'gravelo_tacMSpN_comb1_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(252,[fdir 'gravelo_tacPaud_comb1_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(253,[fdir 'gravelo_tacMSpN_comb1_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(3,:);
      print(224,[fdir 'gravelo_tacPaud_comb1_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(225,[fdir 'gravelo_tacMSpN_comb1_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(254,[fdir 'gravelo_tacPaud_comb1_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(255,[fdir 'gravelo_tacMSpN_comb1_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(1,:);
      print(226,[fdir 'gravehi_tacPaud_comb1_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(227,[fdir 'gravehi_tacMSpN_comb1_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(256,[fdir 'gravehi_tacPaud_comb1_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(257,[fdir 'gravehi_tacMSpN_comb1_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(2,:);
      print(228,[fdir 'gravehi_tacPaud_comb1_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(229,[fdir 'gravehi_tacMSpN_comb1_plvabsavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(258,[fdir 'gravehi_tacPaud_comb1_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(259,[fdir 'gravehi_tacMSpN_comb1_plvangavg_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      
      ylim=ylimlo(1,:);
      print(1220,[fdir 'gravelo_tacPaud_comb1_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1221,[fdir 'gravelo_tacMSpN_comb1_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1250,[fdir 'gravelo_tacPaud_comb1_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1251,[fdir 'gravelo_tacMSpN_comb1_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(2,:);
      print(1222,[fdir 'gravelo_tacPaud_comb1_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1223,[fdir 'gravelo_tacMSpN_comb1_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1252,[fdir 'gravelo_tacPaud_comb1_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1253,[fdir 'gravelo_tacMSpN_comb1_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimlo(3,:);
      print(1224,[fdir 'gravelo_tacPaud_comb1_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1225,[fdir 'gravelo_tacMSpN_comb1_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1254,[fdir 'gravelo_tacPaud_comb1_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1255,[fdir 'gravelo_tacMSpN_comb1_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(1,:);
      print(1226,[fdir 'gravehi_tacPaud_comb1_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1227,[fdir 'gravehi_tacMSpN_comb1_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1256,[fdir 'gravehi_tacPaud_comb1_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1257,[fdir 'gravehi_tacMSpN_comb1_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      ylim=ylimhi(2,:);
      print(1228,[fdir 'gravehi_tacPaud_comb1_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1229,[fdir 'gravehi_tacMSpN_comb1_plvavgabs_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1258,[fdir 'gravehi_tacPaud_comb1_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      print(1259,[fdir 'gravehi_tacMSpN_comb1_plvavgang_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
      
      if comb2flag
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
    end
    close all
    
    
    
    %       if 0
    %       topoplotTFR_highlight(31,gravelo_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       topoplotTFR_highlight(32,gravelo_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       topoplotTFR_highlight(33,gravehi_tacPaud_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       topoplotTFR_highlight(34,gravehi_tacMSpN_comb2,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       topoplotTFR_highlight(41,gravelo_tacPaud_fftaddbn,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       topoplotTFR_highlight(42,gravelo_tacMSpN_fftaddbn,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       topoplotTFR_highlight(43,gravehi_tacPaud_fftaddbn,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       topoplotTFR_highlight(44,gravehi_tacMSpN_fftaddbn,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       topoplotTFR_highlight(51,gravelo_tacPaud_comb1bn,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       topoplotTFR_highlight(52,gravelo_tacMSpN_comb1bn,[tacbasemax(ll); .55-audbasemax(ll)],ylimlo(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       topoplotTFR_highlight(53,gravehi_tacPaud_comb1bn,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       topoplotTFR_highlight(54,gravehi_tacMSpN_comb1bn,[tacbasemax(ll); .55-audbasemax(ll)],ylimhi(2,:),[-1.2 tacbasemax(ll)],[],0.05);
    %       end
    
    if audtacflag
      topoplotTFR_highlight(15,gravelo_audPtac_fftadd,[audbasemax(ll); .55-tacbasemax(ll)],ylimlo(2,:),[-1.2 audbasemax(ll)],[],0.05);
      topoplotTFR_highlight(16,gravelo_audMSpN_fftadd,[audbasemax(ll); .55-tacbasemax(ll)],ylimlo(2,:),[-1.2 audbasemax(ll)],[],0.05);
      topoplotTFR_highlight(17,gravehi_audPtac_fftadd,[audbasemax(ll); .55-tacbasemax(ll)],ylimhi(2,:),[-1.2 audbasemax(ll)],[],0.05);
      topoplotTFR_highlight(18,gravehi_audMSpN_fftadd,[audbasemax(ll); .55-tacbasemax(ll)],ylimhi(2,:),[-1.2 audbasemax(ll)],[],0.05);
      
      topoplotTFR_highlight(25,gravelo_audPtac_comb1,[audbasemax(ll); .55-tacbasemax(ll)],ylimlo(2,:),[-1.2 audbasemax(ll)],[],0.05);
      topoplotTFR_highlight(26,gravelo_audMSpN_comb1,[audbasemax(ll); .55-tacbasemax(ll)],ylimlo(2,:),[-1.2 audbasemax(ll)],[],0.05);
      topoplotTFR_highlight(27,gravehi_audPtac_comb1,[audbasemax(ll); .55-tacbasemax(ll)],ylimhi(2,:),[-1.2 audbasemax(ll)],[],0.05);
      topoplotTFR_highlight(28,gravehi_audMSpN_comb1,[audbasemax(ll); .55-tacbasemax(ll)],ylimhi(2,:),[-1.2 audbasemax(ll)],[],0.05);
      
      if comb2flag
        topoplotTFR_highlight(35,gravelo_audPtac_comb2,[audbasemax(ll); .55-tacbasemax(ll)],ylimlo(2,:),[-1.2 audbasemax(ll)],[],0.05);
        topoplotTFR_highlight(36,gravelo_audMSpN_comb2,[audbasemax(ll); .55-tacbasemax(ll)],ylimlo(2,:),[-1.2 audbasemax(ll)],[],0.05);
        topoplotTFR_highlight(37,gravehi_audPtac_comb2,[audbasemax(ll); .55-tacbasemax(ll)],ylimhi(2,:),[-1.2 audbasemax(ll)],[],0.05);
        topoplotTFR_highlight(38,gravehi_audMSpN_comb2,[audbasemax(ll); .55-tacbasemax(ll)],ylimhi(2,:),[-1.2 audbasemax(ll)],[],0.05);
      end
      
      if 0
        topoplotTFR_highlight(45,gravelo_audPtac_fftaddbn,[audbasemax(ll); .55-tacbasemax(ll)],ylimlo(2,:),[-1.2 audbasemax(ll)],[],0.05);
        topoplotTFR_highlight(46,gravelo_audMSpN_fftaddbn,[audbasemax(ll); .55-tacbasemax(ll)],ylimlo(2,:),[-1.2 audbasemax(ll)],[],0.05);
        topoplotTFR_highlight(47,gravehi_audPtac_fftaddbn,[audbasemax(ll); .55-tacbasemax(ll)],ylimhi(2,:),[-1.2 audbasemax(ll)],[],0.05);
        topoplotTFR_highlight(48,gravehi_audMSpN_fftaddbn,[audbasemax(ll); .55-tacbasemax(ll)],ylimhi(2,:),[-1.2 audbasemax(ll)],[],0.05);
        
        topoplotTFR_highlight(55,gravelo_audPtac_comb1bn,[audbasemax(ll); .55-tacbasemax(ll)],ylimlo(2,:),[-1.2 audbasemax(ll)],[],0.05);
        topoplotTFR_highlight(56,gravelo_audMSpN_comb1bn,[audbasemax(ll); .55-tacbasemax(ll)],ylimlo(2,:),[-1.2 audbasemax(ll)],[],0.05);
        topoplotTFR_highlight(57,gravehi_audPtac_comb1bn,[audbasemax(ll); .55-tacbasemax(ll)],ylimhi(2,:),[-1.2 audbasemax(ll)],[],0.05);
        topoplotTFR_highlight(58,gravehi_audMSpN_comb1bn,[audbasemax(ll); .55-tacbasemax(ll)],ylimhi(2,:),[-1.2 audbasemax(ll)],[],0.05);
      end
    end
    
    
  end % plotflag
  
  if 0
     % get multiplot going as well...want to see TFR and PLV for channel-cluster!
    figure;
    cfg=[];
    cfg.layout='elec1010.lay';
    cfg.baselinetype='relchange';
    cfg.baseline=[-1.2 tacbasemax(ll)];
    cfg.zlim='maxabs';
    cfg.ylim=[4 24];
    ft_multiplotTFR(cfg,gravelo_tacPaud_comb1);
    ft_multiplotTFR(cfg,gravelo_tacMSpN_comb1);
    cfg.parameter='plvavgabs';
    ft_multiplotTFR(cfg,gravelo_tacPaud_comb1);
    ft_multiplotTFR(cfg,gravelo_tacMSpN_comb1);
  end
  
  
  
  % % % Stats begins  % % %
  
  nsub=length(freqloall_tAudAlone_comb);
  nsub_nKD=length(freqloall_tAudAlone_nKD);
  nsub_nSD=length(freqloall_tAudAlone_nSD);
  nsub_nKD_nSD=length(freqloall_tAudAlone_nKD_nSD);
  if statsflag && nsub>1
    
    cfg=[];
    cfg.avgoverfreq='yes';
    cfg.frequency=[4 6.5];
    grindlo_tacPaud_comb1_theta=ft_selectdata(cfg,grindlo_tacPaud_comb1);
    grindlo_tacMSpN_comb1_theta=ft_selectdata(cfg,grindlo_tacMSpN_comb1);
    cfg.frequency=[8 12.5];
    grindlo_tacPaud_comb1_alpha=ft_selectdata(cfg,grindlo_tacPaud_comb1);
    grindlo_tacMSpN_comb1_alpha=ft_selectdata(cfg,grindlo_tacMSpN_comb1);
    cfg.frequency=[14 30];
    grindlo_tacPaud_comb1_beta=ft_selectdata(cfg,grindlo_tacPaud_comb1);
    grindlo_tacMSpN_comb1_beta=ft_selectdata(cfg,grindlo_tacMSpN_comb1);
    if sleep && ss==12
      cfg.frequency=[4 6.5];
      grindlo_tacPaud_nKD_theta=ft_selectdata(cfg,grindlo_tacPaud_nKD);
      grindlo_tacMSpN_nKD_theta=ft_selectdata(cfg,grindlo_tacMSpN_nKD);
      grindlo_tacPaud_nSD_theta=ft_selectdata(cfg,grindlo_tacPaud_nSD);
      grindlo_tacMSpN_nSD_theta=ft_selectdata(cfg,grindlo_tacMSpN_nSD);
      grindlo_tacPaud_nKD_nSD_theta=ft_selectdata(cfg,grindlo_tacPaud_nKD_nSD);
      grindlo_tacMSpN_nKD_nSD_theta=ft_selectdata(cfg,grindlo_tacMSpN_nKD_nSD);
      cfg.frequency=[8 12.5];
      grindlo_tacPaud_nKD_alpha=ft_selectdata(cfg,grindlo_tacPaud_nKD);
      grindlo_tacMSpN_nKD_alpha=ft_selectdata(cfg,grindlo_tacMSpN_nKD);
      grindlo_tacPaud_nSD_alpha=ft_selectdata(cfg,grindlo_tacPaud_nSD);
      grindlo_tacMSpN_nSD_alpha=ft_selectdata(cfg,grindlo_tacMSpN_nSD);
      grindlo_tacPaud_nKD_nSD_alpha=ft_selectdata(cfg,grindlo_tacPaud_nKD_nSD);
      grindlo_tacMSpN_nKD_nSD_alpha=ft_selectdata(cfg,grindlo_tacMSpN_nKD_nSD);
      cfg.frequency=[14 30];
      grindlo_tacPaud_nKD_beta=ft_selectdata(cfg,grindlo_tacPaud_nKD);
      grindlo_tacMSpN_nKD_beta=ft_selectdata(cfg,grindlo_tacMSpN_nKD);
      grindlo_tacPaud_nSD_beta=ft_selectdata(cfg,grindlo_tacPaud_nSD);
      grindlo_tacMSpN_nSD_beta=ft_selectdata(cfg,grindlo_tacMSpN_nSD);
      grindlo_tacPaud_nKD_nSD_beta=ft_selectdata(cfg,grindlo_tacPaud_nKD_nSD);
      grindlo_tacMSpN_nKD_nSD_beta=ft_selectdata(cfg,grindlo_tacMSpN_nKD_nSD);
    end
    if comb2flag
      cfg.frequency=[4 6.5];
      grindlo_tacPaud_comb2_theta=ft_selectdata(cfg,grindlo_tacPaud_comb2);
      grindlo_tacMSpN_comb2_theta=ft_selectdata(cfg,grindlo_tacMSpN_comb2);
      cfg.frequency=[8 12.5];
      grindlo_tacPaud_comb2_alpha=ft_selectdata(cfg,grindlo_tacPaud_comb2);
      grindlo_tacMSpN_comb2_alpha=ft_selectdata(cfg,grindlo_tacMSpN_comb2);
      cfg.frequency=[14 30];
      grindlo_tacPaud_comb2_beta=ft_selectdata(cfg,grindlo_tacPaud_comb2);
      grindlo_tacMSpN_comb2_beta=ft_selectdata(cfg,grindlo_tacMSpN_comb2);
    end
    
    
    
    cfg=[];
    %       cfg.latency=[tacbasemax(ll); .55-audbasemax(ll)];
    cfg.latency=[-.15-audbasemax(ll); 1.05-audbasemax(ll)];
    cfg.neighbours=neighbours;
    cfg.parameter='powspctrm';
    cfg.method='montecarlo';
    cfg.numrandomization=2000;
    % cfg.correctm='holm';
    cfg.correctm='cluster';
    cfg.clusteralpha = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan = 2;
    cfg.statistic='depsamplesT';
    cfg.design=set_cfg_design_depT(nsub);
    cfg.ivar=1;
    cfg.uvar=2;
    cfg.randomseed=mcseed;
    %     if usetr==2
    %       grindlo_TPA_MSPN_comb1_zeros=grindlo_TPA_MSPN_comb1;
    %       grindlo_TPA_MSPN_comb1_zeros.powspctrm=zeros(size(grindlo_TPA_MSPN_comb1.powspctrm));
    %       grindlo_TPA_MSPN_comb1_zeros.plvspctrm=zeros(size(grindlo_TPA_MSPN_comb1.plvspctrm));
    %       grindhi_TPA_MSPN_comb1_zeros=grindhi_TPA_MSPN_comb1;
    %       grindhi_TPA_MSPN_comb1_zeros.powspctrm=zeros(size(grindhi_TPA_MSPN_comb1.powspctrm));
    %       grindhi_TPA_MSPN_comb1_zeros.plvspctrm=zeros(size(grindhi_TPA_MSPN_comb1.plvspctrm));
    %       stattl_mc_comb1=ft_freqstatistics(cfg, grindlo_TPA_MSPN_comb1, grindlo_TPA_MSPN_comb1_zeros);
    %       statth_mc_comb1=ft_freqstatistics(cfg, grindhi_TPA_MSPN_comb1, grindhi_TPA_MSPN_comb1_zeros);
    %       if comb2flag
    %         grindlo_TPA_MSPN_comb2_zeros=grindlo_TPA_MSPN_comb2;
    %         grindlo_TPA_MSPN_comb2_zeros.powspctrm=zeros(size(grindlo_TPA_MSPN_comb2.powspctrm));
    %         grindlo_TPA_MSPN_comb2_zeros.plvspctrm=zeros(size(grindlo_TPA_MSPN_comb2.plvspctrm));
    %         grindhi_TPA_MSPN_comb2_zeros=grindhi_TPA_MSPN_comb2;
    %         grindhi_TPA_MSPN_comb2_zeros.powspctrm=zeros(size(grindhi_TPA_MSPN_comb2.powspctrm));
    %         grindhi_TPA_MSPN_comb2_zeros.plvspctrm=zeros(size(grindhi_TPA_MSPN_comb2.plvspctrm));
    %         stattl_mc_comb2=ft_freqstatistics(cfg, grindlo_TPA_MSPN_comb2, grindlo_TPA_MSPN_comb2_zeros);
    %         statth_mc_comb2=ft_freqstatistics(cfg, grindhi_TPA_MSPN_comb2, grindhi_TPA_MSPN_comb2_zeros);
    %       end
    %     else
    
    %       if fftaddflag
    %         stattl_mc_fftadd=ft_freqstatistics(cfg, grindlo_tacPaud_fftadd, grindlo_tacMSpN_fftadd);
    %         statth_mc_fftadd=ft_freqstatistics(cfg, grindhi_tacPaud_fftadd, grindhi_tacMSpN_fftadd);
    %       end
    stattl_mc_comb1=ft_freqstatistics(cfg, grindlo_tacPaud_comb1, grindlo_tacMSpN_comb1);
    statth_mc_comb1=ft_freqstatistics(cfg, grindhi_tacPaud_comb1, grindhi_tacMSpN_comb1);
    stattl_mc_comb1_theta=ft_freqstatistics(cfg, grindlo_tacPaud_comb1_theta, grindlo_tacMSpN_comb1_theta);
    stattl_mc_comb1_alpha=ft_freqstatistics(cfg, grindlo_tacPaud_comb1_alpha, grindlo_tacMSpN_comb1_alpha);
    stattl_mc_comb1_beta =ft_freqstatistics(cfg, grindlo_tacPaud_comb1_beta,  grindlo_tacMSpN_comb1_beta);
    if sleep && ss==12
      cfg.design=set_cfg_design_depT(nsub_nKD);
      stattl_mc_nKD=ft_freqstatistics(cfg, grindlo_tacPaud_nKD, grindlo_tacMSpN_nKD);
      statth_mc_nKD=ft_freqstatistics(cfg, grindhi_tacPaud_nKD, grindhi_tacMSpN_nKD);
      stattl_mc_nKD_theta=ft_freqstatistics(cfg, grindlo_tacPaud_nKD_theta, grindlo_tacMSpN_nKD_theta);
      stattl_mc_nKD_alpha=ft_freqstatistics(cfg, grindlo_tacPaud_nKD_alpha, grindlo_tacMSpN_nKD_alpha);
      stattl_mc_nKD_beta =ft_freqstatistics(cfg, grindlo_tacPaud_nKD_beta,  grindlo_tacMSpN_nKD_beta);
      cfg.design=set_cfg_design_depT(nsub_nSD);
      stattl_mc_nSD=ft_freqstatistics(cfg, grindlo_tacPaud_nSD, grindlo_tacMSpN_nSD);
      statth_mc_nSD=ft_freqstatistics(cfg, grindhi_tacPaud_nSD, grindhi_tacMSpN_nSD);
      stattl_mc_nSD_theta=ft_freqstatistics(cfg, grindlo_tacPaud_nSD_theta, grindlo_tacMSpN_nSD_theta);
      stattl_mc_nSD_alpha=ft_freqstatistics(cfg, grindlo_tacPaud_nSD_alpha, grindlo_tacMSpN_nSD_alpha);
      stattl_mc_nSD_beta =ft_freqstatistics(cfg, grindlo_tacPaud_nSD_beta,  grindlo_tacMSpN_nSD_beta);
      cfg.design=set_cfg_design_depT(nsub_nKD);
      stattl_mc_nKD_nSD=ft_freqstatistics(cfg, grindlo_tacPaud_nKD_nSD, grindlo_tacMSpN_nKD_nSD);
      statth_mc_nKD_nSD=ft_freqstatistics(cfg, grindhi_tacPaud_nKD_nSD, grindhi_tacMSpN_nKD_nSD);
      stattl_mc_nKD_nSD_theta=ft_freqstatistics(cfg, grindlo_tacPaud_nKD_nSD_theta, grindlo_tacMSpN_nKD_nSD_theta);
      stattl_mc_nKD_nSD_alpha=ft_freqstatistics(cfg, grindlo_tacPaud_nKD_nSD_alpha, grindlo_tacMSpN_nKD_nSD_alpha);
      stattl_mc_nKD_nSD_beta =ft_freqstatistics(cfg, grindlo_tacPaud_nKD_nSD_beta,  grindlo_tacMSpN_nKD_nSD_beta);
      cfg.design=set_cfg_design_depT(nsub);
    end
    if comb2flag
      stattl_mc_comb2=ft_freqstatistics(cfg, grindlo_tacPaud_comb2, grindlo_tacMSpN_comb2);
      statth_mc_comb2=ft_freqstatistics(cfg, grindhi_tacPaud_comb2, grindhi_tacMSpN_comb2);
      stattl_mc_comb2_theta=ft_freqstatistics(cfg, grindlo_tacPaud_comb2_theta, grindlo_tacMSpN_comb2_theta);
      stattl_mc_comb2_alpha=ft_freqstatistics(cfg, grindlo_tacPaud_comb2_alpha, grindlo_tacMSpN_comb2_alpha);
      stattl_mc_comb2_beta =ft_freqstatistics(cfg, grindlo_tacPaud_comb2_beta,  grindlo_tacMSpN_comb2_beta);
    end    
    %       stattl_mc_fftaddbn=ft_freqstatistics(cfg, grindlo_tacPaud_fftaddbn, grindlo_tacMSpN_fftaddbn);
    %       statth_mc_fftaddbn=ft_freqstatistics(cfg, grindhi_tacPaud_fftaddbn, grindhi_tacMSpN_fftaddbn);
    %       stattl_mc_comb1bn=ft_freqstatistics(cfg, grindlo_tacPaud_comb1bn, grindlo_tacMSpN_comb1bn);
    %       statth_mc_comb1bn=ft_freqstatistics(cfg, grindhi_tacPaud_comb1bn, grindhi_tacMSpN_comb1bn);
    %       if synchasynch && ll<5
    %         cfg.latency=[tacbasemax(9); .55-audbasemax(9)];
    %         stattl_synch_comb1=ft_freqstatistics(cfg, grindlo_tMSsynch_comb1, grindlo_tMSasynch_comb1);
    %         statth_synch_comb1=ft_freqstatistics(cfg, grindhi_tMSsynch_comb1, grindhi_tMSasynch_comb1);
    %         stattl_synch_comb2=ft_freqstatistics(cfg, grindlo_tMSsynch_comb2, grindlo_tMSasynch_comb2);
    %         statth_synch_comb2=ft_freqstatistics(cfg, grindhi_tMSsynch_comb2, grindhi_tMSasynch_comb2);
    %       end
    %
    %       if audtacflag
    %         cfg.latency=[audbasemax(ll); .55-tacbasemax(ll)]; % ???
    %         statal_mc_fftadd=ft_freqstatistics(cfg, grindlo_audPtac_fftadd, grindlo_audMSpN_fftadd);
    %         statah_mc_fftadd=ft_freqstatistics(cfg, grindhi_audPtac_fftadd, grindhi_audMSpN_fftadd);
    %         statal_mc_comb1=ft_freqstatistics(cfg, grindlo_audPtac_comb1, grindlo_audMSpN_comb1);
    %         statah_mc_comb1=ft_freqstatistics(cfg, grindhi_audPtac_comb1, grindhi_audMSpN_comb1);
    %         if comb2flag
    %           statal_mc_comb2=ft_freqstatistics(cfg, grindlo_audPtac_comb2, grindlo_audMSpN_comb2);
    %           statah_mc_comb2=ft_freqstatistics(cfg, grindhi_audPtac_comb2, grindhi_audMSpN_comb2);
    %         end
    %         %       statal_mc_fftaddbn=ft_freqstatistics(cfg, grindlo_audPtac_fftaddbn, grindlo_audMSpN_fftaddbn);
    %         %       statah_mc_fftaddbn=ft_freqstatistics(cfg, grindhi_audPtac_fftaddbn, grindhi_audMSpN_fftaddbn);
    %         %       statal_mc_comb1bn=ft_freqstatistics(cfg, grindlo_audPtac_comb1bn, grindlo_audMSpN_comb1bn);
    %         %       statah_mc_comb1bn=ft_freqstatistics(cfg, grindhi_audPtac_comb1bn, grindhi_audMSpN_comb1bn);
    %       end
    %     end % usetr
    
    %       cfg.latency=[tacbasemax(ll); .55-audbasemax(ll)];
    cfg.latency=[-.15-audbasemax(ll); 1.05-audbasemax(ll)];
    cfg.parameter='plvspctrm';
    %       cfg.clusterthreshold = 'nonparametric_common'; % or 'nonparametric_individual' see clusterstat.m
    %       stattl_mc_comb1plvdac=ft_freqstatistics(cfg, grindlo_tacPaud_comb1, grindlo_tacMSpN_comb1);
    %       statth_mc_comb1plvdac=ft_freqstatistics(cfg, grindhi_tacPaud_comb1, grindhi_tacMSpN_comb1);
    %       if comb2flag
    %         stattl_mc_comb2plvdac=ft_freqstatistics(cfg, grindlo_tacPaud_comb2, grindlo_tacMSpN_comb2);
    %         statth_mc_comb2plvdac=ft_freqstatistics(cfg, grindhi_tacPaud_comb2, grindhi_tacMSpN_comb2);
    %       end
    cfg.clusterthreshold = 'nonparametric_individual'; % or 'nonparametric_individual' see clusterstat.m
    %     if usetr==2
    %       cfg.statistic='diff_itc';
    %       cfg.complex='diffabs'; % default; not sensitive to phase differences between conditions
    %       stattl_mc_comb1plvdai=ft_freqstatistics(cfg, grindlo_TPA_MSPN_comb1, grindlo_TPA_MSPN_comb1_zeros);
    %       statth_mc_comb1plvdai=ft_freqstatistics(cfg, grindhi_TPA_MSPN_comb1, grindhi_TPA_MSPN_comb1_zeros);
    %       if comb2flag
    %         stattl_mc_comb2plvdai=ft_freqstatistics(cfg, grindlo_TPA_MSPN_comb2, grindlo_TPA_MSPN_comb2_zeros);
    %         statth_mc_comb2plvdai=ft_freqstatistics(cfg, grindhi_TPA_MSPN_comb2, grindhi_TPA_MSPN_comb2_zeros);
    %       end
    %     else
    if itcdepsampflag
      cfg.statistic='depsamplesT';
      cfg.design=set_cfg_design_depT(nsub);
      %         grindlo_tacPaud_comb1_abs=grindlo_tacPaud_comb1;
      %         grindlo_tacPaud_comb1_abs.plvspctrm=abs(grindlo_tacPaud_comb1_abs.plvspctrm);
      %         grindlo_tacMSpN_comb1_abs=grindlo_tacMSpN_comb1;
      %         grindlo_tacMSpN_comb1_abs.plvspctrm=abs(grindlo_tacMSpN_comb1_abs.plvspctrm);
      %         grindhi_tacPaud_comb1_abs=grindhi_tacPaud_comb1;
      %         grindhi_tacPaud_comb1_abs.plvspctrm=abs(grindhi_tacPaud_comb1_abs.plvspctrm);
      %         grindhi_tacMSpN_comb1_abs=grindhi_tacMSpN_comb1;
      %         grindhi_tacMSpN_comb1_abs.plvspctrm=abs(grindhi_tacMSpN_comb1_abs.plvspctrm);
      %         stattl_mc_comb1plvDepT=ft_freqstatistics(cfg, grindlo_tacPaud_comb1_abs, grindlo_tacMSpN_comb1_abs);
      %         statth_mc_comb1plvDepT=ft_freqstatistics(cfg, grindhi_tacPaud_comb1_abs, grindhi_tacMSpN_comb1_abs);
      
      cfg.parameter='plvabs';
      stattl_mc_comb1plvabsDepT=ft_freqstatistics(cfg, grindlo_tacPaud_comb1, grindlo_tacMSpN_comb1);
      statth_mc_comb1plvabsDepT=ft_freqstatistics(cfg, grindhi_tacPaud_comb1, grindhi_tacMSpN_comb1);
      stattl_mc_comb1plvabsDepT_theta=ft_freqstatistics(cfg, grindlo_tacPaud_comb1_theta, grindlo_tacMSpN_comb1_theta);
      stattl_mc_comb1plvabsDepT_alpha=ft_freqstatistics(cfg, grindlo_tacPaud_comb1_alpha, grindlo_tacMSpN_comb1_alpha);
      stattl_mc_comb1plvabsDepT_beta =ft_freqstatistics(cfg, grindlo_tacPaud_comb1_beta,  grindlo_tacMSpN_comb1_beta);
      if sleep && ss==12
        cfg.design=set_cfg_design_depT(nsub_nKD);
        stattl_mc_nKDplvabsDepT=ft_freqstatistics(cfg, grindlo_tacPaud_nKD, grindlo_tacMSpN_nKD);
        statth_mc_nKDplvabsDepT=ft_freqstatistics(cfg, grindhi_tacPaud_nKD, grindhi_tacMSpN_nKD);
        stattl_mc_nKDplvabsDepT_theta=ft_freqstatistics(cfg, grindlo_tacPaud_nKD_theta, grindlo_tacMSpN_nKD_theta);
        stattl_mc_nKDplvabsDepT_alpha=ft_freqstatistics(cfg, grindlo_tacPaud_nKD_alpha, grindlo_tacMSpN_nKD_alpha);
        stattl_mc_nKDplvabsDepT_beta =ft_freqstatistics(cfg, grindlo_tacPaud_nKD_beta,  grindlo_tacMSpN_nKD_beta);
        cfg.design=set_cfg_design_depT(nsub_nSD);
        stattl_mc_nSDplvabsDepT=ft_freqstatistics(cfg, grindlo_tacPaud_nSD, grindlo_tacMSpN_nSD);
        statth_mc_nSDplvabsDepT=ft_freqstatistics(cfg, grindhi_tacPaud_nSD, grindhi_tacMSpN_nSD);
        stattl_mc_nSDplvabsDepT_theta=ft_freqstatistics(cfg, grindlo_tacPaud_nSD_theta, grindlo_tacMSpN_nSD_theta);
        stattl_mc_nSDplvabsDepT_alpha=ft_freqstatistics(cfg, grindlo_tacPaud_nSD_alpha, grindlo_tacMSpN_nSD_alpha);
        stattl_mc_nSDplvabsDepT_beta =ft_freqstatistics(cfg, grindlo_tacPaud_nSD_beta,  grindlo_tacMSpN_nSD_beta);
        cfg.design=set_cfg_design_depT(nsub_nKD_nSD);
        stattl_mc_nKD_nSDplvabsDepT=ft_freqstatistics(cfg, grindlo_tacPaud_nKD_nSD, grindlo_tacMSpN_nKD_nSD);
        statth_mc_nKD_nSDplvabsDepT=ft_freqstatistics(cfg, grindhi_tacPaud_nKD_nSD, grindhi_tacMSpN_nKD_nSD);
        stattl_mc_nKD_nSDplvabsDepT_theta=ft_freqstatistics(cfg, grindlo_tacPaud_nKD_nSD_theta, grindlo_tacMSpN_nKD_nSD_theta);
        stattl_mc_nKD_nSDplvabsDepT_alpha=ft_freqstatistics(cfg, grindlo_tacPaud_nKD_nSD_alpha, grindlo_tacMSpN_nKD_nSD_alpha);
        stattl_mc_nKD_nSDplvabsDepT_beta =ft_freqstatistics(cfg, grindlo_tacPaud_nKD_nSD_beta,  grindlo_tacMSpN_nKD_nSD_beta);
      end
      
      if comb2flag
        %           grindlo_tacPaud_comb2_abs=grindlo_tacPaud_comb2;
        %           grindlo_tacPaud_comb2_abs.plvspctrm=abs(grindlo_tacPaud_comb2_abs.plvspctrm);
        %           grindlo_tacMSpN_comb2_abs=grindlo_tacMSpN_comb2;
        %           grindlo_tacMSpN_comb2_abs.plvspctrm=abs(grindlo_tacMSpN_comb2_abs.plvspctrm);
        %           grindhi_tacPaud_comb2_abs=grindhi_tacPaud_comb2;
        %           grindhi_tacPaud_comb2_abs.plvspctrm=abs(grindhi_tacPaud_comb2_abs.plvspctrm);
        %           grindhi_tacMSpN_comb2_abs=grindhi_tacMSpN_comb2;
        %           grindhi_tacMSpN_comb2_abs.plvspctrm=abs(grindhi_tacMSpN_comb2_abs.plvspctrm);
        %           stattl_mc_comb2plvDepT=ft_freqstatistics(cfg, grindlo_tacPaud_comb2_abs, grindlo_tacMSpN_comb2_abs);
        %           statth_mc_comb2plvDepT=ft_freqstatistics(cfg, grindhi_tacPaud_comb2_abs, grindhi_tacMSpN_comb2_abs);
        cfg.parameter='plvabs';
        stattl_mc_comb2plvabsDepT=ft_freqstatistics(cfg, grindlo_tacPaud_comb2, grindlo_tacMSpN_comb2);
        statth_mc_comb2plvabsDepT=ft_freqstatistics(cfg, grindhi_tacPaud_comb2, grindhi_tacMSpN_comb2);
        stattl_mc_comb2plvabsDepT_theta=ft_freqstatistics(cfg, grindlo_tacPaud_comb2_theta, grindlo_tacMSpN_comb2_theta);
        stattl_mc_comb2plvabsDepT_alpha=ft_freqstatistics(cfg, grindlo_tacPaud_comb2_alpha, grindlo_tacMSpN_comb2_alpha);
        stattl_mc_comb2plvabsDepT_beta =ft_freqstatistics(cfg, grindlo_tacPaud_comb2_beta,  grindlo_tacMSpN_comb2_beta);
      end
    else
      cfg.statistic='diff_itc';
      cfg.complex='diffabs'; % default; not sensitive to phase differences between conditions
      
      
      stattl_mc_comb1plvdai=ft_freqstatistics(cfg, grindlo_tacPaud_comb1, grindlo_tacMSpN_comb1);
      statth_mc_comb1plvdai=ft_freqstatistics(cfg, grindhi_tacPaud_comb1, grindhi_tacMSpN_comb1);
      if comb2flag
        stattl_mc_comb2plvdai=ft_freqstatistics(cfg, grindlo_tacPaud_comb2, grindlo_tacMSpN_comb2);
        statth_mc_comb2plvdai=ft_freqstatistics(cfg, grindhi_tacPaud_comb2, grindhi_tacMSpN_comb2);
      end
    end
    %       if synchasynch && ll<5
    %         cfg.latency=[tacbasemax(9); .55-audbasemax(9)];
    %         stattl_synch_comb1plvdai=ft_freqstatistics(cfg, grindlo_tMSsynch_comb1, grindlo_tMSasynch_comb1);
    %         statth_synch_comb1plvdai=ft_freqstatistics(cfg, grindhi_tMSsynch_comb1, grindhi_tMSasynch_comb1);
    %         stattl_synch_comb2plvdai=ft_freqstatistics(cfg, grindlo_tMSsynch_comb2, grindlo_tMSasynch_comb2);
    %         statth_synch_comb2plvdai=ft_freqstatistics(cfg, grindhi_tMSsynch_comb2, grindhi_tMSasynch_comb2);
    %       end
    %       if audtacflag
    %         cfg.latency=[audbasemax(ll); .55-tacbasemax(ll)];
    %         %         cfg.clusterthreshold = 'nonparametric_common'; % or 'nonparametric_individual' see clusterstat.m
    %         %         statal_mc_comb1plvdac=ft_freqstatistics(cfg, grindlo_audPtac_comb1, grindlo_audMSpN_comb1);
    %         %         statah_mc_comb1plvdac=ft_freqstatistics(cfg, grindhi_audPtac_comb1, grindhi_audMSpN_comb1);
    %         %         if comb2flag
    %         %           statal_mc_comb2plvdac=ft_freqstatistics(cfg, grindlo_audPtac_comb2, grindlo_audMSpN_comb2);
    %         %           statah_mc_comb2plvdac=ft_freqstatistics(cfg, grindhi_audPtac_comb2, grindhi_audMSpN_comb2);
    %         %         end
    %         cfg.clusterthreshold = 'nonparametric_individual'; % or 'nonparametric_individual' see clusterstat.m
    %         statal_mc_comb1plvdai=ft_freqstatistics(cfg, grindlo_audPtac_comb1, grindlo_audMSpN_comb1);
    %         statah_mc_comb1plvdai=ft_freqstatistics(cfg, grindhi_audPtac_comb1, grindhi_audMSpN_comb1);
    %         if comb2flag
    %           statal_mc_comb2plvdai=ft_freqstatistics(cfg, grindlo_audPtac_comb2, grindlo_audMSpN_comb2);
    %           statah_mc_comb2plvdai=ft_freqstatistics(cfg, grindhi_audPtac_comb2, grindhi_audMSpN_comb2);
    %         end
    %       end
    
    %       if 0  % absdiff
    %         %       cfg.latency=[tacbasemax(ll); .55-audbasemax(ll)];
    %         cfg.latency=[-.15-audbasemax(ll); 1.05-audbasemax(ll)];
    %         cfg.complex='absdiff'; % sensitive to phase difference between conditions
    %         %       cfg.clusterthreshold = 'nonparametric_common'; % or 'nonparametric_individual' see clusterstat.m
    %         %       stattl_mc_comb1plvadc=ft_freqstatistics(cfg, grindlo_tacPaud_comb1, grindlo_tacMSpN_comb1);
    %         %       statth_mc_comb1plvadc=ft_freqstatistics(cfg, grindhi_tacPaud_comb1, grindhi_tacMSpN_comb1);
    %         %       if comb2flag
    %         %         stattl_mc_comb2plvadc=ft_freqstatistics(cfg, grindlo_tacPaud_comb2, grindlo_tacMSpN_comb2);
    %         %         statth_mc_comb2plvadc=ft_freqstatistics(cfg, grindhi_tacPaud_comb2, grindhi_tacMSpN_comb2);
    %         %       end
    %         cfg.clusterthreshold = 'nonparametric_individual'; % or 'nonparametric_individual' see clusterstat.m
    %         stattl_mc_comb1plvadi=ft_freqstatistics(cfg, grindlo_tacPaud_comb1, grindlo_tacMSpN_comb1);
    %         statth_mc_comb1plvadi=ft_freqstatistics(cfg, grindhi_tacPaud_comb1, grindhi_tacMSpN_comb1);
    %         if comb2flag
    %           stattl_mc_comb2plvadi=ft_freqstatistics(cfg, grindlo_tacPaud_comb2, grindlo_tacMSpN_comb2);
    %           statth_mc_comb2plvadi=ft_freqstatistics(cfg, grindhi_tacPaud_comb2, grindhi_tacMSpN_comb2);
    %         end
    %         if synchasynch && ll<5
    %           cfg.latency=[tacbasemax(9); .55-audbasemax(9)];
    %           stattl_synch_comb1plvadi=ft_freqstatistics(cfg, grindlo_tMSsynch_comb1, grindlo_tMSasynch_comb1);
    %           statth_synch_comb1plvadi=ft_freqstatistics(cfg, grindhi_tMSsynch_comb1, grindhi_tMSasynch_comb1);
    %           stattl_synch_comb2plvadi=ft_freqstatistics(cfg, grindlo_tMSsynch_comb2, grindlo_tMSasynch_comb2);
    %           statth_synch_comb2plvadi=ft_freqstatistics(cfg, grindhi_tMSsynch_comb2, grindhi_tMSasynch_comb2);
    %         end
    %
    %         if audtacflag
    %           cfg.latency=[audbasemax(ll); .55-tacbasemax(ll)];
    %           %         cfg.clusterthreshold = 'nonparametric_common'; % or 'nonparametric_individual' see clusterstat.m
    %           %         statal_mc_comb1plvadc=ft_freqstatistics(cfg, grindlo_audPtac_comb1, grindlo_audMSpN_comb1);
    %           %         statah_mc_comb1plvadc=ft_freqstatistics(cfg, grindhi_audPtac_comb1, grindhi_audMSpN_comb1);
    %           %         if comb2flag
    %           %           statal_mc_comb2plvadc=ft_freqstatistics(cfg, grindlo_audPtac_comb2, grindlo_audMSpN_comb2);
    %           %           statah_mc_comb2plvadc=ft_freqstatistics(cfg, grindhi_audPtac_comb2, grindhi_audMSpN_comb2);
    %           %         end
    %           cfg.clusterthreshold = 'nonparametric_individual'; % or 'nonparametric_individual' see clusterstat.m
    %           statal_mc_comb1plvadi=ft_freqstatistics(cfg, grindlo_audPtac_comb1, grindlo_audMSpN_comb1);
    %           statah_mc_comb1plvadi=ft_freqstatistics(cfg, grindhi_audPtac_comb1, grindhi_audMSpN_comb1);
    %           if comb2flag
    %             statal_mc_comb2plvadi=ft_freqstatistics(cfg, grindlo_audPtac_comb2, grindlo_audMSpN_comb2);
    %             statah_mc_comb2plvadi=ft_freqstatistics(cfg, grindhi_audPtac_comb2, grindhi_audMSpN_comb2);
    %           end
    %         end
    %       end
    
    %     end % usetr
    
    
  end % statsflag
  
  if plotflag
    plvabsdiffmaxz=plvabsmaxz/2;
    timwin=[]; % not needed to narrow down, as stat will have already been preselected.
    base=[]; % stats weren't baseline corrected, so viewing of data shouldn't here be either.
    for yy=1:size(ylimlo,1)
      ylim=ylimlo(yy,:);
      
      try
        fig=240+yy;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb1,timwin,ylim,base,stattl_mc_comb1,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb1_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=840+yy;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb1,timwin,ylim,base,stattl_mc_comb1plvdai,.05,'plvavgabs',[0 plvabsdiffmaxz(yy)]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb1plvabsdai_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
        fig=1840+yy;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb1,timwin,ylim,base,stattl_mc_comb1plvdai,.05,'plvavgang',[-pi pi]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb1plvangdai_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=940+yy;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb1,timwin,ylim,base,stattl_mc_comb1plvadi,.05,'plvavgabs',[0 plvabsdiffmaxz(yy)]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb1plvabsadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
        fig=1940+yy;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb1,timwin,ylim,base,stattl_mc_comb1plvadi,.05,'plvavgang',[-pi pi]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb1plvangadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=340+yy;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb2,timwin,ylim,base,stattl_mc_comb2,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb2_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=2840+yy;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb2,timwin,ylim,base,stattl_mc_comb2plvdai,.05,'plvavgabs',[0 plvabsdiffmaxz(yy)]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb2plvabsdai_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
        fig=3840+yy;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb2,timwin,ylim,base,stattl_mc_comb2plvdai,.05,'plvavgang',[-pi pi]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb2plvangdai_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=2940+yy;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb2,timwin,ylim,base,stattl_mc_comb2plvadi,.05,'plvavgabs',[0 plvabsdiffmaxz(yy)]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb2plvabsadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
        fig=3940+yy;
        topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb2,timwin,ylim,base,stattl_mc_comb2plvadi,.05,'plvavgang',[-pi pi]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb2plvangadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=1340+yy;
        topoplotTFR_highlight(fig,gravelo_TMSs_TMSa_comb1,timwin,ylim,base,stattl_synch_comb1,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_synch_ta_comb1_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=2340+yy;
        topoplotTFR_highlight(fig,gravelo_TMSs_TMSa_comb2,timwin,ylim,base,stattl_synch_comb2,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
        if printflag
          print(fig,[fdir 'grdifflo_topoOverTime_synch_ta_comb2_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      
      
      if 0 % obsolete
        try
          fig=140+yy;
          topoplotTFR_highlight(fig,gravelo_TPA_MSPN_fftadd,timwin,ylim,base,stattl_mc_fftadd,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_fftadd_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=640+yy;
          topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb1,timwin,ylim,base,stattl_mc_comb1plvdac,.05,'plvavgabs',[0 plvabsdiffmaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb1plvdac_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=740+yy;
          topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb1,timwin,ylim,base,stattl_mc_comb1plvadc,.05,'plvavgabs',[0 plvabsdiffmaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb1plvabsadc_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
          fig=1740+yy;
          topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb1,timwin,ylim,base,stattl_mc_comb1plvadc,.05,'plvavgang',[-pi pi]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb1plvangadc_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=440+yy;
          topoplotTFR_highlight(fig,gravelo_TPA_MSPN_fftaddbn,timwin,ylim,base,stattl_mc_fftaddbn,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_fftaddbn_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=540+yy;
          topoplotTFR_highlight(fig,gravelo_TPA_MSPN_comb1bn,timwin,ylim,base,stattl_mc_comb1bn,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_ta_comb1bn_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=160+yy;
          topoplotTFR_highlight(fig,gravelo_APT_MSPN_fftadd,timwin,ylim,base,statal_mc_fftadd,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_at_fftadd_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=260+yy;
          topoplotTFR_highlight(fig,gravelo_APT_MSPN_comb1,timwin,ylim,base,statal_mc_comb1,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_at_comb1_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=660+yy;
          topoplotTFR_highlight(fig,gravelo_APT_MSPN_comb1,timwin,ylim,base,statal_mc_comb1plvdac,.05,'plvavgabs',[0 plvabsdiffmaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_at_comb1plvdac_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=760+yy;
          topoplotTFR_highlight(fig,gravelo_APT_MSPN_comb1,timwin,ylim,base,statal_mc_comb1plvadc,.05,'plvavgabs',[0 plvabsdiffmaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_at_comb1plvabsadc_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
          fig=1760+yy;
          topoplotTFR_highlight(fig,gravelo_APT_MSPN_comb1,timwin,ylim,base,statal_mc_comb1plvadc,.05,'plvavgang',[-pi pi]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_at_comb1plvangadc_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=860+yy;
          topoplotTFR_highlight(fig,gravelo_APT_MSPN_comb1,timwin,ylim,base,statal_mc_comb1plvdai,.05,'plvavgabs',[0 plvabsdiffmaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_at_comb1plvdai_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=960+yy;
          topoplotTFR_highlight(fig,gravelo_APT_MSPN_comb1,timwin,ylim,base,statal_mc_comb1plvadi,.05,'plvavgabs',[0 plvabsdiffmaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_at_comb1plvabsadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
          fig=1960+yy;
          topoplotTFR_highlight(fig,gravelo_APT_MSPN_comb1,timwin,ylim,base,statal_mc_comb1plvadi,.05,'plvavgang',[-pi pi]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_at_comb1plvangadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=360+yy;
          topoplotTFR_highlight(fig,gravelo_APT_MSPN_comb2,timwin,ylim,base,statal_mc_comb2,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_at_comb2_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=460+yy;
          topoplotTFR_highlight(fig,gravelo_APT_MSPN_fftaddbn,timwin,ylim,base,statal_mc_fftaddbn,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_at_fftaddbn_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=560+yy;
          topoplotTFR_highlight(fig,gravelo_APT_MSPN_comb1bn,timwin,ylim,base,statal_mc_comb1bn,.05,'powspctrm',[-powalonemaxz(yy) powalonemaxz(yy)]);
          if printflag
            print(fig,[fdir 'grdifflo_topoOverTime_mc_at_comb1bn_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
      end
      
    end
    
    for yy=1:size(ylimhi,1)
      ylim=ylimhi(yy,:);
      
      try
        fig=250+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb1,timwin,ylim,base,statth_mc_comb1,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb1_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=850+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb1,timwin,ylim,base,statth_mc_comb1plvdai,.05,'plvavgabs',[0 plvabsdiffmaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb1plvabsdai_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
        fig=1850+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb1,timwin,ylim,base,statth_mc_comb1plvdai,.05,'plvavgang',[-pi pi]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb1plvangdai_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=950+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb1,timwin,ylim,base,statth_mc_comb1plvadi,.05,'plvavgabs',[0 plvabsdiffmaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb1plvabsadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
        fig=1950+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb1,timwin,ylim,base,statth_mc_comb1plvadi,.05,'plvavgang',[-pi pi]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb1plvangadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=350+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb2,timwin,ylim,base,statth_mc_comb2,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb2_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=2850+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb2,timwin,ylim,base,statth_mc_comb2plvdai,.05,'plvavgabs',[0 plvabsdiffmaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb2plvabsdai_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
        fig=3850+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb2,timwin,ylim,base,statth_mc_comb2plvdai,.05,'plvavgang',[-pi pi]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb2plvangdai_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=2950+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb2,timwin,ylim,base,statth_mc_comb2plvadi,.05,'plvavgabs',[0 plvabsdiffmaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb2plvabsadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
        fig=3950+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb2,timwin,ylim,base,statth_mc_comb2plvadi,.05,'plvavgang',[-pi pi]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb2plvangadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=1350+yy;
        topoplotTFR_highlight(fig,gravehi_TMSs_TMSa_comb1,timwin,ylim,base,statth_synch_comb1,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_synch_ta_comb1_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=2350+yy;
        topoplotTFR_highlight(fig,gravehi_TMSs_TMSa_comb2,timwin,ylim,base,statth_synch_comb2,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_synch_ta_comb2_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      
      try
        fig=150+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_fftadd,timwin,ylim,base,statth_mc_fftadd,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_fftadd_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=650+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb1,timwin,ylim,base,statth_mc_comb1plvdac,.05,'plvavgabs',[0 plvabsdiffmaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb1plvdac_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=750+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb1,timwin,ylim,base,statth_mc_comb1plvadc,.05,'plvavgabs',[0 plvabsdiffmaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb1plvabsadc_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
        fig=1750+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb1,timwin,ylim,base,statth_mc_comb1plvadc,.05,'plvavgang',[-pi pi]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb1plvangadc_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=450+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_fftaddbn,timwin,ylim,base,statth_mc_fftaddbn,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_fftaddbn_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      try
        fig=550+yy;
        topoplotTFR_highlight(fig,gravehi_TPA_MSPN_comb1bn,timwin,ylim,base,statth_mc_comb1bn,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
        if printflag
          print(fig,[fdir 'grdiffhi_topoOverTime_mc_ta_comb1bn_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
        end
      catch
      end
      
      if 0 % obsolete
        try
          fig=170+yy;
          topoplotTFR_highlight(fig,gravehi_APT_MSPN_fftadd,timwin,ylim,base,statah_mc_fftadd,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
          if printflag
            print(fig,[fdir 'grdiffhi_topoOverTime_mc_at_fftadd_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=270+yy;
          topoplotTFR_highlight(fig,gravehi_APT_MSPN_comb1,timwin,ylim,base,statah_mc_comb1,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
          if printflag
            print(fig,[fdir 'grdiffhi_topoOverTime_mc_at_comb1_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=670+yy;
          topoplotTFR_highlight(fig,gravehi_APT_MSPN_comb1,timwin,ylim,base,statah_mc_comb1plvdac,.05,'plvavgabs',[0 plvabsdiffmaxz(yy+3)]);
          if printflag
            print(fig,[fdir 'grdiffhi_topoOverTime_mc_at_comb1plvdac_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=770+yy;
          topoplotTFR_highlight(fig,gravehi_APT_MSPN_comb1,timwin,ylim,base,statah_mc_comb1plvadc,.05,'plvavgabs',[0 plvabsdiffmaxz(yy+3)]);
          if printflag
            print(fig,[fdir 'grdiffhi_topoOverTime_mc_at_comb1plvabsadc_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
          fig=1770+yy;
          topoplotTFR_highlight(fig,gravehi_APT_MSPN_comb1,timwin,ylim,base,statah_mc_comb1plvadc,.05,'plvavgang',[-pi pi]);
          if printflag
            print(fig,[fdir 'grdiffhi_topoOverTime_mc_at_comb1plvangadc_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=870+yy;
          topoplotTFR_highlight(fig,gravehi_APT_MSPN_comb1,timwin,ylim,base,statah_mc_comb1plvdai,.05,'plvavgabs',[0 plvabsdiffmaxz(yy+3)]);
          if printflag
            print(fig,[fdir 'grdiffhi_topoOverTime_mc_at_comb1plvdai_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=970+yy;
          topoplotTFR_highlight(fig,gravehi_APT_MSPN_comb1,timwin,ylim,base,statah_mc_comb1plvadi,.05,'plvavgabs',[0 plvabsdiffmaxz(yy+3)]);
          if printflag
            print(fig,[fdir 'grdiffhi_topoOverTime_mc_at_comb1plvabsadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
          fig=1970+yy;
          topoplotTFR_highlight(fig,gravehi_APT_MSPN_comb1,timwin,ylim,base,statah_mc_comb1plvadi,.05,'plvavgang',[-pi pi]);
          if printflag
            print(fig,[fdir 'grdiffhi_topoOverTime_mc_at_comb1plvangadi_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=370+yy;
          topoplotTFR_highlight(fig,gravehi_APT_MSPN_comb2,timwin,ylim,base,statah_mc_comb2,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
          if printflag
            print(fig,[fdir 'grdiffhi_topoOverTime_mc_at_comb2_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=470+yy;
          topoplotTFR_highlight(fig,gravehi_APT_MSPN_fftaddbn,timwin,ylim,base,statah_mc_fftaddbn,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
          if printflag
            print(fig,[fdir 'grdiffhi_topoOverTime_mc_at_fftaddbn_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
        try
          fig=570+yy;
          topoplotTFR_highlight(fig,gravehi_APT_MSPN_comb1bn,timwin,ylim,base,statah_mc_comb1bn,.05,'powspctrm',[-powalonemaxz(yy+3) powalonemaxz(yy+3)]);
          if printflag
            print(fig,[fdir 'grdiffhi_topoOverTime_mc_at_comb1bn_cond' num2str(ll) num2str(tt) num2str(ss) num2str(ylim(1)) num2str(ylim(2)) '.png'],'-dpng')
          end
        catch
        end
      end
      
    end % yy
  end % plotflag
  
  %   save([edir 'statsgrave_TFR_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '.mat'],'stat*','grave*');
  %   save([edir 'statsgrave_TFR_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '_mcseed' num2str(mcseed) '.mat'],'stat*','grave*');
  %  save([edir 'statsgrave_TFR_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '_iter' num2str(iter) '_usetr' num2str(usetr) '_mcseed' num2str(mcseed) '.mat'],'stat*','grave*');
  if itcdepsampflag
    save([edir 'statsgrave_TFR_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '_usetr' num2str(usetr) '_mcseed' num2str(mcseed) '_itc' '.mat'],'stat*','grave*');
  else
    save([edir 'statsgrave_TFR_cond' num2str(ll) num2str(tt) num2str(ss) num2str(sleep) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '_usetr' num2str(usetr) '_mcseed' num2str(mcseed) '.mat'],'stat*','grave*');
  end
  
  clear stat*mc* gr*
  
  
  
  
  
  
end % ll
%   end % tt
% end % sleep


