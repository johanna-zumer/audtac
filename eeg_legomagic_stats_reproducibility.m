eeg_legomagic_preamble

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

