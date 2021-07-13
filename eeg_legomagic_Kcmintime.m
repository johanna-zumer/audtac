function eeg_legomagic_Kcmintime
% eeg_legomagic_preamble  % often called at beginning to set paths etc.

nsub=19;
load([edir 'Kcmintime.mat']);
anova1(mindATmean([1 3 4 5 6 7 9],:)'); 
anova1([mindATmean([1 3 4 5 6 7 9],:); mindAmean(5,:); mindTmean(5,:)]'); 
meanMinInd1=mean(mean(mindATmean([1 3 4 5 6 7 9],:)')); % %9.7217
meanMinInd2=mean(mean([mindATmean([1 3 4 5 6 7 9],:); mindAmean(5,:); mindTmean(5,:)]')); % 9.667
% conclusion: okay to use same time point for all;
anova1(mindATmean([4 5 6 7 9],:)'); 
anova1([mindATmean([4 5 6 7 9],:); mindAmean(5,:);]'); 
anova1([mindATmean([1 3 ],:);  mindTmean(5,:)]'); 

meanMinInd3=mean(mean([mindATmean([4 5 6 7 9],:); mindAmean(5,:);]')); % 9.667
meanMinInd4=mean(mean([mindATmean([1 3 ],:);  mindTmean(5,:)]')); % 9.667

meanMinTime=kcdist{5,1}.AT.time(round(meanMinInd1));  % interpreted as post 2nd stim

MinIndUse=round(meanMinInd1);
MinIndUse=7;  % literature 550ms

%% Since ANOVA shows difference across conditions, can we use common time of ANOVA across conditions?
load ANOVA_sleep1_ERP.mat
figure;plot(mean(stat_TPA_MSPN_1wayANOVA_long.mask,1))
for cc=1:60,ind{cc}=find(stat_TPA_MSPN_1wayANOVA_long.mask(cc,:));end 
stat_TPA_MSPN_1wayANOVA_long.time(median([ind{:}]))  % gives 684 ms, rounded to 675 (nearest 25ms)
MinIndUse=12;  



%%  
% % not valid, as values in KcDist are not Gaussian.
% 
% for ll=[1 3 4 5 6 7 9]
%   for ii=1:19
%     meandiff(ll,ii)=mean(kcdist{ll,ii}.ATpN.trial(:,1,MinIndUse))-mean(kcdist{ll,ii}.ApT.trial(:,1,MinIndUse));
%   end
% end
% [hh,pp]=ttest(meandiff');

%% plotting kcdist
for ll=[1 3 4 5 6 7 9]
  for ii=1:19
    meankcdist_atpn(ll,:,ii)=squeeze(mean(kcdist{ll,ii}.ATpN.trial,1));
    meankcdist_apt(ll,:,ii)=squeeze(mean(kcdist{ll,ii}.ApT.trial,1));
  end
end


%% 
for ll=[1 3 4 5 6 7 9]
  for ii=1:19
    minatpn(ll,ii)=min(kcdist{ll,ii}.ATpN.trial(:,1,MinIndUse));
    minapt(ll,ii)=min(kcdist{ll,ii}.ApT.trial(:,1,MinIndUse));
    maxatpn(ll,ii)=max(kcdist{ll,ii}.ATpN.trial(:,1,MinIndUse));
    maxapt(ll,ii)=max(kcdist{ll,ii}.ApT.trial(:,1,MinIndUse));    
  end
end
min([minatpn(:) minapt(:)])
max([maxatpn(:) maxapt(:)])
% edges=-218:10:137;  % gleaned from above max & min
% edges=-260:10:171;  % gleaned from above max & min
edges=-283:10:217;  % gleaned from above max & min
for ll=[1 3 4 5 6 7 9]
  for ii=1:19
    figure(ll);
    subplot(4,5,ii);hist(kcdist{ll,ii}.ATpN.trial(:,1,MinIndUse),edges);xlim([edges(1) edges(end)])
    figure(ll+10);
    subplot(4,5,ii);hist(kcdist{ll,ii}.ApT.trial(:,1,MinIndUse),edges);xlim([edges(1) edges(end)])
    figure(ll+20);
    subplot(4,5,ii);hist(kcdist{ll,ii}.AT.trial(:,MinIndUse),edges);xlim([edges(1) edges(end)])
    figure(ll+30);
    subplot(4,5,ii);hist(kcdist{ll,ii}.A.trial(:,MinIndUse),edges);xlim([edges(1) edges(end)])
    figure(ll+40);
    subplot(4,5,ii);hist(kcdist{ll,ii}.T.trial(:,MinIndUse),edges);xlim([edges(1) edges(end)])
    figure(ll+50);
    subplot(4,5,ii);hist(kcdist{ll,ii}.N.trial(:,MinIndUse),edges);xlim([edges(1) edges(end)])
  end
  print(ll,[fdir 'KcDist_ATpN_cond' num2str(ll) '.eps'],'-depsc')
  print(ll+10,[fdir 'KcDist_ApT_cond' num2str(ll) '.eps'],'-depsc')
  print(ll+20,[fdir 'KcDist_AT_cond' num2str(ll) '.eps'],'-depsc')
  print(ll+30,[fdir 'KcDist_A_cond' num2str(ll) '.eps'],'-depsc')
  print(ll+40,[fdir 'KcDist_T_cond' num2str(ll) '.eps'],'-depsc')
  print(ll+50,[fdir 'KcDist_N_cond' num2str(ll) '.eps'],'-depsc')
  print(ll,[fdir 'KcDist_ATpN_cond' num2str(ll) '.png'],'-dpng')
  print(ll+10,[fdir 'KcDist_ApT_cond' num2str(ll) '.png'],'-dpng')
  print(ll+20,[fdir 'KcDist_AT_cond' num2str(ll) '.png'],'-dpng')
  print(ll+30,[fdir 'KcDist_A_cond' num2str(ll) '.png'],'-dpng')
  print(ll+40,[fdir 'KcDist_T_cond' num2str(ll) '.png'],'-dpng')
  print(ll+50,[fdir 'KcDist_N_cond' num2str(ll) '.png'],'-dpng')
end



%% Different types of fitting/distribution assessment
% 1) KLdiv (take KLdiv relative to permuation distribution to 2nd level)
% 2) Kolmogorov-smirnov test (take ks-stat to 2nd level?)
% 3) ex-Gaussian fit: take weight of exponential to 2nd level?
% 4) Bimodal Gaussian: take mean of 2nd Gaussian to 2nd level?  weighte by proportion?
% 5) contrasting medians and/or rank test

%% 1) KLdiv

toadd=eps;
for ll=[1 3 4 5 6 7 9]
  for nn=1:19
    kcdist{ll,nn}.ATpN_ApT_KLdiv=computeAT_KLdiv(kcdist{ll,nn}.ATpN.trial(:,1,MinIndUse),kcdist{ll,nn}.ApT.trial(:,1,MinIndUse),edges,toadd,1);
    print(123,[fdir 'dist_diff_sub' num2str(nn) '.png'],'-dpng');close(123)
    kcdist{ll,nn}.AT_A_KLdiv=computeAT_KLdiv(kcdist{ll,nn}.AT.trial(:,MinIndUse),kcdist{ll,nn}.A.trial(:,MinIndUse),edges,toadd,0);
    kcdist{ll,nn}.T_N_KLdiv=computeAT_KLdiv(kcdist{ll,nn}.T.trial(:,MinIndUse),kcdist{ll,nn}.N.trial(:,MinIndUse),edges,toadd,0);
    
    kldivdist_ATpN_ApT(:,nn,ll)=permKLdiv(kcdist{ll,nn}.ATpN.trial(:,1,MinIndUse),kcdist{ll,nn}.ApT.trial(:,1,MinIndUse),edges,toadd);
    trueKLpercent_ATpN_ApT(ll,nn)=dsearchn(kldivdist_ATpN_ApT(:,nn,ll),kcdist{ll,nn}.ATpN_ApT_KLdiv)/length(kldivdist_ATpN_ApT(:,nn,ll))-.5;
    
    kldivdist_AT_A(:,nn,ll)=permKLdiv(kcdist{ll,nn}.AT.trial(:,MinIndUse),kcdist{ll,nn}.A.trial(:,MinIndUse),edges,toadd);
    trueKLpercent_AT_A(ll,nn)=dsearchn(kldivdist_AT_A(:,nn,ll),kcdist{ll,nn}.AT_A_KLdiv)/length(kldivdist_AT_A(:,nn,ll))-.5;
    
    if nn==1 && ll~=7
      kldivdist_T_N(:,nn,ll)=permKLdiv(kcdist{ll,nn}.T.trial(:,MinIndUse),kcdist{ll,nn}.N.trial(1:end-1,MinIndUse),edges,toadd);
    else
      kldivdist_T_N(:,nn,ll)=permKLdiv(kcdist{ll,nn}.T.trial(:,MinIndUse),kcdist{ll,nn}.N.trial(:,MinIndUse),edges,toadd);
    end
    trueKLpercent_T_N(ll,nn)=dsearchn(kldivdist_T_N(:,nn,ll),kcdist{ll,nn}.T_N_KLdiv)/length(kldivdist_T_N(:,nn,ll))-.5;
  end
end

% testing simulation of code
for np=1:100
  for nn=1:19
    aa(np,nn,:)=30*randn(1,200);
    bb(np,nn,:)=30*randn(1,200);
  end
end
edges_test=[min([min(bb(:)) min(aa(:))])  max([max(bb(:)) max(aa(:))])  ]

edges_test=[-140:10:140];
for np=1:100
  for nn=1:19
    kcdist_test{np,nn}.kldiv=computeAT_KLdiv(squeeze(aa(np,nn,:)),squeeze(bb(np,nn,:)),edges_test,toadd,0);
    kldivdist_test(:,nn,np)=permKLdiv(squeeze(aa(np,nn,:)),squeeze(bb(np,nn,:)),edges_test,toadd);
    trueKLpercent_test(np,nn)=dsearchn(kldivdist_test(:,nn,np),kcdist_test{np,nn}.kldiv)/length(kldivdist_test(:,nn,np))-.5;
  end
end
% length(find(sort(trueKLpercent_test(:))>.45))/1900
%     0.0632
% length(find(sort(trueKLpercent_test(:))<-.45))/1900
%     0.0416

% Taking KLpercent to 2nd level

for ll=[1 3 4 5 6 7 9]
  truemedian(ll)=median(trueKLpercent_ATpN_ApT(ll,:));
  truemean(ll)=mean(trueKLpercent_ATpN_ApT(ll,:));
  [permmedian(ll,:),permmean(ll,:)] = permSign(trueKLpercent_ATpN_ApT(ll,:));
  secLevelpercent_median(ll)=length(find(sort(permmedian(ll,:))<truemedian(ll)))/100;
  secLevelpercent_mean(ll)=length(find(sort(permmean(ll,:))<truemean(ll)))/100;
end

% secLevelpercent_median =
% % 0.7900         0    0.0200    0.0100    0.6400    0.2000    0.2400         0    0.7100
% 0.2000         0    0.4500    0.3300    0.5500    0.6000    0.5700         0    0.6300
% secLevelpercent_mean
% % 0.7400         0    0.1300    0.1000    0.6100    0.4000    0.1300         0    0.7900
% 0.3600         0    0.5200    0.3300    0.3700    0.6000    0.7400         0    0.7900

% trueKLpercent_ATpN_ApT is not Gaussain, but is symmetric.  Thus use 
for ll=[1 3 4 5 6 7 9]
  [p(ll),h]=signrank(trueKLpercent_ATpN_ApT(ll,:));
end
% % p =    0.3979         0    0.2271    0.1589    0.6007    0.6292    0.4458         0    0.4812
% p =    0.5195         0    0.8721    0.7781    1.0000    0.8248    0.5276         0    0.3836
% MinIndUse=12: p =    0.5325         0    0.9358    0.8562    0.3545 0.7781    0.2349         0    0.0268  ******

% % % Across asynchrony
llind=0;
for ll=[1 3 4 5 6 7 9]
  llind=llind+1;
  for nn=1:19
    kldiv_ATvA(llind,nn)=kcdist{ll,nn}.AT_A_KLdiv;
    kldiv_ATpNvApT(llind,nn)=kcdist{ll,nn}.ATpN_ApT_KLdiv;
  end
end

% [p,anovatab,stats]=kruskalwallis(kldiv_ATvA');
[p,anovatab,stats]=kruskalwallis(kldiv_ATpNvApT');
% MinIndUse=12:  p =    0.2581

% llind=0;
% for ll=[1 3 4 5 6 7 9]
%   llind=llind+1;
%   for nn=1:19
%     if ll<5
%       kldiv_ATvs2nd1(llind,nn)=computeAT_KLdiv(kcdist{ll,nn}.AT.trial(:,MinIndUse),kcdist{ll,nn}.T.trial(:,MinIndUse),edges,toadd,0);
%       kldiv_ATvs2nd2(llind,nn)=computeAT_KLdiv(kcdist{ll,nn}.AT.trial(:,MinIndUse),kcdist{ll,nn}.T.trial(:,MinIndUse),edges,toadd,0);
%     elseif ll==5
%       kldiv_ATvs2nd1(llind,nn)=computeAT_KLdiv(kcdist{ll,nn}.AT.trial(:,MinIndUse),kcdist{ll,nn}.A.trial(:,MinIndUse),edges,toadd,0);
%       kldiv_ATvs2nd2(llind,nn)=computeAT_KLdiv(kcdist{ll,nn}.AT.trial(:,MinIndUse),kcdist{ll,nn}.T.trial(:,MinIndUse),edges,toadd,0);
%     elseif ll>5
%       kldiv_ATvs2nd1(llind,nn)=computeAT_KLdiv(kcdist{ll,nn}.AT.trial(:,MinIndUse),kcdist{ll,nn}.A.trial(:,MinIndUse),edges,toadd,0);
%       kldiv_ATvs2nd2(llind,nn)=computeAT_KLdiv(kcdist{ll,nn}.AT.trial(:,MinIndUse),kcdist{ll,nn}.A.trial(:,MinIndUse),edges,toadd,0);
%     end
%   end
% end
% 
% [p,anovatab,stats]=kruskalwallis(kldiv_ATvs2nd1');
% [p,anovatab,stats]=kruskalwallis(kldiv_ATvs2nd2');



%% 2) KS test

for ll=[1 3 4 5 6 7 9]
  for nn=1:19
    [h,p,kcdist{ll,nn}.ATpN_ApT_KStest]=kstest2(kcdist{ll,nn}.ATpN.trial(:,1,MinIndUse),kcdist{ll,nn}.ApT.trial(:,1,MinIndUse));
    [h,p,kcdist{ll,nn}.AT_A_KStest]=kstest2(kcdist{ll,nn}.AT.trial(:,MinIndUse),kcdist{ll,nn}.A.trial(:,MinIndUse));
    [h,p,kcdist{ll,nn}.T_N_KStest]=kstest2(kcdist{ll,nn}.T.trial(:,MinIndUse),kcdist{ll,nn}.N.trial(:,MinIndUse));
    
    ksstatdist_ATpN_ApT(:,nn,ll)=permKSstat(kcdist{ll,nn}.ATpN.trial(:,1,MinIndUse),kcdist{ll,nn}.ApT.trial(:,1,MinIndUse));
    trueKSpercent_ATpN_ApT(ll,nn)=dsearchn(ksstatdist_ATpN_ApT(:,nn,ll),kcdist{ll,nn}.ATpN_ApT_KStest)/length(ksstatdist_ATpN_ApT(:,nn,ll))-.5;
    
    ksstatdist_AT_A(:,nn,ll)=permKSstat(kcdist{ll,nn}.AT.trial(:,MinIndUse),kcdist{ll,nn}.A.trial(:,MinIndUse));
    trueKSpercent_AT_A(ll,nn)=dsearchn(ksstatdist_AT_A(:,nn,ll),kcdist{ll,nn}.AT_A_KStest)/length(ksstatdist_AT_A(:,nn,ll))-.5;
    
    if nn==1 && ll~=7
      ksstatdist_T_N(:,nn,ll)=permKSstat(kcdist{ll,nn}.T.trial(:,MinIndUse),kcdist{ll,nn}.N.trial(1:end-1,MinIndUse));
    else
      ksstatdist_T_N(:,nn,ll)=permKSstat(kcdist{ll,nn}.T.trial(:,MinIndUse),kcdist{ll,nn}.N.trial(:,MinIndUse));
    end
    trueKSpercent_T_N(ll,nn)=dsearchn(ksstatdist_T_N(:,nn,ll),kcdist{ll,nn}.T_N_KStest)/length(ksstatdist_T_N(:,nn,ll))-.5;
  end
end

% testing simulation of code
for np=1:100
  for nn=1:19
    aa(np,nn,:)=30*randn(1,200);
    bb(np,nn,:)=30*randn(1,200);
  end
end

for np=1:100
  for nn=1:19
    [h,p,kcdist_test{np,nn}.ksstat]=kstest2(squeeze(aa(np,nn,:)),squeeze(bb(np,nn,:)));
    ksstatdist_test(:,nn,np)=permKSstat(squeeze(aa(np,nn,:)),squeeze(bb(np,nn,:)));
    trueKSpercent_test(np,nn)=dsearchn(ksstatdist_test(:,nn,np),kcdist_test{np,nn}.ksstat)/length(ksstatdist_test(:,nn,np))-.5;
  end
end
% length(find(sort(trueKSpercent_test(:))>.45))/1900
%     0.0521
% length(find(sort(trueKSpercent_test(:))<-.45))/1900
%     0.0453


% Taking KSpercent to 2nd level

for ll=[1 3 4 5 6 7 9]
  truemedianKS(ll)=median(trueKSpercent_ATpN_ApT(ll,:));
  truemeanKS(ll)=mean(trueKSpercent_ATpN_ApT(ll,:));
  [permmedianKS(ll,:),permmeanKS(ll,:)] = permSign(trueKSpercent_ATpN_ApT(ll,:));
  secLevelpercentKS_median(ll)=length(find(sort(permmedianKS(ll,:))<truemedianKS(ll)))/100;
  secLevelpercentKS_mean(ll)=length(find(sort(permmeanKS(ll,:))<truemeanKS(ll)))/100;
end

% secLevelpercentKS_median =
%  0.2100         0    0.6900    0.9500    0.2200    0.5300    0.4000         0    0.5100
%  0.0100         0    0.4000    0.5800    0.6100    0.2800    0.2600         0    0.7500
% secLevelpercent_mean
%  0.4000         0    0.9300    0.9600    0.3200    0.6300    0.5100         0    0.5400
%  0.3600         0    0.5200    0.3300    0.3700    0.6000    0.7400         0    0.7900
 
% [h,p]=ttest(trueKSpercent_ATpN_ApT');
% p =    0.6964       NaN    0.1722    0.0501    0.6825    0.8692    0.9774       NaN    0.8967
for ll=[1 3 4 5 6 7 9]
 [p(ll),h]=signrank(trueKSpercent_ATpN_ApT(ll,:));
end
% p =  0.7475         0    0.1117    0.0464    0.6872    0.9198    0.8721         0    0.9826
% p =  0.1649         0    0.6580    0.7022    0.7438    0.6874    0.8445         0    0.4093    
% mininduse=12;  p =  0.8248         0    0.4091    0.8106    0.2954    0.7020    0.7781         0    0.3046    

    % % % Across asynchrony
llind=0;
for ll=[1 3 4 5 6 7 9]
  llind=llind+1;
  for nn=1:19
    kstest_ATvA(llind,nn)=kcdist{ll,nn}.AT_A_KStest;
    kstest_ATpNvApT(llind,nn)=kcdist{ll,nn}.ATpN_ApT_KStest;
  end
end

% [p,anovatab,stats]=kruskalwallis(kstest_ATvA');
[p,anovatab,stats]=kruskalwallis(kstest_ATpNvApT');

% llind=0;
% for ll=[1 3 4 5 6 7 9]
%   llind=llind+1;
%   for nn=1:19
%     if ll<5
%       [h,p,ksstat_ATvs2nd1(llind,nn)]=kstest2(kcdist{ll,nn}.AT.trial(:,MinIndUse),kcdist{ll,nn}.T.trial(:,MinIndUse));
%       [h,p,ksstat_ATvs2nd2(llind,nn)]=kstest2(kcdist{ll,nn}.AT.trial(:,MinIndUse),kcdist{ll,nn}.T.trial(:,MinIndUse));
%     elseif ll==5
%       [h,p,ksstat_ATvs2nd1(llind,nn)]=kstest2(kcdist{ll,nn}.AT.trial(:,MinIndUse),kcdist{ll,nn}.A.trial(:,MinIndUse));
%       [h,p,ksstat_ATvs2nd2(llind,nn)]=kstest2(kcdist{ll,nn}.AT.trial(:,MinIndUse),kcdist{ll,nn}.T.trial(:,MinIndUse));
%     elseif ll>5
%       [h,p,ksstat_ATvs2nd1(llind,nn)]=kstest2(kcdist{ll,nn}.AT.trial(:,MinIndUse),kcdist{ll,nn}.A.trial(:,MinIndUse));
%       [h,p,ksstat_ATvs2nd2(llind,nn)]=kstest2(kcdist{ll,nn}.AT.trial(:,MinIndUse),kcdist{ll,nn}.A.trial(:,MinIndUse));
%     end
%   end
% end
% 
% [p,anovatab,stats]=kruskalwallis(ksstat_ATvs2nd1');
% [p,anovatab,stats]=kruskalwallis(ksstat_ATvs2nd2');
% 




%% 3) ex-Gaussian fit

for ll=[1 3 4 5 6 7 9]
  for nn=1:19
    [X_ATpN(:,ll,nn),NLL_ATpN(ll,nn)]=exgauss_fit(-kcdist{ll,nn}.ATpN.trial(:,1,MinIndUse));
    [X_ApT(:,ll,nn),NLL_ApT(ll,nn)]=exgauss_fit(-kcdist{ll,nn}.ApT.trial(:,1,MinIndUse));
    [X_AT(:,ll,nn),NLL_AT(ll,nn)]=exgauss_fit(-kcdist{ll,nn}.AT.trial(:,MinIndUse));
    [X_A(:,ll,nn),NLL_A(ll,nn)]=exgauss_fit(-kcdist{ll,nn}.A.trial(:,MinIndUse));
    [X_T(:,ll,nn),NLL_T(ll,nn)]=exgauss_fit(-kcdist{ll,nn}.T.trial(:,MinIndUse));
    [X_N(:,ll,nn),NLL_N(ll,nn)]=exgauss_fit(-kcdist{ll,nn}.N.trial(:,MinIndUse));
  end
end
% X is [mu,sigma,tau]

% [h,p]=ttest(squeeze(X_A([1 3],5,:))',squeeze(X_N([1 3],5,:))')
% % p =0.0021    0.0000
% [h,p]=ttest(squeeze(X_T([1 3],5,:))',squeeze(X_N([1 3],5,:))')
% p =0.5008    0.2198
% for ll=[1 3 4 5 6 7 9]
%   [h,p_ATvA(ll,:)]=ttest(squeeze(X_AT([1 3],ll,:))',squeeze(X_A([1 3],ll,:))');
%   [h,p_ATpNvApT(ll,:)]=ttest(squeeze(X_ATpN([1 3],ll,:))',squeeze(X_ApT([1 3],ll,:))');
% end
% p_ATpNvApT =
%     0.5101    0.4749
%          0         0
%     0.8268    0.7759
%     0.1095    0.4446
%     0.9521    0.7308
%     0.4684    0.8076
%     0.8838    0.5031
%          0         0
%     0.6502    0.0644

for ll=[1 3 4 5 6 7 9]
  for param=1:3
    [p_ATvA(ll,param)]=signrank(squeeze(X_AT(param,ll,:))',squeeze(X_A(param,ll,:))');
    [p_ATpNvApT(ll,param)]=signrank(squeeze(X_ATpN(param,ll,:))',squeeze(X_ApT(param,ll,:))');
  end
end

% p_ATpNvApT =
%     0.3547    0.8405    0.4209
%          0         0         0
%     0.8721    0.2432    0.6292
%     0.3760    0.8721    0.8405
%     0.6292    0.9679    0.3144
%     0.7172    0.8092    0.7782
%     0.4688    0.5197    0.7172
%          0         0         0
%     0.8405    0.4939    0.0534

% p_ATpNvApT =
%     0.2772    0.1978    0.0990
%          0         0         0
%     0.4209    0.6292    0.2598
%     0.2432    0.8721    1.0000
%     0.6292    0.8405    0.7782
%     0.8092    0.9679    0.8721
%     0.1590    0.8092    0.1590
%          0         0         0
%     0.8405    0.3341    0.2954

% p_ATpNvApT =   % MinIndUse=12
% 
%     0.1165    0.9359    0.4445
%          0         0         0
%     0.1075    0.9359    0.3547
%     0.2273    0.9039    0.7782
%     0.4209    0.1262    0.8405
%     0.3981    0.8721    0.7475
%     0.6874    0.7782    0.3144
%          0         0         0
%     0.4688    0.3144    0.1474

median([squeeze(X_ATpN(param,ll,:))'; squeeze(X_ApT(param,ll,:))']')


kruskalwallis([squeeze(X_ATpN(1,[1 3 4 5 6 7 9],:))-squeeze(X_ApT(1,[1 3 4 5 6 7 9],:))]')
kruskalwallis([squeeze(X_ATpN(2,[1 3 4 5 6 7 9],:))-squeeze(X_ApT(2,[1 3 4 5 6 7 9],:))]')
kruskalwallis([squeeze(X_ATpN(3,[1 3 4 5 6 7 9],:))-squeeze(X_ApT(3,[1 3 4 5 6 7 9],:))]')

% None significant


%% 4) bimodal Gaussian:  we decided distribution doesn't fit this model well enough, so ignoring this option.

% we've decided to only look at TA500
ll=9;
for nn=1:19
  alldata{nn}=[kcdist{ll,nn}.ATpN.trial(:,1,MinIndUse); kcdist{ll,nn}.ApT.trial(:,1,MinIndUse)];
  try
    GM_1{nn}=fitgmdist(alldata{nn},2);
  end
%   GM_ATpN{nn}=fitgmdist(kcdist{ll,nn}.ATpN.trial(:,1,10),2);
%   GM_ApT{nn}=fitgmdist(kcdist{ll,nn}.ApT.trial(:,1,10),2);
end

GM_1{nn}=fitgmdist_JZ(alldata{nn},2);

Sigma(:,:,1)=3;
Sigma(:,:,2)=10;
LogL=computeLLforGMM_JZ(alldata{nn},[0 -5], Sigma, [.75 .25])

f2mg = @computeLLforGMM_JZ;
c1=[0 -5];
c2=Sigma;
c3=[.75 .25];
output=fminsearch( @(X) f2mg(x,c1,c2,c3),alldata{nn});

% 
% fminsearch
% 
% f = @(x,c) x(1).^2+c.*x(2).^2;  % The parameterized function.
% c = 1.5;                        % The parameter.
% X = fminsearch(@(x) f(x,c),[0.3;1])
         


%% 5) Medians:  Wilcoxon Sign Rank test is for paired samples, non-parametric

for ll=[1 3 4 5 6 7 9]
  for nn=1:19
    [p,anovatab{ll,nn},stats]=kruskalwallis([kcdist{ll,nn}.ATpN.trial(:,1,MinIndUse), kcdist{ll,nn}.ApT.trial(:,1,MinIndUse)])
    chisquare(ll,nn)=anovatab{ll,nn}{2,5};
    [p,h,signrank_stats_sort{ll,nn}]=signrank(sort(kcdist{ll,nn}.ATpN.trial(:,1,MinIndUse)),sort(kcdist{ll,nn}.ApT.trial(:,1,MinIndUse)));
    zval_signrank(ll,nn)=signrank_stats_sort{ll,nn}.zval;
    [p,h,ranksum_stats{ll,nn}]=ranksum((kcdist{ll,nn}.ATpN.trial(:,1,MinIndUse)),(kcdist{ll,nn}.ApT.trial(:,1,MinIndUse)));
    zval_ranksum(ll,nn)=ranksum_stats{ll,nn}.zval;
  end
end

% [hh,pp]=ttest(zval_signrank')
% % pp =    0.7977       NaN    0.2490    0.1035    0.4007    0.2674    0.0559       NaN    0.0026
[hh,pp]=ttest(zval_ranksum')
% pp =    0.9813       NaN    0.6008    0.1148    0.4443    0.2688    0.0794       NaN    0.0154
% pp =    0.2734       NaN    0.6084    0.1295    0.8005    0.7168    0.4510       NaN    0.0105
% MinIndUse=12: pp =  0.6459       NaN    0.2750    0.2905    0.6366    0.1143    0.2110       NaN    0.0314    

    
% kruskalwallis(zval_signrank([1 3 4 5 6 7 9],:)')
kruskalwallis(zval_ranksum([1 3 4 5 6 7 9],:)')
kruskalwallis(chisquare([1 3 4 5 6 7 9],:)')


