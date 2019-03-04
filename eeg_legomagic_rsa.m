eeg_legomagic_preamble

%% Computing RSA (Cecere 2017 modified) and RDM (Nili 2014 MRC-CAM RSA toolbox)

clearvars -except sub *dir
iterflag=1;
sleep=0;
tt=3;
statwinorig=0;  % =1 means the latency range of .1 to .45 (first way);  =0 means 0 to .5 (final /better way)
ftver=0;
mcseed=13;
usetr=0;
plotallmask=1; % =0 means each subplot will use significane mask relevant for those sensors & time plotted only
% =1 means each subplot will use mask relevant for all sensors even those not plotted.

if sleep
  iter=11;
  ss=12;
  trialkc=0;
else
  iter=27;  % previous result; final way (easier!)
  %   iter=31;  % smart sampling
  ss=10;
  trialkc=-1;
end

if iterflag
  if sleep
    try
      %       load([edir 'tlock_statmc_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '_statwinorig' num2str(statwinorig) '.mat']);
      load([edir 'tlock_grind_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '_statwinorig' num2str(statwinorig) '.mat']);
    catch
      load([edir 'tlock_statmc_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat']);
      load([edir 'tlock_grind_sleep' num2str(sleep) '_ss' num2str(ss) '_iter' num2str(iter) '_trialkc' num2str(trialkc) '.mat']);
    end
  else
    try
      %       load([edir 'tlock_statmc_sleep' num2str(sleep) '_iter' num2str(iter) '_usetr' num2str(usetr) '_statwinorig' num2str(statwinorig) '_ftver' num2str(ftver) 'mcseed' num2str(mcseed) '.mat']);
      load([edir 'tlock_grind_sleep' num2str(sleep) '_iter' num2str(iter) '_usetr' num2str(usetr) '_statwinorig' num2str(statwinorig) '_ftver' num2str(ftver) '.mat']);
    catch
      try
        %         load([edir 'tlock_statmc_sleep' num2str(sleep) '_iter' num2str(iter) '_statwinorig' num2str(statwinorig) '.mat']);
        load([edir 'tlock_grind_sleep' num2str(sleep) '_iter' num2str(iter) '_statwinorig' num2str(statwinorig) '.mat']);
      catch
        %         load([edir 'tlock_statmc_sleep' num2str(sleep) '_iter' num2str(iter) '.mat']);
        load([edir 'tlock_grind_sleep' num2str(sleep) '_iter' num2str(iter) '.mat']);
      end
    end
  end
else
  %   load([edir 'tlock_statmc_sleep' num2str(sleep) '.mat']);
  load([edir 'tlock_grind_sleep' num2str(sleep) '.mat']);
end

soalist=[1 3 4 5 6 7 9];
timwinstatflag=1; % =1 for using only stat.time, =0 for using [-0.5 1.0];
timwin=[-0.5 1];
stattime=[0 0.5];
soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];

timestart=[0 nan 0 0 0 0.2 .07 nan .5];
timeend=timestart+.5;
startind=[500 nan 500 500 500 520 570 nan 1000];
endind=startind+500;

llcnt=0;lllcnt=0;
rsa_scorr=zeros(7,7,25,22);
rsa_edist=zeros(7,7,25,22);
% rsa_cov=zeros(7,7,25,22);
for ll=soalist
  llcnt=llcnt+1;
  lllcnt=0;
  for lll=soalist
    lllcnt=lllcnt+1;
    for ii=1:22 % number of participants
      val1=squeeze(grind_TPA_MSPN{ll,tt,10}.individual(ii,:,startind(ll):endind(ll)));
      val2=squeeze(grind_TPA_MSPN{lll,tt,10}.individual(ii,:,startind(lll):endind(lll)));
      for jj=1:25 % every 20 ms
        valslide1=mean(val1(:,jj*20-19:jj*20+1),2);
        valslide2=mean(val2(:,jj*20-19:jj*20+1),2);
        rsa_scorr(llcnt,lllcnt,jj,ii)=1-corr(valslide1,valslide2,'type','Spearman'); % 1-r for dissimilarity matrix
        rsa_edist(llcnt,lllcnt,jj,ii)=norm(valslide1-valslide2);
        %         tmpcov=cov(valslide1,valslide2);
        %         rsa_cov(llcnt,lllcnt,jj,ii)=tmpcov(1,2); % how to make this dissimilar?   1./ or 1-?  or - ?
      end
    end
  end
end

figure;for jj=1:25,subplot(5,5,jj);imagesc(mean(rsa_scorr(:,:,jj,:),4));caxis([0 1]);title(num2str(jj*20-10));end;
figure;for jj=1:25,subplot(5,5,jj);imagesc(mean(rsa_edist(:,:,jj,:),4));caxis([0 32]);title(num2str(jj*20-10));end;
figure;for jj=1:25,subplot(5,5,jj);imagesc(mean(rsa_cov(:,:,jj,:),4));caxis([-10 10]);title(num2str(jj*20-10));end;


% Models:

% 1) TIW, fall off as Gaussian with SOA
rsmodel_1a=zeros(7,7);
rsmodel_1a(3:5,3:5)=1;

rsmodel_1b=zeros(7,7);
rsmodel_1b(2:6,2:6)=1;

rsmodel_1c=zeros(7,7);
rsmodel_1c(2:5,2:5)=1;

rsmodel_1d=zeros(7,7);
rsmodel_1d(3:6,3:6)=1;

rsmodel_1e=zeros(7,7);
rsmodel_1e(2:4,2:4)=1;

rsmodel_1f=zeros(7,7);
rsmodel_1f(4:6,4:6)=1;

% this one is linear combo of above two, thus not needed.
% rsmodel_1c=zeros(7,7);
% rsmodel_1c(2,2:6)=1;
% rsmodel_1c(3:5,2)=1;
% rsmodel_1c(3:5,6)=1;
% rsmodel_1c(6,2:6)=1;

figure;imagesc(rsmodel_1a);caxis([-1 1]);
figure;imagesc(rsmodel_1b);caxis([-1 1]);
figure;imagesc(rsmodel_1c);caxis([-1 1]);
figure;imagesc(rsmodel_1d);caxis([-1 1]);
figure;imagesc(rsmodel_1e);caxis([-1 1]);
figure;imagesc(rsmodel_1f);caxis([-1 1]);


% 2) asymmetry: AT will be similar to each other but different from TA (vice versa)
rsmodel_2a=zeros(7,7);
rsmodel_2a(1:4,1:4)=1;
rsmodel_2a(4:7,4:7)=1;
figure;imagesc(rsmodel_2a);caxis([-1 1]);
% rdmodel_2a=-rsmodel_2a;
% figure;imagesc(rdmodel_2a);

% rsamod2b=eye(7,7);
% rsamod2b(1:3,1:3)=1;
% rsamod2b(5:7,5:7)=1;
% rdmcand2b=-rsamod2b;
% figure;imagesc(rsamod2b);
% figure;imagesc(rdmcand2b);
% rsamod2c=zeros(7,7);
% rsamod2c(1:3,1:3)=1;
% rsamod2c(5:7,5:7)=1;
% figure;imagesc(rsamod2c);

% 3) symmetric: TA70 will be more like itself and AT70 than others

% rsamod3=eye(7,7);
% rsamod3(3,5)=1;
% rsamod3(2,6)=1;
% rsamod3(1,7)=1;
% rsamod3(5,3)=1;
% rsamod3(6,2)=1;
% rsamod3(7,1)=1;
% figure;imagesc(rsamod3);
% rsamod3b=zeros(7,7);
% rsamod3b(3,5)=1;
% rsamod3b(2,6)=1;
% rsamod3b(1,7)=1;
% rsamod3b(5,3)=1;
% rsamod3b(6,2)=1;
% rsamod3b(7,1)=1;
% figure;imagesc(rsamod3b);
% rsamod3c=zeros(7,7);
% rsamod3c(2,6)=1;
% rsamod3c(6,2)=1;
% figure;imagesc(rsamod3c);
rsmodel_3a=zeros(7,7);
rsmodel_3a(1,7)=1;
rsmodel_3a(2,6)=1;
rsmodel_3a(3,5)=1;
rsmodel_3a(4,4)=1;
rsmodel_3a(5,3)=1;
rsmodel_3a(6,2)=1;
rsmodel_3a(7,1)=1;
figure;imagesc(rsmodel_3a);caxis([-1 1]);
% rdmcand3b=ones(7,7);
% rdmcand3b(1,7)=0;
% rdmcand3b(7,1)=0;
% figure;imagesc(rdmcand3b);
% rdmcand3c=ones(7,7);
% rdmcand3c(2,6)=0;
% rdmcand3c(6,2)=0;
% figure;imagesc(rdmcand3c);
% rdmcand3d=ones(7,7);
% rdmcand3d(3,5)=0;
% rdmcand3d(5,3)=0;
% figure;imagesc(rdmcand3d);

% 4) Behavioural RT
load([ddir 'rt_allsubj.mat'],'rt*');
rt_allsubjuse=squeeze(rt_msMminshiftuni(soalist,3,:));
rt_avgsubj=nanmean(rt_allsubjuse,2);
for ii=1:size(rt_allsubjuse,2)
  rt_cov_ind(:,:,ii)=rt_allsubjuse(:,ii)*rt_allsubjuse(:,ii)';
end
rt_cov_all=rt_avgsubj*rt_avgsubj';

rsmodel_4a=rt_cov_all/max(abs(rt_cov_all(:)));
figure;imagesc(rsmodel_4a);caxis([-1 1])

% timevec=-.5:.01:.5;
% rtvec=nan(101,1);
% timeind=[dsearchn(timevec',-.5) dsearchn(timevec',-.07) dsearchn(timevec',-.02) dsearchn(timevec',0) dsearchn(timevec',.02) dsearchn(timevec',.07) dsearchn(timevec',.5)]
% rtvec(1)=rtvals(1);
% rtvec(dsearchn(timevec',-.07))=rtvals(2);
% rtvec(dsearchn(timevec',-.02))=rtvals(3);
% rtvec(dsearchn(timevec',0))=rtvals(4);
% rtvec(dsearchn(timevec',.02))=rtvals(5);
% rtvec(dsearchn(timevec',.07))=rtvals(6);
% rtvec(dsearchn(timevec',.5))=rtvals(7);
% w=gausswin(101,5); % adjust 2nd param for width
% figure;plot(timevec,0.037*w);hold on;
% plot(timevec,rtvec,'o');
% rsamod1=sqrt(w(timeind)*w(timeind)');  % take sqrt since we want diagonal to be originall gaussian values
% rdmcand1a=1-rsamod1;
% figure;imagesc(rsamod1);
% figure;imagesc(rdmcand1a);
% uptriind=find(triu(rsamod1));



% % % 4) simple similarity to self.  (don't use: don't fit diagonal anyway)
% % rsamod4=eye(7,7);
% % figure;imagesc(rsamod4);

% Interaction terms:
interact_1a2a=rsmodel_1a.*rsmodel_2a;
interact_1b2a=rsmodel_1b.*rsmodel_2a;
interact_1a3a=rsmodel_1a.*rsmodel_3a;
interact_1b3a=rsmodel_1b.*rsmodel_3a;
% interact_2a3a=rsmodel_2a.*rsmodel_3a;  % this gives only the centre point 4,4
interact_1a4a=rsmodel_1a.*rsmodel_4a;
interact_1b4a=rsmodel_1b.*rsmodel_4a;
interact_2a4a=rsmodel_2a.*rsmodel_4a;
interact_3a4a=rsmodel_3a.*rsmodel_4a;
figure;imagesc(interact_1a2a);caxis([-1 1]);
figure;imagesc(interact_1b2a);caxis([-1 1]);
figure;imagesc(interact_1a3a);caxis([-1 1]);
figure;imagesc(interact_1b3a);caxis([-1 1]); % scale RT
% figure;imagesc(interact_2a3a);caxis([-1 1]);
figure;imagesc(interact_1a4a);caxis([-1 1]);
figure;imagesc(interact_1b4a);caxis([-1 1]);
figure;imagesc(interact_2a4a);caxis([-1 1]);
figure;imagesc(interact_3a4a);caxis([-1 1]);

% rdmcand4_1a2a=rdmcand1a.*rdmcand2a;
% rdmcand4_1a2b=rdmcand1a.*rdmcand2b;
% rdmcand4_1b2a=rdmcand1b.*rdmcand2a;
% rdmcand4_1b2b=rdmcand1b.*rdmcand2b;
% figure;imagesc(rdmcand4_1a2a)
%
% rdmcand5_1a3a=rdmcand1a.*rdmcand3a;
% rdmcand5_1b3a=rdmcand1b.*rdmcand3a;
% figure;imagesc(rdmcand5_1a3a);

% rdmcand6_2a3a=rdmcand2a.*rdmcand3a;
% rdmcand6_2a3b=rdmcand2a.*rdmcand3b;
% rdmcand6_2b3a=rdmcand2b.*rdmcand3a;
% rdmcand6_2b3b=rdmcand2b.*rdmcand3b;

% % Define Models as inclusion/exclusion of regressors (candidates)
% Model 1: RT alone

% Model 2: RT + asymmetry

% Model 3: RT + sym-pairs

% Model 4: RT + asymmetry + sym-pairs

% Model 5: RT + asymmetry + RT_X_asymmetry

% Model 6: RT + sym-pairs + RT_X_sym-pairs

% Model 7: RT + asymmetry + sym-pairs + RT_X_asymmetry

% Model 8: RT + asymmetry + sym-pairs + RT_X_sym-pairs

% Model 9: RT + asymmetry + sym-pairs + RT_X_asymmetry + RT_X_sym-pairs

