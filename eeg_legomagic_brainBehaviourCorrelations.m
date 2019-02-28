% brain-behaviour correlations
clear all
close all
if ispc
  edir='D:\audtac\eeg_data\';
  esdir='D:\audtac\source_data\';
  ddir='D:\audtac\legomagic\diaries\';
  bdir='D:\audtac\behav_data\';
  sdir='D:\audtac\spss_stuff\';
  fdir='D:\audtac\figs\';
  mdir='D:\audtac\structural_MRI\';
  pdir='D:\audtac\polhemus\';
else
  [~,hostname]=system('hostname');
  if ~isempty(strfind(hostname,'les')) | ~isempty(strfind(hostname,'LES')) % either COLLES-151401 or LES-LINUX_FS3
    edir='/home/zumerj/audtac/eeg_data/';
    esdir='/home/zumerj/audtac/source_data/';
    ddir='/home/zumerj/audtac/legomagic/diaries/';
    bdir='/home/zumerj/audtac/behav_data/';
    sdir='/home/zumerj/audtac/spss_stuff/';
    fdir='/home/zumerj/audtac/figs/';
    mdir='/home/zumerj/audtac/structural_MRI/';
    pdir='/home/zumerj/audtac/polhemus/';
  else % assume on VM linux of psychl-132432
    edir='/mnt/hgfs/D/audtac/eeg_data/';
    esdir='/mnt/hgfs/D/audtac/source_data/';
    ddir='/mnt/hgfs/D/audtac/legomagic/diaries/';
    bdir='/mnt/hgfs/D/audtac/behav_data/';
    sdir='/mnt/hgfs/D/audtac/spss_stuff/';
    fdir='/mnt/hgfs/D/audtac/figs/';
    mdir='/mnt/hgfs/D/audtac/structural_MRI/';
    pdir='/mnt/hgfs/D/audtac/polhemus/';
  end
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
  warning off
  rmpath(genpath('D:\matlab\spm8\external\fieldtrip\'))
  rmpath(genpath('D:\fieldtrip_svn\'))
  rmpath(genpath('D:\fieldtrip_git\'))
  warning on
  addpath('D:\fieldtrip_git\')
else
  rmpath(genpath('/mnt/hgfs/D/matlab/spm8/external/fieldtrip/'))
  rmpath(genpath('/mnt/hgfs/D/fieldtrip_svn/'))
  addpath('/mnt/hgfs/D/fieldtrip_svn/')
end
which ft_defaults.m
ft_defaults;
load([edir 'iikeep.mat'])
soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];
soalist= [1 3 4 5 6 7 9 10 11 12]; % MS conditions, Nul, Aud, Tac
condcode=[1 nan 3 4 5 6 7 nan 9 0  -1 -2]; % MS conditions, Nul, Aud, Tac

%%  RT analysis
cd(ddir)

sleep=0;
if sleep
  iiuse=setdiff(iiBuse,3:7);
else
  iiuse=iiSuse;
end

iiuse=[10:18 20:32]; % added Sept 2018; use all possible

% rtsrall=nan(1,10);
% diffms=nan(1,max(iiuse));
% pcb=nan(1,max(iiuse));

bbind=1;
% the 10 columns of 'rtsr' are: [Talone Aalone Nul 1 3 4 5 6 7 9]
for ii=iiuse
  load([sub{ii} '_audtac.mat'],'r*');
  
%   for tt=2:3
  for tt=3
    for ll=fliplr(soalist)
      rt{ll,tt,ii}=rrts(find(rsoareal(tt-1,:)==condcode(ll)));
      numnan(1,ll,bbind)=length(find(isnan(rt{ll,tt,ii})));
      % added 12 July, 2018
%       rt{ll,tt,ii}(1)=nan; % remove first trial of block, not of stim cond
      rt{ll,tt,ii}(find(rt{ll,tt,ii}<.1))=nan; % anticipatory? (or fault with computing timetouch?);   % Throwing out trials < 100ms
      numnan(2,ll,bbind)=length(find(isnan(rt{ll,tt,ii})));
      rt{ll,tt,ii}(find(rt{ll,tt,ii}>1))=nan; % asleep?   % Throwing out trials greater than 1s
      numnan(3,ll,bbind)=length(find(isnan(rt{ll,tt,ii})));
      numtot(ll,bbind)=length(rt{ll,tt,ii});
      % end add
      

      % save out some per-subject metrics of effect
      % 1) actual RT
      rt_med(ll,tt,bbind)=nanmedian(rt{ll,tt,ii});
      if ll<10
        % 2) difference of fastest unisensory minus each MS
        rt_msMminuni(ll,tt,bbind)=min(rt_med(11:12,tt,bbind))-rt_med(ll,tt,bbind);
        if min(rt_med(11:12,tt,bbind))==rt_med(11,tt,bbind)
          [hh,pp]=ttest2(rt{11,tt,ii},rt{ll,tt,ii});
        else
          [hh,pp]=ttest2(rt{12,tt,ii},rt{ll,tt,ii});
        end
        rt_msMminuni_pval(ll,tt,bbind)=pp;
        % 3) difference of mean of unisensory minus each MS
        rt_msMmeanuni(ll,tt,bbind)=mean(rt_med(11:12,tt,bbind))-rt_med(ll,tt,bbind);
        [hh,pp]=ttest2([rt{11,tt,ii} rt{12,tt,ii}],rt{ll,tt,ii});
        rt_msMmeanuni_pval(ll,tt,bbind)=pp;
        % 4) difference of relevant unisensory minus MS
        if ll<5
          rt_msMreluni(ll,tt,bbind)=rt_med(11,tt,bbind)-rt_med(ll,tt,bbind);
          [hh,pp]=ttest2(rt{11,tt,ii},rt{ll,tt,ii});
          rt_msMreluni_pval(ll,tt,bbind)=pp;
        elseif ll>5
          rt_msMreluni(ll,tt,bbind)=rt_med(12,tt,bbind)-rt_med(ll,tt,bbind);
          [hh,pp]=ttest2(rt{12,tt,ii},rt{ll,tt,ii});
          rt_msMreluni_pval(ll,tt,bbind)=pp;
        elseif ll==5
          rt_msMreluni(ll,tt,bbind)=rt_msMminuni(ll,tt,bbind);
          rt_msMreluni_pval(ll,tt,bbind)=rt_msMminuni_pval(ll,tt,bbind);
        end
        % 5) difference of shifted-unisensory minus MS
        %    e.g.  min(A,T+20) - AT20
        if ll<6
          rt_msMminshiftuni(ll,tt,bbind)=min(rt_med(11,tt,bbind),rt_med(12,tt,bbind)-soades(ll))-rt_med(ll,tt,bbind);
          if min(rt_med(11,tt,bbind),rt_med(12,tt,bbind)-soades(ll))==rt_med(11,tt,bbind)
            [hh,pp]=ttest2(rt{11,tt,ii},rt{ll,tt,ii});
          else
            [hh,pp]=ttest2(rt{12,tt,ii}-soades(ll),rt{ll,tt,ii});
          end
        else
          rt_msMminshiftuni(ll,tt,bbind)=min(rt_med(11,tt,bbind)+soades(ll),rt_med(12,tt,bbind))-rt_med(ll,tt,bbind);
          if min(rt_med(11,tt,bbind)+soades(ll),rt_med(12,tt,bbind))==rt_med(12,tt,bbind)
            [hh,pp]=ttest2(rt{12,tt,ii},rt{ll,tt,ii});
          else
            [hh,pp]=ttest2(rt{11,tt,ii}+soades(ll),rt{ll,tt,ii});
          end
        end
        rt_msMminshiftuni_pval(ll,tt,bbind)=pp;
        
        if ll==3 % ll runs backwards, so 3 happens after 5 and 7
          rt_7m5(tt,bbind)=rt_med(7,tt,bbind)-rt_med(5,tt,bbind);
          rt_3m5(tt,bbind)=rt_med(3,tt,bbind)-rt_med(5,tt,bbind);
          rt_37m5(tt,bbind)=mean([rt_med(3,tt,bbind) rt_med(7,tt,bbind)])-rt_med(5,tt,bbind);
        end
        
      end
    end % ll
    % 4)whether +/- 500ms is diff from unisensory
    rt_audMll1(ll,tt,bbind)=rt_med(11,tt,bbind)-rt_med(1,tt,bbind);
    rt_tacMll9(ll,tt,bbind)=rt_med(12,tt,bbind)-rt_med(9,tt,bbind);
    [hh,pp]=ttest2(rt{11,tt,ii},rt{1,tt,ii});
    rt_audMll1_pval(ll,tt,bbind)=pp;
    [hh,pp]=ttest2(rt{12,tt,ii},rt{9,tt,ii});
    rt_tacMll9_pval(ll,tt,bbind)=pp;
    
    if ii>9
      xpts=[-.5 -.5 -.07 -.02 0 .02 .07 .5 .5]';
      ypts=rt_med([11 1 3:7 9 12],tt,bbind);
    else
      xpts=[-.5 -.07 -.02 0 .02 .07 .5 ]';
      ypts=rt_med([11 3:7 12],tt,bbind);
    end
    
    try % in case curve-fitting toolbox licence taken
      [fg{tt,bbind},G{tt,bbind}]=fit(xpts,ypts,'poly2');
      figure(ii)
      subplot(1,2,tt-1);
      x=-.5:.01:.5;
      plot(x,fg{tt,bbind}.p1*x.^2+fg{tt,bbind}.p2*x+fg{tt,bbind}.p3);
      hold on;plot(xpts,ypts,'o')
      
      xpts=[-.07 -.02 0 .02 .07]';
      ypts=rt_med([3:7],tt,bbind);
      [fg5{tt,bbind},G5{tt,bbind}]=fit(xpts,ypts,'poly2');
      figure(ii+100)
      subplot(1,2,tt-1);
      x=-.1:.01:.1;
      plot(x,fg5{tt,bbind}.p1*x.^2+fg5{tt,bbind}.p2*x+fg5{tt,bbind}.p3);
      hold on;plot(xpts,ypts,'o'); axis([-.1 .1 .16 .44])
      fg5p1(tt,bbind)=fg5{tt,bbind}.p1;  % curvature
      
      fg5p2(tt,bbind)=fg5{tt,bbind}.p2;  % bias to A or T
      fg5p3(tt,bbind)=fg5{tt,bbind}.p3;  % overall RT
    catch ME
      disp(ME.message)
    end
    
    
  end % tt
  
  %
  %   pcb(bb)=[nanmean(reshape(rtsr(:,6:8),[size(rtsr,1)*3 1]))/nanmean(reshape(rtsr(:,1:2),[size(rtsr,1)*2 1]))]-1;
  %   diffms(bb)=[nanmean(reshape(rtsr(:,6:8),[size(rtsr,1)*3 1]))-nanmean(reshape(rtsr(:,1:2),[size(rtsr,1)*2 1]))];
  bbind=bbind+1;
  clear rtime_touch
end % ii
% save([ddir 'rt_allsubj.mat'],'rt*','fg*','G*');
save([ddir 'rt_allsubj22.mat'],'rt*','fg*','G*','num*');
return

% added Sept 2018; using rt*22.mat
tt=3;
for ll=soalist(1:7)
  if ll<6
    [h,pval(ll),ci{ll},stats{ll}]=ttest( min(squeeze(rt_med(11,tt,:)),squeeze(rt_med(12,tt,:))-soades(ll)) , squeeze(rt_med(ll,tt,:)))
  else
    [h,pval(ll),ci{ll},stats{ll}]=ttest( min(squeeze(rt_med(11,tt,:))+soades(ll),squeeze(rt_med(12,tt,:))) , squeeze(rt_med(ll,tt,:)));
  end
end

% stats{9}.sd/sqrt(22)  for reporting SEM in paper

% percent no response
ll=9;100*[mean(squeeze(numnan(1,ll,:))./numtot(ll,:)') std(squeeze(numnan(1,ll,:))./numtot(ll,:)')/sqrt(22)]
% percent response <100ms 
ll=12;100*[mean(squeeze(numnan(2,ll,:)-numnan(1,ll,:))./numtot(ll,:)') std(squeeze(numnan(2,ll,:)-numnan(1,ll,:))./numtot(ll,:)')/sqrt(22)]
% percent response >1000ms 
ll=1;100*[mean(squeeze(numnan(3,ll,:)-numnan(2,ll,:))./numtot(ll,:)') std(squeeze(numnan(3,ll,:)-numnan(2,ll,:))./numtot(ll,:)')/sqrt(22)]

% total for all 3 reasons, collapsed over all conditions
mean(reshape(squeeze(numnan(3,[1 3 4 5 6 7 9 11 12],:)),[1 9*22])./reshape(squeeze(numtot([1 3 4 5 6 7 9 11 12],:)),[1 9*22]))
std(reshape(squeeze(numnan(3,[1 3 4 5 6 7 9 11 12],:)),[1 9*22])./reshape(squeeze(numtot([1 3 4 5 6 7 9 11 12],:)),[1 9*22]))/sqrt(22)


%% 
% load([ddir 'rt_allsubj.mat'],'rt_med')
load([ddir 'rt_allsubj22.mat'],'rt_med')
tt=3;
numsubj=size(rt_med,3);
subuse=3:numsubj;
numuse=length(subuse);

figure(1);
errorbar([-.5 .5 ],nanmean(squeeze(rt_med([11 12],tt,subuse)),2),nanstd(squeeze(rt_med([11 12],tt,subuse)),[],2)/sqrt(numuse-1),'bo')
hold on;
errorbar([-.5 -.07 -.02 0 .02 .07 .5 ],nanmean(squeeze(rt_med([1 3 4 5 6 7 9],tt,subuse)),2),nanstd(squeeze(rt_med([1 3 4 5 6 7 9],tt,subuse)),[],2)/sqrt(numuse-1),'ro-')
axis([-.8 .8 .22 .32])

figure(2);
h1=errorbar([-.7 .7 ],nanmean(squeeze(rt_med([11 12],tt,subuse)),2),nanstd(squeeze(rt_med([11 12],tt,subuse)),[],2)/sqrt(numuse-1),'bo');
hold on;
h2=bar([-.7 .7 ],nanmean(squeeze(rt_med([11 12],tt,subuse)),2),.15);
hold on;
h3=errorbar([-.5 -.07 -.02 0 .02 .07 .5 ],nanmean(squeeze(rt_med([1 3 4 5 6 7 9],tt,subuse)),2),nanstd(squeeze(rt_med([1 3 4 5 6 7 9],tt,subuse)),[],2)/sqrt(numuse-1),'ro-');
axis([-.9 .9 .22 .32])
set(h3,'linewidth',3)
set(gca,'FontSize',15)
set(gca,'Xtick',[-.7 -.5 -.07 .07 .5 .7])
set(gca,'XTickLabel',{'A' '-500' '-70' '70' '500' 'T'})

print(1,[fdir 'rt_med.png'],'-dpng')
print(2,[fdir 'rt_med_bar.png'],'-dpng')

colorblindD  =[204 101 0]/256;
colorblindApT=[5 165 255]/256;
colorblindMpN=[0 59 179]/256;
colorblindT  =[255 51 166]/256;
colorblindA  =[0 158 115]/256;
colorblindM  =[0 0 0]/256;
colorblindN  =[128 128 128]/256;
% % 
figure(3);
h1=errorbar([-.8 ],nanmean(squeeze(rt_med([11],tt,subuse)),1),nanstd(squeeze(rt_med([11],tt,subuse)),[],1)/sqrt(numuse-1),'o','Color',colorblindA);
hold on;
h1=errorbar([.8 ],nanmean(squeeze(rt_med([12],tt,subuse)),1),nanstd(squeeze(rt_med([12],tt,subuse)),[],1)/sqrt(numuse-1),'o','Color',colorblindT);
hold on;
h2=bar([-.8 ],nanmean(squeeze(rt_med([11],tt,subuse)),1),.15,'FaceColor',colorblindA);
hold on;
h2=bar([.8 ],nanmean(squeeze(rt_med([12],tt,subuse)),1),.15,'FaceColor',colorblindT);
hold on;
% h3=errorbar([-.6 -.35 -.1 0 .1 .35 .6],nanmean(squeeze(rt_med([1 3 4 5 6 7 9],tt,subuse)),2),nanstd(squeeze(rt_med([1 3 4 5 6 7 9],tt,subuse)),[],2)/sqrt(numuse-1),'mo:');
h3=errorbar([-.6 -.35],nanmean(squeeze(rt_med([1 3],tt,subuse)),2),nanstd(squeeze(rt_med([1 3 ],tt,subuse)),[],2)/sqrt(numuse-1),'o:','Color',colorblindM);
h4=errorbar([-.35 -.1 0 .1 .35],nanmean(squeeze(rt_med([ 3 4 5 6 7 ],tt,subuse)),2),nanstd(squeeze(rt_med([3 4 5 6 7],tt,subuse)),[],2)/sqrt(numuse-1),'o-','Color',colorblindM);
h5=errorbar([.35 .6],nanmean(squeeze(rt_med([7 9],tt,subuse)),2),nanstd(squeeze(rt_med([7 9],tt,subuse)),[],2)/sqrt(numuse-1),'o:','Color',colorblindM);
axis([-1 1 .22 .32])
set(h3,'linewidth',3)
set(h4,'linewidth',3)
set(h5,'linewidth',3)
set(gca,'FontSize',15)
set(gca,'Xtick',[-.8 -.6 -.35 -.1 0 .1 .35 .6 .8])
set(gca,'XTickLabel',{'A' '-500' '-70' '-20' '0' '20' '70' '500' 'T'})
xlabel('Unisensory|--  Multisensory Asynchrony (ms)  --|Unisensory')
ylabel('Reaction time (s)')
print(3,[fdir 'rt_med_bar_colXlog.png'],'-dpng')
print(3,[fdir 'rt_med_bar_colXlog.eps'],'-painters','-depsc')


% % 
figure(4);
h3=errorbar([-.6 -.35 -.1 0 .1 .35 .6],nanmean(squeeze(-rt_msMminshiftuni([1 3 4 5 6 7 9],tt,subuse)),2),nanstd(squeeze(rt_msMminshiftuni([1 3 4 5 6 7 9],tt,subuse)),[],2)/sqrt(numuse-1),'mo-');
axis([-.8 .8 -.05 .03])
set(h3,'linewidth',3)
set(gca,'FontSize',15)
set(gca,'Xtick',[ -.6 -.35 -.1 0 .1 .35 .6 ])
set(gca,'XTickLabel',{'-500' '-70' '-20' '0' '20' '70' '500'})
xlabel('Multisensory Asynchrony (ms)')
ylabel('Redundant Target Effect (s)')

print(4,[fdir 'rt_med_Xlog_RTE.png'],'-dpng')



%% Analyse RT
% load([ddir 'rt_allsubj.mat'],'rt_ms*');
load([ddir 'rt_allsubj22.mat'],'rt_ms*');
rtuse=squeeze(rt_msMminshiftuni(:,3,:));
[pa,ptable]=anova_rm(rtuse([1 3 4 5 6 7 9],:)'); % p=0000
[pa70,ptable70]=anova_rm(rtuse([3 4 5 6 7],:)'); % p=0.0103
[pa20,ptable20]=anova_rm(rtuse([4 5 6],:)'); % p=0.1657

p=nan(7,7);
[h,p(3,4)]=ttest(rtuse(3,:),rtuse(4,:))
[h,p(3,5)]=ttest(rtuse(3,:),rtuse(5,:))
[h,p(3,6)]=ttest(rtuse(3,:),rtuse(6,:))
[h,p(3,7)]=ttest(rtuse(3,:),rtuse(7,:))

[h,p(4,5)]=ttest(rtuse(4,:),rtuse(5,:))
[h,p(4,6)]=ttest(rtuse(4,:),rtuse(6,:))
[h,p(4,7)]=ttest(rtuse(4,:),rtuse(7,:))

[h,p(5,6)]=ttest(rtuse(5,:),rtuse(6,:))
[h,p(5,7)]=ttest(rtuse(5,:),rtuse(7,:))

[h,p(6,7)]=ttest(rtuse(6,:),rtuse(7,:))


% old
load([ddir 'rt_allsubj.mat'],'rt_med');
tt=3;

anova_rm(squeeze(rt_med([1 3 4 5 6 7 9],tt,:))')


p=nan(12)
[h,p(1,3)]=ttest(rt_med(1,tt,:),rt_med(3,tt,:))
[h,p(1,4)]=ttest(rt_med(1,tt,:),rt_med(4,tt,:))
[h,p(1,5)]=ttest(rt_med(1,tt,:),rt_med(5,tt,:))
[h,p(1,6)]=ttest(rt_med(1,tt,:),rt_med(6,tt,:))
[h,p(1,7)]=ttest(rt_med(1,tt,:),rt_med(7,tt,:))
[h,p(1,9)]=ttest(rt_med(1,tt,:),rt_med(9,tt,:))

[h,p(3,4)]=ttest(rt_med(3,tt,:),rt_med(4,tt,:))
[h,p(3,5)]=ttest(rt_med(3,tt,:),rt_med(5,tt,:))
[h,p(3,6)]=ttest(rt_med(3,tt,:),rt_med(6,tt,:))
[h,p(3,7)]=ttest(rt_med(3,tt,:),rt_med(7,tt,:))
[h,p(3,9)]=ttest(rt_med(3,tt,:),rt_med(9,tt,:))

[h,p(4,5)]=ttest(rt_med(4,tt,:),rt_med(5,tt,:))
[h,p(4,6)]=ttest(rt_med(4,tt,:),rt_med(6,tt,:))
[h,p(4,7)]=ttest(rt_med(4,tt,:),rt_med(7,tt,:))
[h,p(4,9)]=ttest(rt_med(4,tt,:),rt_med(9,tt,:))

[h,p(5,6)]=ttest(rt_med(5,tt,:),rt_med(6,tt,:))
[h,p(5,7)]=ttest(rt_med(5,tt,:),rt_med(7,tt,:))
[h,p(5,9)]=ttest(rt_med(5,tt,:),rt_med(9,tt,:))

[h,p(6,7)]=ttest(rt_med(6,tt,:),rt_med(7,tt,:))
[h,p(6,9)]=ttest(rt_med(6,tt,:),rt_med(9,tt,:))

[h,p(7,9)]=ttest(rt_med(7,tt,:),rt_med(9,tt,:))

[h,p(11,1)]=ttest(rt_med(11,tt,:),rt_med(1,tt,:)) % Compare AT500 to Aalone
[h,p(11,3)]=ttest(rt_med(11,tt,:),rt_med(3,tt,:))
[h,p(11,4)]=ttest(rt_med(11,tt,:),rt_med(4,tt,:))
[h,p(11,5)]=ttest(rt_med(11,tt,:),rt_med(5,tt,:))
[h,p(11,6)]=ttest(rt_med(11,tt,:),rt_med(6,tt,:))
[h,p(11,7)]=ttest(rt_med(11,tt,:),rt_med(7,tt,:))
[h,p(11,9)]=ttest(rt_med(11,tt,:),rt_med(9,tt,:))

[h,p(12,1)]=ttest(rt_med(12,tt,:),rt_med(1,tt,:)) % Compare AT500 to Aalone
[h,p(12,3)]=ttest(rt_med(12,tt,:),rt_med(3,tt,:))
[h,p(12,4)]=ttest(rt_med(12,tt,:),rt_med(4,tt,:))
[h,p(12,5)]=ttest(rt_med(12,tt,:),rt_med(5,tt,:))
[h,p(12,6)]=ttest(rt_med(12,tt,:),rt_med(6,tt,:))
[h,p(12,7)]=ttest(rt_med(12,tt,:),rt_med(7,tt,:))
[h,p(12,9)]=ttest(rt_med(12,tt,:),rt_med(9,tt,:))

anova_rm(squeeze(rt_med([3 4 5 6 7],tt,:))')


%% 

% nsub=length(subuse);
% Y=reshape(squeeze(rt_med([1 3 4 6 7 9],3,subuse))',[nsub*6 1]);
% S=repmat([1:nsub]',[6 1]);
% F1=[ones(nsub,1); ones(nsub,1); ones(nsub,1); 2*ones(nsub,1); 2*ones(nsub,1); 2*ones(nsub,1)]; % A leading; T leading
% F2=[1*ones(nsub,1); 2*ones(nsub,1); 3*ones(nsub,1); 3*ones(nsub,1); 2*ones(nsub,1); 1*ones(nsub,1)]; % 500, 70, 20
% NAMES={'Sense leading', 'Asyncrhony'};
% stats = rm_anova2(Y,S,F1,F2,NAMES)
% % Yes, interaction significant (p=0.005) as well as Asychnrony (p<0.0001), not sense-leading (p>0.1)
% % This permits to look at asynchrony with rm-1-way-anova for sense-leading, 
% % AND
% % sense-leading with rm-1-way-anova for asychrony
% 
% % Sense-leading plays a role for 500 only
% p=anova_rm(squeeze(rt_med([1 9],3,subuse))','off'); % p=0.12
% p=anova_rm(squeeze(rt_med([3 7],3,subuse))','off'); % p=0.0052
% p=anova_rm(squeeze(rt_med([4 6],3,subuse))','off'); % p=0.0099
% 
% % Asynchrony affects RT, for each sense-leading separately
% p=anova_rm(squeeze(rt_med([1 3 4],3,subuse))','off'); % p<0.0001
% p=anova_rm(squeeze(rt_med([6 7 9],3,subuse))','off'); % p<0.0001


nsub=length(subuse);
Y=reshape(squeeze(-rt_msMminshiftuni([1 3 4 6 7 9],3,subuse))',[nsub*6 1]);
S=repmat([1:nsub]',[6 1]);
F1=[ones(nsub,1); ones(nsub,1); ones(nsub,1); 2*ones(nsub,1); 2*ones(nsub,1); 2*ones(nsub,1)]; % A leading; T leading
F2=[1*ones(nsub,1); 2*ones(nsub,1); 3*ones(nsub,1); 3*ones(nsub,1); 2*ones(nsub,1); 1*ones(nsub,1)]; % 500, 70, 20
NAMES={'Sense leading', 'Asynchrony'};
stats = rm_anova2(Y,S,F1,F2,NAMES)
% Yes, interaction significant (p=0.0006) as well as Asychnrony (p<0.0001), not sense-leading (p>0.1)
% This permits to look at asynchrony with rm-1-way-anova for sense-leading, 
% AND
% sense-leading with rm-1-way-anova for asychrony

% Sense-leading plays a role for 500 only
p=anova_rm(squeeze(-rt_msMminshiftuni([1 9],3,subuse))','off') % p=0.005
p=anova_rm(squeeze(-rt_msMminshiftuni([3 7],3,subuse))','off') % p=0.24
p=anova_rm(squeeze(-rt_msMminshiftuni([4 6],3,subuse))','off') % p=0.31

% Asynchrony affects RT, for each sense-leading separately
p=anova_rm(squeeze(-rt_msMminshiftuni([1 3 4],3,subuse))','off') % p<0.0001
p=anova_rm(squeeze(-rt_msMminshiftuni([6 7 9],3,subuse))','off') % p<0.0001



%% Corr with ERP sensor

printflag=0;
plotflag=0;
statflag=1;
sleep=0;
ss=10;
iterflag=1;

if sleep
else
  iteruse=27;
end

load([ddir 'rt_allsubj.mat']);
alphaval=.001;
load eeg1010_neighb
for iter=iteruse
  
  if iterflag
    load([edir 'tlock_grind_sleep' num2str(sleep) '_iter' num2str(iter) '.mat'],'grind*');
  else
    if iter==1
      load([edir 'tlock_statmc_sleep' num2str(sleep) '.mat'],'grind*');
    else
      break
    end
  end
  
  for tt=[3]
    for ll=[1 3 4 5 6 7 9]
      if ll==1 || ll==3 || ll==4 || ll==5
        stattimwin=[.07 .45];
      elseif ll==6
        stattimwin=[.09 .47];
      elseif ll==7
        stattimwin=[.14 .52];
      elseif ll==9
        stattimwin=[.57 .95];
      end
      cortimwin=[grind_tacMSpN_save{ll,tt,ss}.time(1) grind_tacMSpN_save{ll,tt,ss}.time(end)];
      %     try
      cor_rtmed_TPA_MSPN{ll,tt,ss}.time=grind_tacMSpN_save{ll,tt,ss}.time;
      cor_rtmed_TPA_MSPN{ll,tt,ss}.label=grind_tacMSpN_save{ll,tt,ss}.label;
      cor_rtmed_tacMSpN{ll,tt,ss}.time=grind_tacMSpN_save{ll,tt,ss}.time;
      cor_rtmed_tacMSpN{ll,tt,ss}.label=grind_tacMSpN_save{ll,tt,ss}.label;
      cor_rtmed_tacPaud{ll,tt,ss}.time=grind_tacMSpN_save{ll,tt,ss}.time;
      cor_rtmed_tacPaud{ll,tt,ss}.label=grind_tacMSpN_save{ll,tt,ss}.label;
      cor_rtrel_TPA_MSPN{ll,tt,ss}.time=grind_tacMSpN_save{ll,tt,ss}.time;
      cor_rtrel_TPA_MSPN{ll,tt,ss}.label=grind_tacMSpN_save{ll,tt,ss}.label;
      cor_rtrel_tacMSpN{ll,tt,ss}.time=grind_tacMSpN_save{ll,tt,ss}.time;
      cor_rtrel_tacMSpN{ll,tt,ss}.label=grind_tacMSpN_save{ll,tt,ss}.label;
      cor_rtrel_tacPaud{ll,tt,ss}.time=grind_tacMSpN_save{ll,tt,ss}.time;
      cor_rtrel_tacPaud{ll,tt,ss}.label=grind_tacMSpN_save{ll,tt,ss}.label;
      %     catch
      %       cor_rtmed_TPA_MSPN{ll,tt,ss}.time=grind_tacMSpN.time;
      %       cor_rtmed_TPA_MSPN{ll,tt,ss}.label=grind_tacMSpN.label;
      %       cor_rtmed_tacMSpN{ll,tt,ss}.time=grind_tacMSpN.time;
      %       cor_rtmed_tacMSpN{ll,tt,ss}.label=grind_tacMSpN.label;
      %       cor_rtmed_tacPaud{ll,tt,ss}.time=grind_tacMSpN.time;
      %       cor_rtmed_tacPaud{ll,tt,ss}.label=grind_tacMSpN.label;
      %     end
      cor_rtmed_TPA_MSPN{ll,tt,ss}.dimord='chan_time';
      cor_rtmed_tacMSpN{ll,tt,ss}.dimord='chan_time';
      cor_rtmed_tacPaud{ll,tt,ss}.dimord='chan_time';
      cor_rtrel_TPA_MSPN{ll,tt,ss}.dimord='chan_time';
      cor_rtrel_tacMSpN{ll,tt,ss}.dimord='chan_time';
      cor_rtrel_tacPaud{ll,tt,ss}.dimord='chan_time';
      timevect=dsearchn(cor_rtmed_TPA_MSPN{ll,tt,ss}.time',cortimwin(1)):dsearchn(cor_rtmed_TPA_MSPN{ll,tt,ss}.time',cortimwin(2));
      %     try
      numchan=size(grind_tacMSpN_save{ll,tt,ss}.individual,2);
      %     catch
      %       numchan=size(grind_tacMSpN.individual,2);
      %     end
      cor_rtmed_TPA_MSPN{ll,tt,ss}.avg=nan(numchan,length(timevect));
      cor_rtmed_tacMSpN{ll,tt,ss}.avg=nan(numchan,length(timevect));
      cor_rtmed_tacPaud{ll,tt,ss}.avg=nan(numchan,length(timevect));
      cor_rtmed_TPA_MSPN{ll,tt,ss}.pval=nan(numchan,length(timevect));
      cor_rtmed_tacMSpN{ll,tt,ss}.pval=nan(numchan,length(timevect));
      cor_rtmed_tacPaud{ll,tt,ss}.pval=nan(numchan,length(timevect));
      cor_rtmed_TPA_MSPN{ll,tt,ss}.time=cor_rtmed_TPA_MSPN{ll,tt,ss}.time(timevect);
      cor_rtmed_tacMSpN{ll,tt,ss}.time=cor_rtmed_tacMSpN{ll,tt,ss}.time(timevect);
      cor_rtmed_tacPaud{ll,tt,ss}.time=cor_rtmed_tacPaud{ll,tt,ss}.time(timevect);
      cor_rtrel_TPA_MSPN{ll,tt,ss}.avg=nan(numchan,length(timevect));
      cor_rtrel_tacMSpN{ll,tt,ss}.avg=nan(numchan,length(timevect));
      cor_rtrel_tacPaud{ll,tt,ss}.avg=nan(numchan,length(timevect));
      cor_rtrel_TPA_MSPN{ll,tt,ss}.pval=nan(numchan,length(timevect));
      cor_rtrel_tacMSpN{ll,tt,ss}.pval=nan(numchan,length(timevect));
      cor_rtrel_tacPaud{ll,tt,ss}.pval=nan(numchan,length(timevect));
      cor_rtrel_TPA_MSPN{ll,tt,ss}.time=cor_rtrel_TPA_MSPN{ll,tt,ss}.time(timevect);
      cor_rtrel_tacMSpN{ll,tt,ss}.time=cor_rtrel_tacMSpN{ll,tt,ss}.time(timevect);
      cor_rtrel_tacPaud{ll,tt,ss}.time=cor_rtrel_tacPaud{ll,tt,ss}.time(timevect);
      if ll==1 || ll==9
        subuse=3:size(rt_med,3);
      else
        subuse=1:size(rt_med,3);
      end
      nsub=length(subuse);
      %     try
      grind_tacMSpN_save{ll,tt,ss}.individual=grind_tacMSpN_save{ll,tt,ss}.individual(subuse,:,:);
      grind_tacPaud_save{ll,tt,ss}.individual=grind_tacPaud_save{ll,tt,ss}.individual(subuse,:,:);
      %     catch
      %       grind_tacMSpN.individual=grind_tacMSpN.individual(subuse,:,:);
      %       grind_tacPaud.individual=grind_tacPaud.individual(subuse,:,:);
      %     end
      
      cfg=[];
      cfg.operation='subtract';
      cfg.parameter='individual';
      grind_TPA_MSPN{ll,tt,ss}=ft_math(cfg,grind_tacPaud_save{ll,tt,ss},grind_tacMSpN_save{ll,tt,ss});
      
      for cc=1:numchan
        for time=timevect
          %         try
          [cor_rtmed_TPA_MSPN{ll,tt,ss}.avg(cc,time-min(timevect)+1),cor_rtmed_TPA_MSPN{ll,tt,ss}.pval(cc,time-min(timevect)+1)]=corr(squeeze(rt_med(ll,tt,subuse)),grind_TPA_MSPN{ll,tt,ss}.individual(:,cc,time));
          [cor_rtmed_tacMSpN{ll,tt,ss}.avg(cc,time-min(timevect)+1),cor_rtmed_tacMSpN{ll,tt,ss}.pval(cc,time-min(timevect)+1)]=corr(squeeze(rt_med(ll,tt,subuse)),grind_tacMSpN_save{ll,tt,ss}.individual(:,cc,time));
          [cor_rtmed_tacPaud{ll,tt,ss}.avg(cc,time-min(timevect)+1),cor_rtmed_tacPaud{ll,tt,ss}.pval(cc,time-min(timevect)+1)]=corr(squeeze(rt_med(ll,tt,subuse)),grind_tacPaud_save{ll,tt,ss}.individual(:,cc,time));
          [cor_rtrel_TPA_MSPN{ll,tt,ss}.avg(cc,time-min(timevect)+1),cor_rtrel_TPA_MSPN{ll,tt,ss}.pval(cc,time-min(timevect)+1)]=corr(squeeze(rt_msMreluni(ll,tt,subuse)),grind_TPA_MSPN{ll,tt,ss}.individual(:,cc,time));
          [cor_rtrel_tacMSpN{ll,tt,ss}.avg(cc,time-min(timevect)+1),cor_rtrel_tacMSpN{ll,tt,ss}.pval(cc,time-min(timevect)+1)]=corr(squeeze(rt_msMreluni(ll,tt,subuse)),grind_tacMSpN_save{ll,tt,ss}.individual(:,cc,time));
          [cor_rtrel_tacPaud{ll,tt,ss}.avg(cc,time-min(timevect)+1),cor_rtrel_tacPaud{ll,tt,ss}.pval(cc,time-min(timevect)+1)]=corr(squeeze(rt_msMreluni(ll,tt,subuse)),grind_tacPaud_save{ll,tt,ss}.individual(:,cc,time));
          %         catch
          %           [cor_rtmed_TPA_MSPN{ll,tt,ss}.avg(cc,time-min(timevect)+1),cor_rtmed_TPA_MSPN{ll,tt,ss}.pval(cc,time-min(timevect)+1)]=corr(squeeze(rt_med(ll,tt,subuse)),grind_TPA_MSPN{ll,tt,ss}.individual(:,cc,time));
          %           [cor_rtmed_tacMSpN{ll,tt,ss}.avg(cc,time-min(timevect)+1),cor_rtmed_tacMSpN{ll,tt,ss}.pval(cc,time-min(timevect)+1)]=corr(squeeze(rt_med(ll,tt,subuse)),grind_tacMSpN.individual(:,cc,time));
          %           [cor_rtmed_tacPaud{ll,tt,ss}.avg(cc,time-min(timevect)+1),cor_rtmed_tacPaud{ll,tt,ss}.pval(cc,time-min(timevect)+1)]=corr(squeeze(rt_med(ll,tt,subuse)),grind_tacPaud.individual(:,cc,time));
          %         end
        end
      end
      cor_rtmed_TPA_MSPN{ll,tt,ss}.mask=cor_rtmed_TPA_MSPN{ll,tt,ss}.pval<alphaval;
      cor_rtmed_tacMSpN{ll,tt,ss}.mask=cor_rtmed_tacMSpN{ll,tt,ss}.pval<alphaval;
      cor_rtmed_tacPaud{ll,tt,ss}.mask=cor_rtmed_tacPaud{ll,tt,ss}.pval<alphaval;
      cor_rtrel_TPA_MSPN{ll,tt,ss}.mask=cor_rtrel_TPA_MSPN{ll,tt,ss}.pval<alphaval;
      cor_rtrel_tacMSpN{ll,tt,ss}.mask=cor_rtrel_tacMSpN{ll,tt,ss}.pval<alphaval;
      cor_rtrel_tacPaud{ll,tt,ss}.mask=cor_rtrel_tacPaud{ll,tt,ss}.pval<alphaval;
      
      if plotflag
        figure(ll);
        cfg=[];
        cfg.maskparameter='mask';
        cfg.layout='elec1010.lay';
        %     cfg.xlim=[cor_rtmed_tacMSpN{ll,tt,ss}.time(1) cor_rtmed_tacMSpN{ll,tt,ss}.time(end)];
        cfg.xlim=[cortimwin(1) cortimwin(2)];
        %     ft_multiplotER(cfg,cor_rtmed_tacMSpN{ll,tt,ss});
        %     ft_multiplotER(cfg,cor_rtmed_tacPaud{ll,tt,ss});
        ft_multiplotER(cfg,cor_rtmed_TPA_MSPN{ll,tt,ss});
        figure(ll+10);
        ft_multiplotER(cfg,cor_rtrel_TPA_MSPN{ll,tt,ss});
      end
      
      % Question: way to find sig-diff of corrs between conditions?  or just if
      % cor_rtmed_TPA_MSPN is sig?
      
      
      
      if statflag
        %     try
        rtmed_grind{ll,tt,ss}=grind_tacMSpN_save{ll,tt,ss};
        rtrel_grind{ll,tt,ss}=grind_tacMSpN_save{ll,tt,ss};
        rtrelshift_grind{ll,tt,ss}=grind_tacMSpN_save{ll,tt,ss};
        %     catch
        %       rtmed_grind{ll,tt,ss}=grind_tacMSpN;
        %     end
        rtmed_grind{ll,tt,ss}.individual=repmat(squeeze(rt_med(ll,tt,subuse)),[1 size(rtmed_grind{ll,tt,ss}.individual,2) size(rtmed_grind{ll,tt,ss}.individual,3)]);
        rtrel_grind{ll,tt,ss}.individual=repmat(squeeze(rt_msMreluni(ll,tt,subuse)),[1 size(rtrel_grind{ll,tt,ss}.individual,2) size(rtrel_grind{ll,tt,ss}.individual,3)]);
        rtrelshift_grind{ll,tt,ss}.individual=repmat(squeeze(rt_msMminshiftuni(ll,tt,subuse)),[1 size(rtrel_grind{ll,tt,ss}.individual,2) size(rtrel_grind{ll,tt,ss}.individual,3)]);
        
        %     rtmed_tacPaud{ll,tt,ss}=grind_tacPaud_save{ll,tt,ss};
        %     rtmed_tacPaud{ll,tt,ss}.individual=repmat(squeeze(rt_med(ll,tt,:)),[1 size(grind_tacPaud_save{ll,tt,ss}.individual,2) size(grind_tacPaud_save{ll,tt,ss}.individual,3)]);
        %     rtmed_TPAMSPN{ll,tt,ss}=grind_tacMSpN_save{ll,tt,ss};
        %     rtmed_TPAMSPN{ll,tt,ss}.individual=repmat(squeeze(rt_med(ll,tt,:)),[1 size(grind_tacMSpN_save{ll,tt,ss}.individual,2) size(grind_tacMSpN_save{ll,tt,ss}.individual,3)]);
        
        cfg=[];
        cfg.latency=stattimwin;
        %         cfg.channel=chanuse;
        cfg.neighbours=neighbours;
        %         cfg.parameter='avg';
        cfg.parameter='individual';
        cfg.method='montecarlo';
        % cfg.method='analytic';
        cfg.numrandomization=100;
        % cfg.correctm='holm';
        cfg.correctm='cluster';
        cfg.clusteralpha = 0.05;
        cfg.clusterstatistic = 'maxsum';
        cfg.minnbchan = 2;
        cfg.statistic='correlationT';
        % cfg.statistic='indepsamplesT';
        cfg.design=zeros(2,2*nsub);
        cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)]; % must be these numbers to code the labels for 'correlationT'
        cfg.design(2,:)=[1:nsub 1:nsub];
        cfg.ivar=1;
        cfg.uvar=2;
        %     try
        stattcormed_tacMSpN{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tacMSpN_save{ll,tt,ss}, rtmed_grind{ll,tt,ss});
        stattcormed_tacPaud{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tacPaud_save{ll,tt,ss}, rtmed_grind{ll,tt,ss});
        stattcormed_TPA_MSPN{ll,tt,ss}=ft_timelockstatistics(cfg, grind_TPA_MSPN{ll,tt,ss}, rtmed_grind{ll,tt,ss});
        stattcorrel_tacMSpN{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tacMSpN_save{ll,tt,ss}, rtrel_grind{ll,tt,ss});
        stattcorrel_tacPaud{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tacPaud_save{ll,tt,ss}, rtrel_grind{ll,tt,ss});
        stattcorrel_TPA_MSPN{ll,tt,ss}=ft_timelockstatistics(cfg, grind_TPA_MSPN{ll,tt,ss}, rtrel_grind{ll,tt,ss});
        stattcorrelshift_TPA_MSPN{ll,tt,ss}=ft_timelockstatistics(cfg, grind_TPA_MSPN{ll,tt,ss}, rtrelshift_grind{ll,tt,ss});
        %     catch
        %       stattcor_tacMSpN{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tacMSpN, rtmed_grind{ll,tt,ss});
        %       stattcor_tacPaud{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tacPaud, rtmed_grind{ll,tt,ss});
        %       stattcor_TPA_MSPN{ll,tt,ss}=ft_timelockstatistics(cfg, grind_TPA_MSPN{ll,tt,ss}, rtmed_grind{ll,tt,ss});
        %     end
        
        %     % "depsamplesregrT does only test for linear increases or decreases in your data, depending on the sequence in which you input your data conditions.
        %     %  ... For example, freqstatistics(cfg, cond1{:}, cond2{:}, cond3{:}) will result in positive clusters for data increasing (quasi-)linearly from cond1 to cond3, and negative clusters for data behaving vice versa."
        % %     cfg.statistic='depsamplesregrT';
        % %     cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)]; % exact numbers here don't matter
        % %     try
        % %     stattIRT_tacMSpN{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tacMSpN_save{ll,tt,ss}, rtmed_tacMSpN{ll,tt,ss});
        % %     catch
        % %     stattIRT_tacMSpN{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tacMSpN, rtmed_tacMSpN{ll,tt,ss});
        % %     end
        %
        % %perhaps for within-subject correlations over trials
        % %     cfg.statistic='indepsamplesregrT';
        % %     cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)]; % exact numbers here don't matter
        % %     try
        % %     stattIRT_tacMSpN{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tacMSpN_save{ll,tt,ss}, rtmed_tacMSpN{ll,tt,ss});
        % %     catch
        % %     stattIRT_tacMSpN{ll,tt,ss}=ft_timelockstatistics(cfg, grind_tacMSpN, rtmed_tacMSpN{ll,tt,ss});
        % %     end
        
        %     try
        grave_tacMSpN{ll,tt,ss}=grind_tacMSpN_save{ll,tt,ss};
        grave_tacMSpN{ll,tt,ss}.avg=squeeze(mean(grind_tacMSpN_save{ll,tt,ss}.individual,1));
        grave_tacMSpN{ll,tt,ss}=rmfield(grave_tacMSpN{ll,tt,ss},'individual');
        grave_tacMSpN{ll,tt,ss}.dimord='chan_time';
        grave_tacPaud{ll,tt,ss}=grind_tacPaud_save{ll,tt,ss};
        grave_tacPaud{ll,tt,ss}.avg=squeeze(mean(grind_tacPaud_save{ll,tt,ss}.individual,1));
        grave_tacPaud{ll,tt,ss}=rmfield(grave_tacPaud{ll,tt,ss},'individual');
        grave_tacPaud{ll,tt,ss}.dimord='chan_time';
        grave_TPA_MSPN{ll,tt,ss}=grind_TPA_MSPN{ll,tt,ss};
        grave_TPA_MSPN{ll,tt,ss}.avg=squeeze(mean(grind_TPA_MSPN{ll,tt,ss}.individual,1));
        grave_TPA_MSPN{ll,tt,ss}=rmfield(grave_TPA_MSPN{ll,tt,ss},'individual');
        grave_TPA_MSPN{ll,tt,ss}.dimord='chan_time';
        %     catch
        %       grave_tacMSpN{ll,tt,ss}=grind_tacMSpN;
        %       grave_tacMSpN{ll,tt,ss}.avg=squeeze(mean(grind_tacMSpN.individual,1));
        %       grave_tacMSpN{ll,tt,ss}=rmfield(grave_tacMSpN{ll,tt,ss},'individual');
        %       grave_tacMSpN{ll,tt,ss}.dimord='chan_time';
        %       grave_tacPaud{ll,tt,ss}=grind_tacPaud;
        %       grave_tacPaud{ll,tt,ss}.avg=squeeze(mean(grind_tacPaud.individual,1));
        %       grave_tacPaud{ll,tt,ss}=rmfield(grave_tacPaud{ll,tt,ss},'individual');
        %       grave_tacPaud{ll,tt,ss}.dimord='chan_time';
        %       grave_TPA_MSPN{ll,tt,ss}=grind_TPA_MSPN{ll,tt,ss};
        %       grave_TPA_MSPN{ll,tt,ss}.avg=squeeze(mean(grind_TPA_MSPN{ll,tt,ss}.individual,1));
        %       grave_TPA_MSPN{ll,tt,ss}=rmfield(grave_TPA_MSPN{ll,tt,ss},'individual');
        %       grave_TPA_MSPN{ll,tt,ss}.dimord='chan_time';
        %     end
        
        if plotflag
          try close 101;end
          topoplot_highlight(101,grave_tacMSpN{ll,tt,ss},stattimwin,stattcormed_tacMSpN{ll,tt,ss});
          try close 102;end
          topoplot_highlight(102,cor_rtmed_tacMSpN{ll,tt,ss},stattimwin,stattcormed_tacMSpN{ll,tt,ss},[-1 1]);
          try close 103;end
          topoplot_highlight(103,grave_tacPaud{ll,tt,ss},stattimwin,stattcormed_tacPaud{ll,tt,ss});
          try close 104;end
          topoplot_highlight(104,cor_rtmed_tacPaud{ll,tt,ss},stattimwin,stattcormed_tacPaud{ll,tt,ss},[-1 1]);
          try close 105;end
          topoplot_highlight(105,grave_TPA_MSPN{ll,tt,ss},stattimwin,stattcormed_TPA_MSPN{ll,tt,ss});
          try close 106;end
          topoplot_highlight(106,cor_rtmed_TPA_MSPN{ll,tt,ss},stattimwin,stattcormed_TPA_MSPN{ll,tt,ss},[-1 1]);
          try close 111;end
          topoplot_highlight(111,grave_tacMSpN{ll,tt,ss},stattimwin,stattcorrel_tacMSpN{ll,tt,ss});
          try close 112;end
          topoplot_highlight(112,cor_rtrel_tacMSpN{ll,tt,ss},stattimwin,stattcorrel_tacMSpN{ll,tt,ss},[-1 1]);
          try close 113;end
          topoplot_highlight(113,grave_tacPaud{ll,tt,ss},stattimwin,stattcorrel_tacPaud{ll,tt,ss});
          try close 114;end
          topoplot_highlight(114,cor_rtrel_tacPaud{ll,tt,ss},stattimwin,stattcorrel_tacPaud{ll,tt,ss},[-1 1]);
          try close 115;end
          topoplot_highlight(115,grave_TPA_MSPN{ll,tt,ss},stattimwin,stattcorrel_TPA_MSPN{ll,tt,ss});
          try close 116;end
          topoplot_highlight(116,cor_rtrel_TPA_MSPN{ll,tt,ss},stattimwin,stattcorrel_TPA_MSPN{ll,tt,ss},[-1 1]);
          
          if printflag
            print(101,[fdir 'corstatmed_grave_tacMSpN_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(103,[fdir 'corstatmed_grave_tacPaud_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(105,[fdir 'corstatmed_grave_TPA_MSPN_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(102,[fdir 'corstatmed_cor_tacMSpN_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(104,[fdir 'corstatmed_cor_tacPaud_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(106,[fdir 'corstatmed_cor_TPA_MSPN_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(111,[fdir 'corstatrel_grave_tacMSpN_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(113,[fdir 'corstatrel_grave_tacPaud_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(115,[fdir 'corstatrel_grave_TPA_MSPN_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(112,[fdir 'corstatrel_cor_tacMSpN_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(114,[fdir 'corstatrel_cor_tacPaud_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(116,[fdir 'corstatrel_cor_TPA_MSPN_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            
          end % printflag
          
        end % plotflag
        save([edir 'cor_stat_sleep' num2str(sleep) '_iter' num2str(iter) '.mat'],'stattcor*');
      end % statflag
      
    end % ll
    save([edir 'cor_cor_sleep' num2str(sleep) '_iter' num2str(iter) '.mat'],'cor*');
  end % tt
end % iter

%% Corr RT (and quad fit) with +/-70ms P2 ERP effect

printflag=0;
plotflag=0;
statflag=1;
sleep=0;
ss=10;
iterflag=1;

load([ddir 'rt_allsubj.mat']);
alphaval=.001;
load eeg1010_neighb
for iter=27
  
  if iterflag
    load([edir 'tlock_grind_sleep' num2str(sleep) '_iter' num2str(iter) '.mat'],'grind*');
  else
    if iter==1
      load([edir 'tlock_statmc_sleep' num2str(sleep) '.mat'],'grind*');
    else
      break
    end
  end
  
  for tt=[3]
    subuse=1:size(rt_med,3);
    nsub=length(subuse);
    for ll=[3 7 37] % only testing edge of window of integration
      if ll==3 || ll==7
        grind_tacMSpN_save{ll,tt,ss}.individual=grind_tacMSpN_save{ll,tt,ss}.individual(subuse,:,:);
        grind_tacPaud_save{ll,tt,ss}.individual=grind_tacPaud_save{ll,tt,ss}.individual(subuse,:,:);
        cfg=[];
        cfg.operation='subtract';
        cfg.parameter='individual';
        grind_TPA_MSPN{ll,tt,ss}=ft_math(cfg,grind_tacPaud_save{ll,tt,ss},grind_tacMSpN_save{ll,tt,ss})
      elseif ll==37
        cfg=[];
        cfg.operation='(x1+x2)/2';
        cfg.parameter='individual';
        grind_TPA_MSPN{ll,tt,ss}=ft_math(cfg,grind_TPA_MSPN{3,tt,ss},grind_TPA_MSPN{7,tt,ss})
      end
      
      % initialising variables
      cortimwin=[grind_TPA_MSPN{ll,tt,ss}.time(1) grind_TPA_MSPN{ll,tt,ss}.time(end)];
      cor_rtrel_TPA_MSPN{ll,tt,ss}.time=grind_TPA_MSPN{ll,tt,ss}.time;
      cor_rtrel_TPA_MSPN{ll,tt,ss}.label=grind_TPA_MSPN{ll,tt,ss}.label;
      cor_rtrel_TPA_MSPN{ll,tt,ss}.dimord='chan_time';
      timevect=dsearchn(cor_rtrel_TPA_MSPN{ll,tt,ss}.time',cortimwin(1)):dsearchn(cor_rtrel_TPA_MSPN{ll,tt,ss}.time',cortimwin(2));
      numchan=size(grind_TPA_MSPN{ll,tt,ss}.individual,2);
      cor_rtrel_TPA_MSPN{ll,tt,ss}.avg=nan(numchan,length(timevect));
      cor_rtrel_TPA_MSPN{ll,tt,ss}.pval=nan(numchan,length(timevect));
      cor_rtrel_TPA_MSPN{ll,tt,ss}.time=cor_rtrel_TPA_MSPN{ll,tt,ss}.time(timevect);
      cor_rtfg5p1_TPA_MSPN{ll,tt,ss}=cor_rtrel_TPA_MSPN{ll,tt,ss};
      cor_rtfg5p2_TPA_MSPN{ll,tt,ss}=cor_rtrel_TPA_MSPN{ll,tt,ss};
      cor_rtfg5p3_TPA_MSPN{ll,tt,ss}=cor_rtrel_TPA_MSPN{ll,tt,ss};
      cor_rtNm5_TPA_MSPN{ll,tt,ss}=cor_rtrel_TPA_MSPN{ll,tt,ss};
      
      
      for cc=1:numchan
        for time=timevect
          [cor_rtfg5p1_TPA_MSPN{ll,tt,ss}.avg(cc,time-min(timevect)+1),cor_rtfg5p1_TPA_MSPN{ll,tt,ss}.pval(cc,time-min(timevect)+1)]=corr(squeeze(fg5p1(tt,subuse))',         grind_TPA_MSPN{ll,tt,ss}.individual(:,cc,time));
          [cor_rtfg5p2_TPA_MSPN{ll,tt,ss}.avg(cc,time-min(timevect)+1),cor_rtfg5p3_TPA_MSPN{ll,tt,ss}.pval(cc,time-min(timevect)+1)]=corr(squeeze(fg5p2(tt,subuse))',         grind_TPA_MSPN{ll,tt,ss}.individual(:,cc,time));
          [cor_rtfg5p3_TPA_MSPN{ll,tt,ss}.avg(cc,time-min(timevect)+1),cor_rtfg5p3_TPA_MSPN{ll,tt,ss}.pval(cc,time-min(timevect)+1)]=corr(squeeze(fg5p3(tt,subuse))',         grind_TPA_MSPN{ll,tt,ss}.individual(:,cc,time));
          if ll==3
            [cor_rtrel_TPA_MSPN{ll,tt,ss}.avg(cc,time-min(timevect)+1),  cor_rtrel_TPA_MSPN{ll,tt,ss}.pval(cc,time-min(timevect)+1)]  =corr(squeeze(rt_msMreluni(ll,tt,subuse)),grind_TPA_MSPN{ll,tt,ss}.individual(:,cc,time));
            [cor_rtNm5_TPA_MSPN{ll,tt,ss}.avg(cc,time-min(timevect)+1),cor_rtNm5_TPA_MSPN{ll,tt,ss}.pval(cc,time-min(timevect)+1)]  =corr(squeeze(rt_3m5(tt,subuse))',        grind_TPA_MSPN{ll,tt,ss}.individual(:,cc,time));
          elseif ll==7
            [cor_rtrel_TPA_MSPN{ll,tt,ss}.avg(cc,time-min(timevect)+1),  cor_rtrel_TPA_MSPN{ll,tt,ss}.pval(cc,time-min(timevect)+1)]  =corr(squeeze(rt_msMreluni(ll,tt,subuse)),grind_TPA_MSPN{ll,tt,ss}.individual(:,cc,time));
            [cor_rtNm5_TPA_MSPN{ll,tt,ss}.avg(cc,time-min(timevect)+1),cor_rtNm5_TPA_MSPN{ll,tt,ss}.pval(cc,time-min(timevect)+1)]  =corr(squeeze(rt_7m5(tt,subuse))',        grind_TPA_MSPN{ll,tt,ss}.individual(:,cc,time));
          elseif ll==37
            [cor_rtNm5_TPA_MSPN{ll,tt,ss}.avg(cc,time-min(timevect)+1),cor_rtNm5_TPA_MSPN{ll,tt,ss}.pval(cc,time-min(timevect)+1)]  =corr(squeeze(rt_37m5(tt,subuse))',       grind_TPA_MSPN{ll,tt,ss}.individual(:,cc,time));
          end
        end
      end
      cor_rtrel_TPA_MSPN{ll,tt,ss}.mask=  cor_rtrel_TPA_MSPN{ll,tt,ss}.pval<alphaval;
      cor_rtfg5p1_TPA_MSPN{ll,tt,ss}.mask=cor_rtfg5p1_TPA_MSPN{ll,tt,ss}.pval<alphaval;
      cor_rtfg5p2_TPA_MSPN{ll,tt,ss}.mask=cor_rtfg5p2_TPA_MSPN{ll,tt,ss}.pval<alphaval;
      cor_rtfg5p3_TPA_MSPN{ll,tt,ss}.mask=cor_rtfg5p3_TPA_MSPN{ll,tt,ss}.pval<alphaval;
      cor_rtNm5_TPA_MSPN{ll,tt,ss}.mask=  cor_rtNm5_TPA_MSPN{ll,tt,ss}.pval<alphaval;
      
      if plotflag
        figure(ll);
        cfg=[];
        cfg.maskparameter='mask';
        cfg.layout='elec1010.lay';
        %     cfg.xlim=[cor_rtmed_tacMSpN{ll,tt,ss}.time(1) cor_rtmed_tacMSpN{ll,tt,ss}.time(end)];
        cfg.xlim=[cortimwin(1) cortimwin(2)];
        %     ft_multiplotER(cfg,cor_rtmed_tacMSpN{ll,tt,ss});
        %     ft_multiplotER(cfg,cor_rtmed_tacPaud{ll,tt,ss});
        ft_multiplotER(cfg,cor_rtrel_TPA_MSPN{ll,tt,ss});
        figure(ll+10);
        ft_multiplotER(cfg,cor_rtNm5_TPA_MSPN{ll,tt,ss});
        figure(ll+20);
        ft_multiplotER(cfg,cor_rtfg5p1_TPA_MSPN{ll,tt,ss});
      end
      
      % Question: way to find sig-diff of corrs between conditions?  or just if
      % cor_rtmed_TPA_MSPN is sig?
      
      
      
      if statflag
        rtrel_grind{ll,tt,ss}=grind_TPA_MSPN{ll,tt,ss};
        rtrelshift_grind{ll,tt,ss}=grind_TPA_MSPN{ll,tt,ss};
        rtfg5p1_grind{ll,tt,ss}=grind_TPA_MSPN{ll,tt,ss};
        rtfg5p2_grind{ll,tt,ss}=grind_TPA_MSPN{ll,tt,ss};
        rtfg5p3_grind{ll,tt,ss}=grind_TPA_MSPN{ll,tt,ss};
        rtNm5_grind{ll,tt,ss}=grind_TPA_MSPN{ll,tt,ss};
        
        if ll<10
          rtrel_grind{ll,tt,ss}.individual     =repmat(squeeze(rt_msMreluni(ll,tt,subuse)),[1 size(rtrel_grind{ll,tt,ss}.individual,2) size(rtrel_grind{ll,tt,ss}.individual,3)]);
          rtrelshift_grind{ll,tt,ss}.individual=repmat(squeeze(rt_msMminshiftuni(ll,tt,subuse)),[1 size(rtrel_grind{ll,tt,ss}.individual,2) size(rtrel_grind{ll,tt,ss}.individual,3)]);
        end
        rtfg5p1_grind{ll,tt,ss}.individual=repmat(squeeze(fg5p1(tt,subuse))',          [1 size(rtfg5p1_grind{ll,tt,ss}.individual,2) size(rtfg5p1_grind{ll,tt,ss}.individual,3)]);
        rtfg5p2_grind{ll,tt,ss}.individual=repmat(squeeze(fg5p2(tt,subuse))',          [1 size(rtfg5p1_grind{ll,tt,ss}.individual,2) size(rtfg5p1_grind{ll,tt,ss}.individual,3)]);
        rtfg5p3_grind{ll,tt,ss}.individual=repmat(squeeze(fg5p3(tt,subuse))',          [1 size(rtfg5p1_grind{ll,tt,ss}.individual,2) size(rtfg5p1_grind{ll,tt,ss}.individual,3)]);
        rtNm5_grind{ll,tt,ss}.individual  =repmat(squeeze(rt_3m5(tt,subuse))',        [1 size(rtrel_grind{ll,tt,ss}.individual,2) size(rtrel_grind{ll,tt,ss}.individual,3)]);
        
        grave_TPA_MSPN{ll,tt,ss}=grind_TPA_MSPN{ll,tt,ss};
        grave_TPA_MSPN{ll,tt,ss}.avg=squeeze(mean(grind_TPA_MSPN{ll,tt,ss}.individual,1));
        grave_TPA_MSPN{ll,tt,ss}=rmfield(grave_TPA_MSPN{ll,tt,ss},'individual');
        grave_TPA_MSPN{ll,tt,ss}.dimord='chan_time';
        
        if ll==3
          stattimwin=[.225 .245];
        elseif ll==7
          stattimwin=[.25 .27];
        elseif ll==37
          stattimwin=[.235 .255];
        end
        
        cfg=[];
        cfg.latency=stattimwin;
        cfg.avgovertime='yes';
        cfg.channel={'F2' 'FC2' 'C2' 'CP2'};  % those significant on AT70 ERP difference
%         cfg.channel={'F2' 'FC2' 'C2' 'Cz' 'C1' 'FC1' 'F1' 'Fz'};
        cfg.avgoverchan='yes';
        %         cfg.neighbours=neighbours;
        cfg.parameter='individual';
%         cfg.method='montecarlo';
        cfg.method='analytic';
%         cfg.numrandomization=1000;
%         cfg.alpha=.05; % .05/9 (3 p1,p2,p3  X  3 (3,7,37)
%                 cfg.correctm='fdr';
        %         cfg.correctm='cluster';
        %         cfg.clusteralpha = 0.05;
        %         cfg.clusterstatistic = 'maxsum';
        %         cfg.minnbchan = 2;
        cfg.statistic='correlationT';  %% DATA MUST BE MEAN-centred prior to using this with Montecarlo, else the permutation doesn't make sense!
        cfg.design=zeros(2,2*nsub);
        cfg.design(1,:)=[ones(1,nsub) 2*ones(1,nsub)]; % must be these numbers to code the labels for 'correlationT'
        cfg.design(2,:)=[1:nsub 1:nsub];
        cfg.ivar=1;
        cfg.uvar=2;
        if ll<10
          stattcorrel_TPA_MSPN{ll,tt,ss}  =ft_timelockstatistics(cfg, grind_TPA_MSPN{ll,tt,ss}, rtrel_grind{ll,tt,ss});
          stattcorrelshift_TPA_MSPN{ll,tt,ss}  =ft_timelockstatistics(cfg, grind_TPA_MSPN{ll,tt,ss}, rtrelshift_grind{ll,tt,ss});
        end
        stattcorfg5p1_TPA_MSPN{ll,tt,ss}=ft_timelockstatistics(cfg, grind_TPA_MSPN{ll,tt,ss}, rtfg5p1_grind{ll,tt,ss});
        stattcorfg5p2_TPA_MSPN{ll,tt,ss}=ft_timelockstatistics(cfg, grind_TPA_MSPN{ll,tt,ss}, rtfg5p2_grind{ll,tt,ss});
        stattcorfg5p3_TPA_MSPN{ll,tt,ss}=ft_timelockstatistics(cfg, grind_TPA_MSPN{ll,tt,ss}, rtfg5p3_grind{ll,tt,ss});
        stattcorNm5_TPA_MSPN{ll,tt,ss}  =ft_timelockstatistics(cfg, grind_TPA_MSPN{ll,tt,ss}, rtNm5_grind{ll,tt,ss});
        
        if plotflag
          try  % if cfg.avgovertime='yes';
            if ll<10
              try close 101;end
              topoplot_highlight(101,grave_TPA_MSPN{ll,tt,ss},stattcorrel_TPA_MSPN{ll,tt,ss}.time,stattcorrel_TPA_MSPN{ll,tt,ss});
              try close 102;end
              topoplot_highlight(102,cor_rtrel_TPA_MSPN{ll,tt,ss},stattcorrel_TPA_MSPN{ll,tt,ss}.time,stattcorrel_TPA_MSPN{ll,tt,ss},[-1 1]);
            end
            try close 105;end
            topoplot_highlight(105,grave_TPA_MSPN{ll,tt,ss},stattcorfg5p1_TPA_MSPN{ll,tt,ss}.time,stattcorNm5_TPA_MSPN{ll,tt,ss});
            try close 106;end
            topoplot_highlight(106,cor_rtNm5_TPA_MSPN{ll,tt,ss},stattcorfg5p1_TPA_MSPN{ll,tt,ss}.time,stattcorNm5_TPA_MSPN{ll,tt,ss},[-1 1]);
            try close 103;end
            topoplot_highlight(103,grave_TPA_MSPN{ll,tt,ss},stattcorfg5p1_TPA_MSPN{ll,tt,ss}.time,stattcorfg5p1_TPA_MSPN{ll,tt,ss});
            try close 104;end
            topoplot_highlight(104,cor_rtfg5p1_TPA_MSPN{ll,tt,ss},stattcorfg5p1_TPA_MSPN{ll,tt,ss}.time,stattcorfg5p1_TPA_MSPN{ll,tt,ss},[-1 1]);
            try close 107;end
            topoplot_highlight(107,grave_TPA_MSPN{ll,tt,ss},stattcorfg5p2_TPA_MSPN{ll,tt,ss}.time,stattcorfg5p2_TPA_MSPN{ll,tt,ss});
            try close 108;end
            topoplot_highlight(108,cor_rtfg5p2_TPA_MSPN{ll,tt,ss},stattcorfg5p2_TPA_MSPN{ll,tt,ss}.time,stattcorfg5p2_TPA_MSPN{ll,tt,ss},[-1 1]);
            try close 109;end
            topoplot_highlight(109,grave_TPA_MSPN{ll,tt,ss},stattcorfg5p3_TPA_MSPN{ll,tt,ss}.time,stattcorfg5p3_TPA_MSPN{ll,tt,ss});
            try close 110;end
            topoplot_highlight(110,cor_rtfg5p3_TPA_MSPN{ll,tt,ss},stattcorfg5p3_TPA_MSPN{ll,tt,ss}.time,stattcorfg5p3_TPA_MSPN{ll,tt,ss},[-1 1]);
          catch  % if cfg.avgovertime='no';
            if ll<10
              try close 101;end
              topoplot_highlight(101,grave_TPA_MSPN{ll,tt,ss},stattimwin,stattcorrel_TPA_MSPN{ll,tt,ss});
              try close 102;end
              topoplot_highlight(102,cor_rtrel_TPA_MSPN{ll,tt,ss},stattimwin,stattcorrel_TPA_MSPN{ll,tt,ss},[-1 1]);
            end
            try close 105;end
            topoplot_highlight(105,grave_TPA_MSPN{ll,tt,ss},stattimwin,stattcorNm5_TPA_MSPN{ll,tt,ss});
            try close 106;end
            topoplot_highlight(106,cor_rtNm5_TPA_MSPN{ll,tt,ss},stattimwin,stattcorNm5_TPA_MSPN{ll,tt,ss},[-1 1]);
            try close 103;end
            topoplot_highlight(103,grave_TPA_MSPN{ll,tt,ss},stattimwin,stattcorfg5p1_TPA_MSPN{ll,tt,ss});
            try close 104;end
            topoplot_highlight(104,cor_rtfg5p1_TPA_MSPN{ll,tt,ss},stattimwin,stattcorfg5p1_TPA_MSPN{ll,tt,ss},[-1 1]);
            try close 107;end
            topoplot_highlight(107,grave_TPA_MSPN{ll,tt,ss},stattimwin,stattcorfg5p2_TPA_MSPN{ll,tt,ss});
            try close 108;end
            topoplot_highlight(108,cor_rtfg5p2_TPA_MSPN{ll,tt,ss},stattimwin,stattcorfg5p2_TPA_MSPN{ll,tt,ss},[-1 1]);
            try close 109;end
            topoplot_highlight(109,grave_TPA_MSPN{ll,tt,ss},stattimwin,stattcorfg5p3_TPA_MSPN{ll,tt,ss});
            try close 110;end
            topoplot_highlight(110,cor_rtfg5p3_TPA_MSPN{ll,tt,ss},stattimwin,stattcorfg5p3_TPA_MSPN{ll,tt,ss},[-1 1]);
          end
          
          
          if printflag
            if ll<10
              print(101,[fdir 'corstatrel2_grave_TPA_MSPN_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
              print(102,[fdir 'corstatrel2_cor_TPA_MSPN_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            end
            print(105,[fdir 'corstatNm5_grave_TPA_MSPN_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(106,[fdir 'corstatNm5_cor_TPA_MSPN_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(103,[fdir 'corstatfg5p1_grave_TPA_MSPN_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(104,[fdir 'corstatfg5p1_cor_TPA_MSPN_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(107,[fdir 'corstatfg5p2_grave_TPA_MSPN_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(108,[fdir 'corstatfg5p2_cor_TPA_MSPN_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(109,[fdir 'corstatfg5p3_grave_TPA_MSPN_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            print(110,[fdir 'corstatfg5p3_cor_TPA_MSPN_topoOverTime_ta_cond' num2str(ll) num2str(tt) num2str(sleep) '.png'],'-dpng')
            
          end % printflag
          
        end % plotflag
        save([edir 'cor_statfg_sleep' num2str(sleep) '_iter' num2str(iter) '.mat'],'stattcor*');
      end % statflag
      
    end % ll
    save([edir 'cor_corfg_sleep' num2str(sleep) '_iter' num2str(iter) '.mat'],'cor*');
  end % tt
end % iter


%%  corr of fg5 to MS raw (not diff)

for lluse=1:7
  ll=soalist(lluse);

  figure(1);
  [ccc,ppp]=corr(fg5p1(3,1:22)',squeeze(mean(grind_MStlock_save{ll,tt,ss}.individual(:,match_str(grind_MStlock_save{ll,tt,ss}.label,{'FC2' 'F2' 'C2' 'CP2'}),:),2)));
  subplot(7,1,lluse);
  plot(grind_MStlock_save{ll,tt,ss}.time,ccc);
  hold on;plot(grind_MStlock_save{ll,tt,ss}.time,ppp<.001)
  
  figure(2);
  [ccc,ppp]=corr(fg5p2(3,1:22)',squeeze(mean(grind_MStlock_save{ll,tt,ss}.individual(:,match_str(grind_MStlock_save{ll,tt,ss}.label,{'FC2' 'F2' 'C2' 'CP2'}),:),2)));
  subplot(7,1,lluse);
  plot(grind_MStlock_save{ll,tt,ss}.time,ccc);
  hold on;plot(grind_MStlock_save{ll,tt,ss}.time,ppp<.001)

  figure(3);
  [ccc,ppp]=corr(fg5p3(3,1:22)',squeeze(mean(grind_MStlock_save{ll,tt,ss}.individual(:,match_str(grind_MStlock_save{ll,tt,ss}.label,{'FC2' 'F2' 'C2' 'CP2'}),:),2)));
  subplot(7,1,lluse);
  plot(grind_MStlock_save{ll,tt,ss}.time,ccc);
  hold on;plot(grind_MStlock_save{ll,tt,ss}.time,ppp<.001)
end




%%

clearvars -except sub *dir
% chanplot{1}={'Fz' 'Cz' 'F1' 'F2' 'FC1' 'FC2' 'C1' 'C2'}; % frontocentral
% chanplot{2}={'CP5' 'P5' 'PO3' 'POz' 'PO4' 'P6' 'CP6' 'Pz' 'P3' 'P4'}; % occipital
% chanplot{3}={'FC6' 'FC4' 'F4' 'F6' 'F8' 'FT8' 'T8' 'C6' 'C4'}; % right frontotemporal
% chanlabel{1}='Frontocentral electrodes';
% chanlabel{2}='Occipital-parietal electrodes';
% chanlabel{3}='Right frontotemporal electrodes';
%
chanplot{1}={'Fz' 'Cz' 'F1' 'F2' 'FC1' 'FC2' 'C1' 'C2'}; % frontocentral
chanplot{2}={'CP5' 'POz' 'Pz' 'P3' 'P4' 'C4' 'O1' 'O2' 'P7' 'PO7'}; % occipital
chanlabel{1}='Frontocentral electrodes';
chanlabel{2}='Occipital-parietal electrodes';
iterflag=1;
iter=27;

sleep=0;
soalist=[1 3 4 5 6 7 9];
tt=3;
ss=10;
timwinstatflag=0; % =1 for using only stat.time, =0 for using [-0.5 1.0];

load([ddir 'rt_allsubj.mat']);

if iterflag
  load([edir 'tlock_statmc_sleep' num2str(sleep) '_iter' num2str(iter) '.mat']);
  load([edir 'tlock_grind_sleep' num2str(sleep) '_iter' num2str(iter) '.mat']);
  load([edir 'cor_stat_sleep' num2str(sleep) '_iter' num2str(iter) '.mat']);
  load([edir 'cor_cor_sleep' num2str(sleep) '_iter' num2str(iter) '.mat']);
else
  load([edir 'tlock_statmc_sleep' num2str(sleep) '.mat']);
  load([edir 'tlock_grind_sleep' num2str(sleep) '.mat']);
  load([edir 'cor_stat_sleep' num2str(sleep) '.mat']);
  load([edir 'cor_cor_sleep' num2str(sleep) '.mat']);
end

colorset=varycolor(size(grind_tacMSpN_save{1,3,ss}.individual,1));
colormapuse=jet(size(grind_tacMSpN_save{1,3,ss}.individual,1));
colorset=colormapuse;

close all
for ll=soalist
  stattimwin=[stattcormed_tacMSpN{ll,tt,ss}.time(1) stattcormed_tacMSpN{ll,tt,ss}.time(end)];
  stattime=stattcormed_tacMSpN{ll,tt,ss}.time;
  cfg=[];
  if timwinstatflag==1
    cfg.latency=stattimwin;
  elseif timwinstatflag==0
    cfg.latency=[-0.55 1];
  end
  tmp=ft_selectdata(cfg,cor_rtmed_tacMSpN{ll,tt,ss});
  tmp.mask=stattcormed_tacMSpN{ll,tt,ss}.mask;
  chanlist=tmp.label(find(any(tmp.mask,2)));
  
  chanuse=[];
  if tt==2
    if any(ll==[3 4 5 9])
      chanuse{1}=chanlist;
    elseif any(ll==[1 6 7])
      chanuse={chanplot{1}};
    end
  elseif tt==3
    if any(ll==[1 3 5 9])
      chanuse{1}=chanlist;
    elseif any(ll==[4 6 7])
      chanuse={chanplot{1}};
    end
  end
  
  cfg=[];
  cfg.avgoverchan='yes';
  cfg.channel=chanuse{1};
  tmp1=ft_selectdata(cfg,tmp);
  if timwinstatflag==0
    tmp1mask=mean(tmp.mask(match_str(tmp.label,cfg.channel),:),1);
    tmp1.mask=zeros(1,length(tmp1.time));
    tmp1.mask(dsearchn(tmp1.time',stattimwin(1)):dsearchn(tmp1.time',stattimwin(end)))=tmp1mask;
  end
  tmp1.mask=logical(ceil(tmp1.mask));
  if length(find(tmp1.mask))==1
    tmp1.mask(find(tmp1.mask)+1)=1; % this may seem like cheating, but otherwise the boxline doesn't appear in figure
  end
  
  
  figure(10*ll);
  cfg=[];
  cfg.parameter='avg';
  cfg.layout='elec1010.lay';
  cfg.maskparameter='mask';
  cfg.maskstyle='box'; % default
  %     cfg.maskstyle='thickness';
  cfg.ylim=[-1 1];
  cfg.linewidth=3;
  ft_singleplotER(cfg,tmp1);
  hold on; % add statwindow lines
  plot([stattimwin(1) stattimwin(1)],cfg.ylim,'k--')
  plot([stattimwin(2) stattimwin(2)],cfg.ylim,'k--')
  print(10*ll,[fdir 'cormed_tacMSpN_final_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  
  if ~isempty(chanlist)
    masktime=find(any(tmp.mask,1));
    cfg=[];
    cfg.parameter='avg';
    cfg.layout='elec1010.lay';
    cfg.maskalpha=0.5;
    cfg.zlim=[-1 1];
    cfg.highlight='on';
    cfg.highlightsize=12;
    cfg.xlim=[stattime(masktime(1)) stattime(masktime(end))];
    %     cfg.highlightchannel=tmp.label(find(ceil(mean(tmp.mask(:,dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2))),2))));
    cfg.highlightchannel=tmp.label(find(ceil(mean(tmp.mask(:,dsearchn(stattime',cfg.xlim(1)):dsearchn(stattime',cfg.xlim(2))),2))));
    figure(100*ll);
    ft_topoplotER(cfg,tmp);
    print(100*ll,[fdir 'cormed_topo_tacMSpN_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    
    cg=1;
    [rtsort,rtord]=sort(squeeze(rt_med(ll,tt,:)));
    figure(1000*ll);
    set(gca, 'ColorOrder', colorset);hold all;
    plot(grind_tacMSpN_save{ll,tt,ss}.time,squeeze(mean(grind_tacMSpN_save{ll,tt,ss}.individual(rtord(~isnan(rtsort)),match_str(grind_tacMSpN_save{ll,tt,ss}.label,chanuse{cg}),:),2)))
    hold on;
    plot(grind_tacMSpN_save{ll,tt,ss}.time,mean(squeeze(mean(grind_tacMSpN_save{ll,tt,ss}.individual(rtord(~isnan(rtsort)),match_str(grind_tacMSpN_save{ll,tt,ss}.label,chanuse{cg}),:),2)),1),'k--','LineWidth',3')
    hold on;
    plot(cor_rtmed_tacMSpN{ll,tt,ss}.time,10*mean(cor_rtmed_tacMSpN{ll,tt,ss}.avg(match_str(cor_rtmed_tacMSpN{ll,tt,ss}.label,chanuse{cg}),:),1),'k','LineWidth',3')
    axis([-0.5 1 -inf inf])
    print(1000*ll,[fdir 'grindtacMSpN_sortbyrtmed_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
  end
  
  % % %%%%%%%%%%%%%%%%%%%%%%%%
  stattimwin=[stattcorrel_TPA_MSPN{ll,tt,ss}.time(1) stattcorrel_TPA_MSPN{ll,tt,ss}.time(end)];
  stattime=stattcorrel_TPA_MSPN{ll,tt,ss}.time;
  cfg=[];
  if timwinstatflag==1
    cfg.latency=stattimwin;
  elseif timwinstatflag==0
    cfg.latency=[-0.55 1];
  end
  tmp=ft_selectdata(cfg,cor_rtrel_TPA_MSPN{ll,tt,ss});
  tmp.mask=stattcorrel_TPA_MSPN{ll,tt,ss}.mask;
  
  chanlist=tmp.label(find(any(tmp.mask,2)));
  lc=1;
  rc=1;
  leftchan=[];
  rightchan=[];
  for cc=1:length(chanlist)
    if rem(str2num(chanlist{cc}(end)),2)
      leftchan{lc}=chanlist{cc};
      lc=lc+1;
    else
      rightchan{rc}=chanlist{cc};
      rc=rc+1;
    end
  end
  
  chanuse=[];
  if tt==2
    if ll==1
      chanuse{1}=leftchan;
      chanuse{2}=rightchan;
    elseif any(ll==[3 4 5 9])
      chanuse{1}=chanlist;
      %     elseif ll==4
      %       chanuse{1}=leftchan;
      %       chanuse{2}=[{'P4'} {'O2'} {'P8'}];
      %       chanuse{3}=[{'F6'} {'C6'} {'AF8'} {'FT8'}];
    elseif any(ll==[6 7])
      chanuse={chanplot{2}};
    end
  elseif tt==3
    if any(ll==[])
      chanuse{1}=chanlist;
    elseif ll==9
      chanuse{1}=leftchan;
      chanuse{2}=rightchan;
    elseif ll==5
      chanuse{1}=setdiff(chanlist,[{'P5'} {'C5'} {'FT7'} ]);
      chanuse{2}=[{'P5'} {'C5'} {'FT7'} ];
    elseif ll==4
      chanuse{1}=[{'F5'} {'C5'} {'AF7'}]
      chanuse{2}=[{'P3'} {'O1'} {'P7'}];
      %     elseif ll==5
      % %       chanuse{1}=setdiff(chanlist,[{'Fz'} {'Fpz'} {'F7'}]);
      % %       chanuse{2}=[{'Fz'} {'Fpz'} {'F7'}];
      %       chanuse{1}=setdiff(chanlist,[{'TP9'}]);
      %       chanuse{2}=[{'TP9'}];
      %     elseif ll==9
      %       chanuse{1}=setdiff(chanlist,[{'AF8'} {'F6'} {'C6'}]);
      %       chanuse{2}=[{'AF8'} {'F6'} {'C6'}];
    elseif any(ll==[1 3 6 7])
      chanuse={chanplot{2}};
    end
  end
  
  for cg=1:length(chanuse)
    
    cfg=[];
    cfg.avgoverchan='yes';
    cfg.channel=chanuse{cg};
    tmp1=ft_selectdata(cfg,tmp);
    if timwinstatflag==0
      tmp1mask=mean(tmp.mask(match_str(tmp.label,cfg.channel),:),1);
      tmp1.mask=zeros(1,length(tmp1.time));
      tmp1.mask(dsearchn(tmp1.time',stattimwin(1)):dsearchn(tmp1.time',stattimwin(end)))=tmp1mask;
    end
    tmp1.mask=logical(ceil(tmp1.mask));
    if length(find(tmp1.mask))==1
      tmp1.mask(find(tmp1.mask)+1)=1; % this may seem like cheating, but otherwise the boxline doesn't appear in figure
    end
    
    figure(10*ll+cg);
    cfg=[];
    cfg.parameter='avg';
    cfg.layout='elec1010.lay';
    cfg.maskparameter='mask';
    cfg.maskstyle='box'; % default
    %     cfg.maskstyle='thickness';
    cfg.ylim=[-1 1];
    cfg.linewidth=3;
    ft_singleplotER(cfg,tmp1);
    hold on; % add statwindow lines
    plot([stattimwin(1) stattimwin(1)],cfg.ylim,'k--')
    plot([stattimwin(2) stattimwin(2)],cfg.ylim,'k--')
    print(10*ll+cg,[fdir 'correl_TPAMSPN_chan' num2str(cg) '_final_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
    
    if ~isempty(chanlist)
      %       masktime=find(any(tmp.mask,1));
      masktime=find(any(tmp.mask(match_str(tmp.label,chanuse{cg}),:),1));
      
      maskbreak=find(diff(masktime)>1);
      clear seg*
      if ~isempty(maskbreak)
        for mm=1:[length(maskbreak)+1]
          if mm==1
            seg(mm)=maskbreak(mm);
            segind(mm,:)=[1 maskbreak(mm)];
          elseif mm==[length(maskbreak)+1]
            seg(mm)=length(masktime)-maskbreak(mm-1);
            segind(mm,:)=[maskbreak(mm-1)+1 length(masktime)];
          else
            seg(mm)=maskbreak(mm)-maskbreak(mm-1);
            segind(mm,:)=[maskbreak(mm-1)+1 maskbreak(mm)];
          end
        end
        [mx,mnd]=max(seg);
        masktime=masktime(segind(mnd,1):segind(mnd,2));
      end
      
      cfg=[];
      cfg.parameter='avg';
      cfg.layout='elec1010.lay';
      cfg.maskalpha=0.5;
      cfg.zlim=[-1 1];
      cfg.highlight='on';
      cfg.highlightsize=12;
      %       cfg.xlim=[tmp.time(masktime(1)) tmp.time(masktime(end))];
      cfg.xlim=[stattime(masktime(1)) stattime(masktime(end))];
      %       cfg.highlightchannel=tmp.label(find(ceil(mean(tmp.mask(:,dsearchn(tmp.time',cfg.xlim(1)):dsearchn(tmp.time',cfg.xlim(2))),2))));
      cfg.highlightchannel=tmp.label(find(ceil(mean(tmp.mask(:,dsearchn(stattime',cfg.xlim(1)):dsearchn(stattime',cfg.xlim(2))),2))));
      figure(100*ll+cg);
      ft_topoplotER(cfg,tmp);
      print(100*ll+cg,[fdir 'correl_topo_TPAMSPN_chan' num2str(cg) '_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
      
      cfg=[];
      cfg.operation='subtract';
      cfg.parameter='individual';
      grind_TPA_MSPN{ll,tt,ss}=ft_math(cfg,grind_tacPaud_save{ll,tt,ss},grind_tacMSpN_save{ll,tt,ss});
      
      %       cfg=[];
      %       cfg.latency=[stattime(masktime(1)) stattime(masktime(end))];
      %       cfg.avgovertime='yes';
      %       cfg.channel=chanuse{cg};
      %       tmperp=ft_selectdata(cfg,grind_TPA_MSPN{ll,tt,ss});
      % %       tmpcor=ft_selectdata(cfg,tmp);
      %       [erpsort,erpord]=sort(mean(tmperp.individual,2));
      
      [rtsort,rtord]=sort(squeeze(rt_msMminshiftuni(ll,tt,:)));
      figure(1000*ll+cg);
      set(gca, 'ColorOrder', colorset);hold all;
      plot(grind_TPA_MSPN{ll,tt,ss}.time,squeeze(mean(grind_TPA_MSPN{ll,tt,ss}.individual(rtord(~isnan(rtsort)),match_str(grind_TPA_MSPN{ll,tt,ss}.label,chanuse{cg}),:),2)))
      hold on;
      plot(grind_TPA_MSPN{ll,tt,ss}.time,mean(squeeze(mean(grind_TPA_MSPN{ll,tt,ss}.individual(rtord(~isnan(rtsort)),match_str(grind_TPA_MSPN{ll,tt,ss}.label,chanuse{cg}),:),2)),1),'k--','LineWidth',3')
      hold on;
      plot(cor_rtrel_TPA_MSPN{ll,tt,ss}.time,10*mean(cor_rtrel_TPA_MSPN{ll,tt,ss}.avg(match_str(cor_rtrel_TPA_MSPN{ll,tt,ss}.label,chanuse{cg}),:),1),'k','LineWidth',3')
      axis([-0.5 1 -inf inf])
      print(1000*ll+cg,[fdir 'grindTPAMSPN_sortbyrtrel_chan' num2str(cg) '_' num2str(ll) num2str(tt) num2str(ss) '.png'],'-dpng')
      
      
    end
  end % cg
  
end %ll


