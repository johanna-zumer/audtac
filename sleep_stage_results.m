% analysis of sleep stage results.


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

%%

iiBkeep=zeros(1,32);
iiSkeep=zeros(1,32);
subcnt=1;
% for ii=[3 5:18 20:32]
for ii=[8:18 20:32]
  for sleep=[0 1]
    clearvars -except ii sub  edir ddir sleep subcnt sleep0023 ii*keep stages*
    cd([edir sub{ii} ])
    if sleep
      try
        load CRC_sleep_withTom.mat
      catch
        load CRC_sleep.mat
      end
    else
      try
        load CRC_wake_withTom.mat
      catch
        load CRC_wake.mat
      end
    end
    
    if sleep==0
      figind=1;
    elseif sleep==1
      figind=2;
    end
    figure(figind);
    subplot(4,7,subcnt);
    % yaxis is minutes of data
    bar([sum(CRC.score{1}==0) sum(CRC.score{1}==1) sum(CRC.score{1}==2) sum(CRC.score{1}==3) sum(CRC.score{1}==5)]/2);
    if sleep==1
      stages_totalmins(subcnt,:,2)=[sum(CRC.score{1}==0) sum(CRC.score{1}==1) sum(CRC.score{1}==2) sum(CRC.score{1}==3) sum(CRC.score{1}==5)]/2;
    elseif sleep==0
      stages_totalmins(subcnt,:,1)=[sum(CRC.score{1}==0) sum(CRC.score{1}==1) sum(CRC.score{1}==2) sum(CRC.score{1}==3) sum(CRC.score{1}==5)]/2;
    end
    set(get(figind,'Children'),'xTickLabel',[0 1 2 3 5])
    if sleep==0
      set(get(figind,'Children'),'ylim',[0 50]);
      sleep0023(subcnt,1)=sum(CRC.score{1}==0)/2; % mins awake during Sitting
    elseif sleep==1
      set(get(figind,'Children'),'ylim',[0 130]);
      sleep0023(subcnt,2)=sum(CRC.score{1}==0)/2; % mins awake during Bed
      sleep0023(subcnt,3)=(sum(CRC.score{1}==2)+sum(CRC.score{1}==3))/2; % mins N2+N3 during Bed
    end
  end
  thresh=31;
  figure(3);
  freezeColors;
  subplot(4,7,subcnt);
  % yaxis is minutes of data
  bar(sleep0023(subcnt,:));
  set(get(3,'Children'),'ylim',[0 170]);
  set(get(3,'Children'),'xTickLabel',[0 0 23])
  if sleep0023(subcnt,1)<thresh || sleep0023(subcnt,3)<thresh
    if sleep0023(subcnt,3)<thresh
      colormap([1 0 0]); % red
    elseif sleep0023(subcnt,1)+sleep0023(subcnt,2) <thresh
      colormap([1 0 0]);% red
    else
      colormap([1 1 0]);% yellow
    end
  else
    colormap([0 1 0]);% green
  end
  
  thresh=27.5; % minutes
  figure(4);
  freezeColors;
  subplot(4,7,subcnt);
  % yaxis is minutes of data
  bar(sleep0023(subcnt,:));
  set(get(4,'Children'),'ylim',[0 170]);
  set(get(4,'Children'),'xTickLabel',[0 0 23])
  if sleep0023(subcnt,1)<thresh || sleep0023(subcnt,3)<thresh
    if sleep0023(subcnt,3)<thresh
      colormap([1 0 0]);
    elseif sleep0023(subcnt,1)+sleep0023(subcnt,2) <thresh
      colormap([1 0 0]);
    else
      colormap([1 1 0]);
    end
  else
    colormap([0 1 0]);
    iikeep(ii)=1;
  end
  if sleep0023(subcnt,1)>=thresh
    iiSkeep(ii)=1;
  end
  if sleep0023(subcnt,3)>=thresh
    iiBkeep(ii)=1;
  end
    
  
  
  subcnt=subcnt+1;
end

mean(stages_totalmins(:,:,2))
std(stages_totalmins(:,:,2))/sqrt(25-1)

cd('D:\audtac\figs\sleep_staging');
print(1,'Figure1_sitting','-dpng')
print(2,'Figure2_bed','-dpng')
print(3,'Figure3_90trials','-dpng')
print(4,'Figure4_80trials','-dpng')

iiSuse=find(iiSkeep);
iiBuse=find(iiBkeep);
save([edir 'iikeep.mat'],'ii*use');


