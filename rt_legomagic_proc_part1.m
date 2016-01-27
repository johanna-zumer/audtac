% processing of RT data from lego-device with chinrest and Notts EarTone
clear a* r* b* f* h* p* d* s*
close all

if ispc
  bdir='D:\audtac\behav_data\';
  ddir='D:\audtac\legomagic\diaries\';
else
  bdir='/mnt/hgfs/D/audtac/behav_data/';
  ddir='/mnt/hgfs/D/audtac/legomagic/diaries/';
end

cd(bdir)

% see rt_legomagic_proc.m for these subjects.
% sub{10}='p10'; % 51520 02/05/14 with condition 3
% sub{11}='p11'; % 53188 02/05/14 with condition 3
% sub{12}='p12'; % 61768 06/05/14 with condition 3
% sub{13}='p13'; % 57811 06/05/14 with condition 4 (but epi in only left ear)
% sub{14}='p14'; % 60811 06/05/14 with condition 4 (but epi in only left ear)
% sub{15}='p15'; % jz
% sub{16}='p16'; % jz
% sub{17}='p17'; % 58054 16/05/14 with condition 2 (epi noise, no ear muffs)
% sub{18}='p18'; % 61660 16/05/14
% sub{19}='p19'; % 54046 19/05/14
% sub{20}='p20'; % 51736 19/05/14
% sub{21}='p21'; % 61762 20/05/14
% sub{22}='p22'; % 46708 20/05/14
% sub{23}='p23'; % 61495 20/05/14
% sub{24}='p24'; % 54142 20/05/14
% sub{25}='p25'; % 61840 20/05/14

% sub{100}='p01'; 
sub{101}='e01'; 
sub{102}='e02'; 
sub{103}='e03'; 
sub{104}='e04'; 
%%%  above is pilot, below is real
sub{105}='e05'; 
sub{106}='e06'; 
sub{107}='e07'; 
sub{108}='e08'; 
sub{109}='e09';
sub{110}='e10';
sub{111}='e11';
sub{112}='e12';
sub{113}='e13';
sub{114}='e14';
sub{115}='e15';
sub{116}='e16';
sub{117}='e17';
sub{118}='e18';
sub{119}='e19';
sub{120}='e20';
sub{121}='e21';
sub{122}='e22';
sub{123}='e23';
sub{124}='e24';
sub{125}='e25';
sub{126}='e26';
sub{127}='e27';
sub{128}='e28';
sub{129}='e29';
sub{130}='e30';
sub{131}='e31';
sub{132}='e32';


ii=103;

%%

% note: check fin.aud delay

files=dir([ddir sub{ii} '*r.mat']);

rtsrall=nan(1,10);
for bb=1:length(files)
  clear resptime rts aa ll rtsshuf rtsr touchdelay auddelay
  load([ddir files(bb).name]);
  
  nt=length(info.lightsensor);
  if ii==17 && bb==2
    nt=148; % she stopped to ask a question
  end
  
  if nt<5
    continue
  end
  
  resptime=info.prsout(1:nt)-info.start_time;
  
  
  %   if ~isfield(fin,'opticSR')
  %     cnt=1;
  %     opticSR=0;
  %     while opticSR==0
  %       try
  %         opticSR=diff(info.lighttime{cnt}(1:2));
  %       end
  %       cnt=cnt+1;
  %     end
  %   else
  opticSR=fin.opticSR/1000;
  %   end
  
  
  for nn=1:nt
    lighttime=0:opticSR:opticSR*[length(info.lightsensor{nn})-1];
    if ~isempty(info.lightsensor{nn})
      % this threshold of 80% is arbitrary...should be tested (70? 90?)
      tmp=find(info.lightsensor{nn}> .8*[max(info.lightsensor{nn})-min(info.lightsensor{nn})]+min(info.lightsensor{nn}) );
      %                 touchdelay(nn)=info.lighttime{nn}(tmp(1))-info.lighttime{nn}(1);
      try
        touchdelay(nn)=lighttime(tmp(1))-lighttime(1);
      catch  % this might happen if block crashes in middle (e.g. ii==116,bb==3)
        touchdelay(nn)=nan;
      end
    else
      touchdelay(nn)=nan;
    end
  end
  tmp=nanmedian(touchdelay);
  touchdelay(find(touchdelay>tmp+.04))=nan; % are these too strict?
  touchdelay(find(touchdelay<tmp-.04))=nan; % ?
  touchdelay(find(touchdelay<.1))=nan;
  touchdelay((touchdelay-info.time_touch)>.005)=nan;
  if ii==111 & bb==3 % tapper got messed up for this subject this run
    touchdelay(173:end)=nan;
  end
  figure;plot(touchdelay);axis([-inf inf .17 .25]);
  notkeep=[];
  for nn=4:length(touchdelay)
    numreal=0;
    cnt=0;
    prevtime=[];
    while numreal<3
      if nn-1-cnt<1
        break
      end
      if ~isnan(touchdelay(nn-1-cnt))
        prevtime=[prevtime touchdelay(nn-1-cnt)];
        numreal=numreal+1;
      end
      cnt=cnt+1;
    end
    
    if (touchdelay(nn)-nanmean(prevtime))>.007
      notkeep=[notkeep nn];
    end
  end
  touchdelay(notkeep)=nan;
  notkeep
  
  figure;plot(touchdelay);axis([-inf inf .17 .25]);
  
  if ~isfield(fin,'audDelay')
    auddelay=.019;
  else
    auddelay=fin.audDelay/1000;
  end
  
  
  % just in case some actual audio or tactile didn't occur when it was supposed to
  soa_seq=info.soa_seq;
  soa_seq(isnan(info.audio_time) & isnan(info.valve_time))=0;
  soa_seq(isnan(info.audio_time) & soa_seq~=-2 & soa_seq~=0)=-1;
  soa_seq(isnan(info.valve_time) & soa_seq~=-1 & soa_seq~=0)=-2;
  soa_seq=soa_seq(1:nt);
  
  
  soaeff=[info.audio_time(1:nt)+auddelay] - [info.valve_time(1:nt)+touchdelay];
  soareal=soa_seq; % soa category based on actual SOA times
  
  %   for nn=1:nt
  %     if soaeff(nn)<-.25 & ssoaeff(nn)>-.52 & soareal(nn)>0
  %       soareal(nn)=1;
  %     elseif soaeff(nn)<-.05 & ssoaeff(nn)>-.52 & soareal(nn)>0
  %       soareal(nn)=3;
  %     elseif soaeff(nn)<-.01 & ssoaeff(nn)>-.52 & soareal(nn)>0
  %       soareal(nn)=4;
  %     elseif soaeff(nn)<.01 & ssoaeff(nn)>-.52 & soareal(nn)>0
  %       soareal(nn)=5;
  %     elseif soaeff(nn)<.05 & ssoaeff(nn)>-.52 & soareal(nn)>0
  %       soareal(nn)=6;
  %     elseif soaeff(nn)<.25 & ssoaeff(nn)>-.52 & soareal(nn)>0
  %       soareal(nn)=7;
  %     elseif soaeff(nn)>=.25 & ssoaeff(nn)>-.52 & soareal(nn)>0
  %       soareal(nn)=9;
  %     end
  %   end
  
  for nn=1:nt
    if soaeff(nn)<-.48 & soaeff(nn)>-.52 & soareal(nn)>0 % -500
      soareal(nn)=1;
    elseif soaeff(nn)<-.06 & soaeff(nn)>-.08 & soareal(nn)>0 % -70
      soareal(nn)=3;
    elseif soaeff(nn)<-.014 & soaeff(nn)>-.026 & soareal(nn)>0 % -20
      soareal(nn)=4;
    elseif soaeff(nn)<.006 & soaeff(nn)>-.006 & soareal(nn)>0 % 0
      soareal(nn)=5;
    elseif soaeff(nn)<.026 & soaeff(nn)>.014 & soareal(nn)>0 % 20
      soareal(nn)=6;
    elseif soaeff(nn)<.08 & soaeff(nn)>.06 & soareal(nn)>0 % 70
      soareal(nn)=7;
    elseif soaeff(nn)<.52 & soaeff(nn)>.48 & soareal(nn)>0 % 500
      soareal(nn)=9;
    elseif any(soareal(nn)==[-2 -1 0])
    else
      soareal(nn)=-3; % not use
    end
  end
  
  
  
  disp(['This number was changed based on actual time delay: ' num2str(length(find(soareal-soa_seq)))])
  soa_seq=soareal;
  %         soareal(isnan(soaeff))=nan;
  
  %         if 0
  %         else % we should use the real times unless we know it's perfect
  %             soa_seq=soareal;
  %         end
  
  %   firststimtime=min(info.audio_time(1:nt)+.019-info.start_time, info.valve_time(1:nt)+touchdelay-info.start_time);
  firststimtime=min(info.audio_time(1:nt)+auddelay-info.start_time, info.valve_time(1:nt)+touchdelay-info.start_time);
  
  rts=resptime-firststimtime;
  
  
  [aa,cc]=sort(soa_seq);
  rtsshuf=[rts(cc)];
  
  soaseq=unique(soa_seq(~isnan(soa_seq)));
  mapcode=[-2 -1 0 1 3 4 5 6 7 9];
  % tac, aud, nul, at500, at70, at20, at0, ta20, ta70, ta500
  
  if any(soaseq==-3)
    startsoaseq=2;
  else
    startsoaseq=1;
  end    
  
  prevmax=1;
  for ll=startsoaseq:length(soaseq) % to get max number of trials for a given condition
    prevmax=max(prevmax,length(find(soa_seq(1:nt)==soaseq(ll))));
  end
%   rtsr=nan(prevmax,max(length(soaseq),size(rtsrall,2)));
  rtsr=nan(prevmax,max(length(mapcode),size(rtsrall,2)));
  prevind=1;
%   for ll=1:length(soaseq)
  for ll=startsoaseq:length(soaseq)
    numtr=length(find(soa_seq==soaseq(ll)));
    rtsr(1:numtr,find(mapcode==soaseq(ll)))=rtsshuf(prevind: prevind+numtr-1);
    prevind=prevind+numtr;
  end
  rtsr(rtsr==rts(1))=nan; % first trial always slow
  rtsr(find(rtsr<.1))=nan; % anticipatory? (or fault with computing timetouch?);   % Throwing out trials < 100ms
  rtsr(find(rtsr>1))=nan; % asleep?   % Throwing out trials greater than 1s
  
  rtsrall=[rtsrall; rtsr];
  
end
clear rtsr
rtsr=rtsrall;
% save(['rtsr_' sub{ii} '_cond' num2str(files(bb).name(6)) '.mat'],'rtsr');
save(['rtsr_' sub{ii} '_cond' num2str(files(bb).name(7)) '.mat'],'rtsr');
clear rtsr
% end

figure(ii);errorbar(nanmean(rtsrall(:,[1 2 4:10])),nanstd(rtsrall(:,[1 2 4:10]))./sqrt(sum(~isnan(rtsrall(:,[1 2 4:10])))))
[h,p]=ttest2(rtsrall(:,1),rtsrall(:,7))
[h,p]=ttest2(rtsrall(:,2),rtsrall(:,7))

print(ii,'-deps','-noui',['RT_sub' num2str(ii)])

return

%% First , crudely collapse over conditions and subjects

clear all
if ispc
  bdir='D:\audtac\behav_data\';
  edir='D:\audtac\eeg_data\';
  ddir='D:\audtac\legomagic\diaries\';
else
  bdir='/mnt/hgfs/D/audtac/behav_data/';
  edir='/mnt/hgfs/D/audtac/eeg_data/';
  ddir='/mnt/hgfs/D/audtac/legomagic/diaries/';
end
% sub{100}='p01'; 
sub{101}='e01'; 
sub{102}='e02'; 
sub{103}='e03'; 
sub{104}='e04'; 
%%%  above is pilot, below is real
sub{105}='e05'; 
sub{106}='e06'; 
sub{107}='e07'; 
sub{108}='e08'; 
sub{109}='e09';
sub{110}='e10';
sub{111}='e11';
sub{112}='e12';
sub{113}='e13';
sub{114}='e14';
sub{115}='e15';
sub{116}='e16';
sub{117}='e17';
sub{118}='e18';
sub{119}='e19';
sub{120}='e20';
sub{121}='e21';
sub{122}='e22';
sub{123}='e23';
sub{124}='e24';
sub{125}='e25';
sub{126}='e26';
sub{127}='e27';
sub{128}='e28';
sub{129}='e29';
sub{130}='e30';
sub{131}='e31';
sub{132}='e32';

cd(bdir)

load([edir 'iikeep.mat'])

%%
rtsrall=nan(1,10);
diffms=nan(1,max(iiuse));
pcb=nan(1,max(iiuse));

bbind=1;
for bb=setdiff(union(iiSuse,iiBuse),3)
  load([bdir 'rtsr_' sub{bb+100} '_condr.mat'])
  rtsrall=[rtsrall; rtsr];
  rtsrmed(bbind,:)=nanmedian(rtsr);
  bbind=bbind+1;
  % save out some per-subject metrics of effect: percent change from baseline and difference in milliseconds
  pcb(bb)=[nanmean(reshape(rtsr(:,6:8),[size(rtsr,1)*3 1]))/nanmean(reshape(rtsr(:,1:2),[size(rtsr,1)*2 1]))]-1;
  diffms(bb)=[nanmean(reshape(rtsr(:,6:8),[size(rtsr,1)*3 1]))-nanmean(reshape(rtsr(:,1:2),[size(rtsr,1)*2 1]))];
end

save('rtgroup_pcb_diffms.mat','pcb','diffms');

% files=dir([bdir 'rtsr_e*condr.mat']);
% 
% for bb=1:length(files)
%   %     if str2num(files(bb).name(7:8))<10 | str2num(files(bb).name(7:8))>14
%   %         continue
%   %     end
%   files(bb).name
%   %     keyboard
%   load(files(bb).name);
%   rtsrall=[rtsrall; rtsr];
%   
%   % save out some per-subject metrics of effect: percent change from baseline and difference in milliseconds
%   pcb(bb)=[nanmean(reshape(rtsr(:,6:8),[size(rtsr,1)*3 1]))/nanmean(reshape(rtsr(:,1:2),[size(rtsr,1)*2 1]))]-1;
%   diffms(bb)=[nanmean(reshape(rtsr(:,6:8),[size(rtsr,1)*3 1]))-nanmean(reshape(rtsr(:,1:2),[size(rtsr,1)*2 1]))];
% end

rt_uni=reshape(rtsrall(:,[1 2 4 10]),size(rtsrall,1)*4,1);
rt_mul=reshape(rtsrall(:,[5 6 7 8 9]),size(rtsrall,1)*5,1);
[h,p]=ttest2(rt_uni,rt_mul)
rt_uni=reshape(rtsrall(:,[1 2 ]),size(rtsrall,1)*2,1);
rt_mul=reshape(rtsrall(:,[6 7 8]),size(rtsrall,1)*3,1);
[h,p]=ttest2(rt_uni,rt_mul)
p=anova1(rtsrall(:,[1 2 4:10]),{'tac alone','aud alone', '  -500ms aud first','  -70ms aud first','  -20ms aud first','simult','20ms tac first','70ms tac first','500ms tac first'},'off')
p=anova1(rtsrall(:,[1 2 5:9]),{'tac alone','aud alone', '  -70ms aud first','  -20ms aud first','simult','20ms tac first','70ms tac first'},'off')

% figure(20);errorbar(nanmean(rtsrall(:,[1 2 4:10])),nanstd(rtsrall(:,[1 2 4:10]))./sqrt(sum(~isnan(rtsrall(:,[1 2 4:10])))))
figure(20);errorbar(nanmean(rtsrall(:,[2 4:10 1])),nanstd(rtsrall(:,[2 4:10 1]))./sqrt(sum(~isnan(rtsrall(:,[2 4:10 1])))))
set(get(20,'Children'),'XTickLabel',{'', 'Aud alone', '  -500ms Aud first','  -70ms Aud first','  -20ms Aud first','Simultaneous','20ms Tac first','70ms Tac first','500ms Tac first','Tac alone',''})
xlabel('Stimulus Condition')
ylabel('Reaction Time (s)')
[h,p]=ttest2(rtsrall(:,1),rtsrall(:,7))
[h,p]=ttest2(rtsrall(:,2),rtsrall(:,7))
print(20,'-deps','-noui',['RT_allsub']);

% Now again, but 'correct' way of using mean of medians
rt_uni=reshape(rtsrmed(:,[1 2 4 10]),size(rtsrmed,1)*4,1);
rt_mul=reshape(rtsrmed(:,[5 6 7 8 9]),size(rtsrmed,1)*5,1);
[h,p]=ttest2(rt_uni,rt_mul)
rt_uni=reshape(rtsrmed(:,[1 2 ]),size(rtsrmed,1)*2,1);
rt_mul=reshape(rtsrmed(:,[6 7 8]),size(rtsrmed,1)*3,1);
[h,p]=ttest2(rt_uni,rt_mul)
p=anova1(rtsrmed(:,[1 2 4:10]),{'tac alone','aud alone', '  -500ms aud first','  -70ms aud first','  -20ms aud first','simult','20ms tac first','60ms tac first','500ms tac first'},'off')
p=anova1(rtsrmed(:,[1 2 5:9]),{'tac alone','aud alone', '  -70ms aud first','  -20ms aud first','simult','20ms tac first','60ms tac first'},'off')

% figure(21);errorbar(nanmean(rtsrmed(:,[1 2 4:10])),nanstd(rtsrmed(:,[1 2 4:10]))./sqrt(sum(~isnan(rtsrmed(:,[1 2 4:10])))))
figure(21);errorbar(nanmean(rtsrmed(:,[2 4:10 1])),nanstd(rtsrmed(:,[2 4:10 1]))./sqrt(sum(~isnan(rtsrmed(:,[2 4:10 1])))))
set(get(21,'Children'),'XTickLabel',{'', 'Aud alone', '  -500ms Aud first','  -70ms Aud first','  -20ms Aud first','Simultaneous','20ms Tac first','70ms Tac first','500ms Tac first','Tac alone',''})
xlabel('Stimulus Condition')
ylabel('Reaction Time (s)')
[h,p]=ttest2(rtsrmed(:,1),rtsrmed(:,7))
[h,p]=ttest2(rtsrmed(:,2),rtsrall(:,7))
print(21,'-deps','-noui',['RT_allsubmed']);

%%   Race model

buse=0;
for bb=setdiff(union(iiSuse,iiBuse),3:9)
  buse=buse+1;
  load([bdir 'rtsr_' sub{bb+100} '_condr.mat'])
  
  tac=rtsr(~isnan(rtsr(:,1)),1);
  aud=rtsr(~isnan(rtsr(:,2)),2);
  ms1=rtsr(~isnan(rtsr(:,4)),4);
  ms3=rtsr(~isnan(rtsr(:,5)),5);
  ms4=rtsr(~isnan(rtsr(:,6)),6);
  ms5=rtsr(~isnan(rtsr(:,7)),7);
  ms6=rtsr(~isnan(rtsr(:,8)),8);
  ms7=rtsr(~isnan(rtsr(:,9)),9);
  ms9=rtsr(~isnan(rtsr(:,10)),10);
  % Modified RaceModel so that plotflag is also the figure handle, and
  % further input for subplot
  [Xp, Yp, Zp1(:,buse), Bp1(:,buse)] = RaceModel(round(1000*tac')+500,round(1000*aud'),    round(1000*ms1'),[.1:.1:1],1,[4 6 buse]);
  [Xp, Yp, Zp3(:,buse), Bp3(:,buse)] = RaceModel(round(1000*tac')+70, round(1000*aud'),    round(1000*ms3'),[.1:.1:1],3,[4 6 buse]);
  [Xp, Yp, Zp4(:,buse), Bp4(:,buse)] = RaceModel(round(1000*tac')+20, round(1000*aud'),    round(1000*ms4'),[.1:.1:1],4,[4 6 buse]);
  [Xp, Yp, Zp5(:,buse), Bp5(:,buse)] = RaceModel(round(1000*tac'),    round(1000*aud'),    round(1000*ms5'),[.1:.1:1],5,[4 6 buse]);
  [Xp, Yp, Zp6(:,buse), Bp6(:,buse)] = RaceModel(round(1000*tac'),    round(1000*aud')+20, round(1000*ms6'),[.1:.1:1],6,[4 6 buse]);
  [Xp, Yp, Zp7(:,buse), Bp7(:,buse)] = RaceModel(round(1000*tac'),    round(1000*aud')+70, round(1000*ms7'),[.1:.1:1],7,[4 6 buse]);
  [Xp, Yp, Zp9(:,buse), Bp9(:,buse)] = RaceModel(round(1000*tac'),    round(1000*aud')+500,round(1000*ms9'),[.1:.1:1],9,[4 6 buse]);
end

hh=nan(10,9);pp=nan(10,9);
[hh(:,1),pp(:,1)]=ttest(Zp1',Bp1','tail','left'); % left-tail means Zp-Bp is negative (that MS Z was faster than Race model B)
[hh(:,3),pp(:,3)]=ttest(Zp3',Bp3','tail','left'); % left-tail means Zp-Bp is negative (that MS Z was faster than Race model B)
[hh(:,4),pp(:,4)]=ttest(Zp4',Bp4','tail','left'); % left-tail means Zp-Bp is negative (that MS Z was faster than Race model B)
[hh(:,5),pp(:,5)]=ttest(Zp5',Bp5','tail','left'); % left-tail means Zp-Bp is negative (that MS Z was faster than Race model B)
[hh(:,6),pp(:,6)]=ttest(Zp6',Bp6','tail','left'); % left-tail means Zp-Bp is negative (that MS Z was faster than Race model B)
[hh(:,7),pp(:,7)]=ttest(Zp7',Bp7','tail','left'); % left-tail means Zp-Bp is negative (that MS Z was faster than Race model B)
[hh(:,9),pp(:,9)]=ttest(Zp9',Bp9','tail','left'); % left-tail means Zp-Bp is negative (that MS Z was faster than Race model B)

% plot of % of participants who violate race model
figure(100);imagesc([mean(Zp1-Bp1<0,2) mean(Zp3-Bp3<0,2) mean(Zp4-Bp4<0,2) mean(Zp5-Bp5<0,2) mean(Zp6-Bp6<0,2) mean(Zp7-Bp7<0,2) mean(Zp9-Bp9<0,2)])
caxis([0 1]);
set(get(100,'Children'),'XTickLabel',{'AT500' 'AT70' 'AT20' 'AT0' 'TA20' 'TA70' 'TA500'})

% [Zp1-Bp1 Zp3-Bp3 Zp4-Bp4 Zp5-Bp5 Zp6-Bp6 Zp7-Bp7 Zp9-Bp9]
sum([mean(Zp1(1:2,:)-Bp1(1:2,:)); mean(Zp3(1:2,:)-Bp3(1:2,:)); mean(Zp4(1:2,:)-Bp4(1:2,:)); mean(Zp5(1:2,:)-Bp5(1:2,:)); mean(Zp6(1:2,:)-Bp6(1:2,:)); mean(Zp7(1:2,:)-Bp7(1:2,:)); mean(Zp9(1:2,:)-Bp9(1:2,:))]<0,2)
sum([mean(Zp1(1:3,:)-Bp1(1:3,:)); mean(Zp3(1:3,:)-Bp3(1:3,:)); mean(Zp4(1:3,:)-Bp4(1:3,:)); mean(Zp5(1:3,:)-Bp5(1:3,:)); mean(Zp6(1:3,:)-Bp6(1:3,:)); mean(Zp7(1:3,:)-Bp7(1:3,:)); mean(Zp9(1:3,:)-Bp9(1:3,:))]<0,2)
sum([mean(Zp1(2:3,:)-Bp1(2:3,:)); mean(Zp3(2:3,:)-Bp3(2:3,:)); mean(Zp4(2:3,:)-Bp4(2:3,:)); mean(Zp5(2:3,:)-Bp5(2:3,:)); mean(Zp6(2:3,:)-Bp6(2:3,:)); mean(Zp7(2:3,:)-Bp7(2:3,:)); mean(Zp9(2:3,:)-Bp9(2:3,:))]<0,2)

numviol=[sum(Zp1-Bp1<0,2) sum(Zp3-Bp3<0,2) sum(Zp4-Bp4<0,2) sum(Zp5-Bp5<0,2) sum(Zp6-Bp6<0,2) sum(Zp7-Bp7<0,2) sum(Zp9-Bp9<0,2)];
% myBinomTest(17,22,.5)
%     0.0169
% myBinomTest(16,22,.5)
%     0.0525
% Thus, we should use 17/22 as the threshold but 16/22 is nearly there.
figure(101);imagesc(numviol>=17)
set(get(101,'Children'),'XTickLabel',{'AT500' 'AT70' 'AT20' 'AT0' 'TA20' 'TA70' 'TA500'})


% % thresholded at 75%
% figure;imagesc([mean(Zp1-Bp1<0,2) mean(Zp3-Bp3<0,2) mean(Zp4-Bp4<0,2) mean(Zp5-Bp5<0,2) mean(Zp6-Bp6<0,2) mean(Zp7-Bp7<0,2) mean(Zp9-Bp9<0,2)]>=.75)

% % Tom reminds me that this is not valid since the data points are dependent on each other; use RM 2-way anova instead
% [hhh(1),ppp(1)]=ttest(reshape(Zp1(1:2,:)',1,44),reshape(Bp1(1:2,:)',1,44),'tail','left'); % left-tail means Zp-Bp is negative (that MS Z was faster than Race model B)
% [hhh(3),ppp(3)]=ttest(reshape(Zp3(1:2,:)',1,44),reshape(Bp3(1:2,:)',1,44),'tail','left'); % left-tail means Zp-Bp is negative (that MS Z was faster than Race model B)
% [hhh(4),ppp(4)]=ttest(reshape(Zp4(1:2,:)',1,44),reshape(Bp4(1:2,:)',1,44),'tail','left'); % left-tail means Zp-Bp is negative (that MS Z was faster than Race model B)
% [hhh(5),ppp(5)]=ttest(reshape(Zp5(1:2,:)',1,44),reshape(Bp5(1:2,:)',1,44),'tail','left'); % left-tail means Zp-Bp is negative (that MS Z was faster than Race model B)
% [hhh(6),ppp(6)]=ttest(reshape(Zp6(1:2,:)',1,44),reshape(Bp6(1:2,:)',1,44),'tail','left'); % left-tail means Zp-Bp is negative (that MS Z was faster than Race model B)
% [hhh(7),ppp(7)]=ttest(reshape(Zp7(1:2,:)',1,44),reshape(Bp7(1:2,:)',1,44),'tail','left'); % left-tail means Zp-Bp is negative (that MS Z was faster than Race model B)
% [hhh(9),ppp(9)]=ttest(reshape(Zp9(1:2,:)',1,44),reshape(Bp9(1:2,:)',1,44),'tail','left'); % left-tail means Zp-Bp is negative (that MS Z was faster than Race model B)

% Help for rmanova2:
% column1 - dependent variable
% column2 - grouping variable for subject
% column3 - grouping variable for factor 1
% column4 - grouping variable for factor 2
% rmanova2([col1 col2 col3 col4])
                
col2=[1:22 1:22 1:22 1:22]';
col3=[ones(1,44) 2*ones(1,44)]';
col4=[ones(1,22) 2*ones(1,22) ones(1,22) 2*ones(1,22)]';

col1=[reshape(Zp1(1:2,:)',1,44) reshape(Bp1(1:2,:)',1,44)]';
stats{1}=rm_anova2(col1,col2,col3,col4,{'RaceModel' 'Percentile'})
col1=[reshape(Zp3(1:2,:)',1,44) reshape(Bp3(1:2,:)',1,44)]';
stats{3}=rm_anova2(col1,col2,col3,col4,{'RaceModel' 'Percentile'})
col1=[reshape(Zp4(1:2,:)',1,44) reshape(Bp4(1:2,:)',1,44)]';
stats{4}=rm_anova2(col1,col2,col3,col4,{'RaceModel' 'Percentile'})
col1=[reshape(Zp5(1:2,:)',1,44) reshape(Bp5(1:2,:)',1,44)]';
stats{5}=rm_anova2(col1,col2,col3,col4,{'RaceModel' 'Percentile'})
col1=[reshape(Zp6(1:2,:)',1,44) reshape(Bp6(1:2,:)',1,44)]';
stats{6}=rm_anova2(col1,col2,col3,col4,{'RaceModel' 'Percentile'})
col1=[reshape(Zp7(1:2,:)',1,44) reshape(Bp7(1:2,:)',1,44)]';
stats{7}=rm_anova2(col1,col2,col3,col4,{'RaceModel' 'Percentile'})
col1=[reshape(Zp9(1:2,:)',1,44) reshape(Bp9(1:2,:)',1,44)]';
stats{9}=rm_anova2(col1,col2,col3,col4,{'RaceModel' 'Percentile'})

% RaceModel factor is significant for 4,5,6 with no interaction.
% Interaction significant for 7; follow up t-test per percentile:
[hhh(71),ppp(71)]=ttest(Zp7(1,:),Bp7(1,:),'tail','left'); % left-tail means Zp-Bp is negative (that MS Z was faster than Race model B)
[hhh(72),ppp(72)]=ttest(Zp7(2,:),Bp7(2,:),'tail','left'); % left-tail means Zp-Bp is negative (that MS Z was faster than Race model B)
% 20th% bin significant but not 10th% bin

% RaceModel factor is significant for 1 with no interaction.
[hhh(1),ppp(1)]=ttest(reshape(Zp1(1:2,:)',1,44),reshape(Bp1(1:2,:)',1,44),'tail','left'); % left-tail means Zp-Bp is negative (that MS Z was faster than Race model B)
% goes in opposite direction

%%  See eeg_legomagic_brainBehaviourCorrelations
% for computation of AT70 versus min(A, T+70)

load([ddir 'rt_allsubj.mat'],'rt*');

[hh,pp]=ttest(squeeze(rt_msMminshiftuni(:,3,:))')
1000*nanmean(squeeze(rt_msMminshiftuni(:,3,:)),2)

%% correlation to spindles
% using race model and subtraction

load([edir 'spindles.mat'])

rmodel=[mean(Zp1(1:3,:)-Bp1(1:3,:)); mean(Zp3(1:3,:)-Bp3(1:3,:)); mean(Zp4(1:3,:)-Bp4(1:3,:)); mean(Zp5(1:3,:)-Bp5(1:3,:)); mean(Zp6(1:3,:)-Bp6(1:3,:)); mean(Zp7(1:3,:)-Bp7(1:3,:)); mean(Zp9(1:3,:)-Bp9(1:3,:))];
iruse=setdiff(union(iiSuse,iiBuse),3:9);
corruse=intersect(find(~isnan(spindles(:,1))),iruse);

[cc,pp]=corr(spindles(corruse,:),rmodel(:,dsearchn(iruse',corruse))');
% nothing significant (after correction)

load([ddir 'rt_allsubj.mat'],'rt_msMminshiftuni');
soalist=[1 3 4 5 6 7 9];
rtdiff=squeeze(rt_msMminshiftuni(soalist,3,:));
iruse=iiSuse;
corruse=setdiff(intersect(find(~isnan(spindles(:,1))),iruse),[8 9]);
[cc,pp]=corr(spindles(corruse,:),rtdiff(:,dsearchn(iruse',corruse))');
% nothing significant (after correction)



