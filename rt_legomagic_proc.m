% processing of RT data from lego-device with chinrest and Notts EarTone
clear a* r* b* f* h* p* d* s*

if ispc
  bdir='D:\audtac\behav_data\';
  ddir='D:\audtac\legomagic\diaries\';
else
  bdir='/mnt/hgfs/D/audtac/behav_data/';
  ddir='/mnt/hgfs/D/audtac/legomagic/diaries/';
end

cd(bdir)

sub{10}='p10'; % 51520 02/05/14 with condition 3
sub{11}='p11'; % 53188 02/05/14 with condition 3
sub{12}='p12'; % 61768 06/05/14 with condition 3
sub{13}='p13'; % 57811 06/05/14 with condition 4 (but epi in only left ear)
sub{14}='p14'; % 60811 06/05/14 with condition 4 (but epi in only left ear)
sub{15}='p15'; % jz
sub{16}='p16'; % jz
sub{17}='p17'; % 58054 16/05/14 with condition 2 (epi noise, no ear muffs)
sub{18}='p18'; % 61660 16/05/14
sub{19}='p19'; % 54046 19/05/14
sub{20}='p20'; % 51736 19/05/14
sub{21}='p21'; % 61762 20/05/14
sub{22}='p22'; % 46708 20/05/14
sub{23}='p23'; % 61495 20/05/14
sub{24}='p24'; % 54142 20/05/14
sub{25}='p25'; % 61840 20/05/14

sub{101}='e01'; % ab.m. 21/05/14
sub{102}='e02'; % ma.a. 04/06/14
sub{103}='e03'; % ag.m. 10/06/14
sub{104}='e04'; % re.g. 17/06/14
%%%  above is pilot, below is real
sub{105}='e05'; % ma.a. 25/06/14
sub{106}='e06'; % e.u.  01/07/14
sub{107}='e07'; % a.s.  07/07/14
sub{108}='e08'; % k.t.  09/07/14
sub{109}='e09';% d.a.  14/07/14
sub{110}='e10';% k.l.  15/07/14
sub{111}='e11';% ab.m.  16/07/14
sub{112}='e12';% b.s.  17/07/14


ii=102;

%%

% note: check fin.aud delay

% for ii=101
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
  
  
  if ~isfield(fin,'opticSR')
    cnt=1;
    opticSR=0;
    while opticSR==0
      try
        opticSR=diff(info.lighttime{cnt}(1:2));
      end
      cnt=cnt+1;
    end
  else
    opticSR=fin.opticSR/1000;
  end
  
  
  for nn=1:nt
    lighttime=0:opticSR:opticSR*[length(info.lightsensor{nn})-1];
    if ~isempty(info.lightsensor{nn})
      % this threshold of 80% is arbitrary...should be tested (70? 90?)
      tmp=find(info.lightsensor{nn}> .8*[max(info.lightsensor{nn})-min(info.lightsensor{nn})]+min(info.lightsensor{nn}) );
      %                 touchdelay(nn)=info.lighttime{nn}(tmp(1))-info.lighttime{nn}(1);
      touchdelay(nn)=lighttime(tmp(1))-lighttime(1);
    else
      touchdelay(nn)=nan;
    end
  end
  tmp=nanmedian(touchdelay);
  touchdelay(find(touchdelay>tmp+.05))=nan; % are these too strict?
  touchdelay(find(touchdelay<tmp-.05))=nan; % ?
  touchdelay(find(touchdelay<.1))=nan;
  
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
  
  for nn=1:nt
    if soaeff(nn)<-.25 & soareal(nn)>0
      soareal(nn)=1;
    elseif soaeff(nn)<-.05 & soareal(nn)>0
      soareal(nn)=3;
    elseif soaeff(nn)<-.01 & soareal(nn)>0
      soareal(nn)=4;
    elseif soaeff(nn)<.01 & soareal(nn)>0
      soareal(nn)=5;
    elseif soaeff(nn)<.05 & soareal(nn)>0
      soareal(nn)=6;
    elseif soaeff(nn)<.25 & soareal(nn)>0
      soareal(nn)=7;
    elseif soaeff(nn)>=.25 & soareal(nn)>0
      soareal(nn)=9;
    end
  end
  soa_seq=soareal;
  %         soareal(isnan(soaeff))=nan;
  
  %         if 0
  %         else % we should use the real times unless we know it's perfect
  %             soa_seq=soareal;
  %         end
  
  firststimtime=min(info.audio_time(1:nt)+.019-info.start_time, info.valve_time(1:nt)+touchdelay-info.start_time);
  
  rts=resptime-firststimtime;
  
  
  [aa,cc]=sort(soa_seq);
  rtsshuf=[rts(cc)];
  
  soaseq=unique(soa_seq(~isnan(soa_seq)));
  mapcode=[-2 -1 0 1 3 4 5 6 7 9];
  prevmax=1;
  for ll=1:length(soaseq)
    prevmax=max(prevmax,length(find(soa_seq(1:nt)==soaseq(ll))));
  end
  rtsr=nan(prevmax,max(length(soaseq),size(rtsrall,2)));
  prevind=1;
  for ll=1:length(soaseq)
    numtr=length(find(soa_seq==soaseq(ll)));
    rtsr(1:numtr,find(mapcode==soaseq(ll)))=rtsshuf(prevind: prevind+numtr-1);
    prevind=prevind+numtr;
  end
  rtsr(rtsr==rts(1))=nan; % first trial always slow
  rtsr(find(rtsr<.1))=nan; % anticipatory? (or fault with computing timetouch?)
  rtsr(find(rtsr>1))=nan; % asleep?
  
  rtsrall=[rtsrall; rtsr];
  
end
clear rtsr
rtsr=rtsrall;
save(['rtsr_' sub{ii} '_cond' num2str(files(bb).name(6)) '.mat'],'rtsr');
clear rtsr
% end

figure;errorbar(nanmean(rtsrall(:,[1 2 4:10])),nanstd(rtsrall(:,[1 2 4:10]))./sqrt(sum(~isnan(rtsrall(:,[1 2 4:10])))))
[h,p]=ttest2(rtsrall(:,1),rtsrall(:,7))
[h,p]=ttest2(rtsrall(:,2),rtsrall(:,7))

%% First , crudely collapse over conditions and subjects

clear all
bdir='/mnt/hgfs/D/audtac/behav_data/';
files=dir([bdir 'rtsr_p*cond*.mat']);

rtsrall=nan(1,10);
for bb=16:25
  %     if str2num(files(bb).name(7:8))<10 | str2num(files(bb).name(7:8))>14
  %         continue
  %     end
  files(bb).name
  %     keyboard
  load(files(bb).name);
  rtsrall=[rtsrall; rtsr];
end

rt_uni=reshape(rtsrall(:,[1 2 4 10]),size(rtsrall,1)*4,1);
rt_mul=reshape(rtsrall(:,[5 6 7 8 9]),size(rtsrall,1)*5,1);
[h,p]=ttest2(rt_uni,rt_mul)
rt_uni=reshape(rtsrall(:,[1 2 ]),size(rtsrall,1)*2,1);
rt_mul=reshape(rtsrall(:,[6 7 8]),size(rtsrall,1)*3,1);
[h,p]=ttest2(rt_uni,rt_mul)
p=anova1(rtsrall(:,[1 2 4:10]),{'tac alone','aud alone', '  -500ms aud first','  -70ms aud first','  -20ms aud first','simult','20ms tac first','60ms tac first','500ms tac first'},'off')
p=anova1(rtsrall(:,[1 2 5:9]),{'tac alone','aud alone', '  -70ms aud first','  -20ms aud first','simult','20ms tac first','60ms tac first'},'off')

figure;errorbar(nanmean(rtsrall(:,[1 2 4:10])),nanstd(rtsrall(:,[1 2 4:10]))./sqrt(sum(~isnan(rtsrall(:,[1 2 4:10])))))
[h,p]=ttest2(rtsrall(:,1),rtsrall(:,7))
[h,p]=ttest2(rtsrall(:,2),rtsrall(:,7))

