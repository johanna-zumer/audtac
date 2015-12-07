
cd C:\Users\jmz\Documents\bham\audtac\magicbox\audiotactile_pneumaticmachine\diaries

for ll=[10:16]
    load(['a' num2str(ll) '.mat']);
    for ii=1:length(info.prsout),
        resptime(ii)=info.prsout{ii}(7)-info.start_time;
    end
    resptime=resptime(2:end); % the first is the one to make it start
    
    audtime=info.audio_time-info.start_time;
    valtime=info.valve_time-info.start_time;
    try
        nultime=info.null_time-info.start_time;
    catch
        nultime=[];
    end
    alltime=sort([nultime audtime valtime]);
    alltime(find(diff(alltime)<.1)+1)=nan;
    alltime=sort(alltime);
    alltime=alltime(~isnan(alltime));
    
    keyboard
    switch ll
        case 10
            rts(1:39)=resptime(1:39)-alltime(1:39);
            rts(40:57)=resptime(41:58)-alltime(40:57);
            rts(58:70)=resptime(60:72)-alltime(58:70);
        case {11 12 13 14 15}
            rts=resptime-alltime;
        case 16
            rts(1:26)=resptime(1:26)-alltime(1:26);
            rts(27:60)=resptime(28:61)-alltime(27:60);
    end
    
    [aa,bb]=sort(info.soa_seq);
    [aa; rts(bb)];
    rtsr=reshape(rts(bb),10,size(rts,2)/10);
    
    rtsr(rtsr==rts(1))=nan;
    keyboard
    switch ll
        case 10
            % aberant outlier
            rtsr(9,1)=nan;
            % also exclude first trial
        case {11 16}
        case 12
            rtsr(5,2)=nan; % <100ms (anticipatory?)
        case 13
            rtsr(4,3)=nan; % <100ms (anticipatory?)
        case 14
            rtsr(4,4)=nan; % <100ms (anticipatory?)
            rtsr(6,6)=nan; % <100ms (anticipatory?)
        case 15
            rtsr(4,4)=nan; % <100ms (anticipatory?)
    end
    
    save(['rtsr' num2str(ll) '.mat'],'rtsr');
    clear all
end

rtsrall=nan(70,6);
for ll=10:16
    load(['rtsr' num2str(ll) '.mat'],'rtsr');
    if ll==10
        rtsrall(1:10,:)=rtsr(:,[1 2 4 5 6 7]);
    else
        rtsrall( [(ll-10)*10 + 1] : [(ll-10)*10 + 10], :)=rtsr;
    end
end
rtsrall(rtsrall<.1)=nan;
p=anova1(rtsrall,{'tac alone','aud alone','  -60ms aud first','simult','20ms tac first','60ms tac first'})

p=anova1(rtsrall(:,[1 3 4 5 6]),{'tac alone',' -60ms aud first','simult','20ms tac first','60ms tac first'})

p=anova1(rtsrall(:,[2 3 4 5 6]),{'aud alone',' -60ms aud first','simult','20ms tac first','60ms tac first'})

p=anova1(rtsrall(:,[3 4 5 6]),{'-60ms aud first','simult','20ms tac first','60ms tac first'})

rt_uni=reshape(rtsrall(:,1:2),70*2,1);
rt_mul=reshape(rtsrall(:,3:6),70*4,1);
[h,p]=ttest2(rt_uni,rt_mul) % uni-modal vs multi-modal different

% however
p=anova1(rtsrall(:,[1 2 4]),{'tac alone','aud alone','simult'})





