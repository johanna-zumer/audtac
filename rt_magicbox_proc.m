clear all

bdir='/mnt/hgfs/D/audtac/behav_data/';
ddir='/mnt/hgfs/D/audtac/magicbox/audiotactile_pneumaticmachine/diaries/';
cd(bdir)

sub{1}='jz';
sub{2}='p01';
sub{3}='p03';
sub{4}='p04';
sub{5}='p05'; % jz 29/04/14, with 4 conditions (epinoise & earmuff)
sub{6}='p06'; % 55695 30/04/14 with conditions 4 and 2
sub{7}='p07'; % 54265 30/04/14 with conditions 2 4 3
sub{8}='p08'; % 61651 30/04/14 with conditions 2 4 1
sub{9}='p09'; % 58687 01/05/14 with conditions 1 2


ii=9;
%%

switch ii
    case 1
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
        clear rtsr;
        rtsr=rtsrall;
        save(['rtsr_' sub{ii} '.mat'],'rtsr');
        rtsrall=rtsr;
        clear rtsr info fin
        
    case {2, 3, 4}
        try
            load([sub{ii} '_p3.mat']);
        catch
            load([sub{ii} '.mat']);
        end
        for ll=1:length(info.prsout),
            resptime(ll)=info.prsout{ll}(7)-info.start_time;
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
        %         alltime(find(diff(alltime)<.1)+1)=nan;
        %         alltime=sort(alltime);
        %         alltime=alltime(~isnan(alltime));
        
        if ii==3
            resptimebak=resptime;
            resptime=resptime([1:6 8:521]);
        end
        clear rts
        allcnt=0;
        for ll=1:length(resptime)
            while resptime(ll)-alltime(ll+allcnt) > .9
                allcnt=allcnt+1;
            end
            if resptime(ll)-alltime(ll+allcnt) < 0
                allcnt=allcnt-1;
            end
            rts(ll)=resptime(ll)-alltime(ll+allcnt);
        end
        figure;plot(info.soa_seq,rts,'o');axis([-3 13 0 8])
        
        [aa,bb]=sort(info.soa_seq);
        rtsshuf=[rts(bb)];
        switch ii
            case 2
                rtsr=reshape(rts(bb),size(rts,2)/length(unique(info.soa_seq)),length(unique(info.soa_seq)));
            case {3, 4}
                soaseq=unique(info.soa_seq);
                prevmax=1;
                for ll=1:length(soaseq)
                    prevmax=max(prevmax,length(find(info.soa_seq==soaseq(ll))));
                end
                rtsr=nan(prevmax,length(soaseq));
                prevind=1;
                for ll=1:length(soaseq)
                    numtr=length(find(info.soa_seq==soaseq(ll)));
                    rtsr(1:numtr,ll)=rtsshuf(prevind: prevind+numtr-1);
                    prevind=prevind+numtr;
                end
                
        end
        
        rtsr(rtsr==rts(1))=nan; % first trial always slow
        rtsr(find(rtsr<.1))=nan; % anticipatory?
        
        save(['rtsr_' sub{ii} '.mat'],'rtsr');
        rtsrall=rtsr;
        clear rtsr info fin
    case {5, 6, 7 , 8, 9}
        % after seeing the data, in terms of possible to analyze or if
        % subject was paying attenion / following instructions, still include:
        % subject 5: cond 1 2 3
        % subject 6: cond 4
        % subject 7: cond 2 4
        % subject 8: cond 1 2
        % subject 9: cond 1 2
        
        files=dir([ddir sub{ii} '*.mat']);
        for bb=1:length(files)
            clear resptime rts aa ll rtsshuf
            load([ddir files(bb).name]);
            
            for ll=1:length(info.prsout),
                resptime(ll)=info.prsout{ll}(15)-info.start_time;
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
            %         alltime(find(diff(alltime)<.1)+1)=nan;
            %         alltime=sort(alltime);
            %         alltime=alltime(~isnan(alltime));
            
            
            % Weird: for p05-p09 it seems some MS trials were only
            % uni-sensory. Change soa_seq to reflect this. aud-alone is -1
            soa_seq=info.soa_seq;
            soa_seq(find(isnan(info.aud_timings) & info.soa_seq>0))=-2; % had been all 4
            soa_seq(find(isnan(info.tac_timings) & info.soa_seq>0))=-1; % had been all 5
            
            
            % establish time of 'alltime' that occurs first for each trial
            acnt=1;
            for ll=1:length(soa_seq)
                allfirst(ll)=alltime(acnt);
                if soa_seq(ll)>0
                    acnt=acnt+2;
                else
                    acnt=acnt+1;
                end
            end
            
            
            % removing double-presses within one trial (if was pressed within 1 second of the last press)
            if ii==6 & bb==1
                [resptime(setdiff(1:306,[60 145 173 206 254 269]))-allfirst; soa_seq];
                resptime=resptime(setdiff(1:306,[60 145 173 206 254 269]));
            end

            % works for ii==7 bb==1 and bb==2
            for ll=1:length(allfirst)
                tmp=find(resptime-allfirst(ll)>0);
                respfirst(ll)=resptime(tmp(1));
            end
            rts=respfirst-allfirst;
            [rts; soa_seq]
            

            if 1
                clear rts
                rcnt=1;
                acnt=1;
                for ll=1:length(soa_seq)
                    if resptime(rcnt)-allfirst(acnt) < -0.5
                        acnt=acnt-1;
                    end
                    %                 if resptime(rcnt)-alltime(acnt) < -0.5
                    %                 if resptime(rcnt)-allfirst(acnt) < -0.5
                    %                     acnt=acnt-1;
                    %                 end
                    %                 rts(ll)=resptime(rcnt)-alltime(acnt);
                    rts(ll)=resptime(rcnt)-allfirst(acnt);
                    %                 keyboard
                    
                    rcnt=rcnt+1;
                    acnt=acnt+1;
                    %                 if soa_seq(ll)>0
                    %                     acnt=acnt+1;
                    %                 end
                end
                [rts; soa_seq]
            else
                
                %             resptime=resptime(setdiff(1:length(resptime),find(diff(resptime)<1.2)));
                resptime=resptime(setdiff(1:length(resptime),find(diff(resptime)<1.2)+1));

                while length(resptime)>length(info.tac_timings)
                    resptime=resptime(setdiff(1:length(resptime), find(diff(resptime)<min(diff(resptime))+.001)+1 ));
                end
                
                clear rts
                allcnt=0;
                for ll=1:length(resptime)
                    while resptime(ll)-alltime(ll+allcnt) > .9
                        allcnt=allcnt+1;
                    end
                    if resptime(ll)-alltime(ll+allcnt) < 0.17
                        allcnt=allcnt-1;
                    end
                    try
                        rts(ll)=resptime(ll)-alltime(ll+allcnt);
                    catch
                        ll
                        allcnt=allcnt+1;
                        rts(ll)=resptime(ll)-alltime(ll+allcnt);
                    end
                end
            end
            
            
            
            figure;plot(soa_seq,rts,'o');axis([-3 13 0 8])
            figure;plot(soa_seq,rts,'o');axis([-3 13 0 20])
            [rts; soa_seq]
            %             figure;plot(info.soa_seq,rts,'o');axis([-3 13 0 8])
            %             figure;plot(info.soa_seq,rts,'o');axis([-3 13 0 20])
            %             [rts; info.soa_seq]
            
            [aa,cc]=sort(soa_seq);
            rtsshuf=[rts(cc)];
            
            soaseq=unique(soa_seq);
            prevmax=1;
            for ll=1:length(soaseq)
                prevmax=max(prevmax,length(find(soa_seq==soaseq(ll))));
            end
            rtsr=nan(prevmax,length(soaseq));
            prevind=1;
            for ll=1:length(soaseq)
                numtr=length(find(soa_seq==soaseq(ll)));
                rtsr(1:numtr,ll)=rtsshuf(prevind: prevind+numtr-1);
                prevind=prevind+numtr;
            end
            rtsr(rtsr==rts(1))=nan; % first trial always slow
            rtsr(find(rtsr<.1))=nan; % anticipatory?
            %             if ii==6 & bb==2;  and ii==8 & bb==2
            %                 rtsr(find(rtsr>8))=nan; % falling asleep?
            %                 rtsr(find(rtsr>8))=nan; % falling asleep?
            % if ii==9 and bb==1;   rtsr(find(rtsr>4))=nan;
            % if ii==9 and bb==2;   rtsr(find(rtsr>2))=nan;
            %             end
            
            
            save(['rtsr_' sub{ii} '_cond' num2str(files(bb).name(6)) '.mat'],'rtsr');
            rtsrall=rtsr;
            clear rtsr info fin
        end
end


%%  First, crudely collapse over conditions and subjects
clear all
bdir='/mnt/hgfs/D/audtac/behav_data/';
files=dir([bdir 'rtsr_p*cond*.mat']);
rtsrall=nan(40*length(files),10);
for ii=1:length(files)
    clear rtsr
    load([bdir files(ii).name]);
    rtsrall((ii-1)*40+1: (ii-1)*40+size(rtsr,1) , :)=rtsr;    
end
rt_uni=reshape(rtsrall(:,[1 2 4 10]),size(rtsrall,1)*4,1);
rt_mul=reshape(rtsrall(:,[5 6 7 8 9]),size(rtsrall,1)*5,1);
[h,p]=ttest2(rt_uni,rt_mul)
rt_uni=reshape(rtsrall(:,[1 2 ]),size(rtsrall,1)*2,1);
rt_mul=reshape(rtsrall(:,[6 7 8]),size(rtsrall,1)*3,1);
[h,p]=ttest2(rt_uni,rt_mul)
figure;errorbar(nanmean(rtsrall(:,[1 2 4:10])),nanstd(rtsrall(:,[1 2 4:10]))/sqrt(300))
p=anova1(rtsrall(:,[1 2 4:10]),{'tac alone','aud alone', '  -500ms aud first','  -70ms aud first','  -20ms aud first','simult','20ms tac first','60ms tac first','500ms tac first'},'off');
p=anova1(rtsrall(:,[1 2 5:9]),{'tac alone','aud alone', '  -70ms aud first','  -20ms aud first','simult','20ms tac first','60ms tac first'},'off')

[h,p]=ttest2(rtsrall(:,1),rtsrall(:,7))
[h,p]=ttest2(rtsrall(:,2),rtsrall(:,7))


% still to do: tests within condition type

%%
% case ii = 5 (different conditions)
switch ii
    case 5
        files=dir([bdir 'rtsr_' sub{ii} '*.mat']);
        rt=nan(30,10,4);
        for bb=1:4
            try
                load(files(bb).name); % this should be in order of conditions now not collection
                rt(:,:,bb)=rtsr;
            end
            clear rtsr
            rtsrall=rt(:,:,bb);
            rtsrall(rtsrall>1.5)=nan; % is this fair?
            
            pallcond(bb,ii)=anova1(rtsrall(:,[1 2 4:10]),{'tac alone','aud alone', '  -500ms aud first','  -70ms aud first','  -20ms aud first','simult','20ms tac first','60ms tac first','500ms tac first'},'off');
            %             p=anova1(rtsrall(:,[1 2 5:9]),{'tac alone','aud alone', '  -70ms aud first','  -20ms aud first','simult','20ms tac first','60ms tac first'},'off');
            pnearsync(bb,ii)=anova1(rtsrall(:,[5:9]),{ '  -70ms aud first','  -20ms aud first','simult','20ms tac first','60ms tac first'},'off');
            rt_uni=reshape(rtsrall(:,[1 2 4 10]),size(rtsrall,1)*4,1);
            rt_mul=reshape(rtsrall(:,[5 6 7 8 9]),size(rtsrall,1)*5,1);
            [h,p4uni5ns(bb,ii)]=ttest2(rt_uni,rt_mul) % uni-modal vs multi-modal different
            p=anova1(rtsrall(:,[1 2 7]),{'tac alone','aud alone','simult'},'off');
            
            rt_uni=reshape(rtsrall(:,[1 2]),size(rtsrall,1)*2,1);
            rt_mul=reshape(rtsrall(:,[5 6 7 8 9]),size(rtsrall,1)*5,1);
            [h,p]=ttest2(rt_uni,rt_mul) % uni-modal vs multi-modal different
            
            % with +/- 500 ms counting as 'uni' and ignoring +/- 70ms
            rt_uni=reshape(rtsrall(:,[1 2 4 10]),size(rtsrall,1)*4,1);
            rt_mul=reshape(rtsrall(:,[6 7 8]),size(rtsrall,1)*3,1);
            rt_con=reshape(rtsrall(:,[5 9]),size(rtsrall,1)*2,1);
            
            p=anova1([rt_uni; rt_mul; rt_con],[repmat({'unisensory'},length(rt_uni),1); repmat({'multisens <20ms'},length(rt_mul),1); repmat({'multisens 70ms'},length(rt_con),1)],'off');
            [h,p]=ttest2(rt_uni,rt_mul) % uni-modal vs multi-modal different
            [h,p]=ttest2(rt_uni,rt_con) % uni-modal vs multi-modal different
            [h,p]=ttest2(rt_con,rt_mul) % uni-modal vs multi-modal different
            
        end
end



% for ii=3 and ii=4 (collected on same day, subject JZ)
load rtsr_p03.mat
rtsr1=rtsr;
load rtsr_p04.mat
rtsr2=[rtsr; rtsr1];
clear rtsr rtsr1
rtsr=rtsr2;clear rtsr2

rtsr(rtsr>2.05)=nan; % is this fair?
rtsrall=rtsr;

% ii=>2
p=anova1(rtsrall,{'tac alone','aud alone', 'null', '  -500ms aud first','  -70ms aud first','  -20ms aud first','simult','20ms tac first','60ms tac first','500ms tac first'})
p=anova1(rtsrall(:,[1 2 4:10]),{'tac alone','aud alone', '  -500ms aud first','  -70ms aud first','  -20ms aud first','simult','20ms tac first','60ms tac first','500ms tac first'})
p=anova1(rtsrall(:,[1 2 5:9]),{'tac alone','aud alone', '  -70ms aud first','  -20ms aud first','simult','20ms tac first','60ms tac first'})

p=anova1(rtsrall(:,[5:9]),{ '  -70ms aud first','  -20ms aud first','simult','20ms tac first','60ms tac first'})

% with +/- 500 ms counting as 'uni'
rt_uni=reshape(rtsrall(:,[1 2 4 10]),size(rtsrall,1)*4,1);
rt_mul=reshape(rtsrall(:,[5 6 7 8 9]),size(rtsrall,1)*5,1);
[h,p]=ttest2(rt_uni,rt_mul) % uni-modal vs multi-modal different

p=anova1(rtsrall(:,[1 2 7]),{'tac alone','aud alone','simult'})

rt_uni=reshape(rtsrall(:,[1 2]),size(rtsrall,1)*2,1);
rt_mul=reshape(rtsrall(:,[5 6 7 8 9]),size(rtsrall,1)*5,1);
[h,p]=ttest2(rt_uni,rt_mul) % uni-modal vs multi-modal different

% with +/- 500 ms counting as 'uni' and ignoring +/- 70ms
rt_uni=reshape(rtsrall(:,[1 2 4 10]),size(rtsrall,1)*4,1);
rt_mul=reshape(rtsrall(:,[6 7 8]),size(rtsrall,1)*3,1);
rt_con=reshape(rtsrall(:,[5 9]),size(rtsrall,1)*2,1);

p=anova1([rt_uni; rt_mul; rt_con],[repmat({'unisensory'},length(rt_uni),1); repmat({'multisens <20ms'},length(rt_mul),1); repmat({'multisens 70ms'},length(rt_con),1)]);

[h,p]=ttest2(rt_uni,rt_mul) % uni-modal vs multi-modal different
[h,p]=ttest2(rt_uni,rt_con) % uni-modal vs multi-modal different
[h,p]=ttest2(rt_con,rt_mul) % uni-modal vs multi-modal different


% ii=1
% A T -60 0 20 60
p=anova1(rtsrall,{'tac alone','aud alone','  -60ms aud first','simult','20ms tac first','60ms tac first'})
p=anova1(rtsrall(:,[1 3 4 5 6]),{'tac alone',' -60ms aud first','simult','20ms tac first','60ms tac first'})
p=anova1(rtsrall(:,[2 3 4 5 6]),{'aud alone',' -60ms aud first','simult','20ms tac first','60ms tac first'})
p=anova1(rtsrall(:,[3 4 5 6]),{'-60ms aud first','simult','20ms tac first','60ms tac first'})

p=anova1(rtsrall(:,[1 2 4 5]),{'tac alone','aud alone','simult','20ms tac first'})

rt_uni=reshape(rtsrall(:,1:2),70*2,1);
rt_mul=reshape(rtsrall(:,3:6),70*4,1);
[h,p]=ttest2(rt_uni,rt_mul) % uni-modal vs multi-modal different

% however
p=anova1(rtsrall(:,[1 2 4]),{'tac alone','aud alone','simult'})





