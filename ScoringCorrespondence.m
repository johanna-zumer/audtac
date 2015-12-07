clear all
%% path to JZ's data
JzPath='E:\JohannaProject\sleep_scores\sleep_scores\';
TwPath='E:\JohannaProject\sleep_scores\sleep_scores\TW\';

subs=dir([JzPath 'e*']);
subsTwWake=dir([TwPath 's*']);
subsTwSleep=dir([TwPath 'b*']);
%% wake data
for subNum=5:28
    
    JzData=load([JzPath subs(subNum).name '\CRC_wake.mat']);
    TwData=load([TwPath subsTwWake(subNum).name]);
    ConsensusData=load([JzPath subs(subNum).name '\CRC_wake_withTom.mat']);
    ctr=0;
    % jz vs. consensus
    for epoch=1:length(JzData.CRC.score{1})
       if JzData.CRC.score{1}(epoch)== ConsensusData.CRC.score{1}(epoch)
           ctr=ctr+1;
       else; end
    end
    JzConsensusCtr(subNum,1)=ctr;
    JzConsensusTotal(subNum,1)=length(JzData.CRC.score{1});
    
    % jz vs. tw
    ctr=0;
    for epoch=1:length(TwData.score)
       if JzData.CRC.score{1}(epoch)== TwData.score(epoch)
           ctr=ctr+1;
       else; end
    end
    JzTwCtr(subNum,1)=ctr;
    JzTwTotal(subNum,1)=length(TwData.score);
    
    % tw vs. consensus
    ctr=0;
    for epoch=1:length(TwData.score)
       if ConsensusData.CRC.score{1}(epoch)== TwData.score(epoch)
           ctr=ctr+1;
       else; end
    end
    TwConsensusCtr(subNum,1)=ctr;
    TwConsensusTotal(subNum,1)=length(TwData.score);
    
    
    clear *Data
end
JzTwOverallWake=sum(JzTwCtr(2:16))/sum(JzTwTotal(2:16));
JzConsensusOverallWake=sum(JzConsensusCtr(2:16))/sum(JzConsensusTotal(2:16));
TwConsensusOverallWake=sum(TwConsensusCtr(2:16))/sum(TwConsensusTotal(2:16));

%% for sleep
ctrS1=0
for subNum=5:28
    TwData=load([TwPath subsTwSleep(subNum).name]);
    JzData=load([JzPath subs(subNum).name '\CRC_sleep.mat']);
    ConsensusData=load([JzPath subs(subNum).name '\CRC_sleep_withTom.mat']);
    ctr=0;
    for epoch=1:length(JzData.CRC.score{1})
       if JzData.CRC.score{1}(epoch)== ConsensusData.CRC.score{1}(epoch)
           ctr=ctr+1;
       else; end
    end
    JzConsensusCtrSleep(subNum,1)=ctr;
    JzConsensusTotalSleep(subNum,1)=length(JzData.CRC.score{1});
    
        % jz vs. tw
    ctr=0;
    for epoch=1:length(TwData.score)
       if JzData.CRC.score{1}(epoch)== TwData.score(epoch)
           ctr=ctr+1;
       else; end
    end
    JzTwCtrSleep(subNum,1)=ctr;
    JzTwTotalSleep(subNum,1)=length(TwData.score);
    
    % tw vs. consensus
    ctr=0;
    for epoch=1:length(TwData.score)
       if ConsensusData.CRC.score{1}(epoch)== TwData.score(epoch)
           ctr=ctr+1;
       else; end
    end
    TwConsensusCtrSleep(subNum,1)=ctr;
    TwConsensusTotalSleep(subNum,1)=length(TwData.score);
end
JzTwOverallSleep=sum(JzTwCtrSleep(2:16))/sum(JzTwTotalSleep(2:16));
JzConsensusOverallSleep=sum(JzConsensusCtrSleep(2:16))/sum(JzConsensusTotalSleep(2:16));
TwConsensusOverallSleep=sum(TwConsensusCtrSleep(2:16))/sum(TwConsensusTotalSleep(2:16));



