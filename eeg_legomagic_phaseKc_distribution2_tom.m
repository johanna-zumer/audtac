ang_bins=4;
sleepStage=[10 11 12 13]; 
condTypes=length(tlock_tac);
ctr(1:23,1:4,1:condTypes)=0;
ctrKc(1:23,1:4,1:condTypes)=0;

TT=3;

%% for tlock_tac
for sleep=sleepStage(1:end)
     for conds=1:condTypes
         if isfield(tlock_tac{conds,TT,sleep},'hilbert')
            tlock_tac{conds,TT,sleep}.alpha.upstateKC.KC=[];
            tlock_tac{conds,TT,sleep}.alpha.upmidstateKC.KC=[];
            tlock_tac{conds,TT,sleep}.alpha.downmidstateKC.KC=[];
            tlock_tac{conds,TT,sleep}.alpha.downstateKC.KC=[];
            
            tlock_tac{conds,TT,sleep}.alpha.upstateKC.KC_all=[];
            tlock_tac{conds,TT,sleep}.alpha.upmidstateKC.KC_all=[];
            tlock_tac{conds,TT,sleep}.alpha.downmidstateKC.KC_all=[];
            tlock_tac{conds,TT,sleep}.alpha.downstateKC.KC_all=[];
            
            tlock_tac{conds,TT,sleep}.alpha.upstateKC.sw=[];
            tlock_tac{conds,TT,sleep}.alpha.upmidstateKC.sw=[];
            tlock_tac{conds,TT,sleep}.alpha.downmidstateKC.sw=[];
            tlock_tac{conds,TT,sleep}.alpha.downstateKC.sw=[];
            
            tlock_tac{conds,TT,sleep}.alpha.upstateKC.delta=[];
            tlock_tac{conds,TT,sleep}.alpha.upmidstateKC.delta=[];
            tlock_tac{conds,TT,sleep}.alpha.downmidstateKC.delta=[];
            tlock_tac{conds,TT,sleep}.alpha.downstateKC.delta=[];
            
            for i=1:size(tlock_tac{conds,TT,sleep}.hilbert.alphaT0.angle,1)
                
                    if [tlock_tac{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)>pi/4 && tlock_tac{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)<3*pi/4] %% up-state
                        ctr(sleep,1,conds)=ctr(sleep,1,conds)+1;
                        % is there a post-stim KC of any kind?
                            if [tlock_tac{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_tac{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_tac{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_tac{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                                    ctrKc(sleep,1,conds)=ctrKc(sleep,1,conds)+1; % if there is, count it!
                            end
                                %% now make a list of onsets to compare distributions - ignore, not used any more!
                                if [tlock_tac{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_tac{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_tac{conds,TT,sleep}.alpha.upstateKC.KC=[tlock_tac{conds,TT,sleep}.alpha.upstateKC.KC tlock_tac{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_tac{conds,TT,sleep}.alpha.upstateKC.KC_all=[tlock_tac{conds,TT,sleep}.alpha.upstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_tac{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_tac{conds,TT,sleep}.alpha.upstateKC.delta=[tlock_tac{conds,TT,sleep}.alpha.upstateKC.delta tlock_tac{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_tac{conds,TT,sleep}.alpha.upstateKC.KC_all=[tlock_tac{conds,TT,sleep}.alpha.upstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_tac{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_tac{conds,TT,sleep}.alpha.upstateKC.sw=[tlock_tac{conds,TT,sleep}.alpha.upstateKC.sw tlock_tac{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_tac{conds,TT,sleep}.alpha.upstateKC.KC_all=[tlock_tac{conds,TT,sleep}.alpha.upstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                                end
                    % same for all the other states!      
                    elseif [tlock_tac{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)>0 && tlock_tac{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)<pi/4] || ... %% up-mid-state
                        [tlock_tac{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)>3*pi/4 && tlock_tac{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)<pi]
                        ctr(sleep,2,conds)=ctr(sleep,2,conds)+1;
                            if [tlock_tac{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_tac{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_tac{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_tac{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                            ctrKc(sleep,2,conds)=ctrKc(sleep,2,conds)+1;
                            end
                                %% now make a list of onsets to compare distributions - ignore, not used any more!
                                if [tlock_tac{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_tac{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_tac{conds,TT,sleep}.alpha.upmidstateKC.KC=[tlock_tac{conds,TT,sleep}.alpha.upmidstateKC.KC tlock_tac{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_tac{conds,TT,sleep}.alpha.upmidstateKC.KC_all=[tlock_tac{conds,TT,sleep}.alpha.upmidstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_tac{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_tac{conds,TT,sleep}.alpha.upmidstateKC.delta=[tlock_tac{conds,TT,sleep}.alpha.upmidstateKC.delta tlock_tac{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_tac{conds,TT,sleep}.alpha.upmidstateKC.KC_all=[tlock_tac{conds,TT,sleep}.alpha.upmidstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_tac{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_tac{conds,TT,sleep}.alpha.upmidstateKC.sw=[tlock_tac{conds,TT,sleep}.alpha.upmidstateKC.sw tlock_tac{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_tac{conds,TT,sleep}.alpha.upmidstateKC.KC_all=[tlock_tac{conds,TT,sleep}.alpha.upmidstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                                end
                            
                    elseif [tlock_tac{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)>-pi && tlock_tac{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)<-3*pi/4] || ... %% down-mid-state
                    [tlock_tac{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)>-pi/4 && tlock_tac{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)<0]    
                    ctr(sleep,3,conds)=ctr(sleep,3,conds)+1;
                            if [tlock_tac{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_tac{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_tac{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_tac{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                            ctrKc(sleep,3,conds)=ctrKc(sleep,3,conds)+1;
                            end
                                %% now make a list of onsets to compare distributions
                                if [tlock_tac{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_tac{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_tac{conds,TT,sleep}.alpha.downmidstateKC.KC=[tlock_tac{conds,TT,sleep}.alpha.downmidstateKC.KC tlock_tac{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_tac{conds,TT,sleep}.alpha.downmidstateKC.KC_all=[tlock_tac{conds,TT,sleep}.alpha.downmidstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_tac{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_tac{conds,TT,sleep}.alpha.downmidstateKC.delta=[tlock_tac{conds,TT,sleep}.alpha.downmidstateKC.delta tlock_tac{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_tac{conds,TT,sleep}.alpha.downmidstateKC.KC_all=[tlock_tac{conds,TT,sleep}.alpha.downmidstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_tac{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_tac{conds,TT,sleep}.alpha.downmidstateKC.sw=[tlock_tac{conds,TT,sleep}.alpha.downmidstateKC.sw tlock_tac{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_tac{conds,TT,sleep}.alpha.downmidstateKC.KC_all=[tlock_tac{conds,TT,sleep}.alpha.downmidstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                                end
                            
                    elseif [tlock_tac{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)>-3*pi/4 && tlock_tac{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)<-pi/4] %% down-state
                            ctr(sleep,4,conds)=ctr(sleep,4,conds)+1;
                            if [tlock_tac{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_tac{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_tac{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_tac{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                            ctrKc(sleep,4,conds)=ctrKc(sleep,4,conds)+1;
                            end
                                %% now make a list of onsets to compare distributions
                                if [tlock_tac{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_tac{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_tac{conds,TT,sleep}.alpha.downstateKC.KC=[tlock_tac{conds,TT,sleep}.alpha.downstateKC.KC tlock_tac{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_tac{conds,TT,sleep}.alpha.downstateKC.KC_all=[tlock_tac{conds,TT,sleep}.alpha.downstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_tac{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_tac{conds,TT,sleep}.alpha.downstateKC.delta=[tlock_tac{conds,TT,sleep}.alpha.downstateKC.delta tlock_tac{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_tac{conds,TT,sleep}.alpha.downstateKC.KC_all=[tlock_tac{conds,TT,sleep}.alpha.downstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_tac{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_tac{conds,TT,sleep}.alpha.downstateKC.sw=[tlock_tac{conds,TT,sleep}.alpha.downstateKC.sw tlock_tac{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_tac{conds,TT,sleep}.alpha.downstateKC.KC_all=[tlock_tac{conds,TT,sleep}.alpha.downstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                                end
                             
                    end
                    
            end
            tlock_tac{conds,TT,sleep}.phaseKC.alphaNtrials=ctr(sleep,:,conds);
            tlock_tac{conds,TT,sleep}.phaseKC.alphaNkcs=ctrKc(sleep,:,conds);
            tlock_tac{conds,TT,sleep}.phaseKC.alpha.KCrate=ctrKc(sleep,:,conds)./ctr(sleep,:,conds);
            
            tlock_tac_to_save{conds,TT,sleep}.phaseKC.alphaNtrials=ctr(sleep,:,conds);
            tlock_tac_to_save{conds,TT,sleep}.phaseKC.alphaNkcs=ctrKc(sleep,:,conds);
            tlock_tac_to_save{conds,TT,sleep}.phaseKC.alpha.KCrate=ctrKc(sleep,:,conds)./ctr(sleep,:,conds);
            
            KC{ii,:,:,:}.alpha.tac.ctr=ctr;
            KC{ii,:,:,:}.alpha.tac.ctrKC=ctrKc;
         end
         
        end
        
    end
%% then for tlock_aud
condTypes=length(tlock_aud);
ctr(1:23,1:4,1:condTypes)=0;
ctrKc(1:23,1:4,1:condTypes)=0;
%% for tlock_aud
for sleep=sleepStage(1:end)
     for conds=1:condTypes
         if isfield(tlock_aud{conds,TT,sleep},'hilbert')
            tlock_aud{conds,TT,sleep}.alpha.upstateKC.KC=[];
            tlock_aud{conds,TT,sleep}.alpha.upmidstateKC.KC=[];
            tlock_aud{conds,TT,sleep}.alpha.downmidstateKC.KC=[];
            tlock_aud{conds,TT,sleep}.alpha.downstateKC.KC=[];
            
            tlock_aud{conds,TT,sleep}.alpha.upstateKC.KC_all=[];
            tlock_aud{conds,TT,sleep}.alpha.upmidstateKC.KC_all=[];
            tlock_aud{conds,TT,sleep}.alpha.downmidstateKC.KC_all=[];
            tlock_aud{conds,TT,sleep}.alpha.downstateKC.KC_all=[];
            
            tlock_aud{conds,TT,sleep}.alpha.upstateKC.sw=[];
            tlock_aud{conds,TT,sleep}.alpha.upmidstateKC.sw=[];
            tlock_aud{conds,TT,sleep}.alpha.downmidstateKC.sw=[];
            tlock_aud{conds,TT,sleep}.alpha.downstateKC.sw=[];
            
            tlock_aud{conds,TT,sleep}.alpha.upstateKC.delta=[];
            tlock_aud{conds,TT,sleep}.alpha.upmidstateKC.delta=[];
            tlock_aud{conds,TT,sleep}.alpha.downmidstateKC.delta=[];
            tlock_aud{conds,TT,sleep}.alpha.downstateKC.delta=[];
            
             for i=1:size(tlock_aud{conds,TT,sleep}.hilbert.alphaT0.angle,1)   
                    if [tlock_aud{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)>pi/4 && tlock_aud{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)<3*pi/4] %% up-state
                        ctr(sleep,1,conds)=ctr(sleep,1,conds)+1;
                            if [tlock_aud{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_aud{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_aud{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_aud{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                                    ctrKc(sleep,1,conds)=ctrKc(sleep,1,conds)+1;
                            end
                                %% now make a list of onsets to compare distributions
                                if [tlock_aud{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_aud{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_aud{conds,TT,sleep}.alpha.upstateKC.KC=[tlock_aud{conds,TT,sleep}.alpha.upstateKC.KC tlock_aud{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_aud{conds,TT,sleep}.alpha.upstateKC.KC_all=[tlock_aud{conds,TT,sleep}.alpha.upstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_aud{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_aud{conds,TT,sleep}.alpha.upstateKC.delta=[tlock_aud{conds,TT,sleep}.alpha.upstateKC.delta tlock_aud{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_aud{conds,TT,sleep}.alpha.upstateKC.KC_all=[tlock_aud{conds,TT,sleep}.alpha.upstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_aud{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_aud{conds,TT,sleep}.alpha.upstateKC.sw=[tlock_aud{conds,TT,sleep}.alpha.upstateKC.sw tlock_aud{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_aud{conds,TT,sleep}.alpha.upstateKC.KC_all=[tlock_aud{conds,TT,sleep}.alpha.upstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                                end
                           
                    elseif  [tlock_aud{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)>0 && tlock_aud{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)<pi/4] || ... %% up-mid-state
                        [tlock_aud{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)>3*pi/4 && tlock_aud{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)<pi]
                        ctr(sleep,2,conds)=ctr(sleep,2,conds)+1;
                            if [tlock_aud{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_aud{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_aud{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_aud{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                            ctrKc(sleep,2,conds)=ctrKc(sleep,2,conds)+1;
                            end
                                %% now make a list of onsets to compare distributions
                                if [tlock_aud{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_aud{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_aud{conds,TT,sleep}.alpha.upmidstateKC.KC=[tlock_aud{conds,TT,sleep}.alpha.upmidstateKC.KC tlock_aud{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_aud{conds,TT,sleep}.alpha.upmidstateKC.KC_all=[tlock_aud{conds,TT,sleep}.alpha.upmidstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_aud{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_aud{conds,TT,sleep}.alpha.upmidstateKC.delta=[tlock_aud{conds,TT,sleep}.alpha.upmidstateKC.delta tlock_aud{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_aud{conds,TT,sleep}.alpha.upmidstateKC.KC_all=[tlock_aud{conds,TT,sleep}.alpha.upmidstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_aud{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_aud{conds,TT,sleep}.alpha.upmidstateKC.sw=[tlock_aud{conds,TT,sleep}.alpha.upmidstateKC.sw tlock_aud{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_aud{conds,TT,sleep}.alpha.upmidstateKC.KC_all=[tlock_aud{conds,TT,sleep}.alpha.upmidstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                                end
                           
                    elseif [tlock_aud{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)>-pi && tlock_aud{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)<-3*pi/4] || ... %% down-mid-state
                    [tlock_aud{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)>-pi/4 && tlock_aud{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)<0]
                        ctr(sleep,3,conds)=ctr(sleep,3,conds)+1;
                            if [tlock_aud{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_aud{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_aud{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_aud{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                            ctrKc(sleep,3,conds)=ctrKc(sleep,3,conds)+1;
                            end
                                %% now make a list of onsets to compare distributions
                                if [tlock_aud{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_aud{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_aud{conds,TT,sleep}.alpha.downmidstateKC.KC=[tlock_aud{conds,TT,sleep}.alpha.downmidstateKC.KC tlock_aud{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_aud{conds,TT,sleep}.alpha.downmidstateKC.KC_all=[tlock_aud{conds,TT,sleep}.alpha.downmidstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_aud{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_aud{conds,TT,sleep}.alpha.downmidstateKC.delta=[tlock_aud{conds,TT,sleep}.alpha.downmidstateKC.delta tlock_aud{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_aud{conds,TT,sleep}.alpha.downmidstateKC.KC_all=[tlock_aud{conds,TT,sleep}.alpha.downmidstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_aud{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_aud{conds,TT,sleep}.alpha.downmidstateKC.sw=[tlock_aud{conds,TT,sleep}.alpha.downmidstateKC.sw tlock_aud{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_aud{conds,TT,sleep}.alpha.downmidstateKC.KC_all=[tlock_aud{conds,TT,sleep}.alpha.downmidstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                                end
                            
                    elseif [tlock_aud{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)>-3*pi/4 && tlock_aud{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)<-pi/4] %% down-state
                            ctr(sleep,4,conds)=ctr(sleep,4,conds)+1;
                            if [tlock_aud{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_aud{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_aud{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_aud{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                            ctrKc(sleep,4,conds)=ctrKc(sleep,4,conds)+1;
                            end
                                %% now make a list of onsets to compare distributions
                                if [tlock_aud{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_aud{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_aud{conds,TT,sleep}.alpha.downstateKC.KC=[tlock_aud{conds,TT,sleep}.alpha.downstateKC.KC tlock_aud{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_aud{conds,TT,sleep}.alpha.downstateKC.KC_all=[tlock_aud{conds,TT,sleep}.alpha.downstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_aud{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_aud{conds,TT,sleep}.alpha.downstateKC.delta=[tlock_aud{conds,TT,sleep}.alpha.downstateKC.delta tlock_aud{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_aud{conds,TT,sleep}.alpha.downstateKC.KC_all=[tlock_aud{conds,TT,sleep}.alpha.downstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_aud{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_aud{conds,TT,sleep}.alpha.downstateKC.sw=[tlock_aud{conds,TT,sleep}.alpha.downstateKC.sw tlock_aud{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_aud{conds,TT,sleep}.alpha.downstateKC.KC_all=[tlock_aud{conds,TT,sleep}.alpha.downstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                                end
                             
                    end
                    
            end
            tlock_aud{conds,TT,sleep}.phaseKC.alphaNtrials=ctr(sleep,:,conds);
            tlock_aud{conds,TT,sleep}.phaseKC.alphaNkcs=ctrKc(sleep,:,conds);
            tlock_aud{conds,TT,sleep}.phaseKC.alpha.KCrate=ctrKc(sleep,:,conds)./ctr(sleep,:,conds);
            
            tlock_aud_to_save{conds,TT,sleep}.phaseKC.alphaNtrials=ctr(sleep,:,conds);
            tlock_aud_to_save{conds,TT,sleep}.phaseKC.alphaNkcs=ctrKc(sleep,:,conds);
            tlock_aud_to_save{conds,TT,sleep}.phaseKC.alpha.KCrate=ctrKc(sleep,:,conds)./ctr(sleep,:,conds);
            
%             KC{ii,:,:,:}.alpha.aud.ctr=ctr;
%             KC{ii,:,:,:}.alpha.aud.ctrKC=ctrKc;
         end
         
        end
        
    end

%% then for tlock_nul
condTypes=length(tlock_nul);
ctr(1:23,1:4,1:condTypes)=0;
ctrKc(1:23,1:4,1:condTypes)=0;
%% for tlock_nul
for sleep=sleepStage(1:end)
     for conds=1:condTypes
         if isfield(tlock_nul{conds,TT,sleep},'hilbert')
            tlock_nul{conds,TT,sleep}.alpha.upstateKC.KC=[];
            tlock_nul{conds,TT,sleep}.alpha.upmidstateKC.KC=[];
            tlock_nul{conds,TT,sleep}.alpha.downmidstateKC.KC=[];
            tlock_nul{conds,TT,sleep}.alpha.downstateKC.KC=[];
            
            tlock_nul{conds,TT,sleep}.alpha.upstateKC.KC_all=[];
            tlock_nul{conds,TT,sleep}.alpha.upmidstateKC.KC_all=[];
            tlock_nul{conds,TT,sleep}.alpha.downmidstateKC.KC_all=[];
            tlock_nul{conds,TT,sleep}.alpha.downstateKC.KC_all=[];
            
            tlock_nul{conds,TT,sleep}.alpha.upstateKC.sw=[];
            tlock_nul{conds,TT,sleep}.alpha.upmidstateKC.sw=[];
            tlock_nul{conds,TT,sleep}.alpha.downmidstateKC.sw=[];
            tlock_nul{conds,TT,sleep}.alpha.downstateKC.sw=[];
            
            tlock_nul{conds,TT,sleep}.alpha.upstateKC.delta=[];
            tlock_nul{conds,TT,sleep}.alpha.upmidstateKC.delta=[];
            tlock_nul{conds,TT,sleep}.alpha.downmidstateKC.delta=[];
            tlock_nul{conds,TT,sleep}.alpha.downstateKC.delta=[];
            for i=1:size(tlock_nul{conds,TT,sleep}.hilbert.alphaT0.angle,1)
                
                    if [tlock_nul{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)>pi/4 && tlock_nul{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)<3*pi/4] %% up-state
                        ctr(sleep,1,conds)=ctr(sleep,1,conds)+1;
                            if [tlock_nul{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_nul{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_nul{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_nul{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                                    ctrKc(sleep,1,conds)=ctrKc(sleep,1,conds)+1;
                            end
                                %% now make a list of onsets to compare distributions
                                if [tlock_nul{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_nul{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_nul{conds,TT,sleep}.alpha.upstateKC.KC=[tlock_nul{conds,TT,sleep}.alpha.upstateKC.KC tlock_nul{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_nul{conds,TT,sleep}.alpha.upstateKC.KC_all=[tlock_nul{conds,TT,sleep}.alpha.upstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_nul{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_nul{conds,TT,sleep}.alpha.upstateKC.delta=[tlock_nul{conds,TT,sleep}.alpha.upstateKC.delta tlock_nul{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_nul{conds,TT,sleep}.alpha.upstateKC.KC_all=[tlock_nul{conds,TT,sleep}.alpha.upstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_nul{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_nul{conds,TT,sleep}.alpha.upstateKC.sw=[tlock_nul{conds,TT,sleep}.alpha.upstateKC.sw tlock_nul{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_nul{conds,TT,sleep}.alpha.upstateKC.KC_all=[tlock_nul{conds,TT,sleep}.alpha.upstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                                end
                            
                    elseif  [tlock_nul{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)>0 && tlock_nul{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)<pi/4] || ... %% up-mid-state
                        [tlock_nul{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)>3*pi/4 && tlock_nul{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)<pi]
                        ctr(sleep,2,conds)=ctr(sleep,2,conds)+1;
                            if [tlock_nul{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_nul{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_nul{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_nul{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                            ctrKc(sleep,2,conds)=ctrKc(sleep,2,conds)+1;
                            end
                                %% now make a list of onsets to compare distributions
                                if [tlock_nul{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_nul{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_nul{conds,TT,sleep}.alpha.upmidstateKC.KC=[tlock_nul{conds,TT,sleep}.alpha.upmidstateKC.KC tlock_nul{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_nul{conds,TT,sleep}.alpha.upmidstateKC.KC_all=[tlock_nul{conds,TT,sleep}.alpha.upmidstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_nul{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_nul{conds,TT,sleep}.alpha.upmidstateKC.delta=[tlock_nul{conds,TT,sleep}.alpha.upmidstateKC.delta tlock_nul{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_nul{conds,TT,sleep}.alpha.upmidstateKC.KC_all=[tlock_nul{conds,TT,sleep}.alpha.upmidstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_nul{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_nul{conds,TT,sleep}.alpha.upmidstateKC.sw=[tlock_nul{conds,TT,sleep}.alpha.upmidstateKC.sw tlock_nul{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_nul{conds,TT,sleep}.alpha.upmidstateKC.KC_all=[tlock_nul{conds,TT,sleep}.alpha.upmidstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                                end
                            
                     elseif [tlock_nul{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)>-pi && tlock_nul{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)<-3*pi/4] || ... %% down-mid-state
                    [tlock_nul{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)>-pi/4 && tlock_nul{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)<0]
                        ctr(sleep,3,conds)=ctr(sleep,3,conds)+1;
                            if [tlock_nul{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_nul{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_nul{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_nul{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                            ctrKc(sleep,3,conds)=ctrKc(sleep,3,conds)+1;
                            end
                                %% now make a list of onsets to compare distributions
                                if [tlock_nul{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_nul{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_nul{conds,TT,sleep}.alpha.downmidstateKC.KC=[tlock_nul{conds,TT,sleep}.alpha.downmidstateKC.KC tlock_nul{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_nul{conds,TT,sleep}.alpha.downmidstateKC.KC_all=[tlock_nul{conds,TT,sleep}.alpha.downmidstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_nul{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_nul{conds,TT,sleep}.alpha.downmidstateKC.delta=[tlock_nul{conds,TT,sleep}.alpha.downmidstateKC.delta tlock_nul{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_nul{conds,TT,sleep}.alpha.downmidstateKC.KC_all=[tlock_nul{conds,TT,sleep}.alpha.downmidstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_nul{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_nul{conds,TT,sleep}.alpha.downmidstateKC.sw=[tlock_nul{conds,TT,sleep}.alpha.downmidstateKC.sw tlock_nul{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_nul{conds,TT,sleep}.alpha.downmidstateKC.KC_all=[tlock_nul{conds,TT,sleep}.alpha.downmidstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                                end
                            
                    elseif [tlock_nul{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)>-3*pi/4 && tlock_nul{conds,TT,sleep}.hilbert.alphaT0.angle(i,1)<-pi/4] %% down-state
                            ctr(sleep,4,conds)=ctr(sleep,4,conds)+1;
                            if [tlock_nul{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_nul{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_nul{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_nul{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                            ctrKc(sleep,4,conds)=ctrKc(sleep,4,conds)+1;
                            end
                                %% now make a list of onsets to compare distributions
                                if [tlock_nul{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_nul{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_nul{conds,TT,sleep}.alpha.downstateKC.KC=[tlock_nul{conds,TT,sleep}.alpha.downstateKC.KC tlock_nul{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_nul{conds,TT,sleep}.alpha.downstateKC.KC_all=[tlock_nul{conds,TT,sleep}.alpha.downstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_nul{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_nul{conds,TT,sleep}.alpha.downstateKC.delta=[tlock_nul{conds,TT,sleep}.alpha.downstateKC.delta tlock_nul{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_nul{conds,TT,sleep}.alpha.downstateKC.KC_all=[tlock_nul{conds,TT,sleep}.alpha.downstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_nul{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_nul{conds,TT,sleep}.alpha.downstateKC.sw=[tlock_nul{conds,TT,sleep}.alpha.downstateKC.sw tlock_nul{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_nul{conds,TT,sleep}.alpha.downstateKC.KC_all=[tlock_nul{conds,TT,sleep}.alpha.downstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                                end
                            
                    end
                    
            end
            tlock_nul{conds,TT,sleep}.phaseKC.alphaNtrials=ctr(sleep,:,conds);
            tlock_nul{conds,TT,sleep}.phaseKC.alphaNkcs=ctrKc(sleep,:,conds);
            tlock_nul{conds,TT,sleep}.phaseKC.alpha.KCrate=ctrKc(sleep,:,conds)./ctr(sleep,:,conds);
            
            tlock_nul_to_save{conds,TT,sleep}.phaseKC.alphaNtrials=ctr(sleep,:,conds);
            tlock_nul_to_save{conds,TT,sleep}.phaseKC.alphaNkcs=ctrKc(sleep,:,conds);
            tlock_nul_to_save{conds,TT,sleep}.phaseKC.alpha.KCrate=ctrKc(sleep,:,conds)./ctr(sleep,:,conds);
            
            KC{ii,:,:,:}.alpha.nul.ctr=ctr;
            KC{ii,:,:,:}.alpha.nul.ctrKC=ctrKc;
         end
         
        end
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AND THE SAME FOR THETA %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ang_bins=4;
sleepStage=[10 11 12 13]; 
condTypes=length(tlock_tac);
ctr(1:23,1:4,1:condTypes)=0;
ctrKc(1:23,1:4,1:condTypes)=0;
%% for tlock_tac
for sleep=sleepStage(1:end)
     for conds=1:condTypes
         if isfield(tlock_tac{conds,TT,sleep},'hilbert')
            tlock_tac{conds,TT,sleep}.theta.upstateKC.KC=[];
            tlock_tac{conds,TT,sleep}.theta.upmidstateKC.KC=[];
            tlock_tac{conds,TT,sleep}.theta.downmidstateKC.KC=[];
            tlock_tac{conds,TT,sleep}.theta.downstateKC.KC=[];
            
            tlock_tac{conds,TT,sleep}.theta.upstateKC.KC_all=[];
            tlock_tac{conds,TT,sleep}.theta.upmidstateKC.KC_all=[];
            tlock_tac{conds,TT,sleep}.theta.downmidstateKC.KC_all=[];
            tlock_tac{conds,TT,sleep}.theta.downstateKC.KC_all=[];
            
            tlock_tac{conds,TT,sleep}.theta.upstateKC.sw=[];
            tlock_tac{conds,TT,sleep}.theta.upmidstateKC.sw=[];
            tlock_tac{conds,TT,sleep}.theta.downmidstateKC.sw=[];
            tlock_tac{conds,TT,sleep}.theta.downstateKC.sw=[];
            
            tlock_tac{conds,TT,sleep}.theta.upstateKC.delta=[];
            tlock_tac{conds,TT,sleep}.theta.upmidstateKC.delta=[];
            tlock_tac{conds,TT,sleep}.theta.downmidstateKC.delta=[];
            tlock_tac{conds,TT,sleep}.theta.downstateKC.delta=[];
            for i=1:size(tlock_tac{conds,TT,sleep}.hilbert.thetaT0.angle,1)
                
                    if [tlock_tac{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)>pi/4 && tlock_tac{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)<3*pi/4] %% up-state
                        ctr(sleep,1,conds)=ctr(sleep,1,conds)+1;
                            if [tlock_tac{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_tac{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_tac{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_tac{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                                    ctrKc(sleep,1,conds)=ctrKc(sleep,1,conds)+1;
                            end
                                %% now make a list of onsets to compare distributions
                                if [tlock_tac{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_tac{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_tac{conds,TT,sleep}.theta.upstateKC.KC=[tlock_tac{conds,TT,sleep}.theta.upstateKC.KC tlock_tac{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_tac{conds,TT,sleep}.theta.upstateKC.KC_all=[tlock_tac{conds,TT,sleep}.theta.upstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_tac{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_tac{conds,TT,sleep}.theta.upstateKC.delta=[tlock_tac{conds,TT,sleep}.theta.upstateKC.delta tlock_tac{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_tac{conds,TT,sleep}.theta.upstateKC.KC_all=[tlock_tac{conds,TT,sleep}.theta.upstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_tac{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_tac{conds,TT,sleep}.theta.upstateKC.sw=[tlock_tac{conds,TT,sleep}.theta.upstateKC.sw tlock_tac{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_tac{conds,TT,sleep}.theta.upstateKC.KC_all=[tlock_tac{conds,TT,sleep}.theta.upstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                               
                            end
                    elseif [tlock_tac{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)>0 && tlock_tac{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)<pi/4] || ... %% up-mid-state
                        [tlock_tac{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)>3*pi/4 && tlock_tac{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)<pi]
                        ctr(sleep,2,conds)=ctr(sleep,2,conds)+1;
                            if [tlock_tac{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_tac{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_tac{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_tac{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                            ctrKc(sleep,2,conds)=ctrKc(sleep,2,conds)+1;
                            end
                                %% now make a list of onsets to compare distributions
                                if [tlock_tac{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_tac{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_tac{conds,TT,sleep}.theta.upmidstateKC.KC=[tlock_tac{conds,TT,sleep}.theta.upmidstateKC.KC tlock_tac{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_tac{conds,TT,sleep}.theta.upmidstateKC.KC_all=[tlock_tac{conds,TT,sleep}.theta.upmidstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_tac{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_tac{conds,TT,sleep}.theta.upmidstateKC.delta=[tlock_tac{conds,TT,sleep}.theta.upmidstateKC.delta tlock_tac{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_tac{conds,TT,sleep}.theta.upmidstateKC.KC_all=[tlock_tac{conds,TT,sleep}.theta.upmidstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_tac{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_tac{conds,TT,sleep}.theta.upmidstateKC.sw=[tlock_tac{conds,TT,sleep}.theta.upmidstateKC.sw tlock_tac{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_tac{conds,TT,sleep}.theta.upmidstateKC.KC_all=[tlock_tac{conds,TT,sleep}.theta.upmidstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                                end
                            
                    elseif [tlock_tac{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)>-pi && tlock_tac{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)<-3*pi/4] || ... %% down-mid-state
                    [tlock_tac{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)>-pi/4 && tlock_tac{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)<0]
                        [tlock_tac{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)>3*pi/4 && tlock_tac{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)<pi]
                        ctr(sleep,3,conds)=ctr(sleep,3,conds)+1;
                            if [tlock_tac{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_tac{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_tac{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_tac{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                            ctrKc(sleep,3,conds)=ctrKc(sleep,3,conds)+1;
                            end
                                %% now make a list of onsets to compare distributions
                                if [tlock_tac{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_tac{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_tac{conds,TT,sleep}.theta.downmidstateKC.KC=[tlock_tac{conds,TT,sleep}.theta.downmidstateKC.KC tlock_tac{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_tac{conds,TT,sleep}.theta.downmidstateKC.KC_all=[tlock_tac{conds,TT,sleep}.theta.downmidstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_tac{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_tac{conds,TT,sleep}.theta.downmidstateKC.delta=[tlock_tac{conds,TT,sleep}.theta.downmidstateKC.delta tlock_tac{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_tac{conds,TT,sleep}.theta.downmidstateKC.KC_all=[tlock_tac{conds,TT,sleep}.theta.downmidstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_tac{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_tac{conds,TT,sleep}.theta.downmidstateKC.sw=[tlock_tac{conds,TT,sleep}.theta.downmidstateKC.sw tlock_tac{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_tac{conds,TT,sleep}.theta.downmidstateKC.KC_all=[tlock_tac{conds,TT,sleep}.theta.downmidstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                                end
                            
                    elseif [tlock_tac{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)>-3*pi/4 && tlock_tac{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)<-pi/4] %% down-state
                            ctr(sleep,4,conds)=ctr(sleep,4,conds)+1;
                            if [tlock_tac{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_tac{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_tac{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_tac{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_tac{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                            ctrKc(sleep,4,conds)=ctrKc(sleep,4,conds)+1;
                            end
                                %% now make a list of onsets to compare distributions
                                if [tlock_tac{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_tac{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_tac{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_tac{conds,TT,sleep}.theta.downstateKC.KC=[tlock_tac{conds,TT,sleep}.theta.downstateKC.KC tlock_tac{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_tac{conds,TT,sleep}.theta.downstateKC.KC_all=[tlock_tac{conds,TT,sleep}.theta.downstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_tac{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_tac{conds,TT,sleep}.theta.downstateKC.delta=[tlock_tac{conds,TT,sleep}.theta.downstateKC.delta tlock_tac{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_tac{conds,TT,sleep}.theta.downstateKC.KC_all=[tlock_tac{conds,TT,sleep}.theta.downstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_tac{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_tac{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_tac{conds,TT,sleep}.theta.downstateKC.sw=[tlock_tac{conds,TT,sleep}.theta.downstateKC.sw tlock_tac{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_tac{conds,TT,sleep}.theta.downstateKC.KC_all=[tlock_tac{conds,TT,sleep}.theta.downstateKC.KC_all tlock_tac{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                                end
                          
                    end
                    
            end
            tlock_tac{conds,TT,sleep}.phaseKC.thetaNtrials=ctr(sleep,:,conds);
            tlock_tac{conds,TT,sleep}.phaseKC.thetaNkcs=ctrKc(sleep,:,conds);
            tlock_tac{conds,TT,sleep}.phaseKC.theta.KCrate=ctrKc(sleep,:,conds)./ctr(sleep,:,conds);
            
            tlock_tac_to_save{conds,TT,sleep}.phaseKC.thetaNtrials=ctr(sleep,:,conds);
            tlock_tac_to_save{conds,TT,sleep}.phaseKC.thetaNkcs=ctrKc(sleep,:,conds);
            tlock_tac_to_save{conds,TT,sleep}.phaseKC.theta.KCrate=ctrKc(sleep,:,conds)./ctr(sleep,:,conds);
            
         end
         
        end
        
end
%% then for tlock_aud
condTypes=length(tlock_aud);
ctr(1:23,1:4,1:condTypes)=0;
ctrKc(1:23,1:4,1:condTypes)=0;
for sleep=sleepStage(1:end)
     for conds=1:condTypes
         if isfield(tlock_aud{conds,TT,sleep},'hilbert')
            tlock_aud{conds,TT,sleep}.theta.upstateKC.KC=[];
            tlock_aud{conds,TT,sleep}.theta.upmidstateKC.KC=[];
            tlock_aud{conds,TT,sleep}.theta.downmidstateKC.KC=[];
            tlock_aud{conds,TT,sleep}.theta.downstateKC.KC=[];
            
            tlock_aud{conds,TT,sleep}.theta.upstateKC.KC_all=[];
            tlock_aud{conds,TT,sleep}.theta.upmidstateKC.KC_all=[];
            tlock_aud{conds,TT,sleep}.theta.downmidstateKC.KC_all=[];
            tlock_aud{conds,TT,sleep}.theta.downstateKC.KC_all=[];
            
            tlock_aud{conds,TT,sleep}.theta.upstateKC.sw=[];
            tlock_aud{conds,TT,sleep}.theta.upmidstateKC.sw=[];
            tlock_aud{conds,TT,sleep}.theta.downmidstateKC.sw=[];
            tlock_aud{conds,TT,sleep}.theta.downstateKC.sw=[];
            
            tlock_aud{conds,TT,sleep}.theta.upstateKC.delta=[];
            tlock_aud{conds,TT,sleep}.theta.upmidstateKC.delta=[];
            tlock_aud{conds,TT,sleep}.theta.downmidstateKC.delta=[];
            tlock_aud{conds,TT,sleep}.theta.downstateKC.delta=[];
            for i=1:size(tlock_aud{conds,TT,sleep}.hilbert.thetaT0.angle,1)
                
                    if [tlock_aud{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)>pi/4 && tlock_aud{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)<3*pi/4] %% up-state
                        ctr(sleep,1,conds)=ctr(sleep,1,conds)+1;
                            if [tlock_aud{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_aud{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_aud{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_aud{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                                    ctrKc(sleep,1,conds)=ctrKc(sleep,1,conds)+1;
                            end
                                %% now make a list of onsets to compare distributions
                                if [tlock_aud{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_aud{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_aud{conds,TT,sleep}.theta.upstateKC.KC=[tlock_aud{conds,TT,sleep}.theta.upstateKC.KC tlock_aud{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_aud{conds,TT,sleep}.theta.upstateKC.KC_all=[tlock_aud{conds,TT,sleep}.theta.upstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_aud{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_aud{conds,TT,sleep}.theta.upstateKC.delta=[tlock_aud{conds,TT,sleep}.theta.upstateKC.delta tlock_aud{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_aud{conds,TT,sleep}.theta.upstateKC.KC_all=[tlock_aud{conds,TT,sleep}.theta.upstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_aud{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_aud{conds,TT,sleep}.theta.upstateKC.sw=[tlock_aud{conds,TT,sleep}.theta.upstateKC.sw tlock_aud{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_aud{conds,TT,sleep}.theta.upstateKC.KC_all=[tlock_aud{conds,TT,sleep}.theta.upstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                                end
                            
                    elseif [tlock_aud{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)>0 && tlock_aud{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)<pi/4] || ... %% up-mid-state
                        [tlock_aud{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)>3*pi/4 && tlock_aud{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)<pi]
                        ctr(sleep,2,conds)=ctr(sleep,2,conds)+1;
                            if [tlock_aud{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_aud{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_aud{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_aud{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                            ctrKc(sleep,2,conds)=ctrKc(sleep,2,conds)+1;
                            end
                                %% now make a list of onsets to compare distributions
                                if [tlock_aud{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_aud{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_aud{conds,TT,sleep}.theta.upmidstateKC.KC=[tlock_aud{conds,TT,sleep}.theta.upmidstateKC.KC tlock_aud{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_aud{conds,TT,sleep}.theta.upmidstateKC.KC_all=[tlock_aud{conds,TT,sleep}.theta.upmidstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_aud{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_aud{conds,TT,sleep}.theta.upmidstateKC.delta=[tlock_aud{conds,TT,sleep}.theta.upmidstateKC.delta tlock_aud{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_aud{conds,TT,sleep}.theta.upmidstateKC.KC_all=[tlock_aud{conds,TT,sleep}.theta.upmidstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_aud{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_aud{conds,TT,sleep}.theta.upmidstateKC.sw=[tlock_aud{conds,TT,sleep}.theta.upmidstateKC.sw tlock_aud{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_aud{conds,TT,sleep}.theta.upmidstateKC.KC_all=[tlock_aud{conds,TT,sleep}.theta.upmidstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                                end
                            
                    elseif [tlock_aud{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)>-pi && tlock_aud{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)<-3*pi/4] || ... %% down-mid-state
                    [tlock_aud{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)>-pi/4 && tlock_aud{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)<0]
                        ctr(sleep,3,conds)=ctr(sleep,3,conds)+1;
                            if [tlock_aud{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_aud{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_aud{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_aud{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                            ctrKc(sleep,3,conds)=ctrKc(sleep,3,conds)+1;
                            end
                                %% now make a list of onsets to compare distributions
                                if [tlock_aud{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_aud{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_aud{conds,TT,sleep}.theta.downmidstateKC.KC=[tlock_aud{conds,TT,sleep}.theta.downmidstateKC.KC tlock_aud{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_aud{conds,TT,sleep}.theta.downmidstateKC.KC_all=[tlock_aud{conds,TT,sleep}.theta.downmidstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_aud{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_aud{conds,TT,sleep}.theta.downmidstateKC.delta=[tlock_aud{conds,TT,sleep}.theta.downmidstateKC.delta tlock_aud{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_aud{conds,TT,sleep}.theta.downmidstateKC.KC_all=[tlock_aud{conds,TT,sleep}.theta.downmidstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_aud{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_aud{conds,TT,sleep}.theta.downmidstateKC.sw=[tlock_aud{conds,TT,sleep}.theta.downmidstateKC.sw tlock_aud{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_aud{conds,TT,sleep}.theta.downmidstateKC.KC_all=[tlock_aud{conds,TT,sleep}.theta.downmidstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                                end
                            
                    elseif [tlock_aud{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)>-3*pi/4 && tlock_aud{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)<-pi/4] %% down-state
                            ctr(sleep,4,conds)=ctr(sleep,4,conds)+1;
                            if [tlock_aud{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_aud{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_aud{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_aud{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_aud{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                            ctrKc(sleep,4,conds)=ctrKc(sleep,4,conds)+1;
                            end
                                %% now make a list of onsets to compare distributions
                                if [tlock_aud{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_aud{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_aud{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_aud{conds,TT,sleep}.theta.downstateKC.KC=[tlock_aud{conds,TT,sleep}.theta.downstateKC.KC tlock_aud{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_aud{conds,TT,sleep}.theta.downstateKC.KC_all=[tlock_aud{conds,TT,sleep}.theta.downstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_aud{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_aud{conds,TT,sleep}.theta.downstateKC.delta=[tlock_aud{conds,TT,sleep}.theta.downstateKC.delta tlock_aud{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_aud{conds,TT,sleep}.theta.downstateKC.KC_all=[tlock_aud{conds,TT,sleep}.theta.downstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_aud{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_aud{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_aud{conds,TT,sleep}.theta.downstateKC.sw=[tlock_aud{conds,TT,sleep}.theta.downstateKC.sw tlock_aud{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_aud{conds,TT,sleep}.theta.downstateKC.KC_all=[tlock_aud{conds,TT,sleep}.theta.downstateKC.KC_all tlock_aud{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                                end
                             
                    end
                    
            end
            tlock_aud{conds,TT,sleep}.phaseKC.thetaNtrials=ctr(sleep,:,conds);
            tlock_aud{conds,TT,sleep}.phaseKC.thetaNkcs=ctrKc(sleep,:,conds);
            tlock_aud{conds,TT,sleep}.phaseKC.theta.KCrate=ctrKc(sleep,:,conds)./ctr(sleep,:,conds);
            
            tlock_aud_to_save{conds,TT,sleep}.phaseKC.thetaNtrials=ctr(sleep,:,conds);
            tlock_aud_to_save{conds,TT,sleep}.phaseKC.thetaNkcs=ctrKc(sleep,:,conds);
            tlock_aud_to_save{conds,TT,sleep}.phaseKC.theta.KCrate=ctrKc(sleep,:,conds)./ctr(sleep,:,conds);
         end
         
        end
        
    end

%% then for tlock_nul
%% then for tlock_nul
condTypes=length(tlock_nul);
ctr(1:23,1:4,1:condTypes)=0;
ctrKc(1:23,1:4,1:condTypes)=0;
for sleep=sleepStage(1:end)
     for conds=1:condTypes
         if isfield(tlock_nul{conds,TT,sleep},'hilbert')
            tlock_nul{conds,TT,sleep}.theta.upstateKC.KC=[];
            tlock_nul{conds,TT,sleep}.theta.upmidstateKC.KC=[];
            tlock_nul{conds,TT,sleep}.theta.downmidstateKC.KC=[];
            tlock_nul{conds,TT,sleep}.theta.downstateKC.KC=[];
            
            tlock_nul{conds,TT,sleep}.theta.upstateKC.KC_all=[];
            tlock_nul{conds,TT,sleep}.theta.upmidstateKC.KC_all=[];
            tlock_nul{conds,TT,sleep}.theta.downmidstateKC.KC_all=[];
            tlock_nul{conds,TT,sleep}.theta.downstateKC.KC_all=[];
            
            tlock_nul{conds,TT,sleep}.theta.upstateKC.sw=[];
            tlock_nul{conds,TT,sleep}.theta.upmidstateKC.sw=[];
            tlock_nul{conds,TT,sleep}.theta.downmidstateKC.sw=[];
            tlock_nul{conds,TT,sleep}.theta.downstateKC.sw=[];
            
            tlock_nul{conds,TT,sleep}.theta.upstateKC.delta=[];
            tlock_nul{conds,TT,sleep}.theta.upmidstateKC.delta=[];
            tlock_nul{conds,TT,sleep}.theta.downmidstateKC.delta=[];
            tlock_nul{conds,TT,sleep}.theta.downstateKC.delta=[];
            for i=1:size(tlock_nul{conds,TT,sleep}.hilbert.thetaT0.angle,1)
                
                    if [tlock_nul{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)>pi/4 && tlock_nul{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)<3*pi/4] %% up-state
                        ctr(sleep,1,conds)=ctr(sleep,1,conds)+1;
                            if [tlock_nul{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_nul{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_nul{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_nul{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                                    ctrKc(sleep,1,conds)=ctrKc(sleep,1,conds)+1;
                            end
                                %% now make a list of onsets to compare distributions
                                if [tlock_nul{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_nul{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_nul{conds,TT,sleep}.theta.upstateKC.KC=[tlock_nul{conds,TT,sleep}.theta.upstateKC.KC tlock_nul{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_nul{conds,TT,sleep}.theta.upstateKC.KC_all=[tlock_nul{conds,TT,sleep}.theta.upstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_nul{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_nul{conds,TT,sleep}.theta.upstateKC.delta=[tlock_nul{conds,TT,sleep}.theta.upstateKC.delta tlock_nul{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_nul{conds,TT,sleep}.theta.upstateKC.KC_all=[tlock_nul{conds,TT,sleep}.theta.upstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_nul{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_nul{conds,TT,sleep}.theta.upstateKC.sw=[tlock_nul{conds,TT,sleep}.theta.upstateKC.sw tlock_nul{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_nul{conds,TT,sleep}.theta.upstateKC.KC_all=[tlock_nul{conds,TT,sleep}.theta.upstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                                end
                           
                    elseif [tlock_nul{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)>0 && tlock_nul{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)<pi/4] || ... %% up-mid-state
                        [tlock_nul{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)>3*pi/4 && tlock_nul{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)<pi]
                        ctr(sleep,2,conds)=ctr(sleep,2,conds)+1;
                            if [tlock_nul{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_nul{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_nul{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_nul{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                            ctrKc(sleep,2,conds)=ctrKc(sleep,2,conds)+1;
                            end
                                %% now make a list of onsets to compare distributions
                                if [tlock_nul{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_nul{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_nul{conds,TT,sleep}.theta.upmidstateKC.KC=[tlock_nul{conds,TT,sleep}.theta.upmidstateKC.KC tlock_nul{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_nul{conds,TT,sleep}.theta.upmidstateKC.KC_all=[tlock_nul{conds,TT,sleep}.theta.upmidstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_nul{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_nul{conds,TT,sleep}.theta.upmidstateKC.delta=[tlock_nul{conds,TT,sleep}.theta.upmidstateKC.delta tlock_nul{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_nul{conds,TT,sleep}.theta.upmidstateKC.KC_all=[tlock_nul{conds,TT,sleep}.theta.upmidstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_nul{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_nul{conds,TT,sleep}.theta.upmidstateKC.sw=[tlock_nul{conds,TT,sleep}.theta.upmidstateKC.sw tlock_nul{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_nul{conds,TT,sleep}.theta.upmidstateKC.KC_all=[tlock_nul{conds,TT,sleep}.theta.upmidstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                                end
                           
                    elseif [tlock_nul{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)>-pi && tlock_nul{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)<-3*pi/4] || ... %% down-mid-state
                    [tlock_nul{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)>-pi/4 && tlock_nul{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)<0]
                        ctr(sleep,3,conds)=ctr(sleep,3,conds)+1;
                            if [tlock_nul{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_nul{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_nul{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_nul{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                            ctrKc(sleep,3,conds)=ctrKc(sleep,3,conds)+1;
                            end
                                %% now make a list of onsets to compare distributions
                                if [tlock_nul{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_nul{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_nul{conds,TT,sleep}.theta.downmidstateKC.KC=[tlock_nul{conds,TT,sleep}.theta.downmidstateKC.KC tlock_nul{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_nul{conds,TT,sleep}.theta.downmidstateKC.KC_all=[tlock_nul{conds,TT,sleep}.theta.downmidstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_nul{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_nul{conds,TT,sleep}.theta.downmidstateKC.delta=[tlock_nul{conds,TT,sleep}.theta.downmidstateKC.delta tlock_nul{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_nul{conds,TT,sleep}.theta.downmidstateKC.KC_all=[tlock_nul{conds,TT,sleep}.theta.downmidstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_nul{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_nul{conds,TT,sleep}.theta.downmidstateKC.sw=[tlock_nul{conds,TT,sleep}.theta.downmidstateKC.sw tlock_nul{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_nul{conds,TT,sleep}.theta.downmidstateKC.KC_all=[tlock_nul{conds,TT,sleep}.theta.downmidstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                                end
                            
                    elseif [tlock_nul{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)>-3*pi/4 && tlock_nul{conds,TT,sleep}.hilbert.thetaT0.angle(i,1)<-pi/4] %% down-state
                            ctr(sleep,4,conds)=ctr(sleep,4,conds)+1;
                            if [tlock_nul{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,150)>0] ...
                                    && [[tlock_nul{conds,TT,sleep}.trialinfo(i,80)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,80)<.8] || ...
                                    [tlock_nul{conds,TT,sleep}.trialinfo(i,120)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,120)<.8] || ...
                                    [tlock_nul{conds,TT,sleep}.trialinfo(i,160)>.4 && tlock_nul{conds,TT,sleep}.trialinfo(i,160)<.8]] 
                            ctrKc(sleep,4,conds)=ctrKc(sleep,4,conds)+1;
                            end
                                %% now make a list of onsets to compare distributions
                                if [tlock_nul{conds,TT,sleep}.trialinfo(i,70)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,110)>0 || tlock_nul{conds,TT,sleep}.trialinfo(i,150)>0] 
                                        if tlock_nul{conds,TT,sleep}.trialinfo(i,80)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,80)<1.5 
                                            tlock_nul{conds,TT,sleep}.theta.downstateKC.KC=[tlock_nul{conds,TT,sleep}.theta.downstateKC.KC tlock_nul{conds,TT,sleep}.trialinfo(i,80)];
                                            tlock_nul{conds,TT,sleep}.theta.downstateKC.KC_all=[tlock_nul{conds,TT,sleep}.theta.downstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,80)];
                                        elseif tlock_nul{conds,TT,sleep}.trialinfo(i,120)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,120)<1.5 
                                            tlock_nul{conds,TT,sleep}.theta.downstateKC.delta=[tlock_nul{conds,TT,sleep}.theta.downstateKC.delta tlock_nul{conds,TT,sleep}.trialinfo(i,120)];
                                            tlock_nul{conds,TT,sleep}.theta.downstateKC.KC_all=[tlock_nul{conds,TT,sleep}.theta.downstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,120)];
                                        elseif tlock_nul{conds,TT,sleep}.trialinfo(i,160)>.0 && tlock_nul{conds,TT,sleep}.trialinfo(i,160)<1.5 
                                            tlock_nul{conds,TT,sleep}.theta.downstateKC.sw=[tlock_nul{conds,TT,sleep}.theta.downstateKC.sw tlock_nul{conds,TT,sleep}.trialinfo(i,160)];
                                            tlock_nul{conds,TT,sleep}.theta.downstateKC.KC_all=[tlock_nul{conds,TT,sleep}.theta.downstateKC.KC_all tlock_nul{conds,TT,sleep}.trialinfo(i,160)];
                                        end
                                end
                            
                    end
                    
            end
            tlock_nul{conds,TT,sleep}.phaseKC.thetaNtrials=ctr(sleep,:,conds);
            tlock_nul{conds,TT,sleep}.phaseKC.thetaNkcs=ctrKc(sleep,:,conds);
            tlock_nul{conds,TT,sleep}.phaseKC.theta.KCrate=ctrKc(sleep,:,conds)./ctr(sleep,:,conds);
            
            tlock_nul_to_save{conds,TT,sleep}.phaseKC.thetaNtrials=ctr(sleep,:,conds);
            tlock_nul_to_save{conds,TT,sleep}.phaseKC.thetaNkcs=ctrKc(sleep,:,conds);
            tlock_nul_to_save{conds,TT,sleep}.phaseKC.theta.KCrate=ctrKc(sleep,:,conds)./ctr(sleep,:,conds);
            
            
         end
         
        end
        
end
    

if tacaud

    %% 20 for unisensory tactile
    KC{ii,20}.alpha.downstateKC.KC_all=[tlock_tac{10,3,12}.alpha.downstateKC.KC_all tlock_tac{10,3,13}.alpha.downstateKC.KC_all];
    KC{ii,20}.alpha.downmidstateKC.KC_all=[tlock_tac{10,3,12}.alpha.downmidstateKC.KC_all tlock_tac{10,3,13}.alpha.downmidstateKC.KC_all];
    KC{ii,20}.alpha.upmidstateKC.KC_all=[tlock_tac{10,3,12}.alpha.upmidstateKC.KC_all tlock_tac{10,3,13}.alpha.upmidstateKC.KC_all];
    KC{ii,20}.alpha.upstateKC.KC_all=[tlock_tac{10,3,12}.alpha.upstateKC.KC_all tlock_tac{10,3,13}.alpha.upstateKC.KC_all];
    
    KC{ii,20}.theta.downstateKC.KC_all=[tlock_tac{10,3,12}.theta.downstateKC.KC_all tlock_tac{10,3,13}.theta.downstateKC.KC_all];
    KC{ii,20}.theta.downmidstateKC.KC_all=[tlock_tac{10,3,12}.theta.downmidstateKC.KC_all tlock_tac{10,3,13}.theta.downmidstateKC.KC_all];
    KC{ii,20}.theta.upmidstateKC.KC_all=[tlock_tac{10,3,12}.theta.upmidstateKC.KC_all tlock_tac{10,3,13}.theta.upmidstateKC.KC_all];
    KC{ii,20}.theta.upstateKC.KC_all=[tlock_tac{10,3,12}.theta.upstateKC.KC_all tlock_tac{10,3,13}.theta.upstateKC.KC_all];
    
    %% 50 for nul
    KC{ii,50}.alpha.downstateKC.KC_all=[tlock_nul{10,3,12}.alpha.downstateKC.KC_all tlock_nul{10,3,13}.alpha.downstateKC.KC_all];
    KC{ii,50}.alpha.downmidstateKC.KC_all=[tlock_nul{10,3,12}.alpha.downmidstateKC.KC_all tlock_nul{10,3,13}.alpha.downmidstateKC.KC_all];
    KC{ii,50}.alpha.upmidstateKC.KC_all=[tlock_nul{10,3,12}.alpha.upmidstateKC.KC_all tlock_nul{10,3,13}.alpha.upmidstateKC.KC_all];
    KC{ii,50}.alpha.upstateKC.KC_all=[tlock_nul{10,3,12}.alpha.upstateKC.KC_all tlock_nul{10,3,13}.alpha.upstateKC.KC_all];
    
    KC{ii,50}.theta.downstateKC.KC_all=[tlock_nul{10,3,12}.theta.downstateKC.KC_all tlock_nul{10,3,13}.theta.downstateKC.KC_all];
    KC{ii,50}.theta.downmidstateKC.KC_all=[tlock_nul{10,3,12}.theta.downmidstateKC.KC_all tlock_nul{10,3,13}.theta.downmidstateKC.KC_all];
    KC{ii,50}.theta.upmidstateKC.KC_all=[tlock_nul{10,3,12}.theta.upmidstateKC.KC_all tlock_nul{10,3,13}.theta.upmidstateKC.KC_all];
    KC{ii,50}.theta.upstateKC.KC_all=[tlock_nul{10,3,12}.theta.upstateKC.KC_all tlock_nul{10,3,13}.theta.upstateKC.KC_all];
    
    for ttype=[20 50]

            [Kc_Cdf(ii,ttype).alpha.downstate.f, Kc_Cdf(ii,ttype).alpha.downstate.x]= ecdf(KC{ii,ttype}.alpha.downstateKC.KC_all); 
            [Kc_Cdf(ii,ttype).alpha.downmidstate.f, Kc_Cdf(ii,ttype).alpha.downmidstate.x]= ecdf(KC{ii,ttype}.alpha.downmidstateKC.KC_all);
            [Kc_Cdf(ii,ttype).alpha.upmidstate.f, Kc_Cdf(ii,ttype).alpha.upmidstate.x]= ecdf(KC{ii,ttype}.alpha.upmidstateKC.KC_all);
            [Kc_Cdf(ii,ttype).alpha.upstate.f, Kc_Cdf(ii,ttype).alpha.upstate.x]= ecdf(KC{ii,ttype}.alpha.upstateKC.KC_all);
            
            [Kc_Cdf(ii,ttype).theta.downstate.f, Kc_Cdf(ii,ttype).theta.downstate.x]= ecdf(KC{ii,ttype}.theta.downstateKC.KC_all); 
            [Kc_Cdf(ii,ttype).theta.downmidstate.f, Kc_Cdf(ii,ttype).theta.downmidstate.x]= ecdf(KC{ii,ttype}.theta.downmidstateKC.KC_all);
            [Kc_Cdf(ii,ttype).theta.upmidstate.f, Kc_Cdf(ii,ttype).theta.upmidstate.x]= ecdf(KC{ii,ttype}.theta.upmidstateKC.KC_all);
try            [Kc_Cdf(ii,ttype).theta.upstate.f, Kc_Cdf(ii,ttype).theta.upstate.x]= ecdf(KC{ii,ttype}.theta.upstateKC.KC_all); catch; end % cos it failed for e18
end
    
else 
    for conds=[1 3 4 5 6 7 9 10]
    %% for conds including Aud
    KC{ii,conds}.alpha.downstateKC.KC_all=[tlock_aud{conds,3,12}.alpha.downstateKC.KC_all tlock_aud{conds,3,13}.alpha.downstateKC.KC_all];
    KC{ii,conds}.alpha.downmidstateKC.KC_all=[tlock_aud{conds,3,12}.alpha.downmidstateKC.KC_all tlock_aud{conds,3,13}.alpha.downmidstateKC.KC_all];
    KC{ii,conds}.alpha.upmidstateKC.KC_all=[tlock_aud{conds,3,12}.alpha.upmidstateKC.KC_all tlock_aud{conds,3,13}.alpha.upmidstateKC.KC_all];
    KC{ii,conds}.alpha.upstateKC.KC_all=[tlock_aud{conds,3,12}.alpha.upstateKC.KC_all tlock_aud{conds,3,13}.alpha.upstateKC.KC_all];
    
    KC{ii,conds}.theta.downstateKC.KC_all=[tlock_aud{conds,3,12}.theta.downstateKC.KC_all tlock_aud{conds,3,13}.theta.downstateKC.KC_all];
    KC{ii,conds}.theta.downmidstateKC.KC_all=[tlock_aud{conds,3,12}.theta.downmidstateKC.KC_all tlock_aud{conds,3,13}.theta.downmidstateKC.KC_all];
    KC{ii,conds}.theta.upmidstateKC.KC_all=[tlock_aud{conds,3,12}.theta.upmidstateKC.KC_all tlock_aud{conds,3,13}.theta.upmidstateKC.KC_all];
    KC{ii,conds}.theta.upstateKC.KC_all=[tlock_aud{conds,3,12}.theta.upstateKC.KC_all tlock_aud{conds,3,13}.theta.upstateKC.KC_all];    
    end
    
    for ttype=[1 3 4 5 6 7 9 10 ]

            [Kc_Cdf(ii,ttype).alpha.downstate.f, Kc_Cdf(ii,ttype).alpha.downstate.x]= ecdf(KC{ii,ttype}.alpha.downstateKC.KC_all); 
            [Kc_Cdf(ii,ttype).alpha.downmidstate.f, Kc_Cdf(ii,ttype).alpha.downmidstate.x]= ecdf(KC{ii,ttype}.alpha.downmidstateKC.KC_all);
            [Kc_Cdf(ii,ttype).alpha.upmidstate.f, Kc_Cdf(ii,ttype).alpha.upmidstate.x]= ecdf(KC{ii,ttype}.alpha.upmidstateKC.KC_all);
            [Kc_Cdf(ii,ttype).alpha.upstate.f, Kc_Cdf(ii,ttype).alpha.upstate.x]= ecdf(KC{ii,ttype}.alpha.upstateKC.KC_all);
            
            [Kc_Cdf(ii,ttype).theta.downstate.f, Kc_Cdf(ii,ttype).theta.downstate.x]= ecdf(KC{ii,ttype}.theta.downstateKC.KC_all); 
            [Kc_Cdf(ii,ttype).theta.downmidstate.f, Kc_Cdf(ii,ttype).theta.downmidstate.x]= ecdf(KC{ii,ttype}.theta.downmidstateKC.KC_all);
            [Kc_Cdf(ii,ttype).theta.upmidstate.f, Kc_Cdf(ii,ttype).theta.upmidstate.x]= ecdf(KC{ii,ttype}.theta.upmidstateKC.KC_all);
            [Kc_Cdf(ii,ttype).theta.upstate.f, Kc_Cdf(ii,ttype).theta.upstate.x]= ecdf(KC{ii,ttype}.theta.upstateKC.KC_all);
end







condNames={ 'AT-500' 'AT-70' 'AT-20' 'AT' 'AT+20' 'AT+70' 'AT+500' 'A alone' 'T alone' 'Null'};
figure(ii)
hold on
ctr=1;
for ttype=[1 3 4 5 6 7 9 10 20 50];

subplot(10,2,ctr);
hold on
if ctr==1
    title('Alpha');
end
ylabel({'CDF'});
xlabel({'Time (s)'});
plot(Kc_Cdf(ii,ttype).alpha.downstate.x,Kc_Cdf(ii,ttype).alpha.downstate.f , 'b')
hold on; plot(Kc_Cdf(ii,ttype).alpha.downmidstate.x, Kc_Cdf(ii,ttype).alpha.downmidstate.f , 'g')
hold on; plot(Kc_Cdf(ii,ttype).alpha.upmidstate.x, Kc_Cdf(ii,ttype).alpha.upmidstate.f , 'r')
hold on; plot(Kc_Cdf(ii,ttype).alpha.upstate.x, Kc_Cdf(ii,ttype).alpha.upstate.f , 'y')
ctr=ctr+1;
subplot(10,2,ctr);
ylabel({'CDF'});
xlabel({'Time (s)'});
if ctr==2
    hold on
    title('Theta');
end
plot(Kc_Cdf(ii,ttype).theta.downstate.x,Kc_Cdf(ii,ttype).theta.downstate.f , 'b')
hold on; plot(Kc_Cdf(ii,ttype).theta.downmidstate.x,Kc_Cdf(ii,ttype).theta.downmidstate.f , 'g')
hold on; plot(Kc_Cdf(ii,ttype).theta.upmidstate.x,Kc_Cdf(ii,ttype).theta.upmidstate.f , 'r')
hold on; plot(Kc_Cdf(ii,ttype).theta.upstate.x,Kc_Cdf(ii,ttype).theta.upstate.f , 'y')
ctr=ctr+1;

end

suptitle('Cumulative Distributions by Phase');
saveas(gca,[num2str(ii) '_cdf.tif'])
end