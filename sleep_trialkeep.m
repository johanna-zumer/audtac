function [trlkeep1,trlkeep2]=sleep_trialkeep(feat1,feat2,fieldname)

cfg=[];
cfg.randomseed=13;
ft_preamble randomseed % sets call to rand in Shuffle same starting point, for reproducibility
feat1ind=find(feat1.(fieldname)==0);
feat2ind=find(feat2.(fieldname)==0);

trlnumKD=min(length(feat1ind),length(feat2ind));
tmp1=Shuffle(1:length(find(feat1.(fieldname)==0)));
trlkeep1=feat1ind(tmp1(1:trlnumKD));
tmp2=Shuffle(1:length(find(feat2.(fieldname)==0)));
trlkeep2=feat2ind(tmp2(1:trlnumKD));
