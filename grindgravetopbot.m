function [grindout,graveout]=grindgravetopbot(grindorig,sublist,subval)
% function [grindout,graveout]=grindgravetopbot(grindorig,sublist,subval)

grindout=grindorig;
grindout.powspctrm=grindorig.powspctrm(sublist==subval,:,:,:);
grindout.plvspctrm=grindorig.plvspctrm(sublist==subval,:,:,:);
grindout.plvabs=grindorig.plvabs(sublist==subval,:,:,:);

graveout=grindout;
graveout.dimord='chan_freq_time';
graveout.powspctrm=squeeze(nanmean(grindout.powspctrm,1));
graveout.plvspctrm=squeeze(nanmean(grindout.plvspctrm,1));
graveout.plvabs=squeeze(nanmean(grindout.plvabs,1));


