function plv=getplv(fourierspctrm);

fourierspctrm=fourierspctrm./abs(fourierspctrm);
plv=squeeze(mean(fourierspctrm,1));