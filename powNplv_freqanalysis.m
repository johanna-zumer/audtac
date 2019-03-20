function freqlo_tacPaud_tmp= powNplv_freqanalysis(cfg,tlock_fake)

freqlo_tacPaud_tmp=ft_freqanalysis(cfg,tlock_fake);
freqlo_tacPaud_tmp.powspctrm=squeeze(mean(abs(freqlo_tacPaud_tmp.fourierspctrm).^2,1));
freqlo_tacPaud_tmp.plvspctrm=getplv(freqlo_tacPaud_tmp.fourierspctrm);
%             freqlo_tacPaud_tmp{cc}.plvmag=abs(getplv(freqlo_tacPaud_tmp{cc}.fourierspctrm));
%             freqlo_tacPaud_tmp{cc}.plvang=angle(getplv(freqlo_tacPaud_tmp{cc}.fourierspctrm));
freqlo_tacPaud_tmp=rmfield(freqlo_tacPaud_tmp,'fourierspctrm');
