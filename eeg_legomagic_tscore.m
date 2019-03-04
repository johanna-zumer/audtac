eeg_legomagic_preamble

%%  Individual t-score

% ignore this for now.

time=-0.1:0.001:0.4;

ccc=nan(length(time),63,9,4,28);

for ii=5:18
  cd([edir sub{ii} ])
  load(['tlock_diffs_' sub{ii} '.mat']);
  if ii>7
    soalist=[1 3 4 5 6 7 9];
  else
    soalist=[3 4 5 6 7];
  end
  
  if ii==5
    labels=tlock_tac{5,1}.label;
  end
  
  %   chanuse=match_str(labels,tlock_tac{5,1}.label);
  %
  %   for ll=soalist
  %     for tt=1:4
  %
  %
  %       for teatime=1:length(time)
  %
  %         audtac=tlock_tac{ll,tt}.trial(:,:,dsearchn(tlock_tac{ll,tt}.time',time(teatime)));
  %         nul=tlock_nul_s0{ll+40,tt}.trial(:,:,dsearchn(tlock_nul_s0{ll+40,tt}.time',time(teatime)));
  %         tac=tlock_tac{10,tt}.trial(ttrialkept{ll,tt},:,dsearchn(tlock_tac{10,tt}.time',time(teatime)));
  %         aud=tlock_aud_s0{ll+40,tt}.trial(:,:,dsearchn(tlock_aud_s0{ll+40,tt}.time',time(teatime)));
  %
  %
  %         data=[tac; aud; audtac; nul];
  %         numtr=size(tac,1);
  %         design=zeros(size(data,1),5);
  %         design(:,1)=[zeros(0*numtr,1); ones(numtr,1); zeros(3*numtr,1)];
  %         design(:,2)=[zeros(1*numtr,1); ones(numtr,1); zeros(2*numtr,1)];
  %         design(:,3)=[zeros(2*numtr,1); ones(numtr,1); zeros(1*numtr,1)];
  %         design(:,4)=[zeros(3*numtr,1); ones(numtr,1); zeros(0*numtr,1)];
  %         design(:,5)=[zeros(0*numtr,1); ones(4*numtr,1); zeros(0*numtr,1)];
  %         cfg=[];
  %         cfg.glm.statistic='beta';
  %         cfg.glm.standardise=0;
  %         stat=ft_statfun_glm(cfg,data',design');
  %         beta=reshape(stat.stat,[size(design,2) size(data,2)]);
  %         con(teatime,chanuse,ll,tt,ii)=[1 1 -1 -1 0]*beta;
  %       end
  %     end
  %   end
  
  
  numtr=size(tlock_tac_s0{ll,tt}.trialinfo,1);
  mu=repmat(nanmean(tlock_tac_s0{ll,tt}.trial,1),[numtr 1 1]);
  sigma=repmat((nansum((tlock_tac_s0{ll,tt}.trial-mu).^2,1)/numtr).^0.5,[numtr 1 1]);
  zvalue=(tlock_tac_s0{ll,tt}.trial-mu)./sigma;
  zvalue=tlock_tac_s0{ll,tt}.trial./sigma;
  zvalue=squeeze(mu(1,:,:)./sigma(1,:,:));
  
  cfg.design=[];
  cfg.design(:,1)=[ ones(1,size(tlock_aud_s0{10,tt,ss}.trial,1)) ones(1,size(tlock_tac_s0{30+ll}.trial,1)) 2*ones(1,size(tlock_tac_s0{ll}.trial,1)) 2*ones(1,size(tlock_nul_s0{10,tt,ss}.trial,1))]';
  %   cfg.design(:,2)=[ ones(1,size(tlock_aud{10,tt,ss}.trial,1)) 2*ones(1,size(tlock_tac{30+ll}.trial,1)) 3*ones(1,size(tlock_tac{ll}.trial,1)) 4*ones(1,size(tlock_nul{10,tt,ss}.trial,1))]';
  stat1_s0{ll} = ft_timelockstatistics(cfg, tlock_aud_s0{10,tt,ss}, tlock_tac_s0{30+ll}, tlock_tac_s0{ll}, tlock_nul_s0{10,tt,ss})
  tlock_tpa_mtamn_s0{ll}.mask=stat1_s0{ll}.mask;
  
  cfg.design=[];
  cfg.design(:,1)=[ ones(1,size(tlock_aud_sall{10,tt,ss}.trial,1)) ones(1,size(tlock_tac_sall{30+ll}.trial,1)) 2*ones(1,size(tlock_tac_sall{ll}.trial,1)) 2*ones(1,size(tlock_nul_sall{10,tt,ss}.trial,1))]';
  %   cfg.design(:,2)=[ ones(1,size(tlock_aud{10,tt,ss}.trial,1)) 2*ones(1,size(tlock_tac{30+ll}.trial,1)) 3*ones(1,size(tlock_tac{ll}.trial,1)) 4*ones(1,size(tlock_nul{10,tt,ss}.trial,1))]';
  stat1_sall{ll} = ft_timelockstatistics(cfg, tlock_aud_sall{10,tt,ss}, tlock_tac_sall{30+ll}, tlock_tac_sall{ll}, tlock_nul_sall{10,tt,ss})
  tlock_tpa_mtamn_sall{ll}.mask=stat1_sall{ll}.mask;
end

% end
% thus, stat1_* is with Aud-alone at time zero, shifted tac, and AT with aud-first for 3, and tac-first for 7

for ll=soalist
  cfg.design=[];
  cfg.design(:,1)=[ ones(1,size(tlock_tac_s0{10,tt,ss}.trial,1)) ones(1,size(tlock_aud_s0{30+ll}.trial,1)) 2*ones(1,size(tlock_aud_s0{ll}.trial,1)) 2*ones(1,size(tlock_nul_s0{10,tt,ss}.trial,1))]';
  stat2_s0{ll} = ft_timelockstatistics(cfg, tlock_tac_s0{10,tt,ss}, tlock_aud_s0{30+ll}, tlock_aud_s0{ll}, tlock_nul_s0{10,tt,ss})
  
  cfg.design=[];
  cfg.design(:,1)=[ ones(1,size(tlock_tac_sall{10,tt,ss}.trial,1)) ones(1,size(tlock_aud_sall{30+ll}.trial,1)) 2*ones(1,size(tlock_aud_sall{ll}.trial,1)) 2*ones(1,size(tlock_nul_sall{10,tt,ss}.trial,1))]';
  stat2_sall{ll} = ft_timelockstatistics(cfg, tlock_tac_sall{10,tt,ss}, tlock_aud_sall{30+ll}, tlock_aud_sall{ll}, tlock_nul_sall{10,tt,ss})
end
% thus, stat2_* is with tac-alone at time zero, shifted aud, and AT with aud-first for 3, and tac-first for 7

save(['stat_erp_' sub{ii} '.mat'],'stat*')

figure;for ll=soalist,subplot(1,7,ll);imagesc(stat1_sall{ll}.time,1:62,stat1_sall{ll}.mask);end
figure;for ll=soalist,subplot(1,7,ll);imagesc(stat2_sall{ll}.time,1:62,stat2_sall{ll}.mask);end
figure;for ll=soalist,subplot(1,7,ll);imagesc(stat1_s0{ll}.time,1:62,stat1_s0{ll}.mask);end
figure;for ll=soalist,subplot(1,7,ll);imagesc(stat2_s0{ll}.time,1:62,stat2_s0{ll}.mask);end

cfg=[];
cfg.interactive='yes';
cfg.layout='EEG1010.lay';
ft_multiplotER(cfg,tlock_tpa_mtamn_sall{7});
ft_multiplotER(cfg,tlock_tpa_mtamn_sall{5});


clear *_s0
% end
