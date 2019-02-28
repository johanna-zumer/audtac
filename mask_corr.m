function [rho,rho_nw,rhocond]=mask_corr(mask,runicaout,datareshape,plotflag)
% function [rho,rho_nw]=mask_corr(mask,runicaout,datareshape,plotflag)

mask_vector=reshape(mask,[1 size(mask,1)*size(mask,2)]);
mask_vector_rep=repmat(mask_vector,[1 7]);
mask_induse=find(mask_vector);  
mask_induse_rep=find(mask_vector_rep);  
for llind=1:7
  rica_masked(llind,:)=runicaout.trial{1}(llind,mask_induse);
  rica_masked_w(llind,:)=[rica_masked(llind,:)*runicaout.topo(1,llind) rica_masked(llind,:)*runicaout.topo(2,llind) rica_masked(llind,:)*runicaout.topo(3,llind) rica_masked(llind,:)*runicaout.topo(4,llind) rica_masked(llind,:)*runicaout.topo(5,llind) rica_masked(llind,:)*runicaout.topo(6,llind) rica_masked(llind,:)*runicaout.topo(7,llind)];
  rica_masked_nw(llind,:)=repmat(rica_masked(llind,:),[1 7]);
end
maskdatavec=reshape(datareshape(:,mask_induse)',[1 7*length(mask_induse)]);
for llind=1:7
  [rho(llind),pval(llind)]=corr(maskdatavec',rica_masked_w(llind,:)');
  [rho_nw(llind),pval_nw(llind)]=corr(maskdatavec',rica_masked_nw(llind,:)');
  for llcond=1:7
    rhocond(llind,llcond)=corr(rica_masked(llind,:)',datareshape(llcond,mask_induse)');
  end
end
if plotflag
  figure;imagesc(maskdatavec);caxis([-max(abs(maskdatavec)) max(abs(maskdatavec))])
  figure;imagesc(rica_masked_w(1,:));caxis([-max(abs(rica_masked_w(1,:))) max(abs(rica_masked_w(1,:)))])
end