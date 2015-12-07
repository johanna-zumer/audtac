% simulation of expected ERP peak height with jitter and SNR

erp=spm_hrf(.1,[20 30 .5 .5 6 0 50])/max(spm_hrf(.1,[20 30 .5 .5 6 0 50])); %amplitude normalised to 1
figure;plot(0:.25:125,erp);  % "sampling rate" is like 4000hz

noiseval=[0 1 5 10 50];
for nn=1:length(noiseval)
  for ii=1:100
    erpn(:,nn,ii)=erp+noiseval(nn)*pinknoise(length(erp))';
  end
end

shiftms=[0 1 2 5 10 20]; % plus/minus that amount
for ss=1:length(shiftms)  
  for ii=1:100
    erpns(:,:,ii,ss)=circshift(erpn(:,:,ii),[round(4*shiftms(ss)*[2*round(rand(1))-1]*rand(1)) 0]);
  end
end

figure;plot(squeeze(max(mean(erpns,3),[],1))')

%%  Simulate adding & averaging order

a=randn(1,100);
b=randn(1,100);

mean(a)+mean(b)

mean(a+b)


