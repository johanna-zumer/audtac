function kldiv_sort=permKLdiv(data1,data2,edges,added)
% gives sorted 100 KLdiv values from permutation distribution

randomseed_ft(13);
trialvals1=[data1, data2]; % ntr x 2
numtr=size(trialvals1,1);
for np=1:100  % number of permutations
  permlabel=round(rand(1,numtr));
  for nt=1:numtr
    if permlabel(nt)
      permtrialvals(nt,:)=circshift(trialvals1(nt,:),[0 1]);
    else
      permtrialvals(nt,:)=trialvals1(nt,:);
    end
  end
  kldiv(np)=computeAT_KLdiv(permtrialvals(:,1),permtrialvals(:,2),edges,added,0);
end
kldiv_sort=sort(kldiv);