function ksstat_sort=permKSstat(data1,data2)
% gives sorted 100 KSstat values from permutation distribution

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
  [h,p,ksstat(np)]=kstest2(permtrialvals(:,1),permtrialvals(:,2));
end
ksstat_sort=sort(ksstat);