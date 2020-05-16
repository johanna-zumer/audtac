function kldiv=computeAT_KLdiv(data1,data2,edges,added,printflag)

[tmp1]=histcounts(data1,edges,'Normalization','probability');tmp1=tmp1+added;
[tmp2]=histcounts(data2,edges,'Normalization','probability');tmp2=tmp2+added;
if printflag
  figure(123);plot([tmp1; tmp2]');
end
kldiv=KLDiv(tmp1,tmp2);
