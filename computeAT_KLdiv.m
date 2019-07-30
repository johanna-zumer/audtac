function kldiv=computeAT_KLdiv(data1,data2,edges,added)

[tmp1]=histcounts(data1,edges,'Normalization','probability');tmp1=tmp1+added;
[tmp2]=histcounts(data2,edges,'Normalization','probability');tmp2=tmp2+added;
kldiv=KLDiv(tmp1,tmp2);
