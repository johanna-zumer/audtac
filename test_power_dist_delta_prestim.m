%%  Testing the power distribution of delta prestimulus

sleep=1;
if sleep==1
  subuse=setdiff(iiBuse,3:7);
  iter=11;
  trialkc=-1; % -1 all, 0 no Kc, 1 only Kc
elseif sleep==0
  subuse=iiSuse;
  iter=27;
end
audbineach=nan(4,9,19);
tacbineach=nan(4,9,19);
nulbineach=nan(4,9,19);
ms1bineach=nan(4,9,19);
soalist=[1 3 4 5 6 7 9];
subuseind=0;
for ii=subuse
  subuseind=subuseind+1;
  cd([edir sub{ii} ])
  %   try
  load(['absbin_' sub{ii} '_sleep' num2str(sleep) '_trialkc' num2str(trialkc) '.mat']);
  %   catch
  %     if trialkc==0
  %       load(['absbin_' sub{ii} '_sleep' num2str(sleep) '.mat']);
  %     else
  %       error('where is the file?')
  %     end
  %   end
  
  poweach4{subuseind}=[aud4mscon1abs{5,3,12}; tac4mscon1abs{5,3,12}; nul4mscon1abs{5,3,12}; ms4mscon1abs{5,3,12}]';
  poweach10{subuseind}=[aud4mscon1abs{5,3,12}; tac4mscon1abs{5,3,12}; nul4mscon1abs{5,3,12}; ms4mscon1abs{1,3,12}; ms4mscon1abs{3,3,12}; ms4mscon1abs{4,3,12}; ms4mscon1abs{5,3,12}; ms4mscon1abs{6,3,12}; ms4mscon1abs{7,3,12}; ms4mscon1abs{9,3,12}]';
  poweachIAF4{subuseind}=[aud4mscon1absIAF{5,3,12}; tac4mscon1absIAF{5,3,12}; nul4mscon1absIAF{5,3,12}; ms4mscon1absIAF{5,3,12}]';
  poweachIAF10{subuseind}=[aud4mscon1absIAF{5,3,12}; tac4mscon1absIAF{5,3,12}; nul4mscon1absIAF{5,3,12}; ms4mscon1absIAF{1,3,12}; ms4mscon1absIAF{3,3,12}; ms4mscon1absIAF{4,3,12}; ms4mscon1absIAF{5,3,12}; ms4mscon1absIAF{6,3,12}; ms4mscon1absIAF{7,3,12}; ms4mscon1absIAF{9,3,12}]';
  
  
  audeach{subuseind}=aud4mscon1abs{5,3,12}';
  taceach{subuseind}=tac4mscon1abs{5,3,12}';
  nuleach{subuseind}=nul4mscon1abs{5,3,12}';
  ms1each{subuseind}=ms4mscon1abs{1,3,12}';
  ms3each{subuseind}=ms4mscon1abs{3,3,12}';
  ms4each{subuseind}=ms4mscon1abs{4,3,12}';
  ms5each{subuseind}=ms4mscon1abs{5,3,12}';
  ms6each{subuseind}=ms4mscon1abs{6,3,12}';
  ms7each{subuseind}=ms4mscon1abs{7,3,12}';
  ms9each{subuseind}=ms4mscon1abs{9,3,12}';
  
  for ll=soalist
    audbineach(:,ll,subuseind)=mean(aud4mscon1absbin{ll,3,12});
    tacbineach(:,ll,subuseind)=mean(tac4mscon1absbin{ll,3,12});
    nulbineach(:,ll,subuseind)=mean(nul4mscon1absbin{ll,3,12});
    ms1bineach(:,ll,subuseind)=mean(ms4mscon1absbin{ll,3,12});
  end
  
end % ii

thresh=.05;
mean(audbineach<thresh,3)
mean(tacbineach<thresh,3)
mean(nulbineach<thresh,3)
mean(ms1bineach<thresh,3)


for ii=1:length(poweach4)
  figure(1);
  subplot(4,5,ii);hist(poweach4{ii},[0:10 20]);
  figure(2);
  subplot(4,5,ii);hist(poweach10{ii},[0:10 20]);
  figure(3);
  subplot(4,5,ii);hist(poweachIAF4{ii},[0:10 20]);
  figure(4);
  subplot(4,5,ii);hist(poweachIAF10{ii},[0:10 20]);
end
figure;hist([poweach10{:}],[0:10 20])
figure;hist([poweachIAF10{:}],[0:10 20])

sortpow=sort([poweach10{:}]);
sortpow(round(1*length(sortpow)/4))
%     0.9973
sortpow(round(2*length(sortpow)/4))
%     1.6479
sortpow(round(3*length(sortpow)/4))
%     2.5349

% % with Kcomplexes (trialkc = -1)
%     1.1966
%     2.0076
%     3.3088

figure;
subplot(2,5,1);hist([nuleach{:}],[0:10 20]);
subplot(2,5,2);hist([taceach{:}],[0:10 20]);
subplot(2,5,3);hist([audeach{:}],[0:10 20]);
subplot(2,5,4);hist([ms1each{:}],[0:10 20]);
subplot(2,5,5);hist([ms3each{:}],[0:10 20]);
subplot(2,5,6);hist([ms4each{:}],[0:10 20]);
subplot(2,5,7);hist([ms5each{:}],[0:10 20]);
subplot(2,5,8);hist([ms6each{:}],[0:10 20]);
subplot(2,5,9);hist([ms7each{:}],[0:10 20]);
subplot(2,5,10);hist([ms9each{:}],[0:10 20]);

