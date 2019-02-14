function featind = feat_pre_during_evoked(data,ll)
% function featind = feat_pre_during_evoked(data,ll)

soades=[-.5 nan -.07 -.02 0 .02 .07 nan .5];
numtr=size(data.trialinfo,1);

featind.KcPre=zeros(1,numtr);
featind.KcDuring=zeros(1,numtr);
featind.KcEvoked=zeros(1,numtr);
featind.SpPre=zeros(1,numtr);
featind.SpDuring=zeros(1,numtr);
featind.SpEvoked=zeros(1,numtr);

for kk=1:numtr
  numKc=length(data.kc{kk});
  for nk=1:numKc
    if data.kc{kk}(nk).upend<0+soades(ll) && data.kc{kk}(nk).upend>-0.5+soades(ll)
      featind.KcPre(kk)=1;
    end
    if data.kc{kk}(nk).down<0+soades(ll) && data.kc{kk}(nk).upend>0+soades(ll)
      featind.KcDuring(kk)=1;
    end
    if data.kc{kk}(nk).negmax<0.8+soades(ll) && data.kc{kk}(nk).negmax>0.1+soades(ll)  %  Dang Vu 2011
      featind.KcEvoked(kk)=1;
    end
  end % nk
  numKc=length(data.sw{kk});
  for nk=1:numKc
    if data.sw{kk}(nk).upend<0+soades(ll) && data.sw{kk}(nk).upend>-0.5+soades(ll)
      featind.KcPre(kk)=1;
    end
    if data.sw{kk}(nk).down<0+soades(ll) && data.sw{kk}(nk).upend>0+soades(ll)
      featind.KcDuring(kk)=1;
    end
    if data.sw{kk}(nk).negmax<0.8+soades(ll) && data.sw{kk}(nk).negmax>0.1+soades(ll)  %  Dang Vu 2011
      featind.KcEvoked(kk)=1;
    end
  end % nk

  numSp=length(data.sp_fast{kk});
  for nk=1:numSp
    if data.sp_fast{kk}(nk).endtime<0+soades(ll) && data.sp_fast{kk}(nk).endtime>-0.5+soades(ll)
      featind.SpPre(kk)=1;
    end
    if data.sp_fast{kk}(nk).starttime<0+soades(ll) && data.sp_fast{kk}(nk).endtime>0+soades(ll)
      featind.SpDuring(kk)=1;
    end
    if data.sp_fast{kk}(nk).starttime<3+soades(ll) && data.sp_fast{kk}(nk).starttime>0.1+soades(ll)  %  Sato et al 2007
      featind.SpEvoked(kk)=1;
    end
  end % nk
  numSp=length(data.sp_slow{kk});
  for nk=1:numSp
    if data.sp_slow{kk}(nk).endtime<0+soades(ll) && data.sp_slow{kk}(nk).endtime>-0.5+soades(ll)
      featind.SpPre(kk)=1;
    end
    if data.sp_slow{kk}(nk).starttime<0+soades(ll) && data.sp_slow{kk}(nk).endtime>0+soades(ll)
      featind.SpDuring(kk)=1;
    end
    if data.sp_slow{kk}(nk).starttime<3+soades(ll) && data.sp_slow{kk}(nk).starttime>0.1+soades(ll)  %  Sato et al 2007
      featind.SpEvoked(kk)=1;
    end
  end % nk
end  % kk



end

