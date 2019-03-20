function featind = feat_pre_during_evoked(data,ll)
% function featind = feat_pre_during_evoked(data,ll)
%
% ll: refers to asynchrony index
% data: structure that has already combined either A+T or AT+N
%
% See litreview_evokedKc_evokedSp.xls for guidance on cutoff
% Kc/SW:  If negmax between 400-750 ms, then is KcEvoked, 
%         If down <0ms and upend>0, then KcDuring, where '0ms' means start of first stimulus
%         If upend between -0.5 to 0ms, then KcPre, where '0ms' means start of first stimulus
%
% Spindles: Sato et al. 2007 found spindles evoked anywhere from 0-3s (median 2s)

% Keep in mind the range should be suitable for either A or T to be relevant/evoked, and shift according to asynchrony.

soades= [-.5 nan -.07 -.02 0 .02 .07 nan .5];
preadd= [-.5 nan -.07 -.02 0   0   0 nan  0];
% postadd=[  0 nan    0    0 0 .02 .07 nan .5];
try
  numtr=size(data.trialinfo,1);
catch
  numtr=size(data.trialinfo1,1);
end

featind.KcPre=zeros(1,numtr);
featind.KcDuring=zeros(1,numtr);
featind.KcEvoked=zeros(1,numtr);
featind.SpPre=zeros(1,numtr);
featind.SpDuring=zeros(1,numtr);
featind.SpEvoked=zeros(1,numtr);

for kk=1:numtr
  numKc=length(data.kc{kk});
  for nk=1:numKc
    if data.kc{kk}(nk).upend<0+preadd(ll) && data.kc{kk}(nk).upend>-0.5+preadd(ll)
      featind.KcPre(kk)=1;
    end
    if data.kc{kk}(nk).down<0+preadd(ll) && data.kc{kk}(nk).upend>0+preadd(ll)
      featind.KcDuring(kk)=1;
    end
    if data.kc{kk}(nk).negmax<0.75+soades(ll) && data.kc{kk}(nk).negmax>0.4+soades(ll)  
      featind.KcEvoked(kk)=1;
    elseif data.kc{kk}(nk).negmax<0.75 && data.kc{kk}(nk).negmax>0.4  
      featind.KcEvoked(kk)=1;
    end
  end % nk
  numKc=length(data.sw{kk});
  for nk=1:numKc
    if data.sw{kk}(nk).upend<0+preadd(ll) && data.sw{kk}(nk).upend>-0.5+preadd(ll)
      featind.KcPre(kk)=1;
    end
    if data.sw{kk}(nk).down<0+preadd(ll) && data.sw{kk}(nk).upend>0+preadd(ll)
      featind.KcDuring(kk)=1;
    end
    if data.sw{kk}(nk).negmax<0.75+soades(ll) && data.sw{kk}(nk).negmax>0.4+soades(ll)  
      featind.KcEvoked(kk)=1;
    elseif data.sw{kk}(nk).negmax<0.75 && data.sw{kk}(nk).negmax>0.4
      featind.KcEvoked(kk)=1;
    end
  end % nk

  numSp=length(data.sp_fast{kk});
  for nk=1:numSp
    if data.sp_fast{kk}(nk).endtime<0+preadd(ll) && data.sp_fast{kk}(nk).endtime>-0.5+preadd(ll)
      featind.SpPre(kk)=1;
    end
    if data.sp_fast{kk}(nk).starttime<0+preadd(ll) && data.sp_fast{kk}(nk).endtime>0+preadd(ll)
      featind.SpDuring(kk)=1;
    end
    if data.sp_fast{kk}(nk).starttime<3+soades(ll) && data.sp_fast{kk}(nk).starttime>0.1+soades(ll)  %  Sato et al 2007
      featind.SpEvoked(kk)=1;
    elseif data.sp_fast{kk}(nk).starttime<3 && data.sp_fast{kk}(nk).starttime>0.1  %  Sato et al 2007
      featind.SpEvoked(kk)=1;
    end
  end % nk
  numSp=length(data.sp_slow{kk});
  for nk=1:numSp
    if data.sp_slow{kk}(nk).endtime<0+preadd(ll) && data.sp_slow{kk}(nk).endtime>-0.5+preadd(ll)
      featind.SpPre(kk)=1;
    end
    if data.sp_slow{kk}(nk).starttime<0+preadd(ll) && data.sp_slow{kk}(nk).endtime>0+preadd(ll)
      featind.SpDuring(kk)=1;
    end
    if data.sp_slow{kk}(nk).starttime<3+soades(ll) && data.sp_slow{kk}(nk).starttime>0.1+soades(ll)  %  Sato et al 2007
      featind.SpEvoked(kk)=1;
    elseif data.sp_slow{kk}(nk).starttime<3 && data.sp_slow{kk}(nk).starttime>0.1  %  Sato et al 2007
      featind.SpEvoked(kk)=1;
    end
  end % nk
end  % kk



end

