function featind = feat_pre_during_evoked_individual(data,ll)
% function featind = feat_pre_during_evoked(data,ll)
%
% ll: refers to asynchrony index, or =0 for always using time=0ms.
% data: structure that has already combined either A+T or AT+N
%
% See litreview_evokedKc_evokedSp.xls for guidance on cutoff
% Kc/SW:  Evoked:  If negmax between 400-750 ms, where time window is after second stimulus
%         During (surrounding): If 'down'< 500ms before 1st stim, and 'upend' > 750ms after 2nd stim.
%
% Spindles: Sato et al. 2007 found spindles evoked anywhere from 0-3s (median 2s)

% Keep in mind the range should be suitable for either A or T to be relevant/evoked, and shift according to asynchrony.

% soades= [-.5 nan -.07 -.02 0 .02 .07 nan .5];
stim1st= [-.5 nan -.07 -.02 0   0   0 nan  0];
stim2nd= [  0 nan    0    0 0 .02 .07 nan .5];
numtr=size(data.trialinfo,1);

featind.KcDuring=zeros(1,numtr);
featind.KcEvoked=zeros(1,numtr);
featind.SpDuring=zeros(1,numtr);
featind.SpEvoked=zeros(1,numtr);

for kk=1:numtr
  numKc=length(data.kc{kk});
  for nk=1:numKc
    if [data.kc{kk}(nk).down<(.75+stim2nd(ll)) && data.kc{kk}(nk).down>(-.5+stim1st(ll))] || [data.kc{kk}(nk).upend<(.75+stim2nd(ll)) && data.kc{kk}(nk).upend>(-.5+stim1st(ll))]
      featind.KcDuring(kk)=1;
    end
    if data.kc{kk}(nk).negmax<0.75+stim2nd(ll) && data.kc{kk}(nk).negmax>0.4+stim2nd(ll)
      featind.KcEvoked(kk)=1;
    end
  end % nk
  numKc=length(data.sw{kk});
  for nk=1:numKc
    if [data.sw{kk}(nk).down<(.75+stim2nd(ll)) && data.sw{kk}(nk).down>(-.5+stim1st(ll))] || [data.sw{kk}(nk).upend<(.75+stim2nd(ll)) && data.sw{kk}(nk).upend>(-.5+stim1st(ll))]
      featind.KcDuring(kk)=1;
    end
    if data.sw{kk}(nk).negmax<0.75+stim2nd(ll) && data.sw{kk}(nk).negmax>0.4+stim2nd(ll)
      featind.KcEvoked(kk)=1;
    end
  end % nk
  
  numSp=length(data.sp_fast{kk});
  for nk=1:numSp
    if [data.sp_fast{kk}(nk).starttime<.1+stim2nd(ll) && data.sp_fast{kk}(nk).starttime<-.5+stim1st(ll)] || [data.sp_fast{kk}(nk).endtime<.1+stim2nd(ll) && data.sp_fast{kk}(nk).endtime<-.5+stim1st(ll)]
      featind.SpDuring(kk)=1;
    end
    if data.sp_fast{kk}(nk).starttime<3+stim2nd(ll) && data.sp_fast{kk}(nk).starttime>0.1+stim2nd(ll)
      featind.SpEvoked(kk)=1;
    end
  end % nk
  numSp=length(data.sp_slow{kk});
  for nk=1:numSp
    if [data.sp_slow{kk}(nk).starttime<.1+stim2nd(ll) && data.sp_slow{kk}(nk).starttime<-.5+stim1st(ll)] || [data.sp_slow{kk}(nk).endtime<.1+stim2nd(ll) && data.sp_slow{kk}(nk).endtime<-.5+stim1st(ll)]
      featind.SpDuring(kk)=1;
    end
    if data.sp_slow{kk}(nk).starttime<3+stim2nd(ll) && data.sp_slow{kk}(nk).starttime>0.1+stim2nd(ll)
      featind.SpEvoked(kk)=1;
    end
  end % nk
end  % kk



end

