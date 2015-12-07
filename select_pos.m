%% obsolete with modification to ft_selectdata for cfg.coordinate
function source = select_pos(source,posind)

fnames=fieldnames(source);
for ff=1:length(fnames)
  if ~isempty(strfind(fnames{ff},'dimord'))
    isdimord(ff)=1;
  else
    isdimord(ff)=0;
  end
end

dimordname={fnames{find(isdimord)}(1:findstr(fnames{find(isdimord)},'dimord')-1)}

dimtok=tokenize(source.(fnames{find(isdimord)}),'_');
posdim=find(strcmp(dimtok,'pos'));

for ff=1:length(fnames)
  tmp=strfind(dimordname,fnames{ff});
  if ~isempty(tmp{1})
    isreldimord(ff)=1;
  else
    isreldimord(ff)=0;
  end
end

if posdim==1
elseif posdim==2
source.(fnames{find(isreldimord)})(:,posind,:)
elseif posdim==3
elseif posdim=4
else
  error('not used to so many dimensions')
end
