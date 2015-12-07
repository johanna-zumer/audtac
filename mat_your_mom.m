function source=mat_your_mom(source)
% creating 'mommat' from .avg.mom cells
% assumes a .time field

source.mommat=nan(length(source.avg.mom),length(source.time));
for mm=1:length(source.avg.mom)
  if ~isempty(source.avg.mom{mm})
    source.mommat(mm,:)=source.avg.mom{mm};
  end
end

try source.avg=rmfield(source.avg,'mom');end
