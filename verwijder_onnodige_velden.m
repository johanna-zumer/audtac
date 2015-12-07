function source=verwijder_onnodige_velden(source,timeflag)
% function source=remove_unnecessary_fields(source)

try source=rmfield(source,'filter');end
if timeflag, try source=rmfield(source,'time');end;end
try source.avg=rmfield(source.avg,'inside');end
try source.avg=rmfield(source.avg,'pos');end

try source=rmfield(source,'cfg');end
% try source.cfg=rmfield(source.cfg,'grid');end
% try source.cfg=rmfield(source.cfg,'vol');end
% try source.cfg=rmfield(source.cfg,'elec');end
