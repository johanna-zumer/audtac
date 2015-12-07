function topoplotTFR_highlight(fig,avg,timwin,ylim,base,stat,timestep,parameter,zlim);
% function topoplotTFR_highlight(fig,avg,timwin,ylim,base,stat,timestep,parameter)
%
% see ft_topoplotTFR for more description
%
% fig        figure handle number
% avg        data structure being plotted
% timwin     [starttime endtime]
% ylim       [startfreq endfreq]
% base       [base(1) base(2)] for baseline time window
% stat       stat structure with mask/clusterstats
% timestep   sampling 'rate' of subplots
% parameter  name of subfield of data to plot
% zlim       [zmin zmax] (colourscale limits)

if isempty(timwin)
  try
    timwin=[stat.time(1) stat.time(end)];
  catch
    timwin=[avg.time(1) avg.time(end)];
  end
end

if ~exist('timestep','var') || isempty(timestep)
  fsample=1/(avg.time(2)-avg.time(1));
  timestep = 1/fsample;      %(in seconds)
end
% sample_count = length(stat.time);
% j = [stat.time(1):timestep:stat.time(end)];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
j = [timwin(1):timestep:timwin(2)];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
% m = [1:timestep*fsample:sample_count];  % temporal endpoints in MEEG samples
numplot=length(j)-1;

if ~exist('parameter','var') || isempty(parameter)
  parameter='powspctrm';
end
if strcmp(parameter,'plvabsavg') && isfield(avg,'plvspctrm')
  avg.(parameter)=abs(avg.plvspctrm);
end
if strcmp(parameter,'plvangavg') && isfield(avg,'plvspctrm')
  avg.(parameter)=wrapToPi(unwrap(angle(avg.plvspctrm),[],3));
end
if strcmp(parameter,'plvavgang')
  avg.(parameter)=wrapToPi(avg.(parameter));
end
if strcmp(parameter,'plvspctrm')
  error('please specify whether abs or angle of plvspctrm');
end
if ~isfield(avg,parameter)
  error('this parameter not present in data');
end

if ~isempty(stat)
  % get relevant (significant) values
  if ~isempty(stat.posclusters)
    pos_cluster_pvals = [stat.posclusters(:).prob];
  else
    pos_cluster_pvals=1;
  end
  % pos_signif_clust = find(pos_cluster_pvals < statt.cfg.alpha);
  pos_signif_clust = find(pos_cluster_pvals < 0.05);
  pos = ismember(stat.posclusterslabelmat, pos_signif_clust);
  
  if ~isempty(stat.negclusters)
    neg_cluster_pvals = [stat.negclusters(:).prob];
  else
    neg_cluster_pvals=1;
  end
  % neg_signif_clust = find(neg_cluster_pvals < statt.cfg.alpha);
  neg_signif_clust = find(neg_cluster_pvals < 0.05);
  neg = ismember(stat.negclusterslabelmat, neg_signif_clust);
  
  yuse=find(stat.freq>=ylim(1) & stat.freq<=ylim(2));
  m=dsearchn(stat.time',j');
else
  pos=zeros(size(avg.(parameter)));
  neg=zeros(size(avg.(parameter)));
  m=dsearchn(avg.time',j');
  yuse=find(avg.freq>=ylim(1) & avg.freq<=ylim(2));
end

if ~exist('zlim','var')
  xuse=[dsearchn(avg.time',timwin(1)) dsearchn(avg.time',timwin(2))];
  
  if ~isempty(base)
    baseuse=[dsearchn(avg.time',base(1)) dsearchn(avg.time',base(2))];
    databn=avg.(parameter)(:,yuse(1):yuse(2),xuse(1):xuse(2))./repmat(nanmean(avg.(parameter)(:,yuse(1):yuse(2),baseuse(1):baseuse(2)),3),[1 1 diff(xuse)+1]);
    % the 0.9 factor is due to the temporal smoothing that occurs later, making the real max not represented in the temporally-averaged data
    maxz=.9*max(max(max(abs(databn))));
    %     zlim='maxabs';
    % this sets it to maxabs over all time steps not just the one plotted per subplot
    zlim=[-maxz maxz];
  else
    maxz=.9*max(max(max(abs(avg.(parameter)(:,yuse(1):yuse(2),xuse(1):xuse(2))))));
    minz=1.05*min(min(min(abs(avg.(parameter)(:,yuse(1):yuse(2),xuse(1):xuse(2))))));
    zlim=[minz maxz];
  end
end

% define parameters for plotting
try
  close(fig)
end
figure(fig);
for k = 1:numplot
  subplot(4,ceil(numplot/4),k);
  cfg = [];
  if ~isempty(base)
    cfg.baseline=[base(1) base(2)];
    cfg.baselinetype='relchange';
  end
  cfg.xlim=[j(k) j(k+1)];
  cfg.ylim=ylim;
  cfg.zlim = zlim;
  %      pos_int = all(pos(:, m(k):m(k+1)), 2);
  %   pos_int = any(pos(:, m(k):m(k+1)), 2);
  pos_int = any(reshape(pos(:, yuse, m(k):m(k+1)),[size(pos,1) length(yuse)*length(m(k):m(k+1))]) ,2);
  neg_int = any(reshape(neg(:, yuse, m(k):m(k+1)),[size(neg,1) length(yuse)*length(m(k):m(k+1))]) ,2);
  cfg.highlight = 'on';
  cfg.highlightchannel = [find(pos_int) find(neg_int)];
  %      keyboard
  cfg.comment = 'jz';
  cfg.commentpos = 'title';
  cfg.layout = 'elec1010.lay';
  cfg.parameter=parameter;
  if any(any(any(isnan(avg.(parameter)(:,yuse(1):yuse(2),dsearchn(avg.time',cfg.xlim(1)):dsearchn(avg.time',cfg.xlim(2)))))))
    continue
  else
    ft_topoplotTFR(cfg, avg);
  end
  if k==1
    colorbar;
  end
end
