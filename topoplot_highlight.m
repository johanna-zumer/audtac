function topoplot_highlight(fig,avg,timwin,stat,zlim)

if ~exist('zlim','var') || isempty(zlim)
  zlim=[-3 3];
end
% sample_count = length(avg.time);
% m = [1:timestep*fsample:sample_count];  % temporal endpoints in MEEG samples
if length(avg.time)>1 && length(timwin)>1
  fsample=1/(avg.time(2)-avg.time(1));
  if timwin(end)-timwin(1)<.025
    timestep=timwin(end)-timwin(1)
    subplotrows=1;
  else
    timestep = 0.025;      %(in seconds)
    subplotrows=4;
  end
  j = [timwin(1):timestep:timwin(2)];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
elseif length(avg.time)==1 || length(timwin)==1
  fsample=1;
  timestep = 0;      %(in seconds)
  j=[timwin(1) timwin(1)+.001];
  %   numplot=1;
  subplotrows=1;
else
  error('avg not have correct .time field')
end
numplot=length(j)-1;

if ~isempty(stat)
  % get relevant (significant) values
  if isfield(stat,'posclusters')
    if ~isempty(stat.posclusters)
      pos_cluster_pvals = [stat.posclusters(:).prob];
    else
      pos_cluster_pvals=1;
    end
    % pos_signif_clust = find(pos_cluster_pvals < statt.cfg.alpha);
    pos_signif_clust = find(pos_cluster_pvals < 0.05);
    pos = ismember(stat.posclusterslabelmat, pos_signif_clust);
  elseif length(stat.mask)==1 && stat.stat>0 && stat.mask
    pos=zeros(size(avg.avg));
    pos(match_str(avg.label,stat.cfg.channel),:)=1;
  else
    pos=zeros(size(avg.avg));
  end
  
  if isfield(stat,'negclusters')
    if ~isempty(stat.negclusters)
      neg_cluster_pvals = [stat.negclusters(:).prob];
    else
      neg_cluster_pvals=1;
    end
    neg_signif_clust = find(neg_cluster_pvals < 0.05);
    neg = ismember(stat.negclusterslabelmat, neg_signif_clust);
  elseif length(stat.mask)==1 && stat.stat<0 && stat.mask
    neg=zeros(size(avg.avg));
    neg(match_str(avg.label,stat.cfg.channel),:)=1;
  else
    neg=zeros(size(avg.avg));
  end
  
  m=dsearchn(stat.time',j');
else
  pos=zeros(size(avg.avg));
  neg=zeros(size(avg.avg));
  m=dsearchn(avg.time',j');
end

% define parameters for plotting
try
  close(fig)
end
figure(fig);
for k = 1:numplot
  subplot(subplotrows,ceil(numplot/4),k);
  cfg = [];
  cfg.xlim=[j(k) j(k+1)];
  cfg.zlim = zlim;
  cfg.marker='off';
  %      pos_int = all(pos(:, m(k):m(k+1)), 2);
  pos_int = any(pos(:, m(k):m(k+1)), 2);
  neg_int = any(neg(:, m(k):m(k+1)), 2);
  cfg.highlight = 'on';
  cfg.highlightchannel = unique([find(pos_int); find(neg_int)]);
  %      keyboard
  cfg.comment = 'xlim';
  cfg.commentpos = 'title';
  cfg.layout = 'elec1010.lay';
  if any(any(isnan(avg.avg(:,dsearchn(avg.time',cfg.xlim(1)):dsearchn(avg.time',cfg.xlim(2))))))
    continue
  else
    ft_topoplotER(cfg, avg);
  end
end
