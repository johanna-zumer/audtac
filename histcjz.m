function N=histcjz(x,edges)
% function N=histcjz(x,edges)
%
% calls Matlab histc and does same thing, except if x is empty, then
% returns 'zeros' in N

if isempty(x)
  N=zeros(length(edges),1);
else
  N=histc(x,edges);
end

