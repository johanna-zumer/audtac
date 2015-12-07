function grind = tidy_up_grind(grind)


grind.mommatdimord='subj_pos_time';
if size(grind.mommat,3)>1
  grind.mommat=grind.mommat(:,find(grind.inside),:);
else
  grind.mommat=grind.mommat(find(grind.inside),:);
end
grind.pos=grind.pos(find(grind.inside),:);
grind.inside=1:length(find(grind.inside));
