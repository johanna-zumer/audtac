function vectout=binarise_artifact_begendpoints(artifact,N)
% function vectout=binarise_artifact_begendpoints(artifact,N);
%
% artifact should be Mx2 matrix with beginning and end samples

vectout=zeros(1,N);
for aa=1:size(artifact,1)
  vectout(artifact(aa,1):artifact(aa,2))=1;
end
vectout=vectout(1:N); % in case artifact extended longer than N


