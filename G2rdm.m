function RDMout=G2rdm(Ginput,plotflag,plottext)
% function Gout=rdm2G(rdminput,plotflag,plottext)
% Johanna Zumer, 2017

if 0
  H=eye(size(Ginput))-ones(size(Ginput))/size(Ginput,1);
  iH=pinv(H);
  RDMout = -2 * iH * Ginput * iH; % inverse of Equation 22; Diedrichsen & Kriegeskorte, 2017
else
  C=pcm_indicatorMatrix('allpairs',1:size(Ginput,1));
  RDMout=squareform(diag(C*Ginput*C'));
end

if plotflag
  figure;subplot(1,2,1);imagesc(Ginput);colorbar;title(['G ' plottext])
  subplot(1,2,2);imagesc(RDMout);colorbar;title(['RDM ' plottext])
end

  
