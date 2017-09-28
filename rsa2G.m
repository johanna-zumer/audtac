function Gout=rsa2G(rsainput,plotflag,plottext)

rdm=1-rsainput;
H=eye(size(rsainput))-ones(size(rsainput))/size(rsainput,1); 
Gout=-0.5*H*rdm*H;

if plotflag
  figure;subplot(1,3,1);imagesc(rsainput);colorbar;title(['RSA ' plottext])
  subplot(1,3,2);imagesc(rdm);colorbar;title(['RDM ' plottext])
  subplot(1,3,3);imagesc(Gout);colorbar;title(['G ' plottext])
end

  
