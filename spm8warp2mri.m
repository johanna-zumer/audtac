function mricoord=spm8warp2mri(t,mnicoord)
% function mricoord=spm8warp2mri(t,mnicoord)
%
% t is mri.params resulting from spm warping to MNI with subfields:
% .VG, .Tr, .Affine, .VF

% this is modified from code from NUTMEG nut_mni2mri, which is 
% modified from code originally written by John Ashburner:
% http://www.sph.umich.edu/~nichols/JG2/get_orig_coord2.m

if size(mnicoord,2)~=3, error('mnicoord must be an N x 3 matrix'); end;

Mat = inv(t.VG.mat);
xyz = Mat(1:3,:)*[mnicoord, ones(size(mnicoord,1))]';
Tr  = t.Tr;
Affine = t.Affine;
d   = t.VG.dim(1:3);
Mult = t.VF.mat*Affine;

if numel(Tr) == 0
    affine_only = 1;
    basX = 0; tx = 0;
    basY = 0; ty = 0;
    basZ = 0; tz = 0;
else
    affine_only = 0;
    basX = spm_dctmtx(d(1),size(Tr,1),xyz(1,:)-1);
    basY = spm_dctmtx(d(2),size(Tr,2),xyz(2,:)-1);
    basZ = spm_dctmtx(d(3),size(Tr,3),xyz(3,:)-1);
end;

if affine_only
    xyz2 = Mult(1:3,:)*[xyz ; ones(1,size(xyz,2))];
else
    Tr1 = reshape(Tr(:,:,:,1),size(Tr,1)*size(Tr,2),size(Tr,3));
    Tr2 = reshape(Tr(:,:,:,2),size(Tr,1)*size(Tr,2),size(Tr,3));
    Tr3 = reshape(Tr(:,:,:,3),size(Tr,1)*size(Tr,2),size(Tr,3));
    
    xyz2 = zeros(size(xyz));
    xyztmp = zeros(size(xyz));
    
    for i=1:size(xyz,2),
        bx = basX(i,:);
        by = basY(i,:);
        bz = basZ(i,:);
        % 		tx = reshape(...
        % 			reshape(Tr(:,:,:,1),size(Tr,1)*size(Tr,2),size(Tr,3))...
        % 			*bz', size(Tr,1), size(Tr,2) );
        % 		ty = reshape(...
        % 			reshape(Tr(:,:,:,2),size(Tr,1)*size(Tr,2),size(Tr,3))...
        % 			*bz', size(Tr,1), size(Tr,2) );
        % 		tz =  reshape(...
        % 			reshape(Tr(:,:,:,3),size(Tr,1)*size(Tr,2),size(Tr,3))...
        % 			*bz', size(Tr,1), size(Tr,2) );
        tx = reshape(Tr1*bz', size(Tr,1), size(Tr,2) );
        ty = reshape(Tr2*bz', size(Tr,1), size(Tr,2) );
        tz = reshape(Tr3*bz', size(Tr,1), size(Tr,2) );
        % 		xyz2(:,i) = Mult(1:3,:)*[xyz(:,i) + [bx*tx*by' ; bx*ty*by' ; bx*tz*by']; 1];
        xyztmp(:,i) = [bx*tx*by' ; bx*ty*by' ; bx*tz*by'];
    end;
    
    xyz2 = Mult(1:3,:)*[ xyz + xyztmp ; ones(1,size(xyz,2))];
    
end;
orig_coord = xyz2';

Vnorminput=t.VG;

mricoord=nut_coordtfm(orig_coord,t.VF.mat*inv(Vnorminput.mat));


