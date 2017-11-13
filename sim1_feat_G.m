% sim1_feat_G

clear Mf*

Mf_1a(:,1)=[0 0 1 1 1 0 0]';
Gm_1a=Mf_1a*Mf_1a';
rdmGfeat_1a=G2rdm(Gm_1a,plotflag,'model 1a');

Mf_1b(:,1)=[0 1 1 1 1 1 0]';
Gm_1b=Mf_1b*Mf_1b';
rdmGfeat_1b=G2rdm(Gm_1b,plotflag,'model 1b');

Mf_2(:,1)=[1 1 1 .5 0 0 0]';
Mf_2(:,2)=[0 0 0 .5 1 1 1]';
Gm_2=Mf_2*Mf_2';
rdmGfeat_2=G2rdm(Gm_2,plotflag,'model 2');

Mf_3(:,1)=[1 0 0 0 0 0 1]';
Mf_3(:,2)=[0 1 0 0 0 1 0]';
Mf_3(:,3)=[0 0 1 0 1 0 0]';
Gm_3=Mf_3*Mf_3';
rdmGfeat_3=G2rdm(Gm_3,plotflag,'model 3');

C=pcm_indicatorMatrix('allpairs',[1:size(G_1a,1)]');
Coord_1a=pcm_classicalMDS(Gm_1a,'contrast',C);
Coord_1b=pcm_classicalMDS(Gm_1b,'contrast',C);
Coord_2 =pcm_classicalMDS(Gm_2,'contrast',C);
Coord_3 =pcm_classicalMDS(Gm_3,'contrast',C);
figure;
subplot(4,1,1);plot(Coord_1a(:,1),Coord_1a(:,2),'o');axis equal;
subplot(4,1,2);plot(Coord_1b(:,1),Coord_1b(:,2),'o');axis equal;
subplot(4,1,3);plot(Coord_2(:,1),Coord_2(:,2),'o');axis equal;
subplot(4,1,4);plot(Coord_3(:,1),Coord_3(:,2),'o');axis equal;

figure;
subplot(1,4,1);imagesc(Mf_1a);colorbar;
subplot(1,4,2);imagesc(Mf_1b);colorbar;
subplot(1,4,3);imagesc(Mf_2);colorbar;
subplot(1,4,4);imagesc(Mf_3);colorbar;
