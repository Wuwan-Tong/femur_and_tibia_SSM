% SyntheticBoneData_pre_rot=cell(1,length(SyntheticBoneData_pre));
% SyntheticBoneData_pre_rot(1)=SyntheticBoneData_pre(1);
%%
k=35;
pc1=pointCloud(cell2mat(SyntheticBoneData_pre(1)));
pck=pointCloud(cell2mat(SyntheticBoneData_pre(k)));
% figure
% pcshowpair(pc1,pck)
% daspect([1 1 1])

[tformOut,pck_new] = pcregistericp(pck,pc1);
figure
pcshowpair(pc1,pck_new)
daspect([1 1 1])
%%
Rmat=[-1 0 0;
      0 -1 0;
      0 0 1];
mirror_img=pck_new.Location;
mirror_img(:,2)=-mirror_img(:,2);
mirror_img=mirror_img*Rmat;
pck_new_1=pointCloud(mirror_img);
% figure
% pcshowpair(pc1,pck_new_1)
% daspect([1 1 1])
[tformOut,pck_new_2] = pcregistericp(pck_new_1,pc1);
figure
pcshowpair(pc1,pck_new_2)
daspect([1 1 1])
%%
rotationAngles = [0 0 150];
translation = [0 0 0];
tform = rigidtform3d(rotationAngles,translation);
pck_new_3 = pctransform(pck_new_2,tform);
[tformOut,pck_new_4] = pcregistericp(pck_new_3,pc1);
figure
pcshowpair(pc1,pck_new_4)
daspect([1 1 1])
%%
LocMatrix=pck_new_2.Location;
[rows,cols]=size(LocMatrix); % temp matrix   
SyntheticBoneData_pre_rot(k)=mat2cell(LocMatrix,[rows],[cols]);

%%
pc1=pointCloud(cell2mat(SyntheticBoneData_pre_rot(1)));
figure
pcshow(pc1)
hold on
for i=2:35
    pc2=pointCloud(cell2mat(SyntheticBoneData_pre_rot(i)));
    pcshow(pc2)
end

