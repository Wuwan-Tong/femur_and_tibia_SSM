figure
for i=1:length(BoneData_r)
    tempbone=cell2mat8SyntheticBoneData_pre(i);
    scatter3(tempbone(:,1),tempbone(:,2),tempbone(:,3),'.');
    hold on
end
legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17')

figure
for i=1:length(BoneData_r)
    tempbone=cell2mat8SyntheticBoneData_pre(i);
    scatter3(tempbone(:,1),tempbone(:,2),tempbone(:,3),'.');
    hold on
end
legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18')