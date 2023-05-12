% ImgPath='U:\Documents\Export (IBS_Fitting)\patient12_femur_l_aussen.ply';
% ptCloud= pcread(ImgPath);      
% LocMatrix=ptCloud.Location;
x1=MU(:,1);
y1=MU(:,2);
z1=MU(:,3);
x1 = x1(:); y1 = y1(:); z1 = z1(:);
P = [x1 y1 z1];
P = unique(P,'rows');
shp = alphaShape(P(:,1),P(:,2),P(:,3),9);
[F,V]=boundaryFacets(shp);
F=F.'; V=V.';
Q=mean(reshape(V(:,F),3,3,[]),2);
Q=num2cell(reshape(Q,3,[]).',1);
[x2,y2,z2]=deal(Q{:});   %"Interpolated" points

P = [x2 y2 z2];
shp = alphaShape(P(:,1),P(:,2),P(:,3),12);
[F,V]=boundaryFacets(shp);
F=F.'; V=V.';
Q=mean(reshape(V(:,F),3,3,[]),2);
Q=num2cell(reshape(Q,3,[]).',1);
[x3,y3,z3]=deal(Q{:});   %"Interpolated" points

P = [x3 y3 z3];
shp = alphaShape(P(:,1),P(:,2),P(:,3),12);
hold on
plot(shp);
light
lighting phong;
lightangle(0,-90)
% hg=scatter3(x1,y1,z1,'SizeData',50,'MarkerFaceColor','red','MarkerEdgeColor','none');
% hi=scatter3(x3,y3,z3,'SizeData',50,'MarkerFaceColor','blue','MarkerEdgeColor','none');
% legend([hg,hi],'Given Points','"Interpolated" Points');
hold off
