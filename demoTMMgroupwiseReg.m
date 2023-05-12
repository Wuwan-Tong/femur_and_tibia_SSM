load('syntheticBunnyData.mat')
load('GroundTruthRotations-Bunny.mat')
% Note: Samples 2, 3 and 4 in bunny data set are transformed and modified
% version of Sample 1.

[MU,Transform,TrainingSet,UP,PP,Mcoeffs,nu,convg,SSM]=TMMgroupwiseReg(bunnySet,30,150,1);

D1 = (Rgt{1}*(Transform(2).R*Transform(1).R')');
theta1 = acos((trace(D1) - 1)./2);
deg1 = (theta1*360)./(2*pi);
D2 = (Rgt{2}*(Transform(3).R*Transform(1).R')');
theta2 = acos((trace(D2) - 1)./2);
deg2 = (theta2*360)./(2*pi);
D3 = (Rgt{3}*(Transform(4).R*Transform(1).R')');
theta3 = acos((trace(D3) - 1)./2);
deg3 = (theta3*360)./(2*pi);

disp(['Angular rotation error for Sample 2 (in degrees) = ' num2str(deg1)])
disp(['Angular rotation error for Sample 3 (in degrees) = ' num2str(deg2)])
disp(['Angular rotation error for Sample 4 (in degrees) = ' num2str(deg3)])