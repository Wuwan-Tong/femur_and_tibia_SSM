function [] = plotSSM(SSM,weightEigVec,plotType)
% SMM visualization on each PC

% input: 
%    plotType: 'points' or 'mesh'
%    weightEigVec: weight of standard deviation


% number of PCs anf Points
numPC=length(SSM.eVals);
numPts=length(SSM.MU)/3;

%calculate plus and minus on each PC
PCPlus=zeros(3*numPts,numPC);
PCMinus=zeros(3*numPts,numPC);
for i=1:numPC
    PCPlus(:,i)=SSM.MU+(weightEigVec*sqrt(SSM.eVals(i))*SSM.eVecs(:,i))';
    PCMinus(:,i)=SSM.MU-(weightEigVec*sqrt(SSM.eVals(i))*SSM.eVecs(:,i))';
end

% reshape 1-by-(3*numPts) array to numPts-by-3 matrix, for plot
PtsPlus=zeros(numPts,3,numPC);
PtsMinus=zeros(numPts,3,numPC);
PtsMiddle=reshape(SSM.MU,3,[])';
for i=1:numPC
    PtsPlus(:,:,i)=reshape(PCPlus(:,i),3,[])';
    PtsMinus(:,:,i)=reshape(PCMinus(:,i),3,[])';
end

% plot middle, plus and minus for each PCs
% [t,tnorm]=MyRobustCrust(PtsMiddle);
[t] = MyCrustOpen(PtsMiddle); %Triangle connectivity, specified as a 3-column matrix where each row contains the point vertices defining a triangle face.
for i=1:3 % plot only first 3 PCs, write for i=1:numPC to plot all PCs
    % set the title of each subplots
    figNamePlus=append ('PC',num2str(i),' Plus, (mean+',num2str(weightEigVec),'\sigma)');
    figNameMinus=append ('PC',num2str(i),' Minus,(mean-',num2str(weightEigVec),'\sigma)');   
    figNameMiddle=append (' Mean');
    
    figure     
    colormap bone
    set(gcf,'Color',[1 1 0.88]) % set color for figure handle   
    view(0,0) % view(az,el) set the viewing angle for a three-dimensional plot, az, is the horizontal rotation 
    % about the z-axis as measured in degrees from the negative y-axis. Positive values indicate counterclockwise 
    % rotation of the viewpoint. el is the vertical elevation of the viewpoint in degrees
    subplot(1,3,1);
    if strcmp(plotType,'mesh')
        trisurf(t,PtsMinus(:,1,i),PtsMinus(:,2,i),PtsMinus(:,3,i),'Edgecolor','none');
        light
        lighting phong;
        lightangle(0,-90)
    elseif strcmp(plotType,'points')
        scatter3(PtsMinus(:,1,i),PtsMinus(:,2,i),PtsMinus(:,3,i),'.');
    end
    title(figNameMinus);
    set(gca, 'visible', 'off') % axis not wisible
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
    axis equal
    

    subplot(1,3,2);
    if strcmp(plotType,'mesh')
        trisurf(t,PtsMiddle(:,1),PtsMiddle(:,2),PtsMiddle(:,3),'Edgecolor','none');
        light
        lighting phong;
        lightangle(0,-90)
    elseif strcmp(plotType,'points')
        scatter3(PtsMiddle(:,1),PtsMiddle(:,2),PtsMiddle(:,3),'.');
    end
    title(figNameMiddle);
    set(gca, 'visible', 'off') % axis not wisible
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);    
    axis equal
    

    subplot(1,3,3);
    if strcmp(plotType,'mesh')
        trisurf(t,PtsPlus(:,1,i),PtsPlus(:,2,i),PtsPlus(:,3,i),'Edgecolor','none');
        light
        lighting phong;
        lightangle(0,-90)
    elseif strcmp(plotType,'points')
        scatter3(PtsPlus(:,1,i),PtsPlus(:,2,i),PtsPlus(:,3,i),'.');
    end
    title(figNamePlus);
    set(gca, 'visible', 'off') % axis not wisible
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);   
    axis equal
    
end

%% plot all pcs in one figure
figure
view(0,0) 
for i=1:3 % plot only first 3 PCs, write for i=1:numPC to plot all PCs        
    % set the title of each subplots
    figNamePCs=append('PC',num2str(i));
    
    colormap bone
    set(gcf,'Color',[1 1 0.88]) % set color for figure handle   
    
    subplot(2,4,i);
    if strcmp(plotType,'mesh')
        trisurf(t,PtsMinus(:,1,i),PtsMinus(:,2,i),PtsMinus(:,3,i),'Edgecolor','none');
        light
        lighting phong;
        lightangle(0,-90)
    elseif strcmp(plotType,'points')
        scatter3(PtsMinus(:,1,i),PtsMinus(:,2,i),PtsMinus(:,3,i),'.');
    end
    title(figNamePCs);
    if i==1
        text(min(xlim),(max(ylim)-min(ylim))*0.6,'-3\sigma')
    end
    set(gca, 'visible', 'off') % axis not wisible
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
    axis equal
   
    subplot(2,4,i+4);
    if strcmp(plotType,'mesh')
        trisurf(t,PtsPlus(:,1,i),PtsPlus(:,2,i),PtsPlus(:,3,i),'Edgecolor','none');
        light
        lighting phong;
        lightangle(0,-90)
    elseif strcmp(plotType,'points')
        scatter3(PtsPlus(:,1,i),PtsPlus(:,2,i),PtsPlus(:,3,i),'.');
    end
    if i==1
        text(min(xlim),(max(ylim)-min(ylim))*0.6,'+3\sigma')
    end
    set(gca, 'visible', 'off') % axis not wisible
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);   
    axis equal
    
end
subplot(2,4,[4,8]);
if strcmp(plotType,'mesh')
        trisurf(t,PtsMiddle(:,1),PtsMiddle(:,2),PtsMiddle(:,3),'Edgecolor','none');
        light
        lighting phong;
        lightangle(0,-90)
elseif strcmp(plotType,'points')
    scatter3(PtsMiddle(:,1),PtsMiddle(:,2),PtsMiddle(:,3),'.');
end
title('mean');
set(gca, 'visible', 'off') % axis not wisible
set(findall(gca, 'type', 'text'), 'visible', 'on')
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);   
axis equal
end