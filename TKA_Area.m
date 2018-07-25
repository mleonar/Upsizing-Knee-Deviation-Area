TKA_Area Script:

function [minx,miny,minval,totarea,area, Cx,Cy]...
    =TKA_Area(x, y, dev, biomechid, txtflg,v,c,side,LOI,size,femsize,pathname)
% %Side 1 is Left, Side 2 is right
% %Determine if Medial or Lateral
cd(pathname)
if side==1 %If Left implant
    if txtflg==1 %if left side of left implant
        direction='Lateral'; %lateral
    else
        direction='Medial';
    end
else %if right implant
    if txtflg==1 %if left side of right implant
        direction='Medial';%medial
    else
        direction='Lateral';
    end
end
if size==1
    size='1-2';
elseif size==3
    size='3-4';
elseif size==5
    size='5-6';
else
    size='7-8';
end
%% Plot the entire deviation map of left or rigth side
figure
scatter(x,y,[],dev,'.')

axis(v)
caxis(c)
colormap(jet)
% colorbar

%%
xnew=[]; %preallocating for speed like zeros()
ynew=[];
devnew=[];
%Threshold for deviation
ind=find(dev<-0.08);
xnew=x(ind);
ynew=y(ind);
devnew=dev(ind);
if length(xnew)<10
    return %exit code or not enough points to run algorithm 
end
[minval, minind]=min(devnew); %most negative/max deviation {min value, which point]
minx=xnew(minind);
miny=ynew(minind);
h=figure;
scatter(xnew,ynew,[],devnew,'.')

axis(v)
caxis(c)
colormap(jet)
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
set(gca,'YTick',[])
set(gca,'XTick',[])

% colorbar
% title('Deviation')
% xlabel('x coordinates [mm]')
% ylabel('y coordinates [mm]')

box on
set(gca,'linewidth',1.5, 'fontweight','bold')
saveas(h,[biomechid '_' direction '_DeviationNoTitle.jpg'])

%%

data=[xnew,ynew]; %Mtraix=array array doesnt have to be matrix
figure
scatter(xnew,ynew,[],devnew,'.')

axis(v)
caxis(c)
colormap(jet)
% colorbar
flag=0;
cluster=0;
temparea=0;
pause(0.5)
tic
numdel=0;
while (flag<4 && cluster<50 && numdel>0) || (flag<4&&cluster<50&&numdel==0) %numdel= # of clusters deleted 
    cluster=cluster+1
    if rem(cluster,5)==0 || cluster==1
        h=figure;
        hold on
        axis(v)
%         title(['Number of Clusters: ' num2str(cluster)]);
%         xlabel('x coordinates [mm]')
%         ylabel('y coordinates [mm]')
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
set(gca,'YTick',[])
set(gca,'XTick',[])

        box on
        scatter(xnew,ynew,[],devnew, '.')
        set(gca,'linewidth',1.5, 'fontweight','bold')
        totarea=0;
        IDX=kmeans(data,cluster);
        allarea=zeros(cluster,1);
        numpix=zeros(cluster,1);
        for i=1:cluster
            index=find(IDX==i); %idx=row/point
            xtemp=xnew(index);
            ytemp=ynew(index);
            [k, area]=boundary(xtemp,ytemp,0.5);
            totarea=totarea+area; %tracking area of clusters together
            plot(xtemp(k),ytemp(k),'k')
            allarea(i)=area;
            numpix(i)=length(index);%number of points in cluster
        end
        S=std(allarea);
        M=mean(allarea);
        delclust=allarea(:,1)<(M-2*S); %probability to find cluster due to noise 
        delclust2=(numpix./allarea)<5; %gets rid of empty spaced clusters 
        delclust=delclust+delclust2;
        if sum(delclust>0) && numdel<3 %number of clusters deleted
            for i=1:cluster
                if delclust(i)>0 && numpix(i)<75 %finds if cluster should be deleted
                    index=find(IDX==i);
                    xnew(index)=[];
                    ynew(index)=[];
                    devnew(index)=[];
                    data=[xnew,ynew];
                    cluster=0;
                    temparea=0;
                    flag=0;
                    numdel=numdel+1;
                end
            end
            
        end
        clust2=(numpix./allarea)<20;
        if sum(clust2)>0
            temparea=0;
            flag=0;
        end
        axis(v)
        caxis(c)
        colormap(jet)
%         colorbar
        
        hold off
        saveas(h,[biomechid '_' direction '_Cluster' num2str(cluster) '.jpg'])
        
    else %if number of clusters isn't multiple of 5
        totarea=0;
        IDX=kmeans(data,cluster);
        allarea=zeros(cluster,1); %initializing/reseting
        numpix=zeros(cluster,1);
        for i=1:cluster
            index=find(IDX==i);
            xtemp=xnew(index);
            ytemp=ynew(index);
            [k, area]=boundary(xtemp,ytemp,0.5);
            totarea=totarea+area;
            allarea(i)=area;
            numpix(i)=length(index);
        end
        S=std(allarea);
        M=mean(allarea);
        delclust=allarea(:,1)<(M-2*S); 
        delclust2=(numpix./allarea)<5;
        delclust=delclust+delclust2;
        if sum(delclust>0) && numdel < 3
            for i=1:cluster
                if delclust(i)>0 && numpix(i)<75 %deleting cluster again
                    index=find(IDX==i);
                    xnew(index)=[];
                    ynew(index)=[];
                    devnew(index)=[];
                    data=[xnew,ynew];
                    cluster=0;
                    temparea=0;
                    flag=0;
                    numdel=numdel+1;
                end
            end
        end
        clust2=(numpix./allarea)<20;
        if sum(clust2)>0
            temparea=0;
            flag=0;
        end
    end
    if abs(temparea-totarea)<5 %flag refers to area based on iteration before and after
        flag=flag+1;
    else
        flag=0;
    end
    temparea=totarea;
    
end
toc %end timer
%% Delete noise

%% Centroid of Pressure
Cx=sum(xnew.*devnew)/sum(devnew); %weighted average
Cy=sum(ynew.*devnew)/sum(devnew);

%% Normalized Center of Pressure
% xL=max(x)-min(x);
% yL=max(y)-min(y);
% %X Distance from edge
% if side==1 %If Left implant
%     if txtflg==1 %if left side of left implant
%         direction='Lateral'; %lateral
%         normalizedx=(Cx-min(x))/xL*100;
%     else
%         direction='Medial';
%         normalizedx=(max(x)-Cx)/xL*100;
%     end
% else %if right implant
%     if txtflg==1 %if left side of right implant
%         direction='Medial';%medial
%         normalizedx=(Cx-min(x))/xL*100;
%     else
%         direction='Lateral';
%         normalizedx=(max(x)-Cx)/xL*100;
%     end
% end
% %Y Distance from maximum anterior point
% normalizedy=(max(y)-Cy)/yL*100;
%% plot cluster groups
h=figure;
hold on
scatter(xnew,ynew,[],devnew,'.')

tempmaxdev=zeros(length(cluster),3);
% caxis(c);
% colormap(jet)
% map=colormap;
% numcol=length(map);
% col=linspace(c(1),c(2),numcol);
for i=1:cluster
    index=find(IDX==i);
    xtemp=xnew(index);
    ytemp=ynew(index);
    devtemp=devnew(index);
    k=boundary(xtemp,ytemp,0.5);
    %         avgdev=mean(devtemp);
    %         ind=find(col>avgdev);
    %         plot(xtemp(k),ytemp(k),'Color',map(ind(1),:))
    plot(xtemp(k),ytemp(k),'k', 'markersize',5)
    [minvaltemp, minindtemp]=min(devtemp); %largest deviation
    tempmaxdev(i,:)=[xtemp(minindtemp) ytemp(minindtemp) minvaltemp]; %tells you coordinate & value of max deviatino of each cluster
end
[~,area]=boundary(x,y,0);
maxdev=sortrows(tempmaxdev,3); %sort by deviation in ascending order keeping coordinates (largest negative to smallest negative)
maxdev=maxdev(1,3);
% str={['Maximum Deviation=' num2str(minval)];[ 'at X=' num2str(minx) ' Y=' num2str(miny)]};

% title({[biomechid '; ' direction '; LOI=' num2str(LOI) 'months'];...
%     ['Insert Size=' size '; Femoral Size=' num2str(femsize)];...
%     ['X Centroid of Pressure: ' num2str(Cx) ...
%     'mm; Y Centroid of Pressure ' num2str(Cy) 'mm']})
% xlabel('x coordinates [mm]')
% ylabel('y coordinates [mm]')
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
set(gca,'YTick',[])
set(gca,'XTick',[])

box on
set(gca,'linewidth',1.5, 'fontweight','bold')
axis(v)
caxis(c)
colormap(jet)
% colorbar
% plot(Cx,Cy, 'rs','markersize',10, 'linewidth', 2)
hold off
saveas(h,[biomechid '_' direction '_clustersonly.jpg'])

%% plot deviation > threshold and cluster groups
h=figure;
hold on
scatter(x,y,[],dev,'.');
caxis(c)
colormap(jet)
% colorbar
for i=1:cluster
    index=find(IDX==i);
    xtemp=xnew(index);
    ytemp=ynew(index);
    devtemp=devnew(index);
    k=boundary(xtemp,ytemp,0.5);
    
    plot(xtemp(k),ytemp(k),'k', 'markersize',5)
    [minvaltemp, minindtemp]=min(devtemp);
    
    tempmaxdev(i,:)=[xtemp(minindtemp) ytemp(minindtemp) minvaltemp];
end
axis(v)
% set(gca,'YTickLabel',[])
% set(gca,'XTickLabel',[])
% xlabel('x coordinates [mm]','fontweight','bold')
% ylabel('y coordinates [mm]','fontweight','bold')

% title({[biomechid '; ' direction '; LOI=' num2str(LOI) 'months'];...
%     ['Insert Size=' size '; Femoral Size=' num2str(femsize)];...
%     ['Normalized Area of Deviation=' num2str(totarea/area*100) '%']});
plot(minx,miny, 'r.', 'MarkerSize', 20, 'LineWidth', 5)
str={['Maximum Deviation=' num2str(minval)];...
    [ 'at X=' num2str(minx) ' Y=' num2str(miny)];...
    ['Centroid of Pressure: ' ];[ 'X=' num2str(Cx)...
    'mm' ];['Y=' num2str(Cy) 'mm']};
% if txtflg ==1
%     text(-52, 21, str);
% else
%     text (33,21, str);
% end
box on
set(gca,'linewidth',1.5, 'fontweight','bold')
set(gca, 'FontSize', 20)

plot(Cx,Cy, 'rs','markersize',10, 'linewidth', 2)

hold off
saveas(h,[biomechid '_' direction '_clusters.jpg'])
%% plot deviation > threshold
h=figure;
hold on
scatter(xnew,ynew,[], devnew, '.')


% title({[biomechid '; ' direction '; LOI=' num2str(LOI) 'months'];...
%     ['Insert Size=' size '; Femoral Size=' num2str(femsize)];...
%     ['Normalized Area of Deviation=' num2str(totarea/area*100) '%']});
% ylabel('y coordinates [mm]','fontweight','bold')
% xlabel('x coordinates [mm]','fontweight','bold')
% set(gca,'YTickLabel',[])
% set(gca,'XTickLabel',[])
% set(gca,'YTick',[])
% set(gca,'XTick',[])
plot(Cx,Cy, 'rs','markersize',10, 'linewidth', 2)%Plot centroid
set(gca, 'FontSize', 14)
axis(v)
caxis(c)
colormap(jet)
% colorbar
box on
set(gca,'linewidth',1.5, 'fontweight','bold')
saveas(h,[biomechid '_' direction '_deviation.jpg'])
hold off
%% plot entire image
h=figure;
scatter(xnew,ynew,[], devnew, '.')
set(gca, 'FontSize', 14)
% title(biomechid)
% xlabel('x coordinates [mm]')
% ylabel('y coordinates [mm]')

% title({[biomechid '; ' direction '; LOI=' num2str(LOI) 'months'];...
%     ['Insert Size=' size '; Femoral Size=' num2str(femsize)];...
%     ['Normalized Area of Deviation=' num2str(totarea/area*100) '%']});
% ylabel('y coordinates [mm]','fontweight','bold')
% xlabel('x coordinates [mm]','fontweight','bold')
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
set(gca,'YTick',[])
set(gca,'XTick',[])
axis(v)
caxis(c)
colormap(jet)
% colorbar
box on
set(gca,'linewidth',1.5, 'fontweight','bold')
saveas(h,[biomechid '_justdev_' direction '.jpg'])
% close all
