# Upsizing-Knee-Deviation-Area
Thu's code that can be used to calculate the deviation area on inserts to study the effects of differential sizing in total knee arthroplasty

Main2 script:

clear all
close all
clc
currentfolder=pwd;
directoryname=uigetdir;
[filename2, pathname2]=uigetfile('*.xlsx', 'Select MasterDemographics File');

cd(directoryname)
D=dir('*.csv');
names={D.name};
clear D
dat=zeros(length(names),21);

cd(pathname2)
[num, txt,raw]=xlsread([pathname2 '\' filename2]);
for i=1:length(raw)
    raw{i,1}=num2str(raw{i,1});
end
v=[ -37.3087   36.8055  -29.7594   37.6188];

for z=1:length(names)
%     close all
    cd(directoryname)
    
    data=dlmread(names{z});
    datasort=sortrows(data); %Sort data
    clear data
    %% Isolate bearing surface
    condition=(datasort(:,3)<5);
    datasort(condition,:)=[];
    
    c=[-.75, 0.75]; %colorbar scale
    figure
    scatter(datasort(:,1),datasort(:,2),[],datasort(:,4),'.')
    caxis(c)
    colormap(jet)
    colorbar
%     axis equal
%     v=axis;
axis(v);
    set(gca,'fontsize',14)
    xlabel('x coordinates [mm]','fontweight','bold')
ylabel('y coordinates [mm]','fontweight','bold')
box on
set(gca,'linewidth',1.5, 'fontweight','bold')
    condition=(datasort(:,1)>-6.8 & datasort(:,1)<6.8);
    datasort(condition,:)=[];
    clear condition
    
    
    biomechid=names{z}(1:8)
    
    index=strfind((raw(1:end,1)),biomechid);
    count=0;
    for i=1:length(index)
        count=count+1;
        if index{i}==1
            index=count;
            break
        end
    end
    %% Find Insert type and size
    if strcmp(raw{index,19},'Gen II HF')
        type=1;
    else
        type=2;
    end
    if strcmp(raw{index,23},'1-2')
        size=1;
    elseif strcmp(raw{index,23},'3-4')
        size=3;
    elseif strcmp(raw{index,23},'5-6')
        size=5;
    else
        size=7;
    end
    femsize=raw{index,25};
    %% Find Left or Right
    if strcmp(raw{index,13},'L')
        side=1;
        side1='Left';
    else
        side=2;
        side1='Right';
    end
    LOI=raw{index,16};
    
    %% Separate Data
    x=datasort(:,1);
    y=datasort(:,2);
    dev=datasort(:,4);
    clear data
    %% Center the x and y
    deltax=max(x)+min(x);
    x=x-deltax/2;
    deltay=max(y)+min(y);
    y=y-deltay/2;
    %% Separate left or right
    count=1;
    while x(count)<0;
        count=count+1;
    end
    ind=count-1;
    clear count
    
    
    xL=x(1:ind);
    xR=x(ind+1:end);
    yL=y(1:ind);
    yR=y(ind+1:end);
    devL=dev(1:ind);
    devR=dev(ind+1:end);
    clear x y dev
    %% Run analysis separately, flag 1 is left side, flag 2 is right side
    cd(currentfolder)
    flag=1;
%     v=[-54.4909, 5.8246, -27.0000, 27.0000]; %axes scale
    [xL,yL, devL]=cropdata2new(xL,yL,devL, flag,v,size,type, biomechid);
%     [minxL,minyL,maxdevL,devareaL,areaL, CxL,CyL]...
%         =TKA_Area(xL, yL, devL, biomechid, flag,v,c,side, LOI,size,femsize,pathname2);
    flag=2;
        cd(currentfolder)

    % -7.5877   55.2187  -21.6935   27.8425
%     v=[-5.8246,54.4909, -27.0000, 27.0000]; %axes scale
    [xR,yR, devR]=cropdata2new(xR,yR,devR, flag,v, size, type,biomechid);
    x=[xL;xR];
    y=[yL;yR];
    dev=[devL;devR];
%     [minxR,minyR,maxdevR,devareaR,areaR, CxR,CyR]...
%         =TKA_Area(xR, yR, devR, biomechid, flag,v,c,side, LOI,size,femsize,pathname2);
    %% Plot the entire deviation map
h=figure;
scatter(x,y,[],dev,'.')

axis(v)
caxis(c)
colormap(jet)
colorbar('southoutside')
title ([num2str(round(LOI)) ' Months'])
% colorbar
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
set(gca,'YTick',[])
set(gca,'XTick',[])

set(gca,'fontsize',14)
box on
set(gca,'linewidth',1.5, 'fontweight','bold')
saveas(h,[biomechid '_ioslatedsurface' side1 '.jpg'])

%%
xnew=[];
ynew=[];
devnew=[];
%Threshold for deviation
ind=find(dev<-0.08);
xnew=x(ind);
ynew=y(ind);
devnew=dev(ind);
if length(xnew)<10
    return
end
[minval, minind]=min(devnew);
minx=xnew(minind);
miny=ynew(minind);
h=figure;
scatter(xnew,ynew,[],devnew,'.')
title ([num2str(round(LOI)) ' Months'])
axis(v)
caxis(c)
colormap(jet)
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
set(gca,'YTick',[])
set(gca,'XTick',[])

% axis equal
caxis(c)
colormap(jet)
% colorbar
set(gca,'fontsize',14)
box on
set(gca,'linewidth',1.5, 'fontweight','bold')
saveas(h,[biomechid '_ioslated_dev' side1  '.jpg'])
    %%
    
%     dat(z,1)=str2num(biomechid);
%     dat(z,2)=LOI;
%     dat(z,3)=side; %1 is left, 2 is right
%     dat(z,4)=size;
%     dat(z,5)=minxL;
%     dat(z,6)=minyL;
%     dat(z,7)=maxdevL;
%     dat(z,8)=devareaL;
%     dat(z,9)=areaL;
%     dat(z,10)=devareaL/areaL*100;
%     dat(z,11)=CxL;
%     dat(z,12)=CyL;
%     dat(z,14)=minxR;
%     dat(z,15)=minyR;
%     dat(z,16)=maxdevR;
%     dat(z,17)=devareaR;
%     dat(z,18)=areaR;
%     dat(z,19)=devareaR/areaR*100;
%     dat(z,20)=CxR;
%     dat(z,21)=CyR;
    
end
xlswrite(['DeviationData' datestr(now,'yyyymmdd') '.xlsx'],dat)



Cropdata2new Script:

 function [xcrop,ycrop,devcrop]=cropdata2new(x,y, dev, txtflg,v, size,...
 type,biomechid)
[k,~]=boundary(x,y);
figure
hold on
plot(x(k),y(k))
axis equal
scatter(x,y,[],dev,'.');

% keyboard 

%% Initialize
[~,yidxmin]=min(y);
[~,yidxmax]=max(y);
L6=0;
L7=0;
xcrop=[];
ycrop=[];
devcrop=[];
L=sqrt((x(yidxmin)-x(yidxmax))^2+(y(yidxmin)-y(yidxmax))^2);
L1=L/2;
if type==2
    x0=(x(yidxmin)+x(yidxmax))/2*1.5;
    y0=(y(yidxmin)+y(yidxmax))/2*1.5;
else
    x0=(x(yidxmin)+x(yidxmax))*1.1;
    y0=(y(yidxmin)+y(yidxmax))*1.1;
end

%% if PS
if type==2
    if txtflg==1; % if left side
        [~,idx]=min(x);
        
        t=-pi/2:0.01:pi/2;
        %     ind=find(x<x0);
        %     xcrop=x(ind);
        %     ycrop=y(ind);
        %     devcrop=dev(ind);
        L2=sqrt((x(idx)-x0)^2+(y(idx)-y0)^2);
        [L3,L4,L5,x0]=cropsize(L2,L1, size,type,txtflg,x0, L6,L7);
        for i=1:length(x)
            if x(i)<=x0
                xcrop=[xcrop; x(i)];
                ycrop=[ycrop; y(i)];
                devcrop=[devcrop; dev(i)];
            elseif ((x(i)-x0)/L3)^2+((y(i)-y0)/L5)^2<1 && (y(i)>=0)
                xcrop=[xcrop; x(i)];
                ycrop=[ycrop; y(i)];
                devcrop=[devcrop; dev(i)];
            elseif ((x(i)-x0)/L3)^2+((y(i)-y0)/L4)^2<1 && (y(i)<0)
                xcrop=[xcrop; x(i)];
                ycrop=[ycrop; y(i)];
                devcrop=[devcrop; dev(i)];
            end
        end
    else
        %% if Right Side
        [~,idx]=max(x);
        
        L2=sqrt((x(idx)-x0)^2+(y(idx)-y0)^2);
        [L3,L4,L5,x0]=cropsize(L2,L1, size,type,txtflg,x0,L6,L7);
        t=pi/2:0.01:3*pi/2;
        for i=1:length(x)
            if x(i)>=x0
                xcrop=[xcrop; x(i)];
                ycrop=[ycrop; y(i)];
                devcrop=[devcrop; dev(i)];
            elseif ((x(i)-x0)/L3)^2+((y(i)-y0)/L5)^2<1 && (y(i)>=0)
                xcrop=[xcrop; x(i)];
                ycrop=[ycrop; y(i)];
                devcrop=[devcrop; dev(i)];
            elseif ((x(i)-x0)/L3)^2+((y(i)-y0)/L4)^2<1 && (y(i)<0)
                xcrop=[xcrop; x(i)];
                ycrop=[ycrop; y(i)];
                devcrop=[devcrop; dev(i)];
            end
        end
    end
else %%HF
    if txtflg==1 % Left side
        [~,idx]=min(x);
        
        t=-pi/2:0.01:pi/2;
        L2=sqrt((x(idx)-x0)^2+(y(idx)-y0)^2);
        [L3,L4,L5,x0,L6,L7]=cropsize(L2,L1, size,type,txtflg,x0, L6,L7);
        for i=1:length(x)
            if x(i)<=x0
                xcrop=[xcrop; x(i)];
                ycrop=[ycrop; y(i)];
                devcrop=[devcrop; dev(i)];
            elseif ((x(i)-x0)/L3)^2+((y(i)-y0)/L5)^2<1 && (y(i)>=0)
                xcrop=[xcrop; x(i)];
                ycrop=[ycrop; y(i)];
                devcrop=[devcrop; dev(i)];
            elseif ((x(i)-x0)/L3)^2+((y(i)-y0)/L4)^2<1 && (y(i)<0) && ((x(i)-x0)/L6)^2+((y(i)-y0)/L7)^2<1
                xcrop=[xcrop; x(i)];
                ycrop=[ycrop; y(i)];
                devcrop=[devcrop; dev(i)];
            end
        end
        figure
        scatter(xcrop,ycrop,[],devcrop,'.')
        axis equal
    else
        [~,idx]=max(x);
        
        L2=sqrt((x(idx)-x0)^2+(y(idx)-y0)^2);
        [L3,L4,L5,x0,L6,L7]=cropsize(L2,L1, size,type,txtflg,x0,L6,L7);
        
        for i=1:length(x)
            if x(i)>=x0
                xcrop=[xcrop; x(i)];
                ycrop=[ycrop; y(i)];
                devcrop=[devcrop; dev(i)];
            elseif ((x(i)-x0)/L3)^2+((y(i)-y0)/L5)^2<1 && (y(i)>=0)
                xcrop=[xcrop; x(i)];
                ycrop=[ycrop; y(i)];
                devcrop=[devcrop; dev(i)];
            elseif  ((x(i)-x0)/L3)^2+((y(i)-y0)/L4)^2<1 && (y(i)<0) && ((x(i)-x0)/L6)^2+((y(i)-y0)/L7)^2<1
                xcrop=[xcrop; x(i)];
                ycrop=[ycrop; y(i)];
                devcrop=[devcrop; dev(i)];
            end
        end
        figure
        scatter(xcrop,ycrop,[],devcrop,'.')
        axis equal
    end
end
% if txtflg==2
%     xM1=x0+L3*cos(t(1:end*0.75));
%     xM2=x0+L6*cos(t(end*0.75:end));
% else
%     xM1=x0+L6*cos(t(1:end*0.25));
%     xM2=x0+L3*cos(t(end*0.25:end));
% end
% if txtflg==2
%     yM1=y0+L5*sin(t(1:end/2));
%     yM2=y0+L4*sin(t(end/2:end*0.75));
%     yM3=y0+L2*sin(t(end*0.75:end));
% else
%     yM1=y0+L4*sin(t(1:end/2));
%     yM2=y0+L5*sin(t(end/2:end));
%     yM3=y0+L2*sin(t(end*0.25:end));
%
%
% end
% xM=[xM1,xM2];
% yM=[yM1, yM2,yM3];
% hold on
% plot(xM,yM)

[k2,~]=boundary(xcrop,ycrop);

% plot(x(k),y(k))
% plot(xcrop(k2)+.5,ycrop(k2))
if txtflg==1
    in=inpolygon(xcrop,ycrop,xcrop(k2)+.5,ycrop(k2));
else
    in=inpolygon(xcrop,ycrop,xcrop(k2)-.5, ycrop(k2));
end
[xcrop,ycrop,devcrop]=crop(xcrop,ycrop,devcrop,in);
[k2,~]=boundary(xcrop,ycrop);
in=inpolygon(xcrop,ycrop,xcrop(k2),ycrop(k2)+.5);
[xcrop,ycrop,devcrop]=crop(xcrop,ycrop,devcrop,in);
[k2,~]=boundary(xcrop,ycrop);
in=inpolygon(xcrop,ycrop,xcrop(k2),ycrop(k2)-.5);
[xcrop,ycrop,devcrop]=crop(xcrop,ycrop,devcrop,in);
hold off
h=figure;
hold on
% plot(xM,yM)
plot(x(k),y(k));
scatter(xcrop,ycrop,[],devcrop,'.')
ylabel('y coordinates [mm]');
xlabel('x coordinates [mm]');
scatter(xcrop,ycrop,[],devcrop, '.');
axis(v)
if txtflg==1
    side='Left';
else
    side='Right';
end
saveas(h,[biomechid '_' side '.jpg'])
% keyboard
end
%%
function [xcrop,ycrop,devcrop]=crop(xcrop,ycrop,devcrop,in)
xcrop=xcrop(in);
ycrop=ycrop(in);
devcrop=devcrop(in);
end
function [L3,L4,L5,x0, L6,L7]=cropsize(L2,L1, size, type,txtflg,x0, L6,L7)
%type 1=high flex
%type 2=ps
if type==1
    if size==1
        if txtflg==1
            L3=L2*1.85;
            L5=L1*0.85;
            L4=L1*2.5;
            x0=x0*0.87;
            L6=L1*1.15;
            L7=L2*3;
        else
            %             L3=L2*1.6;
            %             L5=L1*0.815;
            %             L4=L1*0.925;
            %             x0=x0*0.8;
            %             L6=L1*1.5;
            L3=L2*1.85;
            L5=L1*0.85;
            L4=L1*2.5;
            x0=x0*0.9;
            L6=L1*1.05;
            L7=L2*2.4;
        end
        
    elseif size==3
        if txtflg==1
             L3=L2*1.4;
            L5=L1*0.74;
            L4=L1*2.5;
%             x0=x0;
            L6=L1*1.05;
            L7=L2*2.1;
%             L3=L2*1.4;
%             L5=L1*0.75;
%             L4=L1*1.025;
        else
            L3=L2*1.5;
            L5=L1*0.75;
            L4=L1*2.5;
            L6=L1*1.05;
            L7=L2*2.1;
        end
        
    elseif size==5
        if txtflg==1
            L3=L2*1.78;
            L5=L1*0.7;
%             x0=x0*0.9;
            L4=L1*2.5;
            L6=L1*1.05;
            L7=L2*2.5;
            
        else
            L3=L2*1.3;
            L5=L1*0.77;
            x0=x0*1.05;
            L4=L1*2.5;
            L6=L1*1.5;
            L7=L2*1.9;
            
        end
    else
        L3=L2;
        L5=L1*0.9;
    end
else
    if size==1
        if txtflg==1
            L3=L2/2.3;
            L5=L1*0.9;
        else
            L3=L2/2.4;
            L5=L1*0.9;
        end
        L4=L1*2.5;
    elseif size==3
        if txtflg==1
            L3=L2/2.6;
%             L4=L1*1.2;
            L5=L1*0.9;
        else
            L3=L2/2;
%             L4=L1*1.2;
            L5=L1*0.9;
        end
        L4=L1*2.5;
    elseif size==5
        if txtflg==1
            L3=L2/2.2;
            L5=L1*0.9;
%             L4=L1*1.2;
        else
            L3=L2/2.9;
            L5=L1*0.9;
%             L4=L1*1.2;
        end
        L4=L1*2.5;
        
        
    else
        if txtflg==1
            L3=L2/4;
            L5=L1*0.9;
%             L5=L1;
            L4=L1*1.5;
        else
            L3=L2/2.9;
            L3=L2/5;
            L5=L1*0.9;
            L4=L1*1.05;
            L4=L1*1.3;
        end
    end
end
end



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
