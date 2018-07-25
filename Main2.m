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



