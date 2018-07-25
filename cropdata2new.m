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
