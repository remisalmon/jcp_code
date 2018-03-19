clear all;
%close all;
clc;

load matrice.out;
data=matrice;
clear matrice;

nx=size(data,1)/2; % x size = i
ny=size(data,2); % y size = j
R=1;

al=0:((2*pi)/6):2*pi;
xh=R*cos(al);
yh=R*sin(al);

x0=0;
y0=0;
x=x0;
y=y0;

x0p=3*R/2;
y0p=sqrt(3)*R/2;
xp=x0p;
yp=y0p;

figure;
hold on
for nxc=10:nx-10
    for nyc=10:ny-10
        val = data(nxc*2-1,nyc);
        if(val == 1)
            fill(y+yh,x+xh,'red');%,'LineStyle','none');
        end
        if(val == 0)
            fill(y+yh,x+xh,'white');%,'LineStyle','none');
        end
        
        y=y+sqrt(3)*R;
        val = data(nxc*2,nyc);
        if(val == 1)
            fill(yp+yh,xp+xh,'red');%,'LineStyle','none');
        end
        if(val == 0)
            fill(yp+yh,xp+xh,'white');%,'LineStyle','none');
        end
        
        yp=yp+sqrt(3)*R;
    end
    x=x+3*R;
    y=y0;

    xp=xp+3*R;
    yp=y0p;
end
axis equal;
axis off;
set(gca,'YDir','reverse');
hold off
clear all;
