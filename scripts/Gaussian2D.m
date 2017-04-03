function [cx,cy,sx,sy,PeakOD] = Gaussian2D(m,tol,center)
% m = image
% tol = fitting tolerance

options = optimset('Display','off','TolFun',tol,'LargeScale','off');

[sizey sizex] = size(m);
[cx,cy,sx,sy] = centerofmass(m);
cx = center(1); cy = center(2);

pOD = max(max(m));

mx = m(round(cy),:);
x1D = 1:sizex;
ip1D = [cx,sx,pOD];
fp1D = fminunc(@fitGaussian1D,ip1D,options,mx,x1D);

cx = fp1D(1);
sx = fp1D(2);
PeakOD = fp1D(3);

my = m(:,round(cx))';
y1D = 1:sizey;
ip1D = [cy,sy,pOD];
fp1D = fminunc(@fitGaussian1D,ip1D,options,my,y1D);

cy = fp1D(1);
sy = fp1D(2);
PeakOD = fp1D(3);
[X,Y] = meshgrid(1:sizex,1:sizey);

initpar = [cx,cy,sx,sy,PeakOD];
fp = fminunc(@fitGaussian2D,initpar,options,m,X,Y);
cx = fp(1);
cy = fp(2);
sx = fp(3);
sy = fp(4);
PeakOD = fp(5);


function [cx,cy,sx,sy] = centerofmass(m);
% PURPOSE: find c of m of distribution
[sizey sizex] = size(m);
vx = sum(m);
vy = sum(m');

vx = vx.*(vx>0);
vy = vy.*(vy>0);

x = [1:sizex];
y = [1:sizey];

cx = sum(vx.*x)/sum(vx);
cy = sum(vy.*y)/sum(vy);

sx = sqrt(sum(vx.*(abs(x-cx).^2))/sum(vx));
sy = sqrt(sum(vy.*(abs(y-cy).^2))/sum(vy));

function [z] = fitGaussian1D(p,v,x)

%cx = p(1);
%wx = p(2);
%amp = p(3);

zx = p(3)*exp(-0.5*(x-p(1)).^2./(p(2)^2)) - v;

z = sum(zx.^2);


function [z] = fitGaussian2D(p,m,X,Y)

%cx = p(1);
%cy = p(2);
%wx = p(3);
%wy = p(4);
%amp = p(5);

ztmp = p(5)*(exp(-0.5*(X-p(1)).^2./(p(3)^2)-0.5*(Y-p(2)).^2./(p(4)^2))) - m;

z = sum(sum(ztmp.^2));