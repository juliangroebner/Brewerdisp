function [wl,fwhm]=slitfit(slit,type,pl,fignb)
%function [wl,fwhm]=slitfit(slit,type,pl)
% 2 11 2011 JG gauss fit to slit function
% type=0, default, fit with isoceles triangle
% type=1, gaussfit (array SRM)
%type=[2 N], gaussN fit, N=1 to 8, result is centroid and FWHM
% type=3, centroid for centerwl, 

% slit is two columns,
% 31 7 2013 add young-ohno supergaussian fit
% 20 2 2014 JG, make gaussfit =1.177;
% 18 8 2015 JG, add centroid for center of wl, fwhm is 1 sigma (68.27%)
% 30 6 2016 JG make gauss2 fit work...
% 15 12 2016 JG optimise centroid fit (assumes many points...);
% 23 1 2017 JG error in centroid calc
% 22 6 2020 JG error in calculating the FWHM from gaussfit
% 26 3 2021 JG remove nans from slit function

fwhm=nan;

if nargin<4,fignb=[];end
if isempty(fignb),fignb=0;end
        
if nargin<3,pl=[];end
if isempty(pl),pl=0;end
if nargin<2,type=[];end
if isempty(type),type=0;end

b=isnan(slit(:,2));
slit(b,:)=[];

slit(:,2)=slit(:,2)/max(slit(:,2));  % normalise to 1
switch type(1),
    case 0, % standatd isoceless fit
        [wl,fwhm]=isofit(slit,pl);
    case 1, % do gaussfit
        [wl,fwhm]=gaussfit(slit,pl,fignb);
    case 2, % gaussN fit
        if length(type)<2,type(2)=2;end  % default is gauss2
        [wl,fwhm]=gaussn(slit,type(2),pl);
    case 3, [wl,fwhm]=centroidfit(slit);

end

function [wl,fwhm]=centroidfit(slit)
ind=logical(ones(size(slit(:,1))));
bg=0;
for j=1:2,
    myslit=slit(ind,:);myslit(:,2)=myslit(:,2)-bg;
cs=cumsum(myslit(:,2))./sum(myslit(:,2));
[buf,i]=unique(cs);
if j==1,
    wl=myslit(myslit(:,2)==max(myslit(:,2)),1);
else
wl=interp1(cs(i),myslit(i,1),0.5);
end
[buf,i]=unique(myslit(:,2));
myslit2=myslit(i,:);
[buf,i]=sort(myslit2(:,1));myslit2=myslit2(i,:);
ind=myslit2(:,1)>wl;
%fwhm=spline(myslit2(ind,2),myslit2(ind,1)/max(myslit2(:,2)),0.5)-spline(myslit2(~ind,2),myslit2(~ind,1)/max(myslit2(:,2)),0.5);
fwhm=interp1(myslit2(ind,2)/max(myslit2(:,2)),myslit2(ind,1),0.5)-interp1(myslit2(~ind,2)/max(myslit2(:,2)),myslit2(~ind,1),0.5);
ind=abs(myslit(:,1)-wl)<(2*fwhm);
bg=mean(myslit(~ind,2));
if isnan(bg),bg=0;end
end

function [center,fwhm]=gaussfit(slit,pl,fignb)

xdata=slit(:,1);
ydata=slit(:,2);
halfwidth = (max(xdata)-min(xdata))/2;
intensity = 1;
background = 0;

if 0,
gauss = @(x,xdata) x(4) + x(3)./(sqrt(2*pi)*x(2)).*exp(-(xdata-x(1)).^2./(2*x(2).^2));
x0=[mean(xdata), halfwidth, 1, background];

out = lsqcurvefit(gauss,x0,xdata,ydata, [xdata(1), 0.01, 0.1, -1000], [xdata(end), 100, 10, 1000],optimset('Display','off'));

center=out(1);
if 0,   % geht nicht, da es oft auch aussreisser in den Daten gibt.
f=@(x,c) abs(normpdf(x,c(1),c(2))/normpdf(c(1),c(1),c(2))-0.5);
fwhm=fminbnd(@(x)f(x,out(1:2)),0,out(2)*2)*2;  % is about 17% larger than MU
end
fwhm=out(2)*2*sqrt(2*log(2));   % 20 2 2014 JG
if pl,   % do plot
    xi=min(xdata):mean(diff(xdata))/10:max(xdata);
    if ~fignb,
        figure;
    elseif strcmp(get(fignb,'type'),'figure'),
        figure(fignb);hold on;
    elseif strcmp(get(fignb,'type'),'axes'),
        axes(fignb);
    end
      buf=gauss(out,xi);
       plot(xdata,ydata/max(buf),'.',xi,gauss(out,xi)/max(buf));grid;
end
else
 [center,fwhm]=gaussn(slit,[ 1],pl);   
end

function [center,fwhm]=isofit(slit,pl);
% 2 11 2011 JG use dsp
x=slit(:,1);
y=slit(:,2);
%indmax=find(y==max(y));
%ind=abs(x-x(indmax))<10;
%x=x(ind);
%y=y(ind);


x=[x;x(end:-1:1)];
y=[y;y(end:-1:1)];
[wl1,wl2,fwhm]=dsp([x y],'',[],[],pl);

center=wl1;

function [center,fwhm]=gaussn(slit,n,pl)
fun=sprintf('gauss%d',n);

opt=fitoptions(fun);

opt.Lower=repmat([0 -10 0],1,n);
opt.Display='off'; 
f=fit(slit(:,1), slit(:,2),fun,opt);
 gg=@(x) f(x)-0.001;
 opt=optimoptions('fsolve','Display','off');
 x1=fsolve(@(x) gg(x),slit(1,1),opt);
 x2=fsolve(@(x) gg(x),slit(end,1),opt);
dx=(x2-x1)/1e3;
x=x1:dx:x2;
y=f(x);

if n==1,
    center=f.b1;
    fwhm=f.c1*2*sqrt(log(2)); % 22 6 2020 JG is now correct conversion to FWHM
else
cs=cumsum(y)./sum(y);
center=interp1(cs,x,0.5);
fwhm=(interp1(cs,x,97.5/100)-interp1(cs,x,0.025))/2;  % 1 sigma now  95% coverage divided  by 2
end
if pl,
    x=min(slit(:,1)):mean(diff(slit(:,1)))/10:max(slit(:,1));
    plot(x,f(x),slit(:,1),slit(:,2),'.');
end

