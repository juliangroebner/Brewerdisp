function [o3abs]=ozonecoeff2(fname,ozonepos,dcfname,outfname,O3xsec,fitgauss);
%function [o3abs]=ozonecoeff2(fname,ozonepos,dcfname,outfname,O3xsec,fitgauss);
% 4 2 98 julian
% ozonepos is total steps from uvzero.
% dcfname is dcfile name for brstps.
% fname is data from alldsp
% dcfname is file ontained from savedsp
% outfname is were log should be written
% 23 11 2012 add o3xsec from Bremen_2012
% 30 4 2015 JG add xsection name to outfilename
% 14 2 2019 JG add gauss fitting for slit function
% 27 5 2020 JG add Teff for Xsec calculation and output of o3abs.
% 28 10 2020 JG add x-sections from ESA project
% 29 10 2020, use BDM with linear t fit 
% 11 11 2020 JG make bremen average of 233 and 223 K

%O3=2;   % default
O3XSEC={'paur&bass','brion','Bremen','Bremen-ATMOZ','Bremen Teff','B&P IGC Teff','Bremen-ATMOZ Teff','Brion Teff','ACS'};

if nargin<6,fitgauss=[];end
if isempty(fitgauss),fitgauss=0;end
if nargin<5,O3xsec=[];end
if isempty(O3xsec),O3xsec=1;end
        
if nargin<3,dcfname='';end   % never happens, otherwise error later on
if nargin<2,ozonepos=[];end
if isempty(dcfname),error('Need filename for dispersion "dcfname"');end

if isempty(ozonepos),ozonepos=brstps2(3063.0,1,[],dcfname);end   

if nargin<4,outfname=[];end

if ~isempty(outfname),
    [ppi,fn,ext]=fileparts(outfname);
    if isempty(ppi),ppi=pwd;end
    outfname=[ppi '/' fn  datestr(now,'YYYYmmDD') ext];
    disp(sprintf('Saving ozonecoeffs to %s',[outfname]));
    fid=fopen(outfname,'at');
    fprintf(fid,'Analysis on %s using files:%s,%s\n',date,fname,dcfname);
    if fitgauss,
        fprintf(fid,'Fitting using gaussian slit function\n');
    else
        fprintf(fid,'Fitting using default trapezoid slit function\n');
    end
    fprintf(fid,'Ozone X-section: %s\n',O3XSEC{O3xsec});
    fprintf(fid,'%s              %8s %8s %8s %8s %8s %8s\n','step','slit#0','slit#1','slit#2','slit#3','slit#4','slit#5');
else
    fid=-1;
end

   
eval(['load ' fname ' -mat']); % get data from powfiu6;
if ~exist('wl'),wl=wlstd;end


switch O3xsec(1),
    
    case 1,% here load paur-bass
        load ozxsec2.dat
        f=ozxsec2;
        o3xsec=f(:,[1 4]);
    case 2,% here load brion.
        load brion
        o3xsec=brion(:,1)/10;
        o3xsec(:,2)=brion(:,2);
        
    case 3,
        load bremen;
        o3xsec=wlair;
        o3xsec(:,2)=mean(bremen(:,7:8),2); %-40 degC, -50°C   11 11 2020, no change in o3abs for br 163 (to the last digit...)
       % o3xsec(:,2)=bremen(:,8); %-40 degC, -50°C
        
    case 4,  % 29 9 2017 JG new ATMOZ x-sec from bremen
        load bremen_atmoz;
        o3xsec=wlair;
        o3xsec(:,2)=mean(xsec(:,7:8),2); %-40 degC -45°C
        
    case 5,   % 27 5 2020 JG, here do Teff for Bremen Xsections
        if length(O3xsec)<2,
            O3xsec(2)=-45;
        end
        Teff=O3xsec(2);
        load Serdyuchenko_O3_temp.dat
        s=Serdyuchenko_O3_temp;
        o3xsec=vac2air(s(:,1));
        o3xsec(:,2)=(s(:,2).*(1+s(:,3).*Teff+s(:,4)*Teff.^2))*1e-20;
        O3XSEC{5}=sprintf('%s=%.1f C',O3XSEC{5},Teff);
    case 6, % B&P Igaco x-secs
        if length(O3xsec)<2,
            O3xsec(2)=-45;
        end
        Teff=O3xsec(2);
        load bp.par
        o3xsec=bp(:,1);
        o3xsec(:,2)=(bp(:,2)+bp(:,3).*Teff+bp(:,4).*Teff.^2)*1e-20;
        O3XSEC{6}=sprintf('%s=%.1f C',O3XSEC{6},Teff);
       
    case 7,   % Bremen ATMOZ xsections
        if length(O3xsec)<2,
            O3xsec(2)=-45;
        end
        Teff=O3xsec(2);
        load gorshelev_11T_201709_v3_temp.dat
        s=gorshelev_11T_201709_v3_temp;
        o3xsec=vac2air(s(:,1));
        o3xsec(:,2)=(s(:,2)+s(:,3).*Teff+s(:,4)*Teff.^2)*1e-20;
        O3XSEC{7}=sprintf('%s=%.1f C',O3XSEC{7},Teff);
        
     case 8,   % Brion x-sections
        if length(O3xsec)<2,
            O3xsec(2)=-45;
        end
        Teff=O3xsec(2);
        br=load('DBM'); 
        o3xsec=br.wlair';
        s=br.p1;
        o3xsec(:,2)=(s(:,2)+s(:,1).*Teff);
        O3XSEC{8}=sprintf('%s=%.1f C',O3XSEC{8},Teff);

     case 9,   % ACS ESA x-sections
        if length(O3xsec)<2,
            O3xsec(2)=-45;
        end
        Teff=O3xsec(2);
        br=load('acs'); 
        o3xsec=br.wlair;
        s=br.xsecquad;
        o3xsec(:,2)=(s(:,3)+s(:,2).*Teff+s(:,1)*Teff.^2);
        O3XSEC{9}=sprintf('%s=%.1f C',O3XSEC{9},Teff);

end
%fprintf(fid,'Using %s xsections',O3XSEC{O3xsec});

% load so2;
%load g:\brewer\xsections\so2_jim\so2295sm.dat;  % calculated by jim.
load so2295sm.dat;  % calculated by jim.
fso2(:,1)=so2295sm(:,6)*10;
fso2(:,2)=so2295sm(:,7);

k=1.38062e-23; % bolzmann
o3x1T=o3xsec(:,2)*1.013*1e5/(k*273.1)*1e-6; %in %1/cm
o3x1wl=o3xsec(:,1)*10;

ind=(wl>2950 & wl<3500);
wl=wl(ind);  % 21 2 98 julian new.
fwhm=fwhm(ind,:);
fwhm(isnan(fwhm))=0;
%absrayl=nicolet(slitwl/10,P); % rayl for natural log

%if nargin<2,ozonepos=977+3700;end % for br#119.

ozonepos=round(ozonepos);
if length(ozonepos)==1,oposrange=ozonepos-5:ozonepos+5;
else
   oposrange=ozonepos;
end

if 0
% 13 7 2001 Julian Use atlas3 to calculate ET at ozone wavelengths
at3=load('atlas3_br.dat'); % wl in vacuum (correct by about 0.1nm)
at3(:,1)=wlrefrac(at3(:,1))*10;  % change to air
else  % load QASUMEFTS, 18 3 2019 JG
    et=load('QASUMEFTS');
    at3(:,1)=et.QASUMEFTS(:,1)*10;  % change to angstroem
    at3(:,2)=et.QASUMEFTS(:,2)*1e3;
end

cnt=0;
for opos=oposrange;
   
 for i=1:6,  % slits 0 to 5
    ind=fwhm(:,i)~=0;  % use only without zero!!!!
    if sum(ind)>1,
     thiswl(i)=brstps2(opos,i-1,1,dcfname); % wl for ozone position
     wlstep(i)=diff(brstps2(opos-1:opos,i-1,1,dcfname)); % one step is ... wl at this position
     pp(i,:)=polyfit(wl(ind),fwhm(ind,i),1); % fit fwhm in step space
     fwhmwl(i)=polyval(pp(i,:),thiswl(i))*wlstep(i); % transform to wl
     
     if ~fitgauss,
     buf=brslit(fwhmwl(i),o3x1wl-thiswl(i)); % calculate slitfunction
     else
         mu=thiswl(i);
         sigma=fwhmwl(i)/2/sqrt(2*log(2));
         buf=o3x1wl;
         buf(:,2)=normpdf(buf(:,1),mu,sigma);
     end
     o3coeff(i)=sum(buf(:,2).*o3x1T)/sum(buf(:,2)); % integrate
     o3coeff(i)=o3coeff(i)/log(10); % need in log10 units, as in brewer
  %figure;semilogy(buf(:,1)+thiswl(i),buf(:,2),o3x1wl,o3x1T/max(o3x1T));
  %disp(sprintf('%10.0f steps:slit %d has FW=%10.4f  A at wl=%15.5f A and o3coeff of %15.5f',opos,i,fwhmwl(i),thiswl(i),o3coeff(i)));
  % now rayleigh:
     %wlx=-100:0.1:100; % in nm % 17 4 98 angstroem
     %slitray=brslit(fwhmwl(i),wlx);
     %raycoeff(i)=sum(slitray(:,2).*nicolet(wlx+thiswl(i)/10,1013.25)')/sum(slitray(:,2));
%     raycoeff(i)=nicolet(thiswl(i)/10);
     raycoeff(i)=bodhaine(thiswl(i)/10);
     
     raycoeff(i)=raycoeff(i)/log(10);
     %disp(sprintf('%10.0f steps:slit %d has rayleigh in log10 units of %15.5f',opos,i,raycoeff(i)));
     % calculate so2
     if ~fitgauss,     
     buf=brslit(fwhmwl(i),fso2(:,1)-thiswl(i));
     else
         mu=thiswl(i);
         sigma=fwhmwl(i)/2/sqrt(2*log(2));
         buf=fso2(:,1);
         buf(:,2)=normpdf(buf(:,1),mu,sigma);
     end
     so2coeff(i)=sum(buf(:,2).*fso2(:,2))/sum(buf(:,2));
     so2coeff(i)=so2coeff(i)/log(10);

% 13 7 2001 julian now I0 from atlas3
     if ~fitgauss,     
     buf=brslit(fwhmwl(i),at3(:,1)-thiswl(i));
     else
         mu=thiswl(i);
         sigma=fwhmwl(i)/2/sqrt(2*log(2));
         buf=at3(:,1);
         buf(:,2)=normpdf(buf(:,1),mu,sigma);
     end
     I0(i)=sum(buf(:,2).*at3(:,2))/sum(buf(:,2));    
  else
     o3coeff(i)=nan;
     raycoeff(i)=nan;
     thiswl(i)=nan;
     fwhm(i)=nan;
  end   
end  

O3Coeff=(o3coeff(2+1))-0.5*(o3coeff(3+1))-2.2*(o3coeff(4+1))+1.7*(o3coeff(5+1));
RAYCoeff=(raycoeff(2+1)-0.5*raycoeff(3+1)-2.2*raycoeff(4+1)+1.7*raycoeff(5+1));
So2coeff=(so2coeff(1+1)-4.2*so2coeff(4+1)+3.2*so2coeff(5+1));
O34So2cc=o3coeff(1+1)-4.2*o3coeff(4+1)+3.2*o3coeff(5+1);
%I0Coeff=((log10(I0(2+1)))-0.5*(log10(I0(3+1)))-2.2*(log10(I0(4+1)))+1.7*(log10(I0(5+1))))*1e4;

s1=sprintf('%5.0f WL(A)        %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f',opos,thiswl);
s2=sprintf('      Res.         %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f',fwhmwl*2);
s3=sprintf('      O3abs(1/cm)  %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f O3:%8.4f',o3coeff,O3Coeff);
s3b=sprintf('     So2abs(1/cm)  %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f',so2coeff);    
s4=sprintf('  1e4*Rayabs(1/cm) %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f R:%8.4f',raycoeff*1e4,RAYCoeff);
s4b=sprintf(' I0(mW m^-2nm^-1)  %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f',I0);
s5=sprintf('Ozone offset due to Rayleigh: %3.1f DU',-RAYCoeff./O3Coeff*1e3);
s6=sprintf('Ratio Ozone for So2(A3)=%8.4f, So2/O3(A2)=%8.4f', O34So2cc,So2coeff/O34So2cc);
disp(s1)
disp(s2)
disp(s3)
disp(s3b)
disp(s4)
disp(s4b)
disp(s5)
disp(s6)

if fid>0,  % print to file
   fprintf(fid,'%s\n',s1);
   fprintf(fid,'%s\n',s2);
   fprintf(fid,'%s\n',s3);
   fprintf(fid,'%s\n',s3b);
   fprintf(fid,'%s\n',s4);
   fprintf(fid,'%s\n',s4b);
   fprintf(fid,'%s\n',s5);
   fprintf(fid,'%s\n\n',s6);
end

%disp(sprintf('%10.0f steps has ozonecoeff:%10.4f',opos,(o3coeff(2))-0.5*(o3coeff(3))-2.2*(o3coeff(4))+1.7*(o3coeff(5))));
%disp(sprintf('%10.0f steps has raycoeff:%10.4f',opos,(raycoeff(2)-0.5*raycoeff(3)-2.2*raycoeff(4)+1.7*raycoeff(5))));
cnt=cnt+1;o3abs(cnt)=O3Coeff;
end

if fid>0,fclose(fid);end

