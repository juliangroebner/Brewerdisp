function analysedsp(day,year,brewnb,polynb,path,comment,O3xsec)
%function analysedsp(day,year,brewnb,polynb,path,comment,O3xsec)
% 30 8 99 julian
% combines all functions to make analysis straightforward
% including ozone calculation
% data has to be in d:\brewer\dsp\...
% 27 7 2008 jg adapt for standard lines
% 14 2 2019 JG use gaussian fitting for dispersion
% 18 3 2019 JG move o3-xsec to upper level 
% 25 6 2020 JG, add additional cross-sections and new parameter
%O3XSEC={'paur&bass','brion','Bremen','Bremen-ATMOZ','Bremen Teff','B&P IGC Teff','Bremen-ATMOZ Teff','Brion Teff','ACS'};
% 28 10 2020 JG add ACS x-sections from ESA project

if nargin<7,O3xsec=[];end
if isempty(O3xsec),O3xsec=1;end

if nargin<6,comment='';end
if nargin<5,path=[];end
if isempty(path),path=[''];end

if nargin<4,polynb=[];end
if isempty(polynb),polynb=3;end

FITGAUSS=0;  % 14 2 2019 JG use gauss fitting instead of standard trapezoidal.

%O3xsec=1;  
%O3XSEC={'paur&bass','brion','Bremen','Bremen-ATMOZ'};
O3XSEC={'paur&bass','brion','Bremen','Bremen-ATMOZ','Bremen-Teff','B&P-IGC-Teff','Bremen-ATMOZ-Teff','Brion-Teff','ACS'};

% first run and save alldsp information
[wl,dsp,dspstd,fwhm,fwhmstd,backlash]=alldsp(day,year,brewnb,[],[],[],[],path,FITGAUSS);

fnamealldsp=sprintf('alldsp_%03d%02d_%s.%03d',day(1),year,comment,brewnb);

disp(sprintf('saving alldsp to %s',[path fnamealldsp]));
save([path fnamealldsp],'wl','dsp','dspstd','fwhm','fwhmstd','backlash');


% now do dispersion fit using powfiu7

fname7=sprintf('dsp_%03d%02d_%s.%03d',day(1),year,comment,brewnb);
disp(sprintf('saving powfiu7 to %s',[path fname7]));

pos0=powfiu7;
if nanmean(fwhm(fwhm(:,1)~=0,1))<65,pos0=powfiu7(1);end

[fstd,slitpos,rmsf]=powfiu7(pos0,wl,dsp,polynb,[],[],[],[path fname7]);
disp('Residuals using powfiu7 [RMS]:');
rmsf

% check slitpos residuals.

dsl=slitpos(end,:)-pos0;
if any(abs(dsl)>0.05),  % do again with this one slit corrected
 ind=find(dsl==max(dsl));
 pos0(ind)=pos0(ind)+dsl(ind);; % shift slitpos and recalculate
 disp(sprintf('Too large slitpos deviation: Recalc with slit #%d shifted by %.3f',ind,dsl(ind)));
 [fstd,slitpos,rmsf]=powfiu7(pos0,wl,dsp,polynb,[],[],[],[path fname7]);
 disp('Residuals using powfiu7 [RMS]:');
 rmsf
end


ind=fstd==0;fstd(ind)=nan;
%figure;plot(wl,fstd(:,:,end)/10,'o-');
figure;plot(wl/10,fstd(:,:,end)/10,'o-');
title(['Powfiu7 RMS=' sprintf('%.4f',(min(rmsf/10))) ' ' fname7 '-' datestr(now)]);
ylabel('Residuals [nm]');
xlabel('wavelength [nm]');
legend('slit0','1','2','3','4','5');

disp('Use brstps2 to calculate steps and wavelengths');


% now calculate normaldsp

[fwl,fstps,pwl,pstps]=normaldsp(wl,dsp);  % quad for every slit
fnameN=sprintf('dspnorm_%03d%02d_%s.%03d',day(1),year,comment,brewnb);

figure;plot(wl,fwl/10,'o-');
title(['Normaldsp ' path fnameN '-' datestr(now)]);
ylabel('Residuals [nm]');
xlabel('wavelength [nm]');
legend('sl0','1','2','3','4','5');

disp(sprintf('saving normaldsp to %s as brewer compatible file',[path fnameN]));

%save([path fnameN],'pwl','pstps');
data=pwl(:,end:-1:1)';
mydata=data(:);

savefmt([path fnameN],[mydata(4:end);mydata(1:3)],'','%.7e');

disp('Use polyval(pwl(2,:),wl) for calculating normal wavelengths')

% And now show difference between both methods for slits 1 and 5.


testwl1=3000:10:3500;
testwl5=3500:10:3650;
stps17=brstps2(testwl1,1,[],[path fname7]);  % reference steps to use
stps57=brstps2(testwl5,5,[],[path fname7]);

quad1=polyval(pwl(2,:),stps17);
quad5=polyval(pwl(6,:),stps57);

figure;plot(testwl1,(quad1-testwl1)/10,testwl5,(quad5-testwl5)/10);
title('Normaldsp versus powfiu7 using slit 1 (3000:3500) and 5 (3500:3650)')
ylabel('normaldsp-powfiu7 [nm]');


% now calculate ozonecoeffs

outfname=sprintf('opos%03d%02d_%s.%03d',day(1),year,comment,brewnb);

disp(sprintf('Saving ozonecoeffs to %s',[path outfname]));

ozonecoeff2([path fnamealldsp],[],[path fname7],[path outfname],O3xsec,FITGAUSS);

disp('Finished with analysing data');

