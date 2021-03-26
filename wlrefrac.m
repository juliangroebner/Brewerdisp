function [wlair,dwl]=wlrefrac(wl,p,t)
%function [wlair,dwl]=wlrefrac(wl,p,t)
% 12 5 98 julian
% calculates wavelength shift in nm for P and T
% dwl = wlvac-wlair
% p should be in mbar
% wl is wl vacuum

if nargin<3,t=15;end
if nargin<2,p=[];end % in Pascal

if isempty(p),p=1013.25;end % p in mbar

p=p*100 ;%p in pascal.

corrfac=p*(1+p*(61.3-t)*1e-10)/(96095.4*(1+0.003661*t));

sigma=1./(wl*1e-3); % change to micrometers

nn=8342.13+2406030./(130-sigma.^2)+15997./(38.9-sigma.^2);

nn=nn*corrfac;

wlair=wl./(nn*1e-8+1);

dwl=wl-wlair;

return

% 1 6 2016 JG use the Birch and Downs version for refractive index of air
% from http://emtoolbox.nist.gov/Wavelength/Documentation.asp#EdlenorCiddor

% Humidity:
h=rh/100;
psv=partialvaporpressure(t);
pv=h.*psv;  % partial pressure

% now edlen equation
A=8342.54;
B=2406147;
C=15998;
D=96095.43;
E=0.601;
F=0.00972;
G=0.003661;
wl=wl/1e3;   % convert to micrometers
S=1./wl.^2;
 n5 = 1 + 10^(-8) .* (A + B ./ (130 - S) + C ./ (38.9 - S));
 X=(1+10-8*(E-F.*t).*p)./(1+G.*t);
 ntp=1+p*(n5-1)*X/D;
 n=ntp-10^-10*(292.75)/(t+273.15)*(3.7345-0.0401*S)*p;
 
 wlair=wlvac./n;
 

function psv=partialvaporpressure(t)   % saturation vapor pressure
K=[1.16705214528E+03 -7.24213167032E+05 -1.70738469401E+01 1.20208247025E+04 -3.23255503223E+06...
    1.49151086135E+01 -4.82326573616E+03 4.05113405421E+05 -2.38555575678E-01 6.50175348448E+02];

T=t+273.15;
omega=T+K(9)/(T-K(10));
A=omega^2+K(1)*omega+K(2);
B=K(3)*omega^2+K(4)*omega+K(5);
C=K(6)*omega^2+K(7)*omega+K(8);
X=-B+sqrt(B^2+4*A*C);
psv=1e6*(2*C/X)^4;

