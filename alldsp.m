function [wl,DSP,DSPstd,fwhm,fwhmstd,backlash,resup]=alldsp(day,year,brew,lines,minslit,maxslit,usewl,dspp,fitgauss)

% function [wl,DSP,DSPstd,fwhm,fwhmstd,backlash]=alldsp(day,year,brew,lines,minslit,maxslit,usewl,dspp,fitgauss);
% 3 dec 96 julian
% calculates all dsp files available.
% sorts by slits.
% 4 12 96 correct for deadtime; -> calculate counts per sec.
% day can be several days,
% 19 2 97 julian
% 25 4 97 julian. adapt for new format
% 2 5 97 julian Add filter patch for days 114,121,122
% 4 6 97 julian. if usewl exists, then calculate wl from brstps.
% 1 10 97 julian change saving to make real mean of all days per line
% 21 10 97 julian change dsp format back to original. 4 cycles
% 22 10 97 julian backlash=up-down. positive means backlash!
% 14 11 97 julian remove background from lines.-> changes the slitwidth.
% 30 8 99 julian improve input handling
% 11 7 2001 julian add 163 filename format
% 27 7 2008 jg adapt for standard lines...
% 14 2 2019 JG add gauss fitting option
% 25 11 2019 JG change reading of files using dir and reading all available files


% this are filter function ND1,ND2,ND3 around wl 3200 Angstrom
pp(1,:)=[1.97227e-011  -4.55518e-008   1.17860e-005   0.36615075];
pp(2,:)=[1.97552e-011  -2.77914e-008   1.89161e-006   0.10661427];
pp(3,:)=[7.00540e-012  -1.04491e-008   2.12749e-006   0.02666995];


if nargin<3 brew=[];end
if isempty(brew),brew=163;end
if nargin<4,lines=[];end
if isempty(lines),lines=1:9;end
if nargin<5,minslit=[];end
if isempty(minslit),minslit=0;end
if nargin<6,maxslit=[];end
if isempty(maxslit),maxslit=5;end

if nargin<7,usewl=[];end

if nargin<8,dspp=[];end
if isempty(dspp),
    dspp=[''];
end
if nargin<8,fitgauss=[];end
if isempty(fitgauss),fitgauss=0;end

if ~ischar(brew),
    brew=sprintf('.%03d',brew);
end


maxlines=max(lines);

% dsp_arrU=zeros(maxlines,6);   % contains data for Up scan
% dsp_arrD=zeros(maxlines,6);   % for DOWN scan
% DSP=zeros(maxlines,6);
% wl=zeros(maxlines,1);
% DSPstd=zeros(maxlines,6);
% fwhm=zeros(maxlines,6);
% fwhmstd=zeros(maxlines,6);
% backlash=zeros(maxlines,6);
% resup=zeros(maxlines,6);

%if year<50, year=year+90;end
if year>=2000,y2=year-2000;%year=y2+100;
elseif year>1900,y2=year-1900;year=y2;end

d=dir([dspp 'W*' brew]);
linecnt=0;
for i=1:length(d),
    fn=d(i).name;
    ind=findstr(fn,'.');
    slit=sscanf(fn(ind-6),'%d');
    days=sscanf(fn(ind-5:ind-3),'%d');
    filename=[dspp fn];
    lnes=sscanf(fn(2:ind-7),'%d');
    
    daystr=sprintf('%03d',days);
    
    [a,s]=liesfile( filename,1,2);
    
    disp(['now:' filename]);
    if (rem(linecnt,8)==0),figure;linecnt=0;end
    if linecnt==-1,linecnt=0;end
    linecnt=linecnt+1;
    wlb=str2num(s{1});
    lnes=fix(wlb);
    wl(lnes)=wlb;
    %if rem(linecnt,8)==0,figure;figcnt=1;linecnt=0;end
    subplot(2,4,linecnt);
    if ~isempty(usewl),a(:,1)=brstps2(a(:,1),slits,1,usewl);end % calculate wl usewl is filename
    
    
    a(:,2)=a(:,2)-min(a(:,2)); % remove background 9 3 98 wrong! put into dsp 30 3 98 ok_new
    if ~fitgauss,
        [dsp_arrU(lnes,slit+1,days),dsp_arrD(lnes,slit+1,days),fwhmbuf(lnes,slit+1,days),resup(lnes,slit+1,days)]=dsp(a,s,0.2,0.8,1);
    else
        nn=size(a,1)/2;
        buf=prctile(a(:,2),5);bkd=mean(a(a(:,2)<=buf,2));
        if 1, % use only 20-80% fitting similar to brewer iso fit
            cutmin=0.1; cutmax=0.9;
            indmax=find(a(:,1)==max(a(:,1)));indmax=indmax(1);
          up=a(1:indmax,2);
          down=a(size(a,1):-1:(indmax+1),2);
          steps=a(1:indmax,1);  % get steps
          cnts=[up down];
          [maxcnts]=max(cnts);  % get max value and position.
          ind1=find((cnts(:,1)<(cutmax*maxcnts(1))) & (cnts(:,1)>(cutmin*maxcnts(1))));
          ind2=find((cnts(:,2)<(cutmax*maxcnts(2))) & (cnts(:,2)>(cutmin*maxcnts(2))));
          bufa=[a(ind1,:)];
          bufb=[a(ind2+indmax,:)];
        else  % default
            bufa=a(1:nn,:)-bkd;
            bufb=a(end:-1:nn+1,:)-bkd;
        end
        
        [dsp_arrU(lnes,slit+1,days),fwhm1]=slitfit(bufa,[1],1,gcf);
        [dsp_arrD(lnes,slit+1,days),fwhm2]=slitfit(bufb(end:-1:1,:),[1],1,gcf);
        fwhmbuf(lnes,slit+1,days)=(fwhm1+fwhm2)/2;
    end
    xlabel(filename);
    title(sprintf(['%5.3f [' char(197) ']'],wl(lnes)),'fontsize',9);
    
end  % days
dsp_arrU(dsp_arrU==0)=nan;
dsp_arrD(dsp_arrD==0)=nan;
fwhmbuf(fwhmbuf==0)=nan;

fwhm=nanmean(fwhmbuf,3);
    fwhmstd=nanstd(fwhmbuf,[],3);
   DSP=nanmean((dsp_arrU+dsp_arrD)/2,3);
   DSPstd=nanstd((dsp_arrU+dsp_arrD)/2,[],3);
   backlash=nanmean((dsp_arrU-dsp_arrD),3);


ind=(wl==0);
wl(ind)=[];
DSP(ind,:)=[];
DSPstd(ind,:)=[];
fwhm(ind,:)=[];
fwhmstd(ind,:)=[];
backlash(ind,:)=[];

[wl,i]=sort(wl);
DSP=DSP(i,:);
DSPstd=DSPstd(i,:);
fwhm=fwhm(i,:);
fwhmstd=fwhmstd(i,:);
backlash=backlash(i,:);
wl=wl(:);
