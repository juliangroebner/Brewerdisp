function [fwl,fstps,pwl,pstps]=normaldsp(wl,A)
% function [fwl,fstps,pwl,pstps]=normaldsp(wl,A)
%22 1 98 julian
% 30 8 99 julian modify input variables


%wl=A(:,1);
fwl=zeros(length(wl),6);
fstps=fwl;
pwl=[];
pstps=[];

for i=1:6,%2:7,
   ind=A(:,i)~=0 & ~isnan(A(:,i));
   if sum(ind)>2,
    p=polyfit(A(ind,i),wl(ind),2);
    p2=polyfit(wl(ind),A(ind,i),2);
    fwl(ind,i)=polyval(p,A(ind,i))-wl(ind);  % i-1 bei fwl
    fstps(ind,i)=polyval(p2,wl(ind))-A(ind,i);   % i-1 bei fstps
    pwl=[pwl;p];
    pstps=[pstps;p2];
 else
    fwl(:,i)=nan;  % hier auch i-1
    fstps(:,i)=nan;  % i-1
    pwl=[pwl ;[nan nan nan]];
    pstps=[pstps ;[nan nan nan]];   
  end  
end
