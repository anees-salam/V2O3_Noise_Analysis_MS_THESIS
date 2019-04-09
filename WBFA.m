% test program average and forward + backward
% Modified on July-24 (Forward window propagation)

clc;
clear;
close all;

qmin=-10;
qstp=1;
qmax=10;
waves='db8'; ii=8;
levels=8;

poo=[];

file = load ('149.txt');

name12 = file(:,2); 


b=name12(1:24000);

%name = input('Enter file name with path: ','s');

%b=load(name);

%------------------------------------------------------%
[a,f1,b1]=wpadB(b); % constant padding

temp=b; m=length(b);
k=1;
for jj=m:-1:1
    b(k)=temp(jj);
    k=k+1;
end

[d2,f2,b2]=wpadB(b); % constant padding

%------------------------------------------------------%
[C1,I1]=wavedec(a,levels,waves);
[C2,I2]=wavedec(d2,levels,waves);
for j=1:levels

    rec=wrcoef('a',C1,I1,waves,j);
    rec1=(a-rec);
    rec1=rec1((length(f1)+1):(length(rec)-length(b1)));
    ha1=rec1.*rec1;

    rc=wrcoef('a',C2,I2,waves,j);
    rc1=(d2-rc);
    rc1=rc1((length(f2)+1):(length(rc)-length(b2)));
    ha21=rc1.*rc1;
%------------------------------------------------------------%
    t1=ha21; mm=length(ha21);
    kx=1;
    for k=mm:-1:1
        ha2(kx)=t1(k);
        kx=kx+1;
    end
%------------------------------------------------------------%
    ha3=(ha1+ha2)/2; % newly added
    m=length(ha3);
    
    k=1;
    for jj=m:-1:1
        ha4(k)=ha3(jj);
        k=k+1;
    end

%---------------------------------------------------------%

    b=2^(j-1)*ii;
    
   hr1=chopB(ha3,b);
   hr2=chopB(ha4,b);

   hr=[hr1 hr2];
%-------------------------------------------------------

plo=[];
for q = qmin:qstp:qmax
    if (q~=0)
        he3=hr.^(q/2); 
        av=(mean(he3)).^(1/q);
    else
        he3=log(hr);
        av=exp(mean(he3)/2);
    end
        plo=[plo av'];
end
poo=[poo plo'];
end

poo=log(poo);
bb=[];
for hh=1:levels
bb=[bb log(2^(hh-1)*ii)];
end

res=[bb' poo'];

%--------------------------------------------------------%

[m,n]=size(res);
h1='%12.8f '; h2=[ ]; 
for i=1:n
    h2=[h2 h1];
end
h3='\n';
h4=strcat(h2,h3);
sname = input('Enter file name for F_q Vs. q : ','s');

ff=res'; % important command to be noticed
fid = fopen(sname,'w');
fprintf(fid,h4,ff);
fclose(fid);

%--------------------------------------------------------%

[m,n]=size(res);

for i=2:n
    f=polyfit(res(2:8,1),res(2:8,i),1);
%    f=polyfit(res(:,1),res(:,i),1);
    fin(i-1)=f(:,1);
end

finx= qmin:qstp:qmax;
rest=[finx' fin'];

%--------------------------------------------------------%

[m,n]=size(rest);
h1='%12.8f '; h2=[ ]; 
for i=1:n
    h2=[h2 h1];
end
h3='\n';
h4=strcat(h2,h3);
sname = input('Enter file name for c_q Vs. q : ','s');

rest=rest'; % important command to be noticed
fid = fopen(sname,'w');
fprintf(fid,h4,rest);
fclose(fid);

%--------------------------------------------------------%