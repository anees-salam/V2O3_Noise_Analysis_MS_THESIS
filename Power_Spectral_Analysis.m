clear all
close all
fina=[152,152.5,153,153.5,154,154.5,155,155.5,156,156.5];


for i=1:length(fina)
    name=strcat(num2str(fina(i)),'.txt');
    nameout=strcat(num2str(fina(i)),'-originwithout-a.txt');
    nameout1=strcat(num2str(fina(i)),'-originfftwithout-a.txt');

file = load (name);

fileID = fopen(nameout,'w');
fileID1=fopen(nameout1,'w');

x = file(:,2); 
temp=file(:,1);
ui=size(x);
y=smooth(x,ui(1)/10);
meanr=mean(x);
delx=x-y;
delx1=x-mean(x);
%delx=delx./(x.^2);


Fs=14.5;
t(1)=0;
for i = 1:length(x)-1
    t(i+1) = t(i)+1/Fs;
end

N = length(delx1);
xdft = fft(delx1);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(delx1):Fs/2;



for i = 1:length(freq)
    if log10(freq(i))<0
fprintf(fileID1,'%.20f %.50f \r\n',freq(i),psdx(i));
    else
        psdx(i)=0;
    end
    
  

end




figure(2);

loglog(freq,psdx,'g.-');

title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (V^2/Hz)')
legend('127','142','149','165','175','200');
hold on;


figure(3);

nx = max(size(delx)); na = 55;
w = hanning(floor(nx/na));
[pxx,f] = pwelch(delx,w,0,[],Fs);

clear nx na w 



for i = 1:length(f)
        if log10(f(i))<0
        fprintf(fileID,'%.20f %.50f \r\n',f(i),pxx(i));
       
        else
    pxx(i)=0;
        end
end

loglog(f,pxx,'.-');
title('Averaging=31')
xlabel('Frequency [Hz]');
ylabel('Power Density [Vˆ2/Hz]');
legend('127','142','149','165','175','200');
hold on;
%}


end



