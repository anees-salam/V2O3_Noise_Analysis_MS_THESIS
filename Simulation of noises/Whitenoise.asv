clear all;
close all;

num=10000;


for i = 1:num
temp(i) = num;
end

x= wgn(num,1,0);
meanr=mean(x);


delx=x-mean(x);

Fs=14.5;

figure(1);
subplot(2,2,1);
hist(x);
title('Histogram')
xlabel('Voltage')
ylabel('Frequency')



t(1)=0;
for i = 1:length(x)-1
t(i+1) = t(i)+1/Fs;
end

subplot(2,2,2)
plot(temp,x);
title('V vs T')
xlabel('Temperature')
ylabel('Voltage')

subplot(2,2,3);
plot(t,x);
title('Time vs Voltage')
xlabel('Time(s) ')
ylabel('Voltage')

subplot(2,2,4)
plot(t,temp);
title('Time vs Temperature')
xlabel('Time ')
ylabel('Temp(K)')


N = length(delx);
xdft = fft(delx);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(delx):Fs/2;

figure(2);
%subplot(1,2,1);
%hold on;
plot(freq,10*log10(psdx))
loglog(freq,psdx,'r.-');

title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (V^2/Hz)')
legend('n4')

figure(31);

%hold on;
nx = max(size(delx)); na = 50;
w = hanning(floor(nx/na));
[pxx,f] = pwelch(delx,w,0,[],Fs);


loglog(f,pxx,'r.-');
title('Averaging=100')
xlabel('Frequency [Hz]');
ylabel('Power Density [Vˆ2/Hz]');
legend('80','100','120','130');


