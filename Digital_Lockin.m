clc
clear
format long Eng
file = load('dphi_dt_t(ms).txt'); % loading data file containing
                            % pickup response
dB_dt = file(:,1);           %assuming area = 1 then dB/dt = dphi/dt
t = file(:,2);           %time in milli second
t = t./1000;            %time in second
fs = 1/(t(2)-t(1));     %actual sampled freq before interpolation
Fs = 1e6;               %sampling freq after interpolation
dT = 1/Fs;               %sampling time
ti = 0:dT:t(length(t));         %interpolated time in second
dB_dti = zeros(size(ti));        %array for interpolated dB/dt


%interpolating using linear interpolation formula
i = 1;
for j = 2:length(t)
while i <= length(ti) && ti(i) <= t(j)
dB_dti(i) = dB_dt(j-1) + (dB_dt(j)-dB_dt(j-1))*(ti(i)-t(j-1))/(t(j)-t(j-1));

i = i+1;
  end
end

%calculating magnetic field taking initial condition B=0 at t=0
B = zeros(size(ti)); %array for magnetic field
clear i
%integrating dB/dt

for i = 1:length(ti)-1
B(i+1) = B(i)+dB_dti(i)*dT;
end

%model for sample resistance R(B) = Ro + beta*B(t)
Ro = 0.1;
beta = 5;
R = Ro + beta.*B;

%reference current
Fm = 2e5; %modulation freq in Hz
Io = 1; %amplitude
I = Io*sin(2*pi*Fm*ti); %current


%sample response
alpha = 0.01;
V = I.*R + alpha*dB_dti;
noise = 0.0005*sin(2*pi*50*ti)+0.001*sin(2*pi*10000*ti)+0.0005*sin(2*pi*500*ti).*(rand(1,length(V))-0.5); %random noise
V = V + noise; %adding noise to sample response
figure(1)
subplot(2,2,1)
plot(ti,V)
title('V vs t (modulated signal)')
V1 = V.*I; % multiplying with reference

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = order of the filter
% f = vector of pairs of normalized freq point, specified in the range
% between 0 and Nyquist freq. The freq must be in increasing order.
% a = vector containning the desired amlitudes at the point specified in f.
f = [1000,70000];
a = [1,0];
rp = 0.001;
rs = 170;
dev = [(10^(rp/20)-1)/(10^(rp/20)+1),10^(-rs/20)];


tic
[n,fo,ao,w] = firpmord(f,a,dev,Fs);
b = firpm(n,fo,ao,w);
V_filtered = filter(b,1,V1)*2/Io;
toc


subplot(2,2,2)
V_bf = Io*R+alpha*dB_dti + noise; %sample response(V) before filter
plot(ti,V_bf,ti,V_filtered)
ylim([0.09,0.14])
legend('before filter','after filter')
title('V vs t')
subplot(2,2,3)
plot(B,V_filtered,B,Io*R)
ylim([0.098,0.127])
legend('after filter','theoritical')
title('V vs B')

[h,f] = freqz(b,1,1e4,Fs); % gives the frequency response of filter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculating Fourier Transform of sample response
x = V1;
m = length(x); % no. of samples
FFTX=fft(x);
N_up = ceil((m+1)/2); % no. of unique points.
% FFTX is symmetric about nyquist freq point
MX=abs(FFTX(1:N_up)); % amplitude of unique FFTX points
MX=2*MX/m; % Scale the FFT so that it is not a function
% of the length(x).
Fn = Fs/2; % Nyquist frequency
F=(0:N_up-1)*Fn/(N_up-1); % frequency vector
subplot(2,2,4)
loglog(f,abs(h),F,MX)
title('frequency response of filter')