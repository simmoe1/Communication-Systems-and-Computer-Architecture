%Ebrahim Simmons 400200042, Allen Mei 400202911

clear
hold off
format long e
N = 4096; %No. of FFT samples
sampling_rate = 100.0e3; %unit Hz
tstep = 1/sampling_rate;
tmax = N*tstep/2;
tmin = -tmax;
tt = tmin:tstep:tmax-tstep;
fmax = sampling_rate/2;
fmin = -fmax;
fstep = (fmax-fmin)/N;
freq = fmin:fstep:fmax-fstep;
fc=10e3;
Ac = 1;
ct=Ac*cos(2*pi*fc*tt);
Am=1;
fm=1e3;
mt = Am*cos(2*pi*fm*tt);
st = mt.*ct;
figure(1)
Hp1 = plot(tt,st);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('s(t) (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('DSB-SC modulated wave : Time domain');
axis([-10/fc 10/fc min(st) max(st)])
Sf = fftshift(fft(fftshift(st)))/(2*fmax);
figure(2)
Hp1=plot(freq,abs(Sf));
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Frequency (Hz) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('|S(f)|');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Spectrum of the DSB-SC wave S(f)');
axis ([-30e3 30e3 0 max(abs(Sf))])

theta1=0;
theta2=pi/2;
theta3=pi/4;
lo_0 = cos(2*pi*fc*tt + theta1);
st1_0 = st .* lo_0;
figure(3)
Hp1=plot(tt,st1_0);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel(' s hat(t) (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Resulting signal, s hat(t) : Time domain');
axis([-10/fc 10/fc min(st1_0) max(st1_0)])
Sf1_0 = fftshift(fft(fftshift(st1_0)))/(2*fmax);
figure(4)
Hp1=plot(freq,abs(Sf1_0));
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Frequency (Hz) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('|S hat(f)|');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Spectrum S hat(f)');
axis ([-50e3 50e3 0 max(abs(Sf1_0))])
%Low pass filtering
f_cutoff = fc;
%ideal low pass filter
n=1;
for f = freq
 if abs(f) < f_cutoff
 Hf(n) = 1;
 else
 Hf(n) = 0;
 end
n=n+1;
end
lo_2 = cos(2*pi*fc*tt + theta2);
st1_2 = st .* lo_2;
Sf1_2 = fftshift(fft(fftshift(st1_2)))/(2*fmax);
lo_4 = cos(2*pi*fc*tt + theta3);
st1_4 = st .* lo_4;
Sf1_4 = fftshift(fft(fftshift(st1_4)))/(2*fmax);
Mf1_0 = Sf1_0 .* Hf;
mt1_0 = 2*fmax*fftshift(ifft(fftshift(Mf1_0)));
Mf1_2 = Sf1_2 .* Hf;
mt1_2 = 2*fmax*fftshift(ifft(fftshift(Mf1_2)));
Mf1_4 = Sf1_4 .* Hf;
mt1_4 = 2*fmax*fftshift(ifft(fftshift(Mf1_4)));
figure(5)
Hp1=plot(tt,mt1_0,'r',tt,mt*0.5,'b.');
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('m hat(t)(Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Output of LPF,m hat(t), theta=0 : Time domain');
axis([-0.01 0.01 min(mt*0.5) max(mt*0.5)])
legend('LPF output', 'message sig');
figure(6)
Hp1=plot(tt,mt1_2,'r',tt,mt*0.5,'b.');
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('m hat(t)(Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Output of LPF,m hat(t), theta=pi/2 : Time domain');
axis([-0.01 0.01 min(mt*0.5) max(mt*0.5)])
legend('LPF output', 'message sig');
figure(7)
Hp1=plot(tt,mt1_4,'r',tt,mt*0.5,'b.');
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('m hat(t)(Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Output of LPF,m hat(t), theta=pi/4 : Time domain');
axis([-0.01 0.01 min(mt*0.5) max(mt*0.5)])
legend('LPF output', 'message sig');