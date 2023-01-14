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
ct1=Ac*cos(2*pi*fc*tt);
Tm = 0.001;
mt1 = 2*sinc(tt/Tm);
st1 = mt1.*ct1;
ct2=Ac*sin(2*pi*fc*tt);
mt2 = (sinc(tt/Tm)).^2;
st2 = mt2.*ct2;
st = st1+st2;


figure(1)
Hp1 = plot(tt,mt1);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('message  m1(t)  (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('First message signal : Time domain');
axis([-0.01 0.01 0 1.1])
pause(1)

figure(2)
Hp1 = plot(tt,mt2);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('message m2(t)(Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Second message signal : Time domain');
axis([-0.01 0.01 min(mt2) max(mt2)])
pause(1)

figure(3)
Hp1 = plot(tt,st);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('s(t)  (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('QAM signal : Time domain');
axis([-0.01 0.01 min(st) max(st)])
pause(1)

Mf1 = fftshift(fft(fftshift(mt1)))/(2*fmax);
figure(4)

Hp1=plot(freq,abs(Mf1))
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Frequency (Hz) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('|M1(f)|');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Spectrum of the message signal 1');
axis ([-15e3 15e3 0 200/(2*fmax)])
pause(1)

Mf2 = fftshift(fft(fftshift(mt2)))/(2*fmax);
figure(5)
Hp1=plot(freq,abs(Mf2))
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Frequency (Hz) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('|M2(f)|');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Spectrum of the message signal 2');
axis ([-15e3 15e3 0 200/(2*fmax)])
pause(1)

figure(6)
Sf = fftshift(fft(fftshift(st)))/(2*fmax);
Hp1=plot(freq,abs(Sf));
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Frequency (Hz) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('|S(f)|');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Spectrum of the QAM signal S(f)');
axis ([-30e3 30e3 0 150/(2*fmax)])
pause(1)
%@Rx
%QAM demux and demodulation

%Local oscillator at the receiver perfectly synchronized
thet=0;
lo1 = cos(2*pi*fc*tt + thet);
lo2 = sin(2*pi*fc*tt + thet); 
stt1 = st .* lo1;
Sff1 = fftshift(fft(fftshift(stt1)))/(2*fmax);
stt2 = st .* lo2;
Sff2 = fftshift(fft(fftshift(stt2)))/(2*fmax);

%Low pass filtering
f_cutoff = 15e3;
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
Mff1 = Sff1 .* Hf;
mt_hat1 = 2*fmax*fftshift(ifft(fftshift(Mff1)));
Mff2 = Sff2 .* Hf;
mt_hat2 = 2*fmax*fftshift(ifft(fftshift(Mff2)));
figure(7)
Hp1=plot(tt,mt_hat1,'r',tt,mt1*0.5,'b.');
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('m hat1(t)(Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Output of LPF1,  m hat1(t)  : Time domain');
axis([-0.01 0.01 min(mt1*0.5) max(mt1*0.5)])
legend('LPF1 output', 'message sig1');

figure(8)
Hp1=plot(tt,mt_hat2,'r',tt,mt2*0.5,'b.');
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('m hat2(t)(Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Output of LPF2,  m hat2(t)  : Time domain');
axis([-0.01 0.01 min(mt_hat2) max(mt_hat2)])
legend('LPF2 output', 'message sig2');













     
