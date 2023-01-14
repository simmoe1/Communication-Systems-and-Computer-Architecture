%%SQUARE WAVE GENERATOR
%3TR4 Lab #1
%Ebrahim Simmons, 400200042, simmoe1
%Allen Mei, meia6, 400202911

clc
clear all
hold off

%fundamental frequency of square wave
f0=10000;     
%period
T0 = 1/f0;
tstep = 0.005*T0;

%number of samples withink 3*T0
num_sample = 3*T0/tstep + 1; 
%number of samples within T0
num_sample1 = T0/tstep + 1; 

tt = -1.5*T0:tstep:1.5*T0;

%time vector for perdio -2T0 - 2T0
tt1 = -2*T0:tstep:2*T0;
%input, square wave in the period -2T0 - 2T0
gp1 = square(2*pi*f0.*tt1,50);
figure(1)
Hp1 = plot(tt1,gp1);
set(Hp1,'LineWidth',2.5)
Ha = gca;
set(Ha,'Fontsize',15)
title('Input - Time Domain')
%pause

%Fourier series representation of signal (Amplitude Spectrum)

%number of harmonics 
N=100;
nvec = -N:N;
c_in = zeros(size(nvec));
for n = nvec
    m = n+N+1;
    c_in(m) = 0.5*sinc(n/2);
end
%frequency vector
f = nvec*f0;
figure(2)
Hp1=stem(f,c_in);
axis([-8*f0 8*f0 min(c_in) max(abs(c_in))])
set(Hp1,'LineWidth',2.5)
Ha = gca;
set(Ha,'Fontsize',15)
title('Magnitude Spectrum of Input')
%pause

%% Fourier series representation of signal (Phase Spectrum)
%% Designing the 2nd order Butterworth filter
%cutoff frequency of filter
fc=12500     
%fc = 11500 to 13900;

%filter transfer function 
Hf = 1./((1i*f/fc).^2+1.414*(1i*f/fc)+1) ;
figure(3)
Hp1 = plot(f,20*log10(abs(Hf)));
xlim([-4*f0 4*f0])
set(Hp1,'LineWidth',2.5)
Ha = gca;
set(Ha,'Fontsize',15)
title('Amplitude Response of Butterworth Filter')
%pause

%fourier coefficients of the filter output 
c_out = c_in .* Hf; 

figure(4)
stem(f,c_in,'r','LineWidth',2.5);
hold on
stem(f,abs(c_out),'b','LineWidth',2.5);
hold off
axis([-8*f0 8*f0 min(c_in) max(abs(c_in))])
Ha = gca;
set(Ha,'Fontsize',15)
title('Magnitude Spectrum of Filter Output and Input')
Ha = gca;
set(Ha,'Fontsize',15)
legend('input','output')
%pause


%% Construct the output signal from the Cout Fourier coefficients

A = zeros(2*N+1,ceil(num_sample));
for n = nvec
    m=n+N+1;
    A(m,:) = c_out(m) .* exp(1i*2*pi*n*f0*tt);
end
gp_out = sum(A);
figure(5)
Hp1 = plot(tt,real(gp_out),'b',tt1,gp1,'r');
set(Hp1,'LineWidth',2.5)
Ha = gca;
set(Ha,'Fontsize',15)
title('Filter Input and Output-Time Domain')
set(Ha,'Fontsize',16)
legend('output','input')
%pause

%% Plot transfer function of filter
figure(6)
Hp1 = plot(f,abs(Hf));
set(Hp1,'LineWidth',2.5)
Ha = gca;
set(Ha,'Fontsize',15)
title('Absolute Transfer Function of Filter')