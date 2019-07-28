%LOGPROCESSING - fajl za procesiranje podataka logovanih u fajl
%   Podaci koji su iz C programa logovani u fajl ucitavaju se
%   u Matlab workspace kako bi se prikazali relevantni grafici
%   i izvrsila analiza rezultata.
%
%   Autori:
%   Dragan Bozinovic (bd150211d@student.etf.bg.ac.rs)
%   Kristijan Mitrovic (mk150214d@student.etf.bg.ac.rs)

%% Prikupljanje izlogovanih podataka
clear all; close all;

% Otvaranje fajla u kom su izlogovani podaci iz C programa
data = csvread('DDSlog.csv');
% i njihovo smestanje u odgovarajuce promenljive
f0 = data(1,1);
phaseAcc = data(2:end,1);
phaseAccOut = data(2:end,2);
sin = data(2:end,3);
cos = data(2:end,4);
sinFilt = data(2:end,5);
cosFilt = data(2:end,6);

% Self-explanatory pomocne promenljive
n = numel(sin); fs = 100e6; dt = 1/fs; df = fs/n;
numOfPeriodesToDisplay = 1; range = 1:fs/f0*numOfPeriodesToDisplay;
t = (range-1)*dt; f = (0:n/2-1)*df;

%% Plotovanje signala 
figure;
%subplot(211);
plot(t, phaseAcc(range)/2^40, 'LineWidth', 1.5); grid on; 
xlabel('t [s]', 'Interpreter', 'latex'); 
title('Normalizovana cetrdesetobitna vrednost faznog akumulatora', 'Interpreter', 'latex');
subplot(212);
plot(t, phaseAccOut(range)/2^12, 'LineWidth', 1.5); grid on; 
xlabel('t [s]', 'Interpreter', 'latex'); 
title('Normalizovan dvanaestobitni izlaz faznog akumulatora', 'Interpreter', 'latex');

figure;
% subplot(211);
% plot(t, sin(range), 'LineWidth', 1.5); grid on;
% xlabel('t [s]', 'Interpreter', 'latex'); 
% title('Sinusni signal na izlazu CORDIC algoritma', 'Interpreter', 'latex');
% subplot(212); 
plot(t, cos(range), 'LineWidth', 1.5); grid on;
xlabel('t [s]', 'Interpreter', 'latex'); 
%title('Kosinusni signal na izlazu CORDIC algoritma', 'Interpreter', 'latex');

figure;
subplot(211);
plot(t*dt, sinFilt(range), 'LineWidth', 1.5); grid on;
xlabel('t [s]', 'Interpreter', 'latex'); 
title('Sinusni signal na izlazu $\sin x/x$ bloka', 'Interpreter', 'latex');
subplot(212); 
plot(t, cosFilt(range), 'LineWidth', 1.5); grid on;
xlabel('t [s]', 'Interpreter', 'latex'); 
title('Kosinusni signal na izlazu $\sin x/x$ bloka', 'Interpreter', 'latex');

%% Plotovanje spektra
% DFT sinusa na izlazu CORDIC-a i DFT filtriranog sinusa
sinFFT = abs(fftshift(fft(sin)))/n;
sinFiltFFT = abs(fftshift(fft(sinFilt)))/n;
SQNR = -(6.02*14 + 1.76);
fftNoiseFloor =  SQNR - 10*log10(n/2);

figure;
subplot(211);
plot(f, 20*log10(1e-25+sinFFT(n/2+1:end)), 'LineWidth', 1.5); grid on; hold on;
plot(f, SQNR*ones(n/2,1), 'r'); plot(f, fftNoiseFloor*ones(n/2,1), 'g');
legend({'Spektar signala', 'SQNR', '$FFT_{noise}$'}, 'Interpreter', 'latex');
xlabel('f [Hz]', 'Interpreter', 'latex'); 
title('Jednostrani spektar signala na izlazu CORDIC algoritma', 'Interpreter', 'latex');
subplot(212);
plot(f, 20*log10(1e-25+sinFiltFFT(n/2+1:end)), 'LineWidth', 1.5); grid on;
xlabel('f [Hz]', 'Interpreter', 'latex');
title('Jednostrani spektar signala na izlazu $\sin x/x$ bloka', 'Interpreter', 'latex');

[~,maxInd] = max(sinFFT(n/2+1:end));
fprintf('Sinusoida je na frekvenciji: %f Hz\n', maxInd*df);
