function analyzeCrossSpecandVol(crossspec,vpks,vlocs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


[theta, rho] = cart2pol(real(crossspec), imag(crossspec));

figure('Position',[1000 430 500 300])
scatter(theta(vlocs),vpks,[],vlocs)
title('Phase Delay from Cross Spectrum')
xlim([-pi pi])
xlabel('Phase Delay [rad]')
ylabel('Tidal Volume [L]')

figure('Position',[1000 50 500 300])
scatter(rho(vlocs),vpks,[],vlocs)
title('Coherence from Cross Spectrum Magnitude')
xlabel('Coherence')
ylabel('Tidal Vol [L]')

figure('Position',[1000 50 500 300])
plot((1:length(theta)).*0.001,theta)
title('Phase Delay Along The Peak Frequency Ridge of the Data')
ylim([-pi pi])
ylabel('Phase Delay [rad]')
xlabel('Time [s]')
end

