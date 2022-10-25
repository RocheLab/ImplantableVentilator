function [time,f_interp,coherences,crossspec]=plotDynFrequencyCoherenceOverTime(wcoh,wcs,f,fDesired,t_fDesired)
%plotDynFrequencyCoherenceOverTime Summary of this function goes here
%   wcoh = comes from wcoherence, matrics of coherence values for different
%   frequencies over time
%   f = vector of frequencies from wcoherence
%   fDesired = vector of frequencies
%   t_fDesired = vector of time points corresponding to fDesired


% 
% %find the desired frequency in vector f

time = [1:size(wcoh,2)].*0.001; %time in seconds
% 
% figure(2)
% plot(time,coherences)
% 

%linear interpolation of dynamic frequencies
f_interp = interp1(t_fDesired,fDesired,time);


for i=1:length(time)
   [M,fIndex] = min(abs(f-f_interp(i)));
   coherences(i) = wcoh(fIndex,i);
   crossspec(i) = wcs(fIndex,i);
end


plot(time,coherences)
ylabel('Coherence')
xlabel('Time [s]')
title('Coherence Along The Peak Frequency Ridge of the Flow Data')


end

