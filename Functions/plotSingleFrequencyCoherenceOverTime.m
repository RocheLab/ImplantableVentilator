function [fDesired,coherences,vpks,vlocs] = plotSingleFrequencyCoherenceOverTime(wcoh,f,VolAutoData)
%plotSingleFrequencyCoherenceOverTime Summary of this function goes here
%   wcoh = comes from wcoherence, matrics of coherence values for different
%   frequencies over time
%   f = vector of frequencies from wcoherence
%   fDesired = frequency look for

%find fixed estimate of BPM
[vpks,vlocs]=findpeaks(VolAutoData,'MinPeakHeight',0.02,'MinPeakProminence',0.01);
fDesired=length(vpks)/(length(VolAutoData)*0.001);% in Hz (breaths/s)
BPM = fDesired*60;


%find the desired frequency in vector f
[M,fIndex] = min(abs(f-fDesired));
coherences = wcoh(fIndex,:);
time = [1:length(coherences)].*0.001; %in s


% plot(time,coherences)
% ylabel('Coherence')
% xlabel('Time [s]')
% title(sprintf('Coherence Along The Estimated BPM, %.2f bpm = %.2f Hz',BPM,fDesired))

end

