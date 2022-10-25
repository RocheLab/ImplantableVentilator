function [time,VolAutoData,VolAbsData,newVolAbsData,Correction] = calibrateVolumeRemoveIntegralNoise(time,VolAutoData,VolAbsData)
%calibrateVolumeRemoveIntegralNoise This removes the linear drift that
%results from integrating over a signal that contains noise. This is a
%broad level algorithm that allows more precise smaller scale algorithms
%like spirometry Normalization to work more accurates



%% Find the V minimums in VolAbsData


Vminind = islocalmin(VolAbsData,'FlatSelection','last','MinSeparation',600,'MinProminence',0.01);
Vminbounds = find(Vminind);
% figure
% plot(time,VolAbsData,time(Vminind),VolAbsData(Vminind),'r*')


% Find the gaps between stretches of consistent breaths

VindGaps = Vminbounds(2:end)-Vminbounds(1:end-1);
VGapInd = find(VindGaps>4000); %find gaps that are greater than 4 seconds long
Vsectionbounds = [1 (VGapInd+1);VGapInd length(Vminbounds)]; %this array will contain each index bound of breath bounds together as 1 column
Vgapbounds = [VGapInd;(VGapInd+1)]; %this array will contain each index bound of stretches without breaths together as 1 column
%now we have index bounds for which stretches have breaths to run an
%interpolation vs which ones dont

%% Run regressions and interpolations
%this will ignore data prior to the first Vmin detected    
Correction = [VolAbsData(Vsectionbounds(1,1)).*ones(1,Vminbounds(1)-1)]; %data prior to first Vmin will have a direct value subtraction of this volume
size(Correction)
    for i = 1:size(Vsectionbounds,2)
        if Vsectionbounds(1,i)==Vsectionbounds(2,i) %this is a single point and not a cluster
        Correction = [Correction VolAbsData(Vminbounds(Vsectionbounds(1,i)))];
        size(Correction)
        else
        %%%%%%%%%%%%%% a bound with breaths inside->run a regression
        y = VolAbsData(Vminbounds(Vsectionbounds(1,i):Vsectionbounds(2,i)))'; %V mins
        x = time(Vminbounds(Vsectionbounds(1,i):Vsectionbounds(2,i)))'; %times of Vmins
        X = [ones(length(x),1) x]; %X for the linear regression, ones column is for y intercept
        B = X\y;%B is the regression coefficient (slope)
        vregress = X*B;%vector of y points for the linear regression
        Rsq = 1 - sum((y - vregress).^2)/sum((y - mean(y)).^2)
%         %plot regression
%         figure
%         plot(x,y,x,vregress)
        
        
        vinterp = interp1([x(1) x(end)],[vregress(1) vregress(end)], time(Vminbounds(Vsectionbounds(1,i)):Vminbounds(Vsectionbounds(2,i))));%+/-1 avoids double assigning a subtraction value, preferentially allows the regression data
        Correction = [Correction vinterp];
        size(Correction)
        end
        %%%%%%%%%%%%%%

        
        %%%%%%%%%%%%%% a bound with no breaths-> first fill with zeros,
        %%%%%%%%%%%%%% then run a linear interpolation after finding the
        %%%%%%%%%%%%%% end bound
        if i<size(Vsectionbounds,2) %add filler zeros 
        filler = zeros(size(time(Vminbounds(Vgapbounds(1,i))+1:Vminbounds(Vgapbounds(2,i))-1)));
        Correction = [Correction filler];
        size(Correction)
        elseif i == size(Vsectionbounds,2) %finish out the time series with the last value
            tail = Correction(end).*ones(1,length(time)-Vminbounds(end));
            Correction = [Correction tail];
        end
        
        if i > 1  %replace the last set of filler zeros with interpolated values
            vinterp = interp1([time(Vminbounds(Vgapbounds(1,i-1))) time(Vminbounds(Vgapbounds(2,i-1)))],[Correction(Vminbounds(Vgapbounds(1,i-1))) Correction(Vminbounds(Vgapbounds(2,i-1)))],time(Vminbounds(Vgapbounds(1,i-1))+1:Vminbounds(Vgapbounds(2,i-1))-1));%+/-1 avoids double assigning a subtraction value, preferentially allows the regression data
            Correction((Vminbounds(Vgapbounds(1,i-1))+1):(Vminbounds(Vgapbounds(2,i-1))-1)) = vinterp;
        end
        
        
%         %%%%%%%%%%%%%% a bound with no breaths ->run a simple linear interpolation
%         if i<size(Vsectionbounds,2) %will ignore data after the last Vmin detected, bc Vgabbounds is 1 shorter than Vsectionbounds
% 
%         vinterp = interp1([time(Vminbounds(Vgapbounds(1,i))) time(Vminbounds(Vgapbounds(2,i)))],[Correction(end) VolAbsData(Vminbounds(Vgapbounds(2,i)))],time(Vminbounds(Vgapbounds(1,i))+1:Vminbounds(Vgapbounds(2,i))-1));%+/-1 avoids double assigning a subtraction value, preferentially allows the regression data
%         Correction = [Correction vinterp];
%         size(Correction)
%     %     v_interp = interp1([time(100) time(end-100)],[VolAbsData(100) VolAbsData(end-100)],time);
%     % Vol_interp = VolAbsData-v_interp;
%         elseif i == size(Vsectionbounds,2) %finish out the time series with the last value
%             tail = VolAbsData(Vminbounds(end)).*ones(1,length(time)-Vminbounds(end));
%             Correction = [Correction tail];
%         end



    end

if size(VolAbsData) == size(Correction)
    newVolAbsData = VolAbsData-Correction;
end

plotON = 0;
if plotON == 1 %for visualization
    figure
    plot(time(1:length(Correction)),Correction,time,VolAbsData)
    figure
    plot(time,newVolAbsData,'r',time,VolAutoData,'b')
end

end

