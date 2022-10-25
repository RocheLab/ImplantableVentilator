function [Rsq] = Rsquared(y,yfit)
%Rsquared Calculates given R2 for polyfit
%   y = y value data 
%   yfit = y data of the fitted curve

SStot = sum((y-mean(y)).^2);                    % Total Sum-Of-Squares
SSres = sum((y-yfit).^2);                       % Residual Sum-Of-Squares
Rsq = 1-SSres/SStot;                            % R^2
end

