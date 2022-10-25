%generateFig_IdealInputCurves_Revision

clear all
close all

load('shapes.mat')
% Color5 = hex2rgb('#C9DEDB');
% Color10 = hex2rgb('#83ABA4');
% Color15 = hex2rgb('#428C80');
Color20 = hex2rgb('#146759');
ColorTri = hex2rgb('#DDCC77');
ColorSq = hex2rgb('#CC6677');
% ColorSpont = [0.3 0.3 0.3];


f4=figure('Position',[100 100 250 150]);
plot(curve.PR1Timems.*0.001+0.02,curve.PR1PPSI,'Color',Color20,'LineWidth',1.5)
plottype4 = 'Curve20_P';
xlabel('Time [s]')
xlim([0 1.818])
ylabel('Input Pressure [psi]')
ylim([0 21])
yticks([0 5 10 15 20])

X1 = curve.PR1Timems.*0.001+0.02;
Y1 = curve.PR1PPSI;

f5=figure('Position',[100 100 250 150]);
plot(triangle.PR1Timems.*0.001+0.02,triangle.PR1PPSI,'Color',ColorTri,'LineWidth',1.5)
plottype5 = 'Tri20_P';
xlabel('Time [s]')
xlim([0 1.818])
ylabel('Input Pressure [psi]')
ylim([0 21])
yticks([0 5 10 15 20])

X2=triangle.PR1Timems.*0.001+0.02;
Y2 = triangle.PR1PPSI;

f6=figure('Position',[100 100 250 150]);
plot(square.PR1Timems.*0.001.*(1818/5000)+0.02,square.PR1PPSI,'Color',ColorSq,'LineWidth',1.5)
plottype6 = 'Sq20_P';
xlabel('Time [s]')
xlim([0 1.818])
ylabel('Input Pressure [psi]')
ylim([0 21])
yticks([0 5 10 15 20])

X3 = square.PR1Timems.*0.001.*(1818/5000)+0.02;
Y3 = square.PR1PPSI;
%% Save
pathWithFolderName =  strcat(pwd,'\Figures For Paper\');
Condition = 'IdealCurve';


figName4 = strcat(plottype4,Condition);
figFileName4 = strcat(pathWithFolderName,figName4);
exportgraphics(f4,strcat(figFileName4,'.eps'),'ContentType','vector') %will save figure f 
figName5 = strcat(plottype5,Condition);
figFileName5 = strcat(pathWithFolderName,figName5);
exportgraphics(f5,strcat(figFileName5,'.eps'),'ContentType','vector') %will save figure f 
figName6 = strcat(plottype6,Condition);
figFileName6 = strcat(pathWithFolderName,figName6);
exportgraphics(f6,strcat(figFileName6,'.eps'),'ContentType','vector') %will save figure f 
