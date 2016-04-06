%DXM_matlab.m
%
% Audrey Lustig and Cedric Gaucherel, 2010-2011
%
% computes DXM according to:
%
% 
%   Gaucherel, C., A study of the possible extended influence of the ENSO phenomenon. Comptes Rendus Geoscience, 2004. 336(3): p. 175-185.
%   Gaucherel, C., Use of wavelet transform for temporal characterization of remote watersheds. Journal of Hydrology, 2002. 269(3-4): p. 101-121.
%   Gaucherel, C., Analysis of ENSO interannual oscillations using non-stationary quasi-periodic statistics: a study of ENSO memory. International Journal of Climatology, DOI: 10.1002/joc.1937 (2010).


% With slight modifications as highlited in:
% For more details see      Cyclostationnarity analysis of enso memory
%                           Student research report - 2011
%                           Audrey LUSTIG


%
% inputs:  - name: name of the analysed signal
%		 - pmin,pmax : range of width of the DXM map, x time resolution (integers)
%          - umax,vmax : width and height of the reference sub-map (integers) 
%          - nb_repetition : number of MC repetition forestimating the confidence interval
%
% calls:   - confident_level_with_RN : computes confidence level based on a red noise
%          - confident_level_with_WN : computes confidence level based on a white noise
%		 - Calculate the EAC value for a DXM map of quasi-period d such as pmin<=d<=pmax.

% For more details see      Cyclostationnarity analysis of enso memory
%                           Student research report - 2011
%                           Audrey LUSTIG


% Analyse DXM
clear all
close all
clc


% Read files
name='SOI.txt';
[S,Time]=readFile(name);
Mmean=mean(S);
Variance=var(S);
S=(S-Mmean)/Variance;
N=length(S);


% DXM analysis
umax=1;
vmax=1;
pmin=5;
pmax=10; 
EAC_profile=EAC_value(S,N,umax,vmax,pmin,pmax);



% DXM confident level
nb_repetition=3;
variance=nanstd(S)^2;
VM=find(isnan(S));
level=0.99;

% white noise
[EAC_cl_down,EAC_cl_up]=confident_level_with_WN(variance,VM,nb_repetition,N,umax,vmax,pmin,pmax,level);

% red noise
temp=S(~isnan(S)); % eliminer les valeurs manquantes pour le calcul de l'autocorrelation
[ACF,lags,bounds] = autocorr(temp,1);
foac=ACF(2);    
[EAC_cl_down_RN,EAC_cl_up_RN]=confident_level_with_RN(foac,VM,nb_repetition,N,umax,vmax,pmin,pmax,level);



% plot
MIN=[min(EAC_profile),min(EAC_cl_down),min(EAC_cl_down_RN)];
MAX=[max(EAC_profile),max(EAC_cl_up),max(EAC_cl_up_RN)];
[ACF,Lags,Bounds]=autocorr(S,pmax+1);
ACF=ACF(pmin:pmax);
figure(1);
clf;
subplot(3,1,1)
plot(Time,S,'k-')
set(gca,'Ylim',[min(S),max(S)])
set(gca,'Xlim',[Time(1),Time(end)])
xlabel('time')
title('Time Series')

subplot(3,1,2)
plot([pmin:pmax],ACF,'b','LineWidth',2)
set(gca,'Xlim',[Time(1),Time(end)])
hold on
plot([pmin pmax],[Bounds(1) Bounds(1)],'k--')
plot([pmin pmax],[Bounds(2) Bounds(2)],'k--')
hold off
xlabel('Lags')
title('Autocorrelation Function')
legend('ACF','99% confident level')
set(gca,'Xlim',[pmin,pmax])

subplot(3,1,3)
plot([pmin:pmax],EAC_profile,'r','LineWidth',2)
hold on
plot([pmin:pmax],EAC_cl_down,'k--')
plot([pmin:pmax],EAC_cl_down_RN,'b--')
plot([pmin:pmax],EAC_cl_up,'k--')
plot([pmin:pmax],EAC_cl_up_RN,'b--')
hold off
xlabel('Pseudo-period')
title('EAC profile')
legend('EAC','99% confident level WN','99% confident level RN')
set(gca,'Ylim',[min(MIN),max(MAX)])
set(gca,'Xlim',[pmin,pmax])


% figure(2);
% clf;
% Time=[1:length(S)];
% for i=1:12
%     period=i+4;
%     subplot(3,4,i);
%     plot_map(S,Time,period);
%     titre=strcat('période : ', num2str(period));
%     title(titre)
% end


period=20;
Time=[1:length(S)];
plot_map(S,Time,period);
