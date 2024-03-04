clear;
clc 
%% Reading Data
TNo5 = xlsread("Data.xlsx", "B5:B58");
FNo5 = xlsread("Data.xlsx", "C5:C58");

TNo6 = xlsread("Data.xlsx", "E5:E45");
FNo6 = xlsread("Data.xlsx", "F5:F45");

TNo7 = xlsread("Data.xlsx", "H5:H41");
FNo7 = xlsread("Data.xlsx", "I5:I41");

TNo8 = xlsread("Data.xlsx", "K5:K16");
FNo8 = xlsread("Data.xlsx", "L5:L16");

TNo9 = xlsread("Data.xlsx", "N5:N65");
FNo9 = xlsread("Data.xlsx", "O5:O65");

TNo10 = xlsread("Data.xlsx", "Q5:Q43");
FNo10 = xlsread("Data.xlsx", "R5:R43");

TNo11 = xlsread("Data.xlsx", "T5:T28");
FNo11 = xlsread("Data.xlsx", "U5:U28");

TNo12 = xlsread("Data.xlsx", "W5:W22");
FNo12 = xlsread("Data.xlsx", "X5:X22");
%% Core
g = 9.81; % Gravity
Lx = 800; % Plane Dimensions
Ly = 1000; 
Bc = 20; % Channel Width (m)
dx = 100; % Grid Size (m)
dy = 100;
dtsec = 1; % Temporal Discretization (sec) 
ERD = 30; % Effective Rainfall Duration (min)
TST = 120; %Total Simulation Time (min)
np = 0.016; % Manning
nc = 0.28; %Channel
Sox = 0.005; %Bed Slope (plane)
Soy = 0.01;
Soc = 0.02; %Bed Slope (channel)
Ie = 6.967; %Rainfall Data

% Calculations
xend = Lx/dx; %No of dx and dy
yend = Ly/dy;
dtmin = dtsec/60; %Temporal Discretization (min)
t = (0:dtmin:TST)'; %t array in (min)
Iemps = Ie*1e-3/3600; %effective rainfall (m/s)
TSTsec = TST*60; %TST in sec
ERDsec = ERD*60; %ERD in sec
tend = length(t); % No. of time intervals
TotalCatchArea_m2 = 2*Lx*Ly; %total catchment area in sq.m
TotalCatchArea_Km2=TotalCatchArea_m2*1e-6; % Total Catchment Area in sq.km.
QpeakActual=Ie*1e-3/3600*TotalCatchArea_m2;
actualvol=Ie*1e-3/3600*ERDsec*TotalCatchArea_m2; % Calculating Actual Vol
% DECLARING AND INITIALIZING VARIABLES
Qx=zeros(yend,xend+1,2); % specific discharges in x-direction
Qy=zeros(yend+1,xend,2); % specific discharges in y-direction
Qxout=zeros(yend,tend); % qxout at the d/s of plane in x-direction
Sox = Sox.*ones(yend,xend); % Bed Slope Gridwise Matrix
Soy = Soy.*ones(yend,xend); % 2D --> 1D if Soy=0
% Travel time
Kx=zeros(yend,xend,2);
Ky=zeros(yend,xend,2);
% Weighting Coefficient
Xx=zeros(yend,xend,2);
Xy=zeros(yend,xend,2);
% Other Variables
yoMx=zeros(yend,xend,2);
yoMy=zeros(yend,xend,2);
Q3x=zeros(yend,xend,2);
Q3y=zeros(yend,xend,2);
yMx=zeros(yend,xend,2);
yMy=zeros(yend,xend,2);
QMx=zeros(yend,xend,2);
QMy=zeros(yend,xend,2);
FMx=zeros(yend,xend,2);
FMy=zeros(yend,xend,2);
yx=zeros(yend,xend,2);
yy=zeros(yend,xend,2);
for k=2:tend % temporal index
if k>ERDsec/dtsec+1
Iemps=0;
end
for j=1:xend % spatial index x-direction
for i=1:yend % spatial index y-direction
if i==yend
Soy(i,j) = 0;
end
% UNREFINED
Kx(i,j,2)=Kx(i,j,1);
Xx(i,j,2)=Xx(i,j,1);
CDx=dtsec+2.*Kx(i,j,2).*(1-Xx(i,j,2));
C1x=(dtsec-2.*Kx(i,j,2).*Xx(i,j,2))./CDx;
C2x=(dtsec+2.*Kx(i,j,1).*Xx(i,j,1))./CDx;
C3x=(-dtsec+2.*Kx(i,j,1).*(1-Xx(i,j,1)))./CDx;
C4x=dtsec./CDx;
% UNREFINED Qx
if i==yend
Qx(i,j+1,2)=C1x.*Qx(i,j,2)+C2x.*Qx(i,j,1)+C3x.*Qx(i,j+1,1)+C4x.*2.*((Iemps)*dx*dy) + C4x.*(Qy(i,j,1)+Qy(i,j,2));
else
Qx(i,j+1,2)=C1x.*Qx(i,j,2)+C2x.*Qx(i,j,1)+C3x.*Qx(i,j+1,1)+C4x.*2.*((Iemps)*dx*dy)-C4x.*2.*(Qy(i+1,j,1)-Qy(i,j,1));
end
% UNREFINED
Q3x(i,j,2)=Xx(i,j,2)*Qx(i,j,2)+(1-Xx(i,j,2))*Qx(i,j+1,2);
% UNREFINED
yMx(i,j,2)=(Q3x(i,j,2).*np./dy./Sox(i,j).^(1/2))^(3/5);
voMx=Sox(i,j).^(1/2)./np.*(yMx(i,j,2)).^(2/3);
% REFINED
Kx(i,j,2)=dx./voMx;
QMx(i,j,2)=0.5*(Qx(i,j,2)+Qx(i,j+1,2));
%FMx(i,j,2)=sqrt(qMx(i,j,2).^2./(g.*yMx(i,j,2).^3));
CoMx=5/3.*voMx;
% REFINED
Xx(i,j,2)=1/2-Q3x(i,j,2)./(2.*Sox(i,j).*dy.*CoMx.*dx);
CDx=dtsec+2.*Kx(i,j,2).*(1-Xx(i,j,2));
C1x=(dtsec-2.*Kx(i,j,2).*Xx(i,j,2))./CDx;
C2x=(dtsec+2.*Kx(i,j,1).*Xx(i,j,1))./CDx;
C3x=(-dtsec+2.*Kx(i,j,1).*(1-Xx(i,j,1)))./CDx;
C4x=dtsec./CDx;
% REFINED
if i==yend
Qx(i,j+1,2)=C1x.*Qx(i,j,2)+C2x.*Qx(i,j,1)+C3x.*Qx(i,j+1,1)+C4x.*2.*((Iemps)*dx*dy) + C4x.*(Qy(i,j,1)+Qy(i,j,2));
else
Qx(i,j+1,2)=C1x.*Qx(i,j,2)+C2x.*Qx(i,j,1)+C3x.*Qx(i,j+1,1)+C4x.*2.*((Iemps)*dx*dy)-C4x.*2.*(Qy(i+1,j,1)-Qy(i,j,1));
end
if Soy(i,j) == 0
else
% UNREFINED
Ky(i,j,2)=Ky(i,j,1);
Xy(i,j,2)=Xy(i,j,1);
CDy=dtsec+2.*Ky(i,j,2).*(1-Xy(i,j,2));
C1y=(dtsec-2.*Ky(i,j,2).*Xy(i,j,2))./CDy;
C2y=(dtsec+2.*Ky(i,j,1).*Xy(i,j,1))./CDy;
C3y=(-dtsec+2.*Ky(i,j,1).*(1-Xy(i,j,1)))./CDy;
C4y=dtsec./CDy;
% UNREFINED Qy
Qy(i+1,j,2)=C1y.*Qy(i,j,2)+C2y.*Qy(i,j,1)+C3y.*Qy(i+1,j,1)+C4y.*2.*((Iemps)*dx*dy)-C4y.*2.*(Qx(i,j+1,1)-Qx(i,j,1));
% UNREFINED
Q3y(i,j,2)=Xy(i,j,2)*Qy(i,j,2)+(1-Xy(i,j,2))*Qy(i+1,j,2);
% UNREFINED
yMy(i,j,2)=(Q3y(i,j,2)./dx.*np./Soy(i,j).^(1/2))^(3/5);
voMy=Soy(i,j).^(1/2)./np.*(yMy(i,j,2)).^(2/3);
% REFINED
Ky(i,j,2)=dy./voMy;
QMy(i,j,2)=0.5*(Qy(i,j,2)+Qy(i+1,j,2));
%FMx(i,j,2)=sqrt(qMx(i,j,2).^2./(g.*yMx(i,j,2).^3));
CoMy=5/3.*voMy;
% REFINED
Xy(i,j,2)=1/2-Q3y(i,j,2)./(2.*Soy(i,j).*dx.*CoMy.*dy);
CDy=dtsec+2.*Ky(i,j,2).*(1-Xy(i,j,2));
C1y=(dtsec-2.*Ky(i,j,2).*Xy(i,j,2))./CDy;
C2y=(dtsec+2.*Ky(i,j,1).*Xy(i,j,1))./CDy;
C3y=(-dtsec+2.*Ky(i,j,1).*(1-Xy(i,j,1)))./CDy;
C4y=dtsec./CDy;
% REFINED
Qy(i+1,j,2)=C1y.*Qy(i,j,2)+C2y.*Qy(i,j,1)+C3y.*Qy(i+1,j,1)+C4y.*2.*((Iemps)*dx*dy)-C4y.*2.*(Qx(i,j+1,1)-Qx(i,j,1));
end % if
end % i spatial index y-direction
end % j spatial index x-direction
% UPDATING
% x-direction
Qx(:,:,1)=Qx(:,:,2);
Kx(:,:,1)=Kx(:,:,2);
Xx(:,:,1)=Xx(:,:,2);
% y-direction
Qy(:,:,1)=Qy(:,:,2);
Ky(:,:,1)=Ky(:,:,2);
Xy(:,:,1)=Xy(:,:,2);
% Updating Qxout: Discharge at the d/s of plane in x-direction
Qxout(:,k)=Qx(:,xend+1,2);
end % k temporal index
% Specific discharge at the d/s of x-direction
Qxout=2.*Qxout; % for two symmetrical planes
% Channel Sub-Reach length
dyc=dy;
% Lateral Flow to the Channel
qLc=Qxout./dyc;
qLc=qLc';
% Summation
sumQxout=sum(Qxout);
Qxout=sumQxout;
% Peak Discharge for plane at the Plane outlet in x-direction
Qeq=max(Qxout);
% Converting row vector to --> Column Vector
Qxout=Qxout';
% PLOTTING Plane Outlet
% plot(t,Qxout)
% grid on
% title('Discharge Hydrograph (Plane)')
% xlabel('Time (min)')
% ylabel('Qp (m3/s)')
% legend('VPMM 2D Model')
% Calculating Model Total Plane Volume
modelvolplane=0.5*dtsec*(2*sum(Qxout)-Qxout(1)-Qxout(tend));
% Percentage of actualvolume
pervolplane=modelvolplane/actualvol*100;
%% -------------------CHANNEL ROUTING-----------------------------------
% VPMM Full Channel Routing with Lateral Flow QLc
% Initializing all the variables
Qc=zeros(2,yend+1);
% Discharge at the d/s
Qcout=zeros(tend,1);
% The travel time K
Kc=zeros(2,yend);
% Weighting coefficient theta, i.e. X
Xc=zeros(2,yend);
% normal depth
yoMc=zeros(2,yend);
% Q3 = discharge at section 3 corresponding yM ( = QoM; Normal Discharge)
Q3c=zeros(2,yend);
% The normal flow depth corresponding to Q3
yMc=zeros(2,yend);
% The Discharge QM, at the mid-section of the Muskingum routing reach
QMc=zeros(2,yend);
%Flow depth corresponding to outflow Q(i+1,j+1)
yc=zeros(2,yend);
% Calculation at time > 0, i.e. time level > 1. [Step: 2 To Step: 14]
for i=2:tend % temporal index
for j=1:yend % spatial index
% UNREFINED
Kc(2,j)=Kc(1,j);
Xc(2,j)=Xc(1,j);
CDc=(dtsec+2.*Kc(2,j).*(1-Xc(2,j)));
C1c=(dtsec-2.*Kc(2,j).*Xc(2,j))./CDc;
C2c=(dtsec+2.*Kc(1,j).*Xc(1,j))./CDc;
C3c=(-dtsec+2.*Kc(1,j).*(1-Xc(1,j)))./CDc;
C4c=dtsec./CDc;
% UNREFINED Q
Qc(2,j+1)=C1c*Qc(2,j)+C2c*Qc(1,j)+C3c*Qc(1,j+1)+C4c*dyc*(qLc(i,j)+qLc(i-1,j));
% UNREFINED
Q3c(2,j)=Xc(2,j)*Qc(2,j)+(1-Xc(2,j))*Qc(2,j+1);
% UNREFINED
yMc(2,j)=NR_RectangularChannel(Q3c(2,j),Bc, Soc, nc);
% Calculating the geometric elements
[T,AM,PM]=georect(yMc(2,j),Bc);
voMc=Soc.^(1/2)./nc.*(AM./PM).^(2/3);
% REFINED
Kc(2,j)=dyc./voMc;
QMc(2,j)=0.5*(Qc(2,j)+Qc(2,j+1));
FMc=sqrt(QMc(2,j).^2.*T./(g.*AM.^3));
CoMc=((1+2./3.*PM./T.*(Bc./(Bc + 2.*yMc(2,j)) -(2.*Bc.*yMc(2,j))./(Bc + 2.*yMc(2,j)).^2))).*voMc;
% REFINED Weighting Factor X for VPMM Full equation
Xc(2,j)=1/2-Q3c(2,j).*(1-4./9.*(FMc.^2).*(PM./T.*(Bc./(Bc +2.*yMc(2,j)) - (2.*Bc.*yMc(2,j))./(Bc +2.*yMc(2,j)).^2)).^2)./(2.*Soc.*T.*CoMc.*dy);
CDc=(dtsec+2.*Kc(2,j).*(1-Xc(2,j)));
C1c=(dtsec-2.*Kc(2,j).*Xc(2,j))./CDc;
C2c=(dtsec+2.*Kc(1,j).*Xc(1,j))./CDc;
C3c=(-dtsec+2.*Kc(1,j).*(1-Xc(1,j)))./CDc;
C4c=dtsec./CDc;
% REFINED
Qc(2,j+1)=C1c.*Qc(2,j)+C2c.*Qc(1,j)+C3c.*Qc(1,j+1)+C4c.*dyc.*(qLc(i,j)+qLc(i-1,j));
end % j spatial index
% Updating
Qc(1,:)=Qc(2,:);
Kc(1,:)=Kc(2,:);
Xc(1,:)=Xc(2,:);
% Discharge at the d/s of Lx, i.e 100 Km
Qcout(i)=Qc(2,yend+1);
end % i temporal index
% Peak value for Qc at the channel outlet
Qceq=max(Qcout);
modelvolchannel=0.5*dtsec*(2*sum(Qcout)-Qcout(1)-Qcout(tend));
pervolchannel=modelvolchannel/actualvol*100;

%%
Qcout = (max(FNo12)/max(Qcout))*Qcout;

Blah = max(numel(FNo12),numel(Qcout));
Blah2 = min(numel(FNo12),numel(Qcout));
Blah3 = interp(FNo12,floor(Blah/Blah2));
Blah4 = zeros(Blah-numel(Blah3),1);
FNo12 = [Blah4;Blah3];

blah = max(numel(TNo12),numel(Qcout));
blah2 = min(numel(TNo12),numel(Qcout));
blah3 = interp(TNo12,floor(blah/blah2));
blah4 = zeros(blah-numel(blah3),1);
TNo12 = [blah4;blah3];

DC = DC(FNo12,Qcout);
RMSE = sqrt(mean((FNo12 - Qcout).^2));
% PLOTTING Qc at the Channel Outlet
figure
plot(t,Qcout)
title({'Event No. 12'})
hold on
plot(TNo12,FNo12)
xlabel('Time (sec)')
ylabel('Discharge (lps)')
legend('Modeled','Observed')