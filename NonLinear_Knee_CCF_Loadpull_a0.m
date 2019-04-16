% import Smith Chart
Xsmith = xlsread('Smith_Chart_Excel.xlsx','C2:C14801');
Ysmith = xlsread('Smith_Chart_Excel.xlsx','D2:D14801');


% fft sampling initialisation

fs=10000;
t=0:1/fs:2.0-1/fs;
f=50;
x=2*pi*f*t;

% 10W GaN initialisation

Vdc=28;
Imax=2.33;
Vknee=8;
Idq = 0.02; %add small current from Class-AB bias to avoid >100% DE
Imin = 0;

alpha=0;

% Ideal Current Wavefrom
IgenIdeal= (Imax/2) * ((2/pi) - sin(x) - ((4/pi)*(1/3*cos(2*x) + 1/15*cos(4*x) + 1/35*cos(6*x) + 1/63*cos(8*x) + 1/99*cos(10*x)))) +Idq ;

%Voltage Waveform
VgenCCF_K1_a0=(Vdc-(Vknee/2))*(1+(0.003 - (alpha*cos(x)) - (0.0003*alpha*cos(x)) + (231/200*sin(x)) - (539/800*alpha*sin(2*x)) + (77/400*sin(3*x)) - (77/800*alpha*sin(4*x))))+Vknee;

%FFT Voltage Waveform
FFTVgenCCF_K1_a0 = fft(VgenCCF_K1_a0,[],2)/10000;
VgenCCF_K1_f1 = FFTVgenCCF_K1_a0(:,101);

%set dV
dV=(0.1:0.01:1.5).';  %0.1:0.01:1.5

% set dtheta
dtheta=(-84:2:84);

%Manipulate fundamental component
NewCCF_F1=((abs(VgenCCF_K1_f1))*dV) .* exp(1i*(angle(VgenCCF_K1_f1)+(dtheta/180*pi)));

% flatten matrices above
NewCCF_f1=reshape(NewCCF_F1,[],1);

%Put fundamental component back into full waveform equation
NewVgenCCF_K1_a0= (real(NewCCF_f1)*cos(x)-imag(NewCCF_f1)*sin(x))+(Vdc-(Vknee/2))*(1+(0.003 - (539/800*alpha*sin(2*x)) + (77/400*sin(3*x)) - (77/800*alpha*sin(4*x))))+Vknee;
%take only positive voltage waveform above 0.3V
NewVgenCCF_K1_a0(any(NewVgenCCF_K1_a0<=0.2,2),:)=[];

% scaling for current waveforms
minNewVgenCCF_K1_a0=min(NewVgenCCF_K1_a0,[],2);
Igenscale_K1_a0(1,:)=(minNewVgenCCF_K1_a0<Vknee) .* (abs(minNewVgenCCF_K1_a0-Vknee)) + (minNewVgenCCF_K1_a0>=Vknee) .* (0.01);

%Imax Scale
ImaxScale(1,:)= (minNewVgenCCF_K1_a0<=Vknee) .* ((((Imax-Imin)/Vknee)*minNewVgenCCF_K1_a0)+Imin) + (minNewVgenCCF_K1_a0>=Vknee) .* (Imax) ;

% New current waveforms
NewIgenCCF_K1_a0=((1-exp(-NewVgenCCF_K1_a0 ./ Igenscale_K1_a0.')) .* ((ImaxScale.'/2) * ((2/pi) - sin(x) - ((4/pi)*(1/3*cos(2*x) + 1/15*cos(4*x) + 1/35*cos(6*x) + 1/63*cos(8*x) + 1/99*cos(10*x)))) +Idq ));



% FFT Voltage Waveforms
FFTNewVgenCCF_K1_a0 = fft(NewVgenCCF_K1_a0,[],2)/10000;
NewVgenCCF_K1_a0_f0 = FFTNewVgenCCF_K1_a0(:,1)/2;
NewVgenCCF_K1_a0_f1 = FFTNewVgenCCF_K1_a0(:,101);
NewVgenCCF_K1_a0_f2 = FFTNewVgenCCF_K1_a0(:,201);
NewVgenCCF_K1_a0_f3 = FFTNewVgenCCF_K1_a0(:,301);


% FFT Current Waveforms

FFTNewIgenCCF_K1_a0 = fft(NewIgenCCF_K1_a0,[],2)/10000;
NewIgenCCF_K1_a0_f0 = FFTNewIgenCCF_K1_a0(:,1)/2;
NewIgenCCF_K1_a0_f1 = FFTNewIgenCCF_K1_a0(:,101);
NewIgenCCF_K1_a0_f2 = FFTNewIgenCCF_K1_a0(:,201);
NewIgenCCF_K1_a0_f3 = FFTNewIgenCCF_K1_a0(:,301);


% Calculate performances
PoutCCF_K1_a0_W = real(NewIgenCCF_K1_a0_f1 .* NewVgenCCF_K1_a0_f1)/2;
PoutCCF_K1_a0_dBm = real((10*log10(PoutCCF_K1_a0_W))+30);
PdcCCF_K1_a0 = real(NewIgenCCF_K1_a0_f0 .* NewVgenCCF_K1_a0_f0);
EfficiencyCCF_K1_a0 = PoutCCF_K1_a0_W ./ PdcCCF_K1_a0 * 100;


% Calculate Impedances

Zref=50;
%fundamental
ZNewCCF_K1_a0_f1= -(NewVgenCCF_K1_a0_f1./NewIgenCCF_K1_a0_f1);
Mag_NewCCF_Kl_a0_f1 = abs((ZNewCCF_K1_a0_f1-Zref)./(ZNewCCF_K1_a0_f1+Zref));
Phase_Mag_NewCCF_Kl_a0_f1 = angle((ZNewCCF_K1_a0_f1-Zref)./(ZNewCCF_K1_a0_f1+Zref));

%2nd harmonic
ZNewCCF_K1_a0_f2= -(NewVgenCCF_K1_a0_f2./NewIgenCCF_K1_a0_f2);
Mag_NewCCF_Kl_a0_f2 = abs((ZNewCCF_K1_a0_f2-Zref)./(ZNewCCF_K1_a0_f2+Zref));
Phase_Mag_NewCCF_Kl_a0_f2 = angle((ZNewCCF_K1_a0_f2-Zref)./(ZNewCCF_K1_a0_f2+Zref));

%3rd harmonic
ZNewCCF_K1_a0_f3= -(NewVgenCCF_K1_a0_f3./NewIgenCCF_K1_a0_f3);
Mag_NewCCF_Kl_a0_f3 = abs((ZNewCCF_K1_a0_f3-Zref)./(ZNewCCF_K1_a0_f3+Zref));
Phase_Mag_NewCCF_Kl_a0_f3 = angle((ZNewCCF_K1_a0_f3-Zref)./(ZNewCCF_K1_a0_f3+Zref));

% Convert to cartesian 
XNewCCF_Kl_a0_f1=Mag_NewCCF_Kl_a0_f1.*cos(Phase_Mag_NewCCF_Kl_a0_f1);
YNewCCF_Kl_a0_f1=Mag_NewCCF_Kl_a0_f1.*sin(Phase_Mag_NewCCF_Kl_a0_f1);

XNewCCF_Kl_a0_f2=Mag_NewCCF_Kl_a0_f2.*cos(Phase_Mag_NewCCF_Kl_a0_f2);
YNewCCF_Kl_a0_f2=Mag_NewCCF_Kl_a0_f2.*sin(Phase_Mag_NewCCF_Kl_a0_f2);

XNewCCF_Kl_a0_f3=Mag_NewCCF_Kl_a0_f3.*cos(Phase_Mag_NewCCF_Kl_a0_f3);
YNewCCF_Kl_a0_f3=Mag_NewCCF_Kl_a0_f3.*sin(Phase_Mag_NewCCF_Kl_a0_f3);



% Plotting Contour
[xi, yi] = meshgrid(linspace(min(XNewCCF_Kl_a0_f1),max(XNewCCF_Kl_a0_f1)),linspace(min(YNewCCF_Kl_a0_f1),max(YNewCCF_Kl_a0_f1)));
PoutContour_W = griddata(XNewCCF_Kl_a0_f1,YNewCCF_Kl_a0_f1,PoutCCF_K1_a0_W,xi,yi);
PoutContour_dBm = griddata(XNewCCF_Kl_a0_f1,YNewCCF_Kl_a0_f1,PoutCCF_K1_a0_dBm,xi,yi);
EfficiencyContour = griddata(XNewCCF_Kl_a0_f1,YNewCCF_Kl_a0_f1,EfficiencyCCF_K1_a0,xi,yi);





figure
plot(Xsmith,Ysmith,'k')
hold on
[C,h]=contour(xi,yi,PoutContour_dBm,max(PoutCCF_K1_a0_dBm)-10:1:max(PoutCCF_K1_a0_dBm),'blue','ShowText','on');
clabel(C,h,'FontSize',10,'Color','blue','LabelSpacing',300)
hold on
[C,h]=contour(xi,yi,EfficiencyContour,max(EfficiencyCCF_K1_a0)-80:4:max(EfficiencyCCF_K1_a0),'red','ShowText','on');
clabel(C,h,'FontSize',10,'Color','red','LabelSpacing',300)
hold off
title('Pout and DE Contours')

figure
plot(Xsmith,Ysmith,'k')
hold on
[C,h]=contour(xi,yi,PoutContour_dBm,max(PoutCCF_K1_a0_dBm)-10:1:max(PoutCCF_K1_a0_dBm),'blue','ShowText','on');
clabel(C,h,'FontSize',10,'Color','blue','LabelSpacing',300)
hold off
title('Pout Contours')

figure
plot(Xsmith,Ysmith,'k')
hold on
[C,h]=contour(xi,yi,EfficiencyContour,max(EfficiencyCCF_K1_a0)-30:2:max(EfficiencyCCF_K1_a0),'red','ShowText','on');
clabel(C,h,'FontSize',10,'Color','red','LabelSpacing',300)
hold off
title('DE Contours')

%Plot calculated f1, f2 and f3 impedances
figure
plot(Xsmith,Ysmith,'k')
hold on
scatter(XNewCCF_Kl_a0_f1,YNewCCF_Kl_a0_f1,'blue','filled')
hold off
title('Fundamental Impedances')

figure
plot(Xsmith,Ysmith,'k')
hold on
scatter(XNewCCF_Kl_a0_f2,YNewCCF_Kl_a0_f2,'red','filled')
hold on
axis([-1.5 1.5 -1.5 1.5])
hold off
title('2nd Harmonic Impedances')

figure
plot(Xsmith,Ysmith,'k')
hold on
scatter(XNewCCF_Kl_a0_f3,YNewCCF_Kl_a0_f3,'green','filled')
hold on
axis([-1.5 1.5 -1.5 1.5])
hold off
title('3rd Harmonic Impedances')

figure
scatter(PoutCCF_K1_a0_dBm,EfficiencyCCF_K1_a0,'k','filled',...
    'SizeData',20)
axis([30 43 0 100])
xticks(30:1:43)
hold off
title('OPBO')




BO_0dBm = max(PoutCCF_K1_a0_dBm);
DE_BO_0dBm = max(EfficiencyCCF_K1_a0(PoutCCF_K1_a0_dBm < (BO_0dBm+0.05) & PoutCCF_K1_a0_dBm > (BO_0dBm-0.05)));

BO_1dBm = max(PoutCCF_K1_a0_dBm)-1;
DE_BO_1dBm = max(EfficiencyCCF_K1_a0(PoutCCF_K1_a0_dBm < (BO_1dBm+0.1) & PoutCCF_K1_a0_dBm > (BO_1dBm-0.1)));

BO_2dBm = max(PoutCCF_K1_a0_dBm)-2;
DE_BO_2dBm = max(EfficiencyCCF_K1_a0(PoutCCF_K1_a0_dBm < (BO_2dBm+0.1) & PoutCCF_K1_a0_dBm > (BO_2dBm-0.1)));

BO_3dBm = max(PoutCCF_K1_a0_dBm)-3;
DE_BO_3dBm = max(EfficiencyCCF_K1_a0(PoutCCF_K1_a0_dBm < (BO_3dBm+0.1) & PoutCCF_K1_a0_dBm > (BO_3dBm-0.1)));

BO_4dBm = max(PoutCCF_K1_a0_dBm)-4;
DE_BO_4dBm = max(EfficiencyCCF_K1_a0(PoutCCF_K1_a0_dBm < (BO_4dBm+0.1) & PoutCCF_K1_a0_dBm > (BO_4dBm-0.1)));

BO_5dBm = max(PoutCCF_K1_a0_dBm)-5;
DE_BO_5dBm = max(EfficiencyCCF_K1_a0(PoutCCF_K1_a0_dBm < (BO_5dBm+0.1) & PoutCCF_K1_a0_dBm > (BO_5dBm-0.1)));

BO_6dBm = max(PoutCCF_K1_a0_dBm)-6;
DE_BO_6dBm = max(EfficiencyCCF_K1_a0(PoutCCF_K1_a0_dBm < (BO_6dBm+0.1) & PoutCCF_K1_a0_dBm > (BO_6dBm-0.1)));

BO_7dBm = max(PoutCCF_K1_a0_dBm)-7;
DE_BO_7dBm = max(EfficiencyCCF_K1_a0(PoutCCF_K1_a0_dBm < (BO_7dBm+0.1) & PoutCCF_K1_a0_dBm > (BO_7dBm-0.1)));

BO_9dBm = max(PoutCCF_K1_a0_dBm)-9;
DE_BO_9dBm = max(EfficiencyCCF_K1_a0(PoutCCF_K1_a0_dBm < (BO_9dBm+0.1) & PoutCCF_K1_a0_dBm > (BO_9dBm-0.1)));

BO_10dBm = max(PoutCCF_K1_a0_dBm)-10;
DE_BO_10dBm = max(EfficiencyCCF_K1_a0(PoutCCF_K1_a0_dBm < (BO_10dBm+0.1) & PoutCCF_K1_a0_dBm > (BO_10dBm-0.1)));

BO_11dBm = max(PoutCCF_K1_a0_dBm)-11;
DE_BO_11dBm = max(EfficiencyCCF_K1_a0(PoutCCF_K1_a0_dBm < (BO_11dBm+0.1) & PoutCCF_K1_a0_dBm > (BO_11dBm-0.1)));

BO_8dBm = max(PoutCCF_K1_a0_dBm)-8;
DE_BO_8dBm = max(EfficiencyCCF_K1_a0(PoutCCF_K1_a0_dBm < (BO_8dBm+0.1) & PoutCCF_K1_a0_dBm > (BO_8dBm-0.1)));

DE_BO_8dBm_index = find(abs(EfficiencyCCF_K1_a0-DE_BO_8dBm) < 0.001);
X_DE_BO_8dBm = XNewCCF_Kl_a0_f1(DE_BO_8dBm_index);
Y_DE_BO_8dBm = YNewCCF_Kl_a0_f1(DE_BO_8dBm_index);

DE_BO_6dBm_index = find(abs(EfficiencyCCF_K1_a0-DE_BO_6dBm) < 0.001);
X_DE_BO_6dBm = XNewCCF_Kl_a0_f1(DE_BO_6dBm_index);
Y_DE_BO_6dBm = YNewCCF_Kl_a0_f1(DE_BO_6dBm_index);

DE_BO_4dBm_index = find(abs(EfficiencyCCF_K1_a0-DE_BO_4dBm) < 0.001);
X_DE_BO_4dBm = XNewCCF_Kl_a0_f1(DE_BO_4dBm_index);
Y_DE_BO_4dBm = YNewCCF_Kl_a0_f1(DE_BO_4dBm_index);

DE_BO_2dBm_index = find(abs(EfficiencyCCF_K1_a0-DE_BO_2dBm) < 0.001);
X_DE_BO_2dBm = XNewCCF_Kl_a0_f1(DE_BO_2dBm_index);
Y_DE_BO_2dBm = YNewCCF_Kl_a0_f1(DE_BO_2dBm_index);

Max_DE_index = find(EfficiencyCCF_K1_a0 == max(EfficiencyCCF_K1_a0));
Max_Pout_index = find(PoutCCF_K1_a0_dBm == max(PoutCCF_K1_a0_dBm));


%Plot calculated f1, f2 and f3 impedances
figure
plot(Xsmith,Ysmith,'k')
hold on
scatter(XNewCCF_Kl_a0_f1(Max_Pout_index),YNewCCF_Kl_a0_f1(Max_Pout_index),'blue','filled')
hold on
scatter(XNewCCF_Kl_a0_f1(Max_DE_index),YNewCCF_Kl_a0_f1(Max_DE_index),'red','filled')
hold on
scatter(X_DE_BO_8dBm,Y_DE_BO_8dBm,'green','filled')
hold on
scatter(X_DE_BO_6dBm,Y_DE_BO_6dBm,'green','filled')
hold on
scatter(X_DE_BO_4dBm,Y_DE_BO_4dBm,'green','filled')
hold on
scatter(X_DE_BO_2dBm,Y_DE_BO_2dBm,'green','filled')
hold off
title('8dB Back-off DE Impedances')


Pout_BO = [BO_11dBm;BO_10dBm;BO_9dBm;BO_8dBm;BO_7dBm;BO_6dBm;BO_5dBm;BO_4dBm;BO_3dBm;BO_2dBm;BO_1dBm;BO_0dBm];
DE_BO = [DE_BO_11dBm;DE_BO_10dBm;DE_BO_9dBm;DE_BO_8dBm;DE_BO_7dBm;DE_BO_6dBm;DE_BO_5dBm;DE_BO_4dBm;DE_BO_3dBm;DE_BO_2dBm;DE_BO_1dBm;DE_BO_0dBm];

NewPout_BO = linspace(min(Pout_BO),max(Pout_BO),100);
F_DE_BO = griddedInterpolant(Pout_BO,DE_BO,'spline');
DE_BO_Interp = F_DE_BO(NewPout_BO);


%Interp 8dB BO DE
DE_BO_8dBm_interp = max(DE_BO_Interp(NewPout_BO < (BO_8dBm+0.05) & NewPout_BO > (BO_8dBm-0.05)));


figure
scatter(PoutCCF_K1_a0_dBm,EfficiencyCCF_K1_a0,'k','filled',...
    'SizeData',20)
hold on
plot(NewPout_BO,DE_BO_Interp)
hold on
scatter(BO_8dBm,DE_BO_8dBm_interp,'green')
axis([30 43 0 100])
xticks(30:1:43)
hold off
title('OPBO')



figure
plot(Xsmith,Ysmith,'k','HandleVisibility','off')
hold on
contour(xi,yi,PoutContour_dBm,max(PoutCCF_K1_a0_dBm)-8:2:max(PoutCCF_K1_a0_dBm),'b','LineWidth',1)
hold on
contour(xi,yi,EfficiencyContour,max(EfficiencyCCF_K1_a0)-40:10:max(EfficiencyCCF_K1_a0),'r:','LineWidth',2)
hold on
scatter(XNewCCF_Kl_a0_f1(Max_Pout_index),YNewCCF_Kl_a0_f1(Max_Pout_index),'b','filled')
hold on
scatter(XNewCCF_Kl_a0_f1(Max_DE_index),YNewCCF_Kl_a0_f1(Max_DE_index),'r','filled')
hold on
scatter(X_DE_BO_8dBm(1),0,'green','filled')
hold off
legend('Pout Contours','DE Contours','Max Pout (42 dBm)','Max DE (74%)')



%ADDING QUARTER WAVELENGTH 


New_Phase_Mag_NewCCF_Kl_a0_f1 = deg2rad(rad2deg(Phase_Mag_NewCCF_Kl_a0_f1)+180);

New_Phase_Mag_NewCCF_Kl_a0_f2 = deg2rad(rad2deg(Phase_Mag_NewCCF_Kl_a0_f2)+180);

New_Phase_Mag_NewCCF_Kl_a0_f3 = deg2rad(rad2deg(Phase_Mag_NewCCF_Kl_a0_f3)+180);


New_XNewCCF_Kl_a0_f1=Mag_NewCCF_Kl_a0_f1.*cos(New_Phase_Mag_NewCCF_Kl_a0_f1);
New_YNewCCF_Kl_a0_f1=Mag_NewCCF_Kl_a0_f1.*sin(New_Phase_Mag_NewCCF_Kl_a0_f1);

New_XNewCCF_Kl_a0_f2=Mag_NewCCF_Kl_a0_f2.*cos(New_Phase_Mag_NewCCF_Kl_a0_f2);
New_YNewCCF_Kl_a0_f2=Mag_NewCCF_Kl_a0_f2.*sin(New_Phase_Mag_NewCCF_Kl_a0_f2);

New_XNewCCF_Kl_a0_f3=Mag_NewCCF_Kl_a0_f3.*cos(New_Phase_Mag_NewCCF_Kl_a0_f3);
New_YNewCCF_Kl_a0_f3=Mag_NewCCF_Kl_a0_f3.*sin(New_Phase_Mag_NewCCF_Kl_a0_f3);


[New_xi, New_yi] = meshgrid(linspace(min(New_XNewCCF_Kl_a0_f1),max(New_XNewCCF_Kl_a0_f1)),linspace(min(New_YNewCCF_Kl_a0_f1),max(New_YNewCCF_Kl_a0_f1)));
New_PoutContour_W = griddata(New_XNewCCF_Kl_a0_f1,New_YNewCCF_Kl_a0_f1,PoutCCF_K1_a0_W,New_xi,New_yi);
New_PoutContour_dBm = griddata(New_XNewCCF_Kl_a0_f1,New_YNewCCF_Kl_a0_f1,PoutCCF_K1_a0_dBm,New_xi,New_yi);
New_EfficiencyContour = griddata(New_XNewCCF_Kl_a0_f1,New_YNewCCF_Kl_a0_f1,EfficiencyCCF_K1_a0,New_xi,New_yi);

%8dB OPBO
New_X_DE_BO_8dBm = New_XNewCCF_Kl_a0_f1(DE_BO_8dBm_index);
New_Y_DE_BO_8dBm = New_YNewCCF_Kl_a0_f1(DE_BO_8dBm_index);

figure
plot(Xsmith,Ysmith,'k','HandleVisibility','off')
hold on
contour(New_xi,New_yi,New_PoutContour_dBm,max(PoutCCF_K1_a0_dBm)-8:2:max(PoutCCF_K1_a0_dBm),'b','LineWidth',1)
hold on
contour(New_xi,New_yi,New_EfficiencyContour,max(EfficiencyCCF_K1_a0)-40:10:max(EfficiencyCCF_K1_a0),'r:','LineWidth',2)
hold on
scatter(New_XNewCCF_Kl_a0_f1(Max_Pout_index),New_YNewCCF_Kl_a0_f1(Max_Pout_index),'bd','filled')
hold on
scatter(New_XNewCCF_Kl_a0_f1(Max_DE_index),New_YNewCCF_Kl_a0_f1(Max_DE_index),'rd','filled')
hold on
scatter(New_X_DE_BO_8dBm(1),0,'gd','filled')
hold off
legend('Pout Contours','DE Contours','Max Pout (42 dBm)','Max DE (74%)')




% Combine normal PA with inverted PA
figure
plot(Xsmith,Ysmith,'k','HandleVisibility','off')
hold on
contour(xi,yi,PoutContour_dBm,max(PoutCCF_K1_a0_dBm)-8:2:max(PoutCCF_K1_a0_dBm),'b','LineWidth',2)
hold on
contour(New_xi,New_yi,New_PoutContour_dBm,max(PoutCCF_K1_a0_dBm)-8:2:max(PoutCCF_K1_a0_dBm),'b--','LineWidth',2)
hold on
contour(xi,yi,EfficiencyContour,max(EfficiencyCCF_K1_a0)-40:10:max(EfficiencyCCF_K1_a0),'r','LineWidth',2)
hold on
contour(New_xi,New_yi,New_EfficiencyContour,max(EfficiencyCCF_K1_a0)-40:10:max(EfficiencyCCF_K1_a0),'r--','LineWidth',2)
hold on
scatter(XNewCCF_Kl_a0_f1(Max_Pout_index),YNewCCF_Kl_a0_f1(Max_Pout_index),'b','filled','SizeData',60)
hold on
scatter(XNewCCF_Kl_a0_f1(Max_DE_index),YNewCCF_Kl_a0_f1(Max_DE_index),'r','filled','SizeData',60)
hold on
scatter(X_DE_BO_8dBm(1),0,'green','filled','SizeData',60)
hold on
scatter(New_XNewCCF_Kl_a0_f1(Max_Pout_index),New_YNewCCF_Kl_a0_f1(Max_Pout_index),'bd','filled','SizeData',70)
hold on
scatter(New_XNewCCF_Kl_a0_f1(Max_DE_index),New_YNewCCF_Kl_a0_f1(Max_DE_index),'rd','filled','SizeData',70)
hold on
scatter(New_X_DE_BO_8dBm(1),0,'gd','filled','SizeData',70)
hold off
axis([-1 1 -1 1])
