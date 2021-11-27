%----------------------------------------------%
%           *** SSSP - HOMEWORK #3 ***         %
%----------------------------------------------%
%                WDF modeling                  %
%----------------------------------------------%
% Giovanni Affatato, Roberto Alessandri        %
%----------------------------------------------%

clear all
close all
clc

%% Import Input Audio Signal
[Vin,~] = audioread('ExpSweep.wav');

%% LTSpice Files for Ground-Truth
[OutLowSpice,~]=audioread('outlowsweep.wav');
[OutMidSpice,~]=audioread('outmidsweep.wav');
[OutHighSpice,FsLTSpice]=audioread('outhighsweep.wav');
TsLTSpice=1/FsLTSpice;

%% Sampling frequency (to be varied: FsLTSpice/downSampFact, with downSampFact={4,3,2})
downSampFact=2;
Fs =FsLTSpice/downSampFact; 

%% Downsample Input Signal
Vin=Vin([1:downSampFact:end]);

%% Sampling Period
Ts=1/Fs;
%% Number of Samples
Nsamp=length(Vin);
%% Simulated time
tstop=Nsamp*Ts;
%% Parameters of Dynamic Element
L1 = 0.5*10^(-3);
% Other Inductive parameters...
L2 = 0.2*10^(-3);
L3 = 10^(-3);
L4 = 0.5*10^(-3);

C1 = 270*10^(-6);
% Other Capacitive parameters...
C2 = 15*10^(-6);
C3 = 10*10^(-6);
C4 = 0.8*10^(-6);
C5 = 10^(-6);

%% Resistive Parameters
R1=4;
RspkLow=6;
% Other Resistive parameters...
R2 = 10;
R3 = 1.2;
RspkMid = 8;
RspkHigh = 8;


%% WDF setting of free parameters (adaptation conditions)
%Low Pass
ZLP1 = 2*L1 / Ts;
ZLP9 = Ts / (2*C1);
ZLP8 = R1;
ZLP6 = RspkLow;
ZLP7 = ZLP8 + ZLP9;
ZLP5 = ZLP7;
ZLP4 = ZLP5*ZLP6/(ZLP5+ZLP6);
ZLP3 = ZLP4;
ZLP2 = ZLP1 + ZLP3;

%Mid Pass
ZBP10 = RspkMid;
ZBP9 = Ts / (2*C4);
ZBP8 = Ts / (2*C3);
ZBP7 = 2*L3 / Ts;
ZBP4 = 2*L2 / Ts;
ZBP3 = Ts / (2*C2);
ZBP2 = R2;
ZBP6 = (1/ZBP7 + 1/ZBP8 + 1/ZBP9 + 1/ZBP10)^-1;
ZBP5 = ZBP6;
ZBP1 = ZBP2 + ZBP3 + ZBP4 + ZBP5;

%High Pass
ZHP1 = R3;
ZHP6 = RspkHigh;
ZHP4 = Ts / (2*C5);
ZHP7 = 2*L4 / Ts;
ZHP5 = (ZHP6*ZHP7) / (ZHP6 + ZHP7);
ZHP3 = ZHP5;
ZHP2 = ZHP1 + ZHP4 + ZHP3;

%% Computation of Scattering Matrices
NportLP = 3;

%LOW PASS
%Junction 1-2-3
ZJLP1 = diag([ZLP1, ZLP2, ZLP3]);
BLP1 = [1; 1; 1]';
SLP1 = eye(NportLP) - ((2*ZJLP1*BLP1')/(BLP1*ZJLP1*BLP1'))*BLP1;


%Junction 4-5-6
ZJLP2 = diag([ZLP4, ZLP5, ZLP6]);
QLP2 = [1; 1; 1]';
SLP2 = ((2*QLP2')/(QLP2/ZJLP2*QLP2'))*QLP2/ZJLP2 - eye(NportLP);

%Junction 7-8-9
ZJLP3 = diag([ZLP7, ZLP8, ZLP9]);
BLP3 = [1;1;1]';
SLP3 = eye(NportLP) - ((2*ZJLP3*BLP3')/(BLP3*ZJLP3*BLP3'))*BLP3;

%MID PASS
NportBP = 5;
%Junction 1-2-3-4-5
ZJBP1 = diag([ZBP1, ZBP2, ZBP3, ZBP4, ZBP5]);
BBP1 = [1; 1; 1; 1; 1]';
SBP1 = eye(NportBP) - ((2*ZJBP1*BBP1')/(BBP1*ZJBP1*BBP1'))*BBP1;


%Junction 1-2-3-4-5
ZJBP2 = diag([ZBP6, ZBP7, ZBP8, ZBP9, ZBP10]);
QBP2 = [1; 1; 1; 1; 1]';
SBP2 = ((2*QBP2')/(QBP2/ZJBP2*QBP2'))*QBP2/ZJBP2 - eye(NportBP);

%HIGH PASS
%Juction 1-2-3-4
ZJHP1 = diag([ZHP1, ZHP2, ZHP3, ZHP4]);
BHP1 = [1; 1; 1; 1]';
SHP1 = eye(4) - (2*ZJHP1*BHP1')*((BHP1*ZJHP1*BHP1')\BHP1);

%Juction 5-6-7
ZJHP2 = diag([ZHP5, ZHP6, ZHP7]);
QHP2 = [1; 1; 1]';
SHP2 = 2*QHP2'*((QHP2*(ZJHP2\QHP2'))\(QHP2/(ZJHP2))) - eye(3);



%% Initialization of Waves
aLp = zeros(9,1);
bLp = zeros(9,1);

aBp = zeros(10,1);
bBp = zeros(10,1);

aHp = zeros(7,1);
bHp = zeros(7,1);


%% Initialize Output Signals
% Low
VoutLow=zeros(size(Vin));
% Mid
VoutMid=zeros(size(Vin));
% High
VoutHigh=zeros(size(Vin));

ii=0;
while (ii<Nsamp)
    ii=ii+1;

    %% Manage Dynamic Elements
    %LOW PASS
    aLp(1) = -bLp(1); %Inductance
    aLp(9) = bLp(9); %Capacitance
    
    %MID PASS
    aBp(4) = -bBp(4);
    aBp(7) = -bBp(7);
    aBp(9) = bBp(9);
    aBp(8) = bBp(8);
    aBp(3) = bBp(3);
    
    %HIGH PASS
    aHp(7) = -bHp(7);
    aHp(4) = bHp(4);
    
    
    %% Forward Scan
    %LOW PASS
    %First Layer
    bLp(7:9) = SLP3*aLp(7:9);
    aLp(5) = bLp(7); 
    bLp(4:6) = SLP2*aLp(4:6);
    %Secon Layer
    aLp(3) = bLp(4);
    bLp(1:3) = SLP1*aLp(1:3);
    
    %MID PASS
    bBp(6:10) = SBP2*aBp(6:10);
    aBp(5) = bBp(6);
    bBp(1:5) = SBP1*aBp(1:5);
    
    %HIGH PASS
    bHp(5:7) = SHP2*aHp(5:7);
    aHp(3) = bHp(5);
    bHp(1:4) = SHP1*aHp(1:4);

    %% Local Root Scattering
    %LOW PASS
    aLp(2) = 2 * Vin(ii) - bLp(2);
    
    %MID PASS
    aBp(1) = 2 * Vin(ii) - bBp(1);
    
    %HIGH PASS
    aHp(2) = 2 * Vin(ii) - bHp(2);

    %% Backward Scan
    %LOW PASS
    bLp(1:3) = SLP1*aLp(1:3);
    aLp(4) = bLp(3);
    bLp(4:6) = SLP2*aLp(4:6);
    aLp(7) = bLp(5);
    bLp(7:9) = SLP3*aLp(7:9);
    
    %MID PASS
    bBp(1:5) = SBP1*aBp(1:5);
    aBp(6) = bBp(5);
    bBp(6:10) = SBP2*aBp(6:10);
    
    %HIGH PASS
    bHp(1:4) = SHP1*aHp(1:4);
    aHp(5) = bHp(3);
    bHp(5:7) = SHP2*aHp(5:7);

    %% Read Output
    %LOW PASS
    VoutLow(ii) = -(bLp(6) + aLp(6)) / 2;
    
    %MID PASS
    VoutMid(ii) = -(bBp(10) + aBp(10)) / 2;
    
    %HIGH PASS
    VoutHigh(ii) = -(bHp(6) + aHp(6)) / 2;

end


%% Output Plots
figure
set(gcf, 'Color', 'w');
subplot(311)
plot(TsLTSpice*[1:length(OutLowSpice)],OutLowSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutLow,'b--','Linewidth',1); grid on; xlim([0,tstop]); 
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outLow}}$ [V]','Fontsize',16,'interpreter','latex');
legend('LTspice','WDF','Fontsize',16,'interpreter','latex');
title('Output Signals','Fontsize',18,'interpreter','latex');
subplot(312)
plot(TsLTSpice*[1:length(OutMidSpice)],OutMidSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutMid,'b--','Linewidth',1); grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outMid}}$ [V]','Fontsize',16,'interpreter','latex');
subplot(313)
plot(TsLTSpice*[1:length(OutHighSpice)],OutHighSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutHigh,'b--','Linewidth',1); grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outHigh}}$ [V]','Fontsize',16,'interpreter','latex');

%% Error Plots
figure
set(gcf, 'Color', 'w');
subplot(311)
plot(Ts*[1:Nsamp],OutLowSpice([1:downSampFact:end])-VoutLow,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outLow}}$ [V]','Fontsize',16,'interpreter','latex');
title(['Error Signals. $F_{\mathrm{s}}=$ ',num2str(Fs),' Hz'],'Fontsize',18,'interpreter','latex');
subplot(312)
plot(Ts*[1:Nsamp],OutMidSpice([1:downSampFact:end])-VoutMid,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outMid}}$ [V]','Fontsize',16,'interpreter','latex');
subplot(313)
plot(Ts*[1:Nsamp],OutHighSpice([1:downSampFact:end])-VoutHigh,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outHigh}}$ [V]','Fontsize',16,'interpreter','latex');