close all

dragFile = "/media/frederk/Dump drive/Simulations/Cylinder extrap Re200 Ma0.25 800x402x3/output/drag.dat";
intervalStart = 17000;
dt = 0.0068;

drag_all = readmatrix(dragFile);
drag = drag_all(intervalStart:end);
t=[1:length(drag)]*dt;
Fs=1/dt;
drag_dft=fft(drag);
freq = 0:Fs/length(drag):Fs/2;
drag_dft = drag_dft(1:length(drag)/2+1);
figure(1);
semilogx(freq, abs(drag_dft))
title("Drag");

liftFile = "/media/frederk/Dump drive/Simulations/Cylinder extrap Re200 Ma0.25 800x402x3/output/lift.dat";

lift_all = readmatrix(liftFile);
lift = lift_all(intervalStart:end);
t=[1:length(lift)]*dt;
Fs=1/dt;
lift_dft=fft(lift);
freq = 0:Fs/length(lift):Fs/2;
lift_dft = lift_dft(1:length(lift)/2+1);
figure(2);
semilogx(freq, abs(lift_dft))
title("Lift");
