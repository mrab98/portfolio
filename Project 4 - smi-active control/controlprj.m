
%---------------------------------------------------------------------
%  Program name:Control of SDOF system using magnetorheological damper
%      written by Mehrab Zamanain 
%      Email:mehrab.zamanian@ut.ac.ir
%      July,2021
%---------------------------------------------------------------------
clc
clear

%----structure specification(mass,damping coefficient & stiffnes)-----
m=5873;%kg
c=0;
k=29127; %N/cm
initial_disp=0.01;%m
initial_velocity=0.1;%m/s

%----Solver options(ODE 45)---------------------------------------------------
start_time=0;
stop_time=40;
min_step=0.001;
%----Excitation type(1 for earthquke load & -1 for Chirp signal load)----------------
switchF=1;

%----determine gravitational acceleration(g)--------------------------
g=9.80665; %m/s2

%----Chirp signal specification---------------------------------------
target_t=30;
initial_freq=0;
freq_target=2;
amp=2000;

%----time steps (for plotting quantities graph &variable size) --------------------------
delta_t=0.001;
time=0:delta_t:stop_time;

%----MR damper specifications-----------------------------------------
k0=46.9e2;    %N/m
gama=363e-4;   %m-2
beta=363e-4;   %m-2
noo=190;    %s-1
n=2;
A=301;
k1=5e2;        %N/m
alpha_a=140e2; %N/m
alpha_b=695e2; %N/m.V
c0_a=21e2;     %Ns/m
c0_b=3.5e2;    %Ns/m
c1_a=283e2;    %Ns/m
c1_b=2.95e2;   %Ns/m

%----Actuator specifications-----------------------------------------
alpha=5000;

%----Run simulink models----------------------------------------
 sim('mrdamper.slx');
 sim('uncontroled.slx');
 sim('actuator.slx');
 
%----plotting quantities graph---------------------------------------
subplot(3,1,1)
 plot(time,displacement,'k')
 hold on
 plot(time,displacement1,'r')
 plot(time,displacement2,'b')
 xlabel('time(sec)')
 ylabel('displacement(m)')
 legend({'controled using MR damper(black)','uncontroled(red)','controled using actuator(blue)'},'Location','best')
 
subplot(3,1,2)
 plot(time,velocity,'k')
 hold on
 plot(time,velocity1,'r')
 plot(time,velocity2,'b')
 xlabel('time(sec)')
 ylabel('velocity(m/s)')
 
subplot(3,1,3)
 plot(time,acceleration,'k')
 hold on
 plot(time,acceleration1,'r')
 plot(time,acceleration2,'b')
 xlabel('time(sec)')
 ylabel('acceleration(m/s2)')
 
 %----print RMS and maximum of quantities---------------------------

 
 formatSpec1 = 'The RMS of acceleration of structure using MR damper is %g m/s2 and the maximum is %g m/s2\n';
 formatSpec2 = 'The RMS of velocity of structure using MR damper is %g m/s and the maximum is %g m/s \n';
 formatSpec3 = 'The RMS of displacement of structure using MR damper is %g m and the maximum is %g m\n';
  formatSpec4 = 'The RMS of acceleration of uncontroled structure is %g m/s2 and the maximum is %g m/s2\n';
 formatSpec5 = 'The RMS of velocity of uncontroled structure is %g m/s and the maximum is %g m/s \n';
 formatSpec6 = 'The RMS of displacement of uncontroled structure is %g m and the maximum is %g m\n';
  formatSpec7 = 'The RMS of acceleration of structure using actuator is %g m/s2 and the maximum is %g m/s2 \n\n';
 formatSpec8 = 'The RMS of velocity of structure using actuator is %g m/s and the maximum is %g m/s \n\n';
 formatSpec9 = 'The RMS of displacement of structure using actuator is %g m and the maximum is %g m\n\n';

   
 accrms=rms(acceleration);
 maxacc=max(abs(acceleration));
 velrms=rms(velocity);
 maxvel=max(abs(velocity));
 disprms=rms(displacement);
 maxdisp=max(abs(displacement)); 
 
 accrms1=rms(acceleration1);
 maxacc1=max(abs(acceleration1));
 velrms1=rms(velocity1);
 maxvel1=max(abs(velocity1));
 disprms1=rms(displacement1);
 maxdisp1=max(abs(displacement1));
 
 accrms2=rms(acceleration2);
 maxacc2=max(abs(acceleration2));
 velrms2=rms(velocity2);
 maxvel2=max(abs(velocity2));
 disprms2=rms(displacement2);
 maxdisp2=max(abs(displacement2));
 
    fprintf(formatSpec1,accrms, maxacc)
   fprintf(formatSpec4,accrms1, maxacc1)
   fprintf(formatSpec7,accrms2, maxacc2)
   
   fprintf(formatSpec2,velrms,maxvel)
   fprintf(formatSpec5,velrms1,maxvel1)
   fprintf(formatSpec8,velrms2,maxvel2)
   
   fprintf(formatSpec3,disprms,maxdisp)
   fprintf(formatSpec6,disprms1,maxdisp1)
   fprintf(formatSpec9,disprms2,maxdisp2)
 