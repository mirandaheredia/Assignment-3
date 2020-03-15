%% Part 1  A3
% Extension of Assignment 1
% Miranda Heredia 100996160
clear 
close all

%Constants given
length = 200e-9;
height = 100e-9;
% iterations = 50; not used
Noelectron = 300;
K = 1.38064e-23;
T = 300;        %temp in Kelvin
mo = 9.1093856e-31;
m = 0.26*mo;
vth = sqrt(2*K*T/m);
Time_step = .002e-12;
q_0 = 1.60217653e-19;             % electron charge


% Voltage applied in x-direction of electric field
Voltage = 0.1;

Vel_x = zeros(Noelectron,1);
Vel_y = zeros(Noelectron,1);

Pos_x = zeros(Noelectron,1);
Pos_y = zeros(Noelectron,1);

Pos_x(:,1) = length*rand(Noelectron,1);
Pos_y(:,1) = height*rand(Noelectron,1);

e_x = Voltage/(length);     %value of electric field across x (PART A)

fprintf('Part a) The value of the electric field on the electrons is %i\n',e_x);

force = e_x*q_0;         %creates a vector containing forces of all electrons
                            %(PART B)
fprintf('Part b) The value of force on each electron is %i\n',force);  

a_elec = force/m;        %creates a vector containing all acceleration of electrons
                            %will be used to modify the plot. stay tuned
                            %for more (PART C)
fprintf('Part c) The value of acceletation on each electron is %i\n',a_elec);                                  


for n=1:Noelectron
    Vel_x(n,1)= randn(1,1)*vth;
    Vel_y(n,1) = randn(1,1)*vth;
end
AvgV = sqrt(Vel_x.^2 + Vel_y.^2);


%initialize temperature vector
timesteps = 100;
T_avg_V = zeros(timesteps,1);


%PScattering
colorarray = rand(Noelectron,1);
p_scatter = 1 - exp(-Time_step/0.2e-12);
time = 1:timesteps;
for n= 1:timesteps
    
    % velocity gets updated at each time step
     Vel_x = Vel_x + a_elec*(timesteps*4e-17);
    
    Pos_x_old = Pos_x;
    Pos_y_old = Pos_y;

    random = rand(Noelectron,1);
    
    %all electrons with higher probabilities
    new = random < p_scatter;
    
    %all electrons with lower probabilities
    new2 = random >= p_scatter;
    

    rand_v_x = zeros(Noelectron,1);
    rand_v_y = zeros(Noelectron,1);
    
   for i = 1:1:Noelectron
     r1 = randi([1 Noelectron], 1,1);
     r2 = randi([1 Noelectron], 1,1);
        rand_v_x(i,1) = Vel_x(r1,1);
        rand_v_y(i,1) = Vel_y(r2,1);
   end
   %all electrons with lower probabilities will stay the same
   Vel_x = Vel_x.*new2;
   Vel_y = Vel_y.*new2;
   
   rand_v_x=rand_v_x.*new;
   rand_v_y=rand_v_y.*new;
    
   Vel_x = Vel_x+rand_v_x;
   Vel_y = Vel_y+rand_v_y;
    
   Pos_x = Pos_x + Vel_x*Time_step;
   Pos_y = Pos_y +Vel_y*Time_step;
    
   %checking for boundary positions
    idLong = Pos_x>=length;
    Pos_x(idLong) = Pos_x(idLong) - length;
    Pos_x_old(idLong) = 0;
    
    idShort = Pos_x<=0;
    Pos_x(idShort) = Pos_x(idShort) + length;
    Pos_x_old(idShort) = length;
    
    %Check for y boundary and correct
       
    Vel_y(Pos_y>=height) = -1*Vel_y(Pos_y>=height);
    Vel_y(Pos_y<=0) = -1*Vel_y(Pos_y<=0);
    
    Pos_y(Pos_y>height) = height - (Pos_y(Pos_y>height)-height);
   
    %average thermal velocity
    v_avg = mean(sqrt(((Vel_x.^2)+(Vel_y.^2))));
    T_avg_V(n) = (1/sqrt(2))*(m*(v_avg))/K;
    
    
    %%
  
    
    mfp = (0.2e-15)*(v_avg);
    %%
    meantime = (Noelectron./sum(v_avg))*mfp;
    
   
  
    % 2D Trajectory as in Figure 1 in Manual
    figure(1)
    scatter(Pos_x,Pos_y,3,colorarray);
    axis([0 200e-9 0 100e-9])
    title(['The mean free path is ', num2str(mfp)]);
    hold on
    
    % Part D Drift Current of electron
    
    e_conc = 10^15;
    Id = v_avg*e_conc*e_x*q_0;
    
    %Plot for current density of electrons
    
    figure (2)
    scatter(n,Id,'r.')
    title('Current Density of Electrons');
    %% As you can see from the plot, as time goes on, the current density in the electrons increase
    hold on
    
   
end

% Part E - Desnity and Temp Map

figure(3)
hist3([Pos_x Pos_y],'CdataMode','auto','Nbins', [20 20]);
title('Electron Density Map of Final Positions')
colorbar
view(2)

% Still havent figured out the temp map :) 






