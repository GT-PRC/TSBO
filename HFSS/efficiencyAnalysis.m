%Copyright (c) 2020 
%3D Systems Packaging Research Center (PRC), Georgia Tech.

function EfficiencyParams = efficiencyAnalysis(I_Load,V_OUT, VIN, FSW, L, ESR_C, numPhases, R_DC)
if(nargin<8)
    disp('Incorrect Number of Arguments');
    disp('Correct Usage: EfficiencyParams = efficiencyAnalysis(I_Load,V_OUT, VIN, FSW, L, ESR_C, numPhases, R_DC)');
    disp('I_Load: Load current in A');
    disp('V_OUT : Output Voltage in Volts. Use either 1V or 1.05V');
    disp('VIN : Input Voltage in Volts');
    disp('FSW : Switching Frequency in Hz');
    disp('L : Total Inductor Value in H');
    disp('ESR_C : ESR for Output Capacitor in ohm');
else
    EfficiencyParams = struct('efficiency',0,'Iripple',0,'LossBreakdown',[]);
    R_m = csvread('fRL_Soleniodal.csv');

    r=R_m(:,2); %Frequency Dependent Inductor ESR Values
    f=R_m(:,1); %Frequency
    R_ac = @(freq) interp1(f,r,freq);
    
    L_phi=L/numPhases; ESR_C=ESR_C/numPhases;
    
    I_Load = I_Load/numPhases;
    Iopt = I_Load/2; % 50% of peak load
    
    % V_OUT = 1.05;
    RL_phi = 0.02;N = 1; L = L_phi;
    tSW=1/FSW;
    

    tech_r_scale=1;
    tech_c_scale=1;
    
    if(VIN==5)
        % Comment/Uncomment for simulations with 3.3V/1.5V devices in 130nm CMOS
        % 3.3V devices
        Cgsn_f=1.17e-12/tech_c_scale; %per mm
        Cgdn_f=1.184e-12/tech_c_scale; %per mm
        Cgsp_f=1.01e-12/tech_c_scale; %per mm
        Cgdp_f=1.132e-12/tech_c_scale; %per mm
        
        r_scale=1;
        Rp_f = 4.08/r_scale/tech_r_scale; % for 1mm
        Rn_f = 1.108/r_scale/tech_r_scale; % for 1mm
    elseif(VIN==3)
        % 1.2-1.5V Devices. The below numbers are for 1.2V devices, have to scale
        % them up for 1.5V
        Cgsn_f=0.813e-12/tech_c_scale; %per mm
        Cgdn_f=0.583e-12/tech_c_scale; %per mm
        Cgsp_f=0.71e-12/tech_c_scale; %per mm
        Cgdp_f=0.59e-12/tech_c_scale; %per mm
        
        r_scale=1.08/0.78; %~1.39 for for 1.2V to 1.5V scaling
        % r_scale=0.43/0.78; %~0.55 for 1.2V to 0.85V scaling
        Rp_f = 1.94/r_scale/tech_r_scale; % for 1mm
        Rn_f = 0.505/r_scale/tech_r_scale; % for 1mm
        
    elseif(VIN==1.7)
        % 1.2-1.5V Devices. The below numbers are for 1.2V devices, have to scale
        % them up for 1.5V
        Cgsn_f=0.813e-12/tech_c_scale; %per mm
        Cgdn_f=0.583e-12/tech_c_scale; %per mm
        Cgsp_f=0.71e-12/tech_c_scale; %per mm
        Cgdp_f=0.59e-12/tech_c_scale; %per mm
        
        %r_scale=1.08/0.78; %~1.39 for for 1.2V to 1.5V scaling
        r_scale=0.43/0.78; %~0.55 for 1.2V to 0.85V scaling
        Rp_f = 1.94/r_scale/tech_r_scale; % for 1mm
        Rn_f = 0.505/r_scale/tech_r_scale; % for 1mm
    end
    
    
    
    
    D_CCM = V_OUT/VIN;
    I_Ripple=(VIN-V_OUT)*D_CCM/2/L/FSW;
    
    % Optimal Sizing of FETs
    PG_MP_factor = ((Cgsp_f+ Cgdp_f)/2)*3*(VIN/2)^2/tSW;
    PG_MN_factor = ((Cgsn_f+Cgdn_f)/2)*3*(VIN/2)^2/tSW;
    PR_MP_factor =(Iopt^2+(I_Ripple^2/3))*D_CCM*Rp_f;
    PR_MN_factor =(Iopt^2+(I_Ripple^2/3))*(1-D_CCM)*Rn_f;
    
    w_P=sqrt(PR_MP_factor/PG_MP_factor);
    w_N=sqrt(PR_MN_factor/PG_MN_factor);
    
    scale=1;
    w_P=scale*w_P;
    w_N=scale*w_N;
    Rds_on_P = Rp_f/w_P;
    Rds_on_N = Rn_f/w_N;
    Cpfet = (Cgsp_f+Cgdp_f)*w_P;
    Cnfet = (Cgsn_f+Cgdn_f)*w_N;
    

%   P_INDUCTOR_AC=2*I_Ripple^2*(0.405^2*R_ac(FSW) + 0.045^2*R_ac(3*FSW) + 0.016^2*(R_ac(5*FSW)));  
    P_INDUCTOR_AC = (R_ac(FSW) * (2*I_Ripple)^2)/12;
    P_INDUCTOR_DC=I_Load^2*R_DC; 
    P_FET_RESISTIVE=2*((I_Load^2+(I_Ripple^2)/3)*D_CCM*Rp_f/w_P + (I_Load^2+(I_Ripple^2)/3)*(1-D_CCM)*Rn_f/w_N); % factor of 2 for stacked topology
    P_FET_SWITCHING=2*(PG_MP_factor*w_P+ PG_MN_factor*w_N);
    P_CAP_LOSS_AC=I_Ripple^2/3*ESR_C;
    
    P_TOTAL = P_FET_RESISTIVE + P_FET_SWITCHING + P_INDUCTOR_DC + P_INDUCTOR_AC +  P_CAP_LOSS_AC;
    
    EfficiencyParams.efficiency = 100*(V_OUT*I_Load)/(P_TOTAL+V_OUT*I_Load);
    EfficiencyParams.Iripple = I_Ripple;
    EfficiencyParams.LossBreakdown = [P_FET_RESISTIVE,P_FET_SWITCHING,P_INDUCTOR_DC,P_INDUCTOR_AC,P_CAP_LOSS_AC];
    
end
end
