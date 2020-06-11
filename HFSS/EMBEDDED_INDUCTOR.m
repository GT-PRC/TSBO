%Copyright (c) 2020 
%3D Systems Packaging Research Center (PRC), Georgia Tech.
function objective = EMBEDDED_INDUCTOR(s_v,t_a,t_c,w_a,g,N,w_v,delta_wm,t_ratio, t_bottom)
N = round(N);
t_offset = 0;
g_v = 0;
t_m_offset = 0;
w_m = w_a*0.0254-s_v/1000-0.0254*delta_wm;
l_m = N*w_v + (N-1)*g;
t_m = t_a-t_m_offset;
alpha = atan((0.0254*w_v+0.0254*g)/(0.0254*w_a-s_v/1000-0.0254*g_v));
w_bottom = w_v*cos(alpha);
delta_diel_x = 50;
delta_diel_y = 50;
delta_diel_z = 30;
delta_x_gnd = 100;
delta_y_gnd = 100;

A = regexp( fileread('HFSS/inductor_simulation_HFSS.py'), '\n', 'split');
A{24}  = sprintf('					"Value:="		, "%0.15fum"',t_bottom);%t_ratio 
A{482}  = sprintf('					"Value:="		, "%0.15fum"',t_bottom);%t_ratio 

A{24+18}  = sprintf('					"Value:="		, "%0.15f"',t_ratio);%t_ratio 
A{42+18}  = sprintf('					"Value:="		, "%0.15fum"',s_v);%s_v
A{320+18} = sprintf('					"Value:="		, "%0.15fum"',s_v);%s_v

A{60+18} = sprintf('					"Value:="		, "%0.15fum"',t_a); %t_a
A{356+18} = sprintf('					"Value:="		, "%0.15fum"',t_a); %t_a


A{78+18} = sprintf('					"Value:="		, "%0.15fmil"',g);%g
A{392+18} = sprintf('					"Value:="		, "%0.15fmil"',g);%g

A{96+18} = sprintf('					"Value:="		, "%0.15f"',N);%N
A{410+18} = sprintf('					"Value:="		, "%0.15f"',N);%N

A{114+18} = sprintf('					"Value:="		, "%0.15fmil"',w_v);%w_v
A{428+18} = sprintf('					"Value:="		, "%0.15fmil"',w_v);%w_v

A{132+18} = sprintf('					"Value:="		, "%0.15fmil"',delta_wm);%delta_wm
A{222+18} = sprintf('					"Value:="		, "%0.15fum"',t_offset);%t_offset
A{240+18} = sprintf('					"Value:="		, "%0.15fum"',t_m_offset);%tm_offset

A{258+18} = sprintf('					"Value:="		, "%0.15fmil"',w_a);%w_a
A{374+18} = sprintf('					"Value:="		, "%0.15fmil"',w_a);%w_a

A{276+18} = sprintf('					"Value:="		, "%0.15fmil"',g_v);%g_v
A{338+18} = sprintf('					"Value:="		, "%0.15fmil"',g_v);%g_v

A{294+18} = sprintf('					"Value:="		, "%0.15fum"',t_c);%t_c
A{446+18} = sprintf('					"Value:="		, "%0.15fum"',t_c);%t_c
lateral_area = (w_a*0.0254+s_v/1000+g_v*0.0254)*l_m*0.0254;

% 
% A{240} = sprintf('					"Value:="		, "%0.15fmil"',delta_x_gnd);%delta_x_gnd
A{150+18} = sprintf('					"Value:="		, "%0.15fmil"',delta_y_gnd);%delta_y_gnd
A{168+18} = sprintf('					"Value:="		, "%0.15fmil"',delta_diel_x);%delta_diel_x
A{186+18} = sprintf('					"Value:="		, "%0.15fmil"',delta_diel_y);%delta_diel_y
A{204+18} = sprintf('					"Value:="		, "%0.15fmil"',delta_diel_z);%delta_diel_z




A{302+18} = sprintf('oModule.ExportNetworkData("delta_diel_x=\\''%0.15fmil\\'' delta_diel_y=\\''%0.15fmil\\'' delta_diel_z=\\''%0.15fmil\\'' delta_wm=\\''%0.15fmil\\'' delta_x_gnd=\\''%0.15fmil\\'' delta_y_gnd=\\''%0.15fmil\\'' g=\\''%0.15fmil\\'' g_v=\\''%0.15fmil\\'' l_m=\\''%0.15fmil\\'' N=\\''%0.15f\\'' s_v=\\''%0.15fum\\'' t_a=\\''%0.15fum\\'' t_bottom=\\''%0.15fum\\'' t_c=\\''%0.15fum\\'' t_diel=\\''0.16mm\\'' t_gnd=\\''0.07mm\\'' t_m=\\''%0.15fum\\'' t_m_offset=\\''%0.15fum\\'' t_offset=\\''%0.15fum\\'' t_ratio=\\''%0.15f\\'' w_a=\\''%0.15fmil\\'' w_bottom=\\''0\\'' w_m=\\''%0.15fmm\\'' w_v=\\''%0.15fmil\\''", ["Setup1:Sweep"], 3, "C:\\Users\\htorun3\\OneDrive - Georgia Institute of Technology\\CPMT_Paper_Codes\\TSBO_11_28_2018\\HFSS\\Project1_HFSSDesign6.s2p", ',...
    delta_diel_x,    delta_diel_y,    delta_diel_z,...
    delta_wm,    delta_x_gnd,     delta_y_gnd,...
    g,  g_v,   l_m,   N,   s_v,     t_a, t_bottom,  t_c,     t_m,...
    t_m_offset,t_offset, t_ratio,    w_a,    w_m,      w_v);
fid = fopen('HFSS/Embedded_Inductor_HFSS_3.py','w');
fprintf(fid, '%s\n', A{:});
fclose(fid);

system('"C:\Program Files\AnsysEM\AnsysEM19.4\Win64\ansysedt.exe" -RunScriptAndExit "C:\Users\htorun3\OneDrive - Georgia Institute of Technology\CPMT_Paper_Codes\TSBO_11_28_2018\HFSS\Embedded_Inductor_HFSS_3.py"')

data = read(rfdata.data,'HFSS/Project1_HFSSDesign6.s2p');
dc_data = csvread('HFSS/Rdc.csv',1,1);
Rdc = dc_data(1)*1;
freq = data.Freq;
omega = 2*pi*freq;
Y = extract(data,'Y_PARAMETERS');
Z = extract(data,'Z_PARAMETERS');
Z12 = squeeze(Z(1,2,:));
Z1 = squeeze(Z(1,1,:)-Z(1,2,:));

C = imag(1./(omega.*Z12));
L=2*imag(Z1)./omega;
R=2*real(Z1);

L_fromY = -1./(squeeze(imag(Y(1,1,:))).*omega);
R_fromY = real(1./squeeze(Y(1,1,:)));

R_ac = R;
R_ac(1:11) = R_fromY(1:11);
R_ac(12:end) = R(12:end);

fid = fopen('HFSS/fRL_Soleniodal.csv', 'w') ;
for LoopVar = 1:1:length(freq)
 fprintf(fid, '%f, %2.12f, %1.20f\n', freq(LoopVar), R_ac(LoopVar),L(LoopVar));
end
fclose(fid);
fclose('all');

ESR_C = 10e-03;  VIN = 5;   Vout = 1; I_Load = 10; FSW = 10;numPhases = 4;    
index = find(freq == 100e6,1);
Efficiency = efficiencyAnalysis(I_Load, Vout, VIN, FSW*1e6, numPhases*L(index), ESR_C, numPhases, Rdc);
efficiency = Efficiency.efficiency;
objective = (5*efficiency-2*lateral_area);


delete('HFSS/fRL_Soleniodal.csv');
delete('HFSS/Project1_HFSSDesign6.s2p')
delete('HFSS/Rdc.csv')


end

