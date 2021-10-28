tic
clear all
close all

%Parameter

p = parameter_lf22_paper;
options = odeset('RelTol',1e-7,'AbsTol',1e-9,'NonNegative',1);

%% read experimental data

%Files
expdata = 'LF22_withvolume.xlsx';
pumpread = 'Pump_characteristics.xls';

%Read files
CDW_read = xlsread(expdata,2,'P8:P65');
CDW_SD_read = xlsread(expdata,2,'Q8:Q65');
time_read = xlsread(expdata,9,'D7:D25');
CDW = CDW_read(4:3:end);
CDW_sd = CDW_SD_read(4:3:end); 
p.feedingswitch = xlsread(expdata,1,'C15:C17');
Glyc = xlsread(expdata,5,'L38:L56');
Byp = xlsread(expdata,5,'J38:J56');
Phe_read = xlsread(expdata,4,'G31:G49');
Tyr_read = xlsread(expdata,4,'H31:H49');
VR_read = xlsread(expdata,3,'M5:M540')./1000;
tVR_read = xlsread(expdata,3,'A5:A540');
mu_read = xlsread(expdata,8,'I11:I28');
t1 = 0:0.01:p.feedingswitch(1);
t2 = ((p.feedingswitch(1)+0.01):0.01:p.feedingswitch(3));
t3 = ((p.feedingswitch(3)+0.01):0.01:time_read(end));
t = [t1,t2,t3];
Byp(Byp<0) = 0;

%Feed read
f1_setp =  xlsread(pumpread,1,'D12:D15');
f1_mfr =  xlsread(pumpread,1,'I12:I15');                                %Volumetric flow in l/h
f1_gfr =  xlsread(pumpread,1,'J12:J15');                                % flow in g/h
f1_func = polyfit(f1_mfr,f1_setp,1);
f1_func_g = polyfit(f1_gfr,f1_setp,1);
f2_setp =  xlsread(pumpread,2,'D12:D15');
f2_mfr =  xlsread(pumpread,2,'I12:I15');
f2_gfr =  xlsread(pumpread,2,'J12:J15');
f2_func = polyfit(f2_mfr,f2_setp,1);
f2_func_g = polyfit(f2_gfr,f2_setp,1);
f3_setp =  xlsread(pumpread,3,'D12:D15');
f3_mfr =  xlsread(pumpread,3,'I12:I15');
f3_func = polyfit(f3_mfr,f3_setp,1);
feedp_read = [0;xlsread(expdata,3,'L5:L540')];             %inflow in %
qzu.time = [0;xlsread(expdata,3,'A5:A540')];              %time points for feed data
qzu.feed = zeros(length(feedp_read),1);                %feed l/h
qzu.Sin = [120;400;800];

%Tyrosine
Tyr_read(2) = 0.001;
p.Tyr = pchip(time_read,Tyr_read,t);
Tyrs = ones(length(t),1);
Tyrs(p.Tyr <0) = 0.1;
p.time = t;

%calculate volumetric inflow from % inflow data by using profil of pump
%characteristic (points for 1,2,5,10% inflow) - seperately for all 3
%different feeding profiles 

for i=1:length(feedp_read)
    if feedp_read(i) >0
        if qzu.time(i) < p.feedingswitch(2)
            if qzu.time(i) >= p.feedingswitch(1)
                qzu.feed(i) = (feedp_read(i)-f1_func(2))/f1_func(1); %feed l/h
            end    
        elseif qzu.time(i) < p.feedingswitch(3)    
            qzu.feed(i) = (feedp_read(i)-f2_func(2))/f2_func(1);   %feed l/h
        else 
            qzu.feed(i) = (feedp_read(i)-f3_func(2))/f3_func(1);   %feed l/h
        end
    end
end


%% simulation biomass production phase

%initial conditions
p.C0 = 3.5*10^(-4);
p.alpha = 0.5;
y01 = [1;Glyc(1);CDW(1);p.C0;(p.alpha-p.C0*p.mgC)/p.mgP;p.alpha/p.mgM;Phe_read(1);Byp(1)];

%simulation for batch phase
[time1,y1] = ode15s(@ode_paper,t1,y01,options,p,qzu);

y02 = y1(end,:);
% Simulation fed-Batch

[time2,y2] = ode15s(@ode_paper,t2,y02,options,p,qzu);

y03 = y2(end,:);

%% simulation production phase

%different parameters for respiration and overflow metabolism
k6_b = p.k6;
k3_b = p.rmax3;
p.k6 = 8*10^7;
p.rmax3 = 5*10^5;
k6_p = p.k6;
k3_p = p.rmax3;

%simulation with fitted parameters
[time3,y3] = ode15s(@ode_paper,t3,y03,options,p,qzu);

%%
%get all Simulation data

time = [time1;time2;time3];
time_1 = time;
y = [y1;y2;y3];

qzuf = zeros(length(time_1),1);
phit = zeros(length(time_1),1);
Szu = zeros(length(time_1),1);
for i = 1:length(time_1)
    %calculate volumetric flow rate (feed)
    timeim = find(qzu.time<=time_1(i),1,'last');
    timeip = find(qzu.time>=time_1(i),1);
    if time_1(i) < p.tol
        qzuf(i) = 0;
    else
        qzuf(i) = qzu.feed(timeim)+(time_1(i)-qzu.time(timeim))/(qzu.time(timeip)-qzu.time(timeim))*(qzu.feed(timeip)-qzu.feed(timeim));
    end
    if time_1(i) < p.feedingswitch(2)
        if time_1(i) >= p.feedingswitch(1)
            Szu(i) = qzu.Sin(1);           %Feed substrate concentration g/l
        else
            Szu(i) = 0;
        end    
    elseif time_1(i) < p.feedingswitch(3)    
        Szu(i) = qzu.Sin(2);           %Feed substrate concentration g/l
    else 
        Szu(i) = qzu.Sin(3);           %Feed substrate concentration g/l
        phit(i) = p.phi*(time_1(i)-p.feedingswitch(3)).^4./(p.Pheth^4+(time_1(i)-p.feedingswitch(3)).^4);
    end
end

VR = y(:,1);                %reactor volume
S = y(:,2);                 %substrate
B = y(:,3);                 %biomass
C = y(:,4);                 %Metabolites
P = y(:,5);                 %proteins
M = y(:,6);                 %residual biomass
Phe = y(:,7);               %Phenylalanine
O = y(:,8);                 %Acetate

%protein fractions
Q = (1-p.beta)*P;
T = (p.T_Mm*C+p.T_Ml).*P;
R = (p.R_Mm*C+p.R_Ml).*P;

pPhe = zeros(length(time_1),1);
r3 = zeros(length(time_1),1);                   %overflow metabolism
r6 = zeros(length(time_1),1);                   %respiration
for i = 1:length(time_1)
    if time_1(i) >= p.feedingswitch(3)
        pPhe(i) = phit(i).*P(i);
        r6(i) = k6_p.*C(i).*T(i);
        r3(i) = k3_p.*C(i)*T(i);
    else
        r6(i) = k6_b.*C(i).*T(i);
        r3(i) = k3_b.*C(i)*T(i);
    end
end

%reaction rates
r1 = p.rmax1.*T.*S./(p.KS+S);                   %Substrate uptake 
r2 = p.rmax2*R.*C./(p.KC+C).*Tyrs;              %Protein synthesis
r4 = p.rmax4*C./(p.K4+C).*R.*Tyrs;              %residual biomass synthesis
r5 = zeros(length(time_1),1);                   %Phenylalanine production rate
for i = 1:length(time_1)
    if time_1(i) >= p.feedingswitch(3)
        r5(i) = pPhe(i).*(p.rmax5.*C(i));
    end
end

%calculation for mass fractions 
r = [r1';r2';r3';r4';r5';r6'];
mu = p.mg'*p.N*r;

Tmass = T./P;
Rmass = R./P;
Qmass = Q./P;
pPhemass = pPhe./P;
Mmass = M*p.mgM;
Cmass = C*p.mgC;
Pmass = P*p.mgP;

toc

