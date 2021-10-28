function p = parameter_lf22_paper
% Stöchiometrie

p.gammaS = 1;
p.gammaCP = 300;
p.gammaCM = 300;
p.gammaP = 1;
p.gammaM = 1;
p.gammaPhe = 2;

p.beta = 0.5;
p.phi = 0.05;
p.X = 1;
p.tol = 10^(-6);

p.mgC = 88.06;                  %[g/mol]
p.mgM = 102.3*p.gammaCM;          %[g/mol]
p.mgR = p.X*102.3*p.gammaCP;          %[g/mol]
p.mgT = p.X*102.3*p.gammaCP;          %[g/mol]
p.mgQ = p.X*102.3*p.gammaCP;          %[g/mol]
p.mgP = p.X*102.3*p.gammaCP;          %[g/mol]
p.mgM = p.mgC*p.gammaCM;          %[g/mol]
p.mgR = p.X*p.mgC*p.gammaCP;          %[g/mol]
p.mgT = p.X*p.mgC*p.gammaCP;          %[g/mol]
p.mgQ = p.X*p.mgC*p.gammaCP;          %[g/mol]
p.mgP = p.X*p.mgC*p.gammaCP;          %[g/mol]
p.mgS = 92.09;                    %Glycerin [g/mol]
p.mgO = 60;                    %Acetat [g/mol]                     %[g/mol]
p.mgPhe = 165.19;                 %[g/mol]

p.mg = [p.mgC; p.mgP; p.mgM];       %[g/mol]
p.N = [p.gammaS, -p.gammaCP, -1, -p.gammaCM, -p.gammaPhe, -1; 0, p.gammaP, 0,0,0,0;0,0,0,p.gammaM,0,0];

p.rmax1 = 3.5*10^3;                     %Substrate uptake          
p.rmax2 = 1.2;                          %protein production
p.rmax3 = 1*10^4;                       %energy pathway
p.rmax4 = 1.4;                          %residual biomass synthesis
p.rmax5 = 7*10^7;                       %L-Phe production
p.k6 = 8*10^6;                          %respiration

p.KC = 0.1*10^(-4);                     %Michaelis Menten constant for C       [mol/l]
p.KE= 2*10^(-6);                        %Michaelis Menten energy pathway        [mol/l]
p.KS = 0.05;                            %substrate uptake Michaelis Menten      [g/l]
p.K4 = 0.1*10^(-4);                     %Michaelis Menten other cell party      [mol\l]
p.Pheth = 5;


%relation from data
p.T_Mm = -47.23;
p.T_Ml = 0.223;
p.R_Mm = 91.44;
p.R_Ml = 0.131;

end