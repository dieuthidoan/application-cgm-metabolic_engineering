function out = ode_paper(t,y,p,qzu)

% states
VR = y(1);            %reactor volume   [l]
S = y(2);             %substrate        [g/l]
B = y(3);             %biomass          [g/l]
C = y(4);             %metabolites  [mol/g]
P = y(5);              %proteins (T+R+Q) 
M = y(6);             %residual biomass [mol/g]
Phe = y(7);             %Phenylalanin     [g/l]
O = y(8);             %By-Products     [g/l]

%Check for L-Phe production

if t < p.feedingswitch(3)  
    phit = 0;
else
    phit = p.phi*(t-p.feedingswitch(3))./(p.Pheth+(t-p.feedingswitch(3)));
end
pPhe = phit*P;

T = (p.T_Mm*C+p.T_Ml)*P;
R = (p.R_Mm*C+p.R_Ml)*P;


%check for Tyrosine availability 
Tyrtimei = find(p.time>t,1);
Tyr = p.Tyr(Tyrtimei-1);
if Tyr > 0
    Tyrs = 1;
else
    Tyrs = 0.1;
end

%calculate reaction rates
r1 = p.rmax1*T.*S./(p.KS+S);
r2 = p.rmax2*R.*C./(p.KC+C)*Tyrs;
r3 = p.rmax3*C.*T;
r4 = p.rmax4*C./(p.K4+C)*R.*Tyrs;
if t < p.feedingswitch(3)  
    r5 = 0;
else
    r5 = (p.rmax5.*C)*pPhe;
end
r6 = p.k6*C*T;

%calculate volumetric flow rate (feed)
timeim = find(qzu.time<=t,1,'last');
timeip = find(qzu.time>t,1);

if isempty(timeip) == 1
    qzuf = qzu.feed(timeim);
else
    qzuf = qzu.feed(timeim)+(t-qzu.time(timeim))/(qzu.time(timeip)-qzu.time(timeim))*(qzu.feed(timeip)-qzu.feed(timeim));
end
rV   = qzuf;                        % VR
if t < p.feedingswitch(2)
    if t >= p.feedingswitch(1)
        Szu = qzu.Sin(1);           %Feed substrate concentration g/l
    else
        Szu = 0;
    end    
elseif t < p.feedingswitch(3)    
    Szu = qzu.Sin(2);           %Feed substrate concentration g/l
else 
    Szu = qzu.Sin(3);           %Feed substrate concentration g/l
end

%feed
rSin = (qzuf./VR).*(Szu - S);         % Sin

r = [r1;r2;r3;r4;r5;r6];

out = [ rV                              %VR
        rSin - r1*B*p.mgS              %S
        p.mg'*p.N*r.*B - qzuf./VR.*B     %B
        p.N*r-p.mg'*p.N*r.*y(4:6)     %C,M,P
        r5*B*p.mgPhe - Phe*qzuf./VR   %Phe
        r3*B*p.mgO - O*qzuf./VR       % O By-Products 
];
end