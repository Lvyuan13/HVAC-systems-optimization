%{
anthor:Lyu yuan
2018-12-15
shanghaijiaotong University
%}

function HVAC_OPT_pump
disp('start optimization...')
%% test curvefit

disp('test curve...')
A=xlsread('pump_constant.xlsx','H');
Q=A(:,1);
H=A(:,2);
H_new=ones(length(H),1);

for i=1:1:length(H)
    Q_constant=Q(i);
    %H_constant=H(i);   
    H_new(i)=lift_constant(Q_constant);
end
figure(1)
subplot(3,1,1)
plot(Q,H,'O')
hold on
plot(Q,H_new)
title('constant frequency pump H-Q')

B=xlsread('pump_constant.xlsx','N');
Q=B(:,1);
N=B(:,2);
N_new=ones(length(N),1); 
eta=ones(length(N),1); %initialize N_new and eta
for i=1:1:length(N)
    Q_constant=Q(i);
    %N_constant=N(i);
    N_new(i)=shaft_constant(Q_constant);
    eta(i)=eff_constant(Q_constant);
end
subplot(3,1,2)
plot(Q,N,'*')
hold on
plot(Q,N_new)
title('constant frequency pump N-Q')
subplot(3,1,3)
plot(Q,eta)
title('constant frequency pump Efficiency-Q')
suptitle('curvefit of  Constant frequency pump')
C=xlsread('pump_var.xlsx','H');
D=xlsread('pump_var.xlsx','N');
Q_VAR1=C(:,1);
Q_VAR2=D(:,1);
H_VAR=C(:,2);
N_VAR=D(:,2);
H_VAR_Cal=ones(length(H_VAR),1);
N_VAR_Cal=ones(length(N_VAR),1);
COF=ones(length(N_VAR),1);
for i=1:1:length(Q_VAR1)
    H_VAR_Cal(i)=lift_var(Q_VAR1(i),1);  
end
figure(2)
subplot(3,1,1)
plot(Q_VAR1,H_VAR_Cal)
hold on
plot(Q_VAR1,H_VAR,'o')
title('varible frequency pump H-Q')
for i=1:1:length(N_VAR)
    N_VAR_Cal(i)=shaft_var(Q_VAR2(i),1);
    COF(i)=eff_var(Q_VAR2(i),1);
end
subplot(3,1,2)
plot(Q_VAR2,N_VAR_Cal)
hold on
plot(Q_VAR2,N_VAR,'+')
title('varible frequency shaft N-Q')
subplot(3,1,3)
plot(Q_VAR2,COF)
title('varible frequency Efficiency-Q')
suptitle('curvefit of Varible frequency Pump')
% test varible frequency pump's coefficiency
Q_test=ones(100,1);
eta_test=ones(100,1);
kj=1;
for i=1:1:100
Q_test(i)=6*i+50;
eta_test(i)=eff_var(Q_test(i),kj);
end
figure(3)
plot(Q_test,eta_test)
title(['varible frequency Efficiency-Q(set k= ',num2str(kj),')'])
Q=ones(50,1);
eta=ones(50,1);
Hcon=ones(50,1);
for i=1:1:50
    Q(i)=i*3+50;
    eta(i)=eff_constant(Q(i));
    Hcon(i)=lift_constant(Q(i));
end
figure(4)
subplot(2,1,1)
plot(Q,eta)
title('Efficiency curve of Constant frequency pump')
subplot(2,1,2)
plot(Q,Hcon)
title('H-Q curve of Constant frequency pump')
suptitle('performance of Constant frequency pump in wider working field')
%% test algrithom
%just for test
%}
HH=lift_constant(152.1431)
%% test a set specific work point
%{
Qe_test=570.4;
He_test=23.25;
fun=@Epump
% set bounds
x0=initial_guess(Qe_test,He_test)
nlcons=@constraints
tolr=0.01
deltaH=5
kmin=0.7
Qmin_var1=Qmin_var(He_test)
Hel=He_test-tolr
Qmin_c1=70
Qmin_c2=70
Qmax_var1=Qmax_var(He_test)
Qmax_c1=160
Qmax_c2=160
Heu=He_test+tolr
cl=[Qe_test-tolr,kmin,Qmin_var1,Hel,Qmin_c1,Qmin_c2,-tolr,He_test,He_test]'
cu=[Qe_test+tolr,1,Qmax_var1,Heu,Qmax_c1,Qmax_c2,tolr,Inf,Inf]'
% bounds end
%x0=[w_var1,k_var1,Q_var1,w_c1,Q_c1,w_c2,Q_c2]'
O_init=Epump(x0)
opts = optiset('display','iter');
xtype='BCCBCBC'
Opt = opti('fun',fun,'nl',nlcons,cl,cu,'x0',x0,'options',opts,'xtype',xtype)
%Opt = opti('fun',fun,'nl',nlcon,cl,cu,'bounds',lb,ub,'options',opts,'xtype',xtype)
[x,fval,exitflag,info] = solve(Opt)
%% test end
k=H2k_var(He_test,Qe_test)
%}
%% main loop
% import working conditions
disp('importing data...')
WORK=xlsread('demand.xlsx');
Q=WORK(:,1);
H=WORK(:,2);
E_init=ones(length(Q),1);
E_new=ones(length(Q),1);
Time=ones(length(Q),1);
x_store=ones(length(Q),7);
x0_store=ones(length(Q),7);
tolr=0.1
for i=1:1:length(Q)
    disp(['the ',num2str(i),' th data point'])
    Qe=Q(i);
    He=H(i);
    x0=initial_guess(Qe,He);
    opts = optiset('display','iter');
    xtype='BCCBCBC';
    fun=@Epump
    nlcons=@constraints
    kmin=0.7;
    Qmin_var1=Qmin_var(He);
    Qmax_var1=Qmax_var(He);
    Qmin_c1=70;
    Qmin_c2=70;
    Qmax_c1=160;
    Qmax_c2=160;
    cl=[Qe-tolr,kmin,Qmin_var1,He-tolr,Qmin_c1,Qmin_c2,-tolr,He,He]';
    cu=[Qe+tolr,1,Qmax_var1,He+tolr,Qmax_c1,Qmax_c2,tolr,Inf,Inf]';
    Opt=opti('fun',fun,'nl',nlcons,cl,cu,'x0',x0,'options',opts,'xtype',xtype);
    [x,fval,exitflag,info] = solve(Opt)
    
    Time(i)=i;
    % adjust the initial setting for at first,when wi was set as 0,PLR was
    % not 0(to suit high efficency region constraints, it is valid for program 
    % but don't correspond to reality.so we
    % set PLR=0 when the corresponding w was 0
    if x(1)==0
        x(2)=0;
        x(3)=0;
    end
    if x(4)==0
        x(5)=0;
    end
    if x(6)==0
        x(7)=0;
    end
   if x0(1)==0
        x0(2)=0;
        x0(3)=0;
    end
    if x0(4)==0
        x0(5)=0;
    end
    if x0(6)==0
        x0(7)=0;
    end
   
    x_store(i,:)=x;
    x0_store(i,:)=x0;
    E_init(i)=Epump(x0);
    E_new(i)=Epump(x);
end

figure(5)
plot(Time,E_init)
hold on
plot(Time,E_new)
legend('init','optimized')
title('Power consumption')
figure(6)
plot(Time,x_store(:,1))
hold on
plot(Time,x_store(:,4))
hold on 
plot(Time,x_store(:,6))
ylim([0,1.2])
title('stop/start')
legend('Varible frequency pump','constant frequency pump 1','constant frequency pump2')

figure(7)
plot(Time,x_store(:,2))
title('frequency of Varible frequency pump')
ylim([0,1.2])
figure(8)
plot(Time,x_store(:,3))
hold on 
plot(Time,x_store(:,5))
hold on 
plot(Time,x_store(:,7))
title('flow of pumps')
legend('Varible frequency pump','constant frequency pump 1','constant frequency pump2')
end


function E=Epump(x)
% object function
% setting
% m:Dc variable frequency pump m=1
%  costant frequency pump n=2
%set coefficent for varible and constant frequency pump
% pi_varj means: j th pump's i th coefficent 
p1_var1= -6.753151610294475e-04;
p2_var1=0.558138861138861;
p3_var1=20.233939393939380;
p1_c1=-8.411258951347560e-04;
p2_c1=0.343335855781084;
p3_c1=15.418916091414905;
p1_c2=-8.411258951347560e-04;
p2_c2=0.343335855781084;
p3_c2=15.418916091414905;
w_var1=x(1);%stop or start for varible frequency pump NO1
k_var1=x(2); %Frequency modulation ratio of varible frequency pump NO1
Q_var1=x(3);  %massflow of varible frequency pump NO1
w_c1=x(4);  %stop or start of constant pump NO1
Q_c1=x(5);
w_c2=x(6);
Q_c2=x(7);
N_var1=k_var1*p1_var1*Q_var1^2+k_var1^2*p2_var1*Q_var1+k_var1^3*p3_var1;
N_c1=p1_c1*Q_c1^2+p2_c1*Q_c1+p3_c1;
N_c2=p1_c2*Q_c2^2+p2_c2*Q_c2+p3_c2;
E=w_var1*N_var1+w_c1*N_c1+w_c2*N_c2;
end

function c=constraints(x)
% subject function
% 8 constraints
w_var1=x(1);%stop or start for varible frequency pump NO1
k_var1=x(2); %Frequency modulation ratio of varible frequency pump NO1
Q_var1=x(3);  %massflow of varible frequency pump NO1
w_c1=x(4);  %stop or start of constant pump NO1
Q_c1=x(5);
w_c2=x(6);
Q_c2=x(7);
c(1)=w_var1*Q_var1+w_c1*Q_c1+w_c2*Q_c2;%[Qe,Qe]/[Qe,inf]
c(2)=k_var1;% limitation of frequency of varible frequency pump NO1 [kmin,1]
c(3)=Q_var1;% limitation of mass flow of varible frequency pump NO1
%c(3)=eff_var(Q_var1,k_var1);% the same as Q limited in [Qmin,Qmax] [0.7,inf] eta_min=0.7
c(4)=lift_var(Q_var1,k_var1); % only one varible frequency pump [He,inf]/[He,He]
% the H provided by varible frequency pump must be equal as He
c(5)=Q_c1; 
c(6)=Q_c2;%[Qmin,Qmax]
%c(5)=eff_constant(Q_c1); % the same as Q limited in [Qmin,Qmax] [0.7,inf]
%c(6)=eff_constant(Q_c2); % the same as Q limited in [Qmin,Qmax] [0.7,inf]
%c(7)=lift_constant(Q_c1)-lift_constant(Q_c2); %[0,0]
c(7)=Q_c1-Q_c2;   %same as lift_constant(Q_c1)-lift_constant(Q_c2)@[0,0]
c(8)=lift_constant(Q_c1); %[He,inf]
c(9)=lift_constant(Q_c2); %[He,inf]
end





function H=lift_var(Q,k)
h1_var=-9.047783436272624e-05;
h2_var=0.017885073491081;
h3_var=35.884545794186700;
H=h1_var*Q^2+k*h2_var*Q+k^2*h3_var;
end

function N=shaft_var(Q_var1,k_var1)
p1_var1= -6.753151610294475e-04;
p2_var1=0.558138861138861;
p3_var1=20.233939393939380;
N=k_var1*p1_var1*Q_var1^2+k_var1^2*p2_var1*Q_var1+k_var1^3*p3_var1;
end


function eta_var=eff_var(Q,k)
H=lift_var(Q,k);
N=shaft_var(Q,k);
rho=1000;
g=9.8;
eta_var=g*H*Q/N/rho;
end

function H=lift_constant(Q)
h1_c=-0.001076015088502;
h2_c=0.077175596073539;
h3_c=36.415360796790438;
H=h1_c*Q^2+h2_c*Q+h3_c;
end


function N=shaft_constant(Q_c1)
p1_c1=-8.411258951347560e-04;
p2_c1=0.343335855781084;
p3_c1=15.418916091414905;
N=p1_c1*Q_c1^2+p2_c1*Q_c1+p3_c1;
end


function eta_constant=eff_constant(Q)
H=lift_constant(Q);
N=shaft_constant(Q);
rho=1000;
g=9.8;
eta_constant=g*H*Q/N/rho;
end

function Qmax_var=Qmax_var(He)
HA=35.2118;
HB=12.1936;
QA=230;
QB=620;
h1_var=-9.047783436272624e-05;
h2_var=0.017885073491081;
h3_var=35.884545794186700;
if He>HB
    Qmax_var=(-h2_var-sqrt(h2_var^2-4*h1_var*(h3_var-He)))/(2*h1_var);
else
    Qmax_var=sqrt(He)/sqrt(HB)*QB
end

end

function Qmin_var=Qmin_var(He)
HA=35.2118;
HB=12.1936;
QA=230;
QB=620;
HC=17.2539;
h1_var=-9.047783436272624e-05;
h2_var=0.017885073491081;
h3_var=35.884545794186700;
kmin=0.7;
if He>HC
    Qmin_var=sqrt(He)/sqrt(HA)*QA;
else
    Qmin_var=(-kmin*h2_var-sqrt(h2_var^2*kmin^2-4*h1_var*(kmin^2*h3_var-He)))/(2*h1_var);
end
end

function Q=H2Q_var(He,k)
h1_var=-9.047783436272624e-05;
h2_var=0.017885073491081;
h3_var=35.884545794186700;
Q=(-k*h2_var-sqrt(h2_var^2*k^2-4*h1_var*(k^2*h3_var-He)))/(2*h1_var); 
end

function k=H2k_var(H,Q)
h1_var=-9.047783436272624e-05;
h2_var=0.017885073491081;
h3_var=35.884545794186700;
p=[h3_var,h2_var*Q,h1_var*Q^2-H];
root=roots(p);
for i=1:1:length(root)
    if (root(i)>0)&&(root(i)<=1)
        k=root(i);
    else
        k=1;
    end
end
end


function Q=H2Q_constant(He)
h1_c=-0.001076015088502;
h2_c=0.077175596073539;
h3_c=36.415360796790438;
Q=(-h2_c-sqrt(h2_c^2-4*h1_c*(h3_c-He)))/(2*h1_c);
end

function x=initial_guess(Qe,He)
% aim to set a valid initial guess straetgy
% x has seven dimensions
Qmin_var1=Qmin_var(He);
Qmax_var1=Qmax_var(He);
if Qe<Qmin_var1
    w_var1=0; 
    %close varible frequency pump
    Q_var1=(Qmin_var1+Qmax_var1)/2;
    k_var1=H2k_var(He,Q_var1);
    if Qe<=70
        w_c1=1;
        Q_c1=Qe;
        w_c2=0;
        Q_c2=115;
    elseif (Qe>70)&&(Qe<=160)
        w_c1=1;
        Q_c1=Qe;
        w_c2=0;
        Q_c2=115;
    else
        w_c1=1;
        Q_c1=Qe/2;
        w_c2=1;
        Q_c2=Qe/2;
    end
elseif Qe>Qmax_var1
    w_var1=1;
    Q_var1=Qmax_var1;
    Qexcess=Qe-Q_var1;
    k_var1=H2k_var(He,Q_var1);
    if Qexcess<160
        w_c1=1;
        Q_c1=Qexcess;
        w_c2=0;
        Q_c2=115;
    else
        w_c1=1;
        Q_c1=Qexcess/2;
        w_c2=1;
        Q_c2=Qexcess/2;   
    end
else
    w_var1=1;
    Q_var1=Qe;
    k_var1=H2k_var(He,Q_var1);
    w_c1=0;
    Q_c1=115;
    w_c2=0;
    Q_c2=115;
end
x(1)=w_var1;%stop or start for varible frequency pump NO1
x(2)=k_var1; %Frequency modulation ratio of varible frequency pump NO1
x(3)=Q_var1;  %massflow of varible frequency pump NO1
x(4)=w_c1;  %stop or start of constant pump NO1
x(5)=Q_c1;
x(6)=w_c2;
x(7)=Q_c2;
end
