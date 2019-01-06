%{
anthor:Lyu yuan
2018-12-15
%}
function HVAC_OPT_chiller
global Qload CAP1 CAP2 CAP3
%% test Qload-PLR curvefit
CAP1=350; %chiller1:350KW
CAP2=500; %chiller2:500KW
CAP3=200; %chiller3:200KW 
COP1=ones(100,1);
PLR1=ones(100,1);
COP2=ones(100,1);
PLR2=ones(100,1);
COP3=ones(100,1);
PLR3=ones(100,1);% initialize matrix to promote speed 
for i=1:1:100
    PLR1(i)=0.01*i;
    COP1(i)=COP_1(PLR1(i));
    PLR2(i)=0.01*i;
    COP2(i)=COP_2(PLR2(i));
    PLR3(i)=0.01*i;
    COP3(i)=COP_3(PLR3(i));
end
figure(1)
subplot(3,1,1)
plot(PLR1,COP1)
title('chiller NO.1')
subplot(3,1,2)
plot(PLR2,COP2)
title('chiller NO.2')
subplot(3,1,3)
plot(PLR3,COP3)
title('chiller NO.3')
suptitle('COP-PLR curve')
%}

%A=xlsread('coldload.xlsx');
A=xlsread('coldload_huge.xlsx')
Qe=A(:,2);
x0_store=ones(length(Qe),6); %store each optimized results
x_store=ones(length(Qe),6);
Time=ones(length(Qe),1); %time points
Echiller_init=ones(length(Qe),1); % inital straetgy without optimization
tolr=0.01
Echiller_opt=ones(length(Qe),1); % straetgy with optimization
for i=1:1:length(Qe)
    Qload=Qe(i);
    Time(i)=i;
    Q=Qe(i);
    x0=inital(Q)';
    init=Echiller_init(i)
    fun=@Echiller;
    nlcons=@constraints;
    cl=[0 0 0 0 0 0 1-tolr];%lower bounds
    cu=[1 1 1 1 1 1 1+tolr];%upper bounds
    xtype='BCBCBC'; %varibles type B:{0,1} C:float  
    Opt=opti('fun',fun,'nl',nlcons,cl,cu,'x0',x0,'xtype',xtype)
    [x,fval,exitflag,info] = solve(Opt)
    % adjust the initial setting for at first,when wi was set as 0,PLR was
    % not 0(to suit high efficency region constraints, it is valid for program 
    % but don't correspond to reality.so we
    % set PLR=0 when the corresponding w was 0
    if x(1)==0
        x(2)=0;
    end
    if x(3)==0
        x(4)=0;
    end
    if x(5)==0
        x(6)=0;
    end   
    if x0(1)==0
        x0(2)=0;
    end    
    if x0(3)==0
        x0(4)=0;
    end
    if x0(5)==0
        x0(6)=0;
    end    
    
    x0_store(i,:)=x0;
    x_store(i,:)=x;       
    Echiller_init(i)=Echiller(x0);
    Echiller_opt(i)=fval;%save data
end
% plot
figure(2)
plot(Time,Echiller_init,'r')
hold on 
plot(Time,Echiller_opt,'b')
xlabel('Time Point')
title('Power consumption')
legend('inital','after optimization');
figure(3)
plot(Time,x0_store(:,1))
hold on
plot(Time,x_store(:,1))
hold on
plot(Time,x0_store(:,3))
hold on
plot(Time,x_store(:,3))
hold on 
plot(Time,x0_store(:,5))
hold on
plot(Time,x_store(:,5))
xlabel('Time Point')
ylabel('start/stop')
title('start or stop of chillers')
ylim([0,1.2])
legend('chiller1','chiller1opt','chiller2','chiller2opt','chiller3','chilleropt')
figure(4)
plot(Time,x0_store(:,2))
hold on
plot(Time,x_store(:,2))
hold on
plot(Time,x0_store(:,4))
hold on
plot(Time,x_store(:,4))
hold on 
plot(Time,x0_store(:,6))
hold on
plot(Time,x_store(:,6))
xlabel('Time Point')
ylabel('PLR')
ylim([0,1.2]);
title('PLR')
legend('chiller1','chiller1opt','chiller2','chiller2opt','chiller3','chiller3opt')
end



function obj=Echiller(x)
% object of optimization
global Qload CAP1 CAP2 CAP3
w1=x(1);
PLR1=x(2);
w2=x(3);
PLR2=x(4);
w3=x(5);
PLR3=x(6);
L1=PLR1*CAP1/Qload;
L2=PLR2*CAP2/Qload;
L3=PLR3*CAP3/Qload;
COP1=COP_1(PLR1);
COP2=COP_2(PLR2);
COP3=COP_3(PLR3);
obj=(L1*w1/COP1+L2*w2/COP2+L3*w3/COP3)*Qload;
end

function cons=constraints(x)
% constraints for optimization problems
global Qload CAP1 CAP2 CAP3
w1=x(1);
PLR1=x(2);
w2=x(3);
PLR2=x(4);
w3=x(5);
PLR3=x(6);
L1=PLR1*CAP1/Qload;
L2=PLR2*CAP2/Qload;
L3=PLR3*CAP3/Qload;
COP1=COP_1(PLR1);
COP2=COP_2(PLR2);
COP3=COP_3(PLR3);
cons(1)=PLR1;
cons(2)=PLR2;
cons(3)=PLR3;
cons(4)=L1;
cons(5)=L2;
cons(6)=L3;
cons(7)=L1*w1+L2*w2+L3*w3;
%cons(8)=Qload*L1*w1
%cons(9)=Qload*L2*w2
%cons(10)=Qload*L3*w3
end

function COP=COP_1(PLR)
a=-1.0510;
b=2.8079;
c=2.1644;
COP=a*PLR^2+b*PLR+c;
end
function COP=COP_2(PLR)
a=-1.0970;
b=3.0584;
c=2.3589;
COP=a*PLR^2+b*PLR+c;
end
function COP=COP_3(PLR)
a=-0.7852;
b=2.4275;
c=2.0460;
COP=a*PLR^2+b*PLR+c;
end

function x=inital(Q)
% initial strategy without optimization
global Qload CAP1 CAP2 CAP3
if (0<Q)&&(Q<=200)
    w1=0;
    PLR1=1;
    w2=0;
    PLR2=1;
    w3=1;
    PLR3=Q/CAP3;  
elseif (200<Q)&&(Q<=350)
    w1=1;
    PLR1=Q/CAP1;
    w2=0;
    PLR2=1;
    w3=0;
    PLR3=1;
elseif (350<Q)&&(Q<=500)
    w1=0;
    PLR1=1;
    w2=1;
    PLR2=Q/CAP2;
    w3=0;
    PLR3=1;  
elseif (Q>500)&&(Q<=550)
    w1=1;
    PLR1=7*Q/11/CAP1;
    w2=0;
    PLR2=1;
    w3=1;
    PLR3=4*Q/11/CAP3;
elseif (Q>550)&&(Q<=700)
    w1=0;
    PLR1=1;
    w2=1;
    PLR2=10*Q/14/CAP2;
    w3=1;
    PLR3=4*Q/14/CAP3;
elseif (Q>700)&&(Q<=850)
    w1=1;
    PLR1=7*Q/17/CAP1;
    w2=1;
    PLR2=10*Q/17/CAP2;
    w3=0;
    PLR3=1;
elseif (Q>850)&&(Q<=1050)
    w1=1;
    PLR1=7*Q/21/CAP1;
    w2=1;
    PLR2=10*Q/21/CAP2;
    w3=1;
    PLR3=4*Q/21/CAP3;
end
x(1)=w1;
x(2)=PLR1;
x(3)=w2;
x(4)=PLR2;
x(5)=w3;
x(6)=PLR3;
end