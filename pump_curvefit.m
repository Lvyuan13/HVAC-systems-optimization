A=xlsread('pump_var.xlsx','H');
disp('varible frequency pump H=h1*Q^2+h2*Q+h3, fitting...')
Q=A(:,1);
H=A(:,2);
Qt=[ones(size(Q)),Q,Q.^2];
h=Qt\H;
format long

h3=h(1)
h2=h(2)
h1=h(3)

B=xlsread('pump_var.xlsx','N');
Q=B(:,1);
N=B(:,2);
Qt=[ones(size(Q)),Q,Q.^2];
p=Qt\N;
disp('varible frequency pump N=p1*Q^2+p2*Q+p3,  fitting...')
p3=p(1)
p2=p(2)
p1=p(3)

C=xlsread('pump_constant.xlsx','H');
Q=C(:,1);
H=C(:,2);
Qt=[ones(size(Q)),Q,Q.^2];
h=Qt\H;
disp('constant frequency pump H=h1*Q^2+h2*Q+h3,  fitting...')
h3=h(1)
h2=h(2)
h1=h(3)

D=xlsread('pump_constant.xlsx','N');
Q=D(:,1);
N=D(:,2);
Qt=[ones(size(Q)),Q,Q.^2];
p=Qt\N;
disp('constant frequency pump N=p1*Q^2+p2*Q+p3,  fitting...')
p3=p(1)
p2=p(2)
p1=p(3)