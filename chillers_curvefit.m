A=xlsread('chillers.xlsx','350KW');
B=xlsread('chillers.xlsx','500KW');
C=xlsread('chillers.xlsx','200KW');
PLRA=A(:,1);
COPA=A(:,2);
PLRB=B(:,1);
COPB=B(:,2);
PLRC=C(:,1);
COPC=C(:,2);
cofA=[];
cofB=[];
cofC=[];
XA=[PLRA.^2,PLRA,ones(size(PLRA))];
XB=[PLRB.^2,PLRB,ones(size(PLRB))];
XC=[PLRC.^2,PLRC,ones(size(PLRC))];
disp('COP=aPLR^2+bPLR+c')
cofA=XA\COPA
disp('COP=aPLR^2+bPLR+c')
cofB=XB\COPB
disp('COP=aPLR^2+bPLR+c')
cofC=XC\COPC