clear all;
clc,
L=input('Enter Length in [m]=  ');
% insert 1 (i.e.)
R=input('Enter Radius in [m]=');
% insert 0.02 (i.e.)
Ix=(1/4)*pi*R^4;
A=pi*R^2;
E=input('Enter Young"s modulus in [Pa]=  ');
% insert 70000000000 Pa i.e. for aluminium
Ro=input('Enter material density in [kg/m^3]=  ');
% insert 27000 Kg/m3 i.e. for aluminium
% initial values for fzero
aO = [4.7 7.2 10.8 14 17];
zO = [1.8 4.7 7.2 10.8 14];
uO = [1.5 3.7 4.7 7 7.7];
pO = [3 6.2 9.3 12.5 15.6];
for i=1:1:5
% fzero to obtain the spatial frecuencies
an(i)=fzero('cos(x)*cosh(x)-1',aO(i));
zn(i)=fzero('cos(x)*cosh(x)+1',zO(i));
un(i)=fzero('tan(x)-tanh(x)',uO(i));
pn(i)=fzero('sin(x)',pO(i));
% coeficients
F(i)=(cosh(an(i))-cos(an(i)))./(sinh(an(i))-sin(an(i)));
K(i)=(sinh(zn(i))-sin(zn(i)))./(cosh(zn(i))-cos(zn(i)));
Q(i)=(cosh(un(i))-cos(un(i)))./(sinh(un(i))-sin(un(i)));
ax(i)=(an(i)/L)';
zx(i)=(zn(i)/L)';
ux(i)=(un(i)/L)';
px(i)=(pn(i)/L)';
% naturel frecuencies
fna(i)=(ax(i)^2*sqrt((E*Ix)/(Ro*A)))/(2*pi);
fnz(i)=(zx(i)^2*sqrt((E*Ix)/(Ro*A)))/(2*pi);
fnu(i)=(ux(i)^2*sqrt((E*Ix)/(Ro*A)))/(2*pi);
fnp(i)=(px(i)^2*sqrt((E*Ix)/(Ro*A)))/(2*pi);
end
x=linspace(0, L, 180);
xl=x./L; %important
Tc0='(cosh(ax(i).*x(j))+cos(ax(i).*x(j)))+F(i).*(sinh(ax(i).*x(j))+sin(ax(i)*x(j)))';
Tc1='(cosh(zx(i).*x(j))-cos(zx(i).*x(j)))-K(i).*(sinh(zx(i).*x(j))-sin(zx(i)*x(j)))';
Tc2='(cosh(ux(i).*x(j))-cos(ux(i).*x(j)))-Q(i).*(sinh(ux(i).*x(j))-sin(ux(i)*x(j)))';
Tc3='(cosh(ax(i).*x(j))-cos(ax(i).*x(j)))-F(i).*(sinh(ax(i).*x(j))-sin(ax(i)*x(j)))';
Tc4='(sin(px(i)*x(j)))';
% matrices for the results
Xnx1=zeros(5,length(x));
Xnx2=zeros(5,length(x));
Xnx3=zeros(5,length(x));
Xnx4=zeros(5,length(x));
% loop to obtain all the values of each boundary condition mode shape
for i=1:1:5
    for j=1:length(x)
        Xnx0(i,j)=eval(Tc0);
        Xnx1(i,j)=eval(Tc1);
        Xnx2(i,j)=eval(Tc2);
        Xnx3(i,j)=eval(Tc3);
        Xnx4(i,j)=eval(Tc4);
    end
end
% Ploting
figure(1);
    plot(xl,Xnx0(1,:), 'b-');hold on;
    plot(xl,Xnx0(2,:), 'r-');
    plot(xl,Xnx0(3,:), 'g-');
    plot(xl,Xnx0(4,:), 'm-');
    plot(xl,Xnx0(5,:), 'k-');
    title('Mode shape of the Free-free');
    legend('Mode #1','Mode #2','Mode #3','Mode #4','Mode #5', 0); xlabel('x/L'); ylabel('Mode shape'); grid;
figure(2);
    plot(xl,Xnx1(1,:), 'b-');hold on;
    plot(xl,Xnx1(2,:), 'r-');
    plot(xl,Xnx1(3,:), 'g-');
    plot(xl,Xnx1(4,:), 'm-');
    plot(xl,Xnx1(5,:), 'k-');
    title('Mode shape of the Fixed-free');
    legend('Mode #1','Mode #2','Mode #3','Mode #4','Mode #5', 0); xlabel('x/L'); ylabel('Mode shape'); grid;  
figure(3);
    plot(xl,Xnx2(1,:), 'b-');hold on;
    plot(xl,Xnx2(2,:), 'r-');
    plot(xl,Xnx2(3,:), 'g-');
    plot(xl,Xnx2(4,:), 'm-');
    plot(xl,Xnx2(5,:), 'k-');
    title('Mode shape of the Fixed-pinned');
    legend('Mode #1','Mode #2','Mode #3','Mode #4','Mode #5', 0); xlabel('x/L'); ylabel('Mode shape'); grid;    
figure(4);
    plot(xl,Xnx3(1,:), 'b-');hold on;
    plot(xl,Xnx3(2,:), 'r-');
    plot(xl,Xnx3(3,:), 'g-');
    plot(xl,Xnx3(4,:), 'm-');
    plot(xl,Xnx3(5,:), 'k-');
    title('Mode shape of the Fixed-fixed');
    legend('Mode #1','Mode #2','Mode #3','Mode #4','Mode #5', 0); xlabel('x/L'); ylabel('Mode shape'); grid;
figure(5);
    plot(xl,Xnx4(1,:), 'b-');hold on;
    plot(xl,Xnx4(2,:), 'r-');
    plot(xl,Xnx4(3,:), 'g-');
    plot(xl,Xnx4(4,:), 'm-');
    plot(xl,Xnx4(5,:), 'k-');
    title('Mode shape of the Pinned-pinned');
    legend('Mode #1','Mode #2','Mode #3','Mode #4','Mode #5', 0); xlabel('x/L'); ylabel('Mode shape'); grid;
  
  
    

    
    
    
    
    
 
    
    