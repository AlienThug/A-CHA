clear;
clc;
format long

data = [0	0
1.87642	7.2206
6.25411	11.08126
11.46328	20.05961
18.35849	25.82725
27.07567	32.13845
36.86963	37.78315
47.10743	43.73476
57.14421	48.22835
66.72987	51.92034
76.01699	58.21984
85.11598	61.69597
94.96533	63.75993
104.89791	67.74635
115.25174	70.05555
126.37758	75.77379
138.05354	77.51748
149.9564	80.02024
161.68191	83.7643
];

global k1 k2 N s count tspan x0 ii jj k1a k2a a
k1a = 2.09e-05;
k2a = 2.50e-05;
k1 = k1a;
k2 = k2a;
N = 200;
count = 0;
a = 5;
tspan = [0:10:180];
x0 = [0 0 0 200 300 300];
kp = [];
Yp = [];

for i=1:N
    k1 = k1a+i*1e-7;
    for j=1:N
        k2 = k2a+j*1e-6;
        [t,Y] = ode45(@kinetics,tspan,x0);
        count = count+1;
        Y3(:,count) = a*Y(:,3);
    end
end

Y4 = data(:,1);
[m,n] = size(Y3);

for i=1:n
    R(i)=1-(sum((Y3(:,i)-Y4).^2)/sum((Y4-mean(Y4)).^2));
end

[R,s]=max(R);

c = 0;
for i=1:N
    for j=1:N
        c = c+1;
        if c==s 
           ii = i;
           jj = j;
        end
    end
end

k1 = k1a+ii*1e-7;
k2 = k2a+jj*1e-6;
Yp = Y3(:,s);

plot(t,Y3,'k','linewidth',2)
hold on
plot(t,Y4,'r','linewidth',2)
hold on

function [dx]=kinetics(t,x)
global k1 k2;
dx = zeros(6,1); %Initiate dx
dx(1)=-k2*x(1)*x(5)+k2*x(3)*x(4); % 1 refers H1-(H1-H2-H3)
dx(2)=-k2*x(2)*x(6)+k2*x(1)*x(5)-k1*x(2)*x(6)+k1*x(5)*x(4); % 2 refers H1-H2-(H1-H2-H3)/H1-H2
dx(3)=-k2*x(3)*x(4)+k1*x(2)*x(6); % 3 refers H1-H2-H3
dx(4)=-k1*x(4)*x(5)-k2*x(3)*x(4); % 4 refers H1
dx(5)=-k1*x(4)*x(5)-k2*x(1)*x(5); % 5 refers H2
dx(6)=-k2*x(2)*x(6)-k1*x(2)*x(6); % 6 refers H3
end