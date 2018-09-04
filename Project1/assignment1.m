%IIT KGP Computational Neuroscience
%Project-1 
%Author- Aditya Rathore, 16IE10002
u = 100
F1 = @(t1,y1) [u*y1(2);u*(1-y1(1)^2)*y1(2)-y1(1)/u]
y0 = [1; 0]
opts = odeset('stats','on');
[t1,y1]=ode15s(F1, [0 600], y0, opts);
u2 = 1
F2 = @(t2,y2) [u2*y2(2);u2*(1-y2(1)^2)*y2(2)-y2(1)/u2]
[t2,y2]=ode45(F2, [0 100], y0, opts);
u3 = 0.1
F3 = @(t3,y3) [u3*y3(2);u3*(1-y3(1)^2)*y3(2)-y3(1)/u3]
[t3,y3]=ode45(F3, [0 200], y0, opts);
subplot(3,2,1)
plot(y1(:,1),y1(:,2),'.-'),title('Phase Plane'),legend('mu=100')
subplot(3,2,2)
plot(t1,y1(:,1)),title('Oscillatory Motion'),legend('mu=100')
subplot(3,2,3)
plot(y2(:,1),y2(:,2),'.-'),legend('mu=1')
subplot(3,2,4)
plot(t2,y2(:,1)),legend('mu=1')
subplot(3,2,5)
plot(y3(:,1),y3(:,2),'.-'),legend('mu=0.1')
subplot(3,2,6)
plot(t3,y3(:,1)),legend('mu=0.1')