figure
subplot(3,1,1)
hold on
plot(x,Splot(:,1),"LineWidth",1.5,LineStyle="--");
plot(x,Iplot(:,1),'r');
plot(x,Rplot(:,1),'k');
title('time = 0')
xlabel('x')
ylabel('Population')
legend('Susceptible', 'Infected', 'Recovered')
subplot(3,1,2)
hold on
plot(x,Splot(:,50));
plot(x,Iplot(:,50),'r');
plot(x,Rplot(:,50),'k');
title('time = 25')
xlabel('x')
ylabel('Population')
legend('Susceptible', 'Infected', 'Recovered')
subplot(3,1,3)
hold on
plot(x,Splot(:,100));    
plot(x,Iplot(:,100),'r');
plot(x,Rplot(:,100),'k');
title('time = 100')
xlabel('x')
ylabel('Population')
legend('Susceptible', 'Infected', 'Recovered')
sgtitle('Time Evolved SIRDDFT Solution');


STotalPop = sum(Splot(:,1));
ITotalPop = sum(Iplot(:,1));
RTotalPop = sum(Rplot(:,1));
Splot(:,1) = Splot(:,1)/STotalPop;
Iplot(:,1) = Iplot(:,1)/ITotalPop;
Rplot(:,1) = Rplot(:,1)/RTotalPop;

STotalPop = sum(Splot(:,50));
ITotalPop = sum(Iplot(:,50));
RTotalPop = sum(Rplot(:,50));
Splot(:,50) = Splot(:,50)/STotalPop;
Iplot(:,50) = Iplot(:,50)/ITotalPop;
Rplot(:,50) = Rplot(:,50)/RTotalPop;
STotalPop = sum(Splot(:,100));
ITotalPop = sum(Iplot(:,100));
RTotalPop = sum(Rplot(:,100));
Splot(:,100) = Splot(:,100)/STotalPop;
Iplot(:,100) = Iplot(:,100)/ITotalPop;
Rplot(:,100) = Rplot(:,100)/RTotalPop;
figure
subplot(3,1,1)
hold on
plot(x,Splot(:,1),"LineWidth",1.5,LineStyle="--");
plot(x,Iplot(:,1),'r');
plot(x,Rplot(:,1),'k');
title('time = 0')
xlabel('x')
ylabel('population/Respective pop')
legend('Susceptible', 'Infected', 'Recovered')
subplot(3,1,2)
hold on
plot(x,Splot(:,50));
plot(x,Iplot(:,50),'r');
plot(x,Rplot(:,50),'k');
title('time = 25')
xlabel('x')
ylabel('population/Respective pop')
legend('Susceptible', 'Infected', 'Recovered')
subplot(3,1,3)
hold on
plot(x,Splot(:,100));    
plot(x,Iplot(:,100),'r');
plot(x,Rplot(:,100),'k');
title('time = 100')
xlabel('x')
ylabel('population/Respective pop')
legend('Susceptible', 'Infected', 'Recovered')
sgtitle('Time Evolved SIRDDFT Solution');
