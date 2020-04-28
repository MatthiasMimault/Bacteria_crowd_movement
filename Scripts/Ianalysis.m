%% plot of various functional of I
clear
close all

% Parameters
n = 100;
kappa = 0.1;
Xmin = -5;
Xmax = 5;

% Functions
I = @(b) kappa*b./sqrt(1+b.^2);
Ip = @(b) kappa./sqrt(1+b.^2).^3;
IL = @(b) (b<=0).*max(-kappa,kappa*b)+(b>0).*min(kappa*b,kappa);

% Variables
xx = linspace(Xmin,Xmax,n);

%% Plot
plot(xx,I(xx),xx,IL(xx))
axis([Xmin, Xmax, -1.2*kappa, 1.2*kappa])
legend('I','Il')