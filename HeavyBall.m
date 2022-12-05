%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: HeavyBall.m
%--------------------------------------------------------------------------
% Project: Testing out parameters lambda and gamma for fast, oscillatory
% convergence globally and slow, smooth convergence locally.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00

clear all

set(0,'defaultTextInterpreter','tex');

% global variables
global gamma lambda gamma_0 gamma_1 lambda_0 lambda_1 c_0 c_10 V0 V1 delta
%%%%%%%%%%%%%%%%%%%%%
% setting the globals
%%%%%%%%%%%%%%%%%%%%%

% Heavy-ball constants: 
lambda = 10.5; % Gravity. 
            % For gamma fixed, "large values of  lambda are seen to give rise to slowly converging 
            % solutions resembling the steepest descent’s while smaller values give 
            % rise to fast solutions with oscillations getting wilder as lambda decreases."
gamma = 1/2; % Viscous friction to mass ratio.

lambda_0 = 10.5;
lambda_1 = 1/2;

gamma_0 = 1/2;
gamma_1 = 1/2;

delta = 0.1;

CalculateLStar();

c_0 = 1200; % \mathcal{U}_0 
c_10 = 638.3701; % \mathcal{T}_{1,0} 

%%%%%%%%%%%%%%%%%%%%%
% setting the locals
%%%%%%%%%%%%%%%%%%%%%
timeToDeltaSlow = 0;
timeToDeltaOscillate = 0;
timeToDeltaUniting = 0;


timeToDeltaIdxUniting = 1;
timeToDeltaIdxOscillate = 1;
timeToDeltaIdxUniting = 1;

z1deltaSlow = 0;
z1deltaOscillate = 0;
z1deltaUniting = 0;

z2deltaSlow = 0;
z2deltaOscillate = 0;
z2deltaUniting = 0;

lDeltaSlow = 0;
lDeltaOscillate = 0;
lDeltaUniting = 0;

jumpsx11 = [];
jumpsx12 = [];
jumpsx21 = [];
jumpsx22 = [];
jumpst = [];
jumpsL = [];
jumpIndex = 1;

% initial conditions
z1_0 = zeros(100,1);
for i=1:length(z1_0)
    z1_0(i) = -10;  % Just making all elements the same for now. 
end

z2_0 = zeros(100,1);

q_0 = 0;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0];

% simulation horizon
TSPAN=[0 200];
JSPAN = [0 20];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',.1);

% Simulate the slow individual heavy ball controller
[tSlow,jSlow,xSlow] = HyEQsolver(@f,@g,@C,@D,...
    x0,TSPAN,JSPAN,rule,options);

% Update lambda, so we can run the fast, oscillating individual heavy ball controller 
lambda = 1/2;

% Then simulate the fast, oscillating individual heavy ball controller
[tOscillate,jOscillate,xOscillate] = HyEQsolver(@f,@g,@C,@D,...
    x0,TSPAN,JSPAN,rule,options);

% Finally, simulate the hybrid closed-loop heavy ball system
[tUniting,jUniting,xUniting] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options);

% Finding time of convergence for the slow algorithm:
    for i=2:length(xSlow(:,1))
        if (distance(xSlow(i,1)) <= delta) && (distance(xSlow(i-1,1)) > delta)
            timeToDeltaIdxSlow = i;
            z1deltaSlow = xSlow(i,1);
        end
    end
    z2deltaSlow = xSlow(timeToDeltaIdxSlow,2);
    timeToDeltaSlow = tSlow(timeToDeltaIdxSlow,1);
    
%Find the L values for the slow algorithm:
lSlow = zeros(1,length(tSlow));
for i=1:length(tSlow(:,1))
    lSlow(i) = (CalculateL(xSlow(i,1)) - CalculateLStar());
end

lDeltaSlow = lSlow(timeToDeltaIdxSlow);
    
% Finding time of convergence for the fast, oscillatory algorithm
    for i=2:length(xOscillate(:,1))
        if (distance(xOscillate(i,1)) <= delta) && (distance(xOscillate(i-1,1)) > delta)
            timeToDeltaIdxOscillate = i;
            z1deltaOscillate = xOscillate(i,1);
        end
    end
    z2deltaOscillate = xOscillate(timeToDeltaIdxOscillate,2);
    timeToDeltaOscillate = tOscillate(timeToDeltaIdxOscillate,1);
    
%Find the L values for the fast, oscillatory algorithm:
lOscillate = zeros(1,length(tOscillate));
lOscillateAvg = zeros(length(tOscillate),1);
for i=1:length(tOscillate(:,1))
    lOscillate(i) = (CalculateL(xOscillate(i,1)) - CalculateLStar());
    lOscillateAvg(i,1) = (CalculateL(xOscillate(i,1)) - CalculateLStar()); 
end

lDeltaOscillate = lOscillate(timeToDeltaIdxOscillate);

% This is the dotted line indicating the average value:
PO = polyfit(tOscillate,log10(lOscillateAvg(:,1)),1);
yO = polyval(PO,tOscillate);
    
% Finding time of convergence for the hybrid closed-loop system
    for i=2:length(xUniting(:,1))
        if (distance(xUniting(i,1)) <= delta) && (distance(xUniting(i-1,1)) > delta)  
            timeToDeltaIdxUniting = i;
            z1deltaUniting = xUniting(i,1);
        end
    end
    z2deltaUniting = xUniting(timeToDeltaIdxUniting,2);
    timeToDeltaUniting = tUniting(timeToDeltaIdxUniting,1);
    
%Find the L values:
lUniting = zeros(1,length(tUniting));
for i=1:length(tUniting(:,1))
    lUniting(i) = (CalculateL(xUniting(i,1)) - CalculateLStar());
end

lDeltaUniting = lUniting(timeToDeltaIdxUniting);
    
% Find the jumps for plotting:
for i=2:length(jUniting)
    if(jUniting(i,1) ~= jUniting(i-1,1))
        jumpsx11(jumpIndex) = xUniting(i,1);
        jumpsx12(jumpIndex) = xUniting(i,2);
        jumpsx21(jumpIndex) = xUniting(i,3);
        jumpsx22(jumpIndex) = xUniting(i,4);
        jumpsL(jumpIndex) = lUniting(i);
        jumpst(jumpIndex) = tUniting(i,1);
        jumpIndex = jumpIndex + 1;
    end
end

figure(1)
clf
semilogy(tUniting,lUniting,'LineWidth',1.5);
hold on
semilogy(tSlow,lSlow,'Color',[0.4660 0.6740 0.1880],'LineWidth',1.5);
semilogy(tOscillate,lOscillate,'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5,'LineStyle','--');
semilogy(tOscillate,10.^(yO(:,1)),'k--','LineWidth',1.5);

semilogy(timeToDeltaUniting,lDeltaUniting,'k.','MarkerSize', 20)
strDeltaUniting = [num2str(timeToDeltaUniting),' s'];
text(timeToDeltaUniting,lDeltaUniting,strDeltaUniting,'HorizontalAlignment','left','VerticalAlignment','top','FontSize',12);

semilogy(timeToDeltaSlow,lDeltaSlow,'k.','MarkerSize', 20)
strDeltaSlow = [num2str(timeToDeltaSlow),' s'];
text(timeToDeltaSlow,lDeltaSlow,strDeltaSlow,'HorizontalAlignment','right','VerticalAlignment','top','FontSize',12);

semilogy(timeToDeltaOscillate,lDeltaOscillate,'k.','MarkerSize', 20)
strDeltaOscillate = [num2str(timeToDeltaOscillate),' s'];
text(timeToDeltaOscillate,lDeltaOscillate,strDeltaOscillate,'HorizontalAlignment','left','VerticalAlignment','top','FontSize',12);

hold off
axis([0 50 10^(-17) 10^(5)]);
ylabel('L(z_1)-L^*','FontSize',12)
xlabel('t','FontSize',12)
legend({'Uniting','$\lambda = 10.5$','$\lambda = 0.5$','$\lambda = 0.5$, average'},'Location','southwest','Interpreter','latex')

axes('Position',[0.7 0.6 0.15 0.1])
box on
hold on
semilogy(tSlow,lSlow,'Color',[0.4660 0.6740 0.1880],'LineWidth',3);
semilogy(timeToDeltaSlow,lDeltaSlow,'k.','MarkerSize', 20)
strDeltaSlow = [num2str(timeToDeltaSlow),' s'];
text(timeToDeltaSlow,lDeltaSlow,strDeltaSlow,'HorizontalAlignment','right','VerticalAlignment','top','FontSize',12);
hold off
set(gca,'xtick',[0 50 100])
set(gca,'ytick',[10^(-2) 10^(2) 10^(5)])
axis([0 100 10^(-2) 10^(5)])
hold off

saveas(gcf,'Plots\Semilog','epsc')


figure(2)
clf
plot(tUniting,lUniting,'LineWidth',1.5);
hold on
plot(tSlow,lSlow,'Color',[0.4660 0.6740 0.1880],'LineWidth',1.5);
plot(tOscillate,lOscillate,'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5,'LineStyle',':');

plot(timeToDeltaUniting,lDeltaUniting,'k.','MarkerSize', 20)
strDeltaUniting = [num2str(timeToDeltaUniting),' s'];
text(timeToDeltaUniting,lDeltaUniting,strDeltaUniting,'HorizontalAlignment','left','VerticalAlignment','top','FontSize',12);

plot(timeToDeltaSlow,lDeltaSlow,'k.','MarkerSize', 20)
strDeltaSlow = [num2str(timeToDeltaSlow),' s'];
text(timeToDeltaSlow,lDeltaSlow,strDeltaSlow,'HorizontalAlignment','right','VerticalAlignment','top','FontSize',12);

plot(timeToDeltaOscillate,lDeltaOscillate,'k.','MarkerSize', 20)
strDeltaOscillate = [num2str(timeToDeltaOscillate),' s'];
text(timeToDeltaOscillate,lDeltaOscillate,strDeltaOscillate,'HorizontalAlignment','left','VerticalAlignment','top','FontSize',12);

hold off
% axis([0 50 -500 5000]);
ylabel('L(z_1)-L^*','FontSize',12)
xlabel('t','FontSize',12)
legend('hybrid','lambda = 10.5','lambda = 0.5')

saveas(gcf,'Plots\Distance','epsc')