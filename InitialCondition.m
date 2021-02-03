%% Getting the Initial Condition for the bubble bursting
% This is used to generate initial condition for the simulation of bubble
% bursting, save (X,Y) of the interface in a data file then use Basilisk
% code InitialCondition.c to dump Basilisk file when can be used as initial
% condition for any simulation.
% Author: Vatsal Sanjay
% vatsalsanjay@gmail.com
% Physics of Fluids
clc
clear
close all
tic
%% Output Folder
folder = 'BasiliskInitialCondition'; % output folder
opFolder = fullfile(cd, folder);
if ~exist(opFolder, 'dir')
mkdir(opFolder);
end 
%%
% Number of elements: Start with n1 = 1e5, n2 = 1e3, TailxMax = 8.0 and RbGuess = sqrt(10)
nParams.n1 = 1e5; nParams.n11 = 5e2; nParams.n2 = 1e3; nParams.nTail = 2e4; nParams.nfillet = 500; % 500 \equi 0.22
nParams.RbGuess = sqrt(2.0);
% Tolerances
nParams.TailTol = 1e-6; nParams.VolTol = 1e-8; 
nParams.TailMaxCounter = 50; nParams.TailxMax = 8.0; nParams.VolMaxCounter = 50; 
%% Solving the set of equations to get all the three interfaces
Bond = 1e-3;
Interface = shape(Bond, nParams);
%% Save interface coordinates to be used in Basilisk
X = Interface.Bas;
filename = sprintf('Bo%5.4f.dat',Bond);
if exist(filename, 'file')
    delete filename;
end
save(filename,'X','-ascii')
X1 = [Interface.x1 Interface.y1];
X2 = [Interface.x2' Interface.y2'];
X1 = [X1;X2];
filename = sprintf('Bo%5.4f-Air.dat',Bond);
% if exist(filename, 'file')
%     delete filename;
% end
save(filename,'X1','-ascii')

%% Running Basilisk to dump file for initial condition
Bostr = num2str(Bond); ll=evalc(['!./InitialCondition ' Bostr]); % executing the command
% ll=evalc('!./getFacet snapshot-0.0000'); % executing the command
% Getting grid & VOF information
facets=textscan(ll,'%f %f\n');
Bas.Xf = facets{1}; Bas.Yf = facets{2}; 
Bas.Xf = reshape(Bas.Xf,[2, int32(length(Bas.Xf)/2.)]); 
Bas.Yf = reshape(Bas.Yf,[2, int32(length(Bas.Yf)/2.)]);
%% Optional Plot the VOF field contour
% ll=evalc('!./getData Initialdump');
% bolo=textscan(ll,'%f %f %f\n','HeaderLines',1);
% X=bolo{1}; Y=bolo{2}; f = bolo{3};
% n = int32(sqrt(length(X)));
% X = reshape(X,[n,n]); Y = reshape(Y,[n,n]); f = reshape(f,[n,n]);
%% Plotting Time
%% subplot 1
figure1 = figure('visible','on','WindowState','fullscreen','Color',[1 1 1]);
subplot1 = subplot(2,1,1,'Parent',figure1);
hold(subplot1,'on');
if (Bond >= 1)
    TextAnnotation = sprintf('\\boldmath{$Bo = %3.2f$}',Bond);
else 
    TextAnnotation = sprintf('\\boldmath{$Bo = %3.2e$}',Bond);
end
title(TextAnnotation, 'Interpreter','latex');
% Interfaces from MATLAB
plot(Interface.Bas(:,2),Interface.Bas(:,1),'c.','MarkerSize',10,'LineWidth',3,'Parent',subplot1);
% plot(Interface.y2,Interface.x2,'b-','MarkerSize',20,'LineWidth',3,'Parent',subplot1);
plot(-Interface.y1,Interface.x1,'r-','MarkerSize',20,'LineWidth',3,'Parent',subplot1);
plot(-Interface.y2,Interface.x2,'b-','MarkerSize',20,'LineWidth',3,'Parent',subplot1);
plot(-Interface.yTail,Interface.xTail,'g-','MarkerSize',20,'LineWidth',3,'Parent',subplot1);
% Interface from Basilisk
plot(Bas.Yf,Bas.Xf,'k-','MarkerSize',20,'LineWidth',3,'Parent',subplot1);
% Contour from Basilisk
% contourf(Y,X,f,[0.5 0.5],'facecolor', 'r','Parent',subplot1)
% Box for next plot
xWindow1 = Interface.xTail(1)-0.5; xWindow2 = Interface.xTail(1); 
yWindow1 = Interface.yTail(1); yWindow2 = Interface.yTail(1)+1.0; 
plot([yWindow1 yWindow2],[xWindow1 xWindow1],'--','Color',[0.5 0.5 0.5],'MarkerSize',20,'LineWidth',2,'Parent',subplot1);
plot([yWindow2 yWindow2],[xWindow1 xWindow2],'--','Color',[0.5 0.5 0.5],'MarkerSize',20,'LineWidth',2,'Parent',subplot1);
plot([yWindow1 yWindow2],[xWindow2 xWindow2],'--','Color',[0.5 0.5 0.5],'MarkerSize',20,'LineWidth',2,'Parent',subplot1);
plot([yWindow1 yWindow1],[xWindow1 xWindow2],'--','Color',[0.5 0.5 0.5],'MarkerSize',20,'LineWidth',2,'Parent',subplot1);
% Plot features
axis equal
xlabel('\boldmath{$X$}','LineWidth',2,'FontWeight','bold','FontSize',50,...
            'FontName','times new roman',...
            'Interpreter','latex');
ylabel('\boldmath{$Y$}','LineWidth',2,'FontWeight','bold','FontSize',50,...
    'FontName','times new roman',...
    'Interpreter','latex');
box(subplot1,'on');
set(subplot1,'FontName','times new roman','FontSize',30,'FontWeight','bold',...
    'LineWidth',3);
xlim(subplot1,[-8.0 8.0]);
ylim(subplot1,[-2 2]);
box(subplot1,'on');
set(subplot1,'FontName','times new roman','FontSize',30,'FontWeight','bold',...
    'LineWidth',3);
%% subplot 2
subplot2 = subplot(2,1,2,'Parent',figure1);
hold(subplot2,'on');
% Interface from MATLAB
plot(Interface.Bas(:,2),Interface.Bas(:,1),'c.','MarkerSize',10,'LineWidth',3,'Parent',subplot2);
% Interface from Basilisk
plot(Bas.Yf,Bas.Xf,'k-','MarkerSize',20,'LineWidth',3,'Parent',subplot2);
% Countor from Basilisk
% contourf(Y,X,f,[0.5 0.5],'facecolor', 'r','Parent',subplot2)
% Box
plot([yWindow1 yWindow2],[xWindow1 xWindow1],'--','Color',[0.5 0.5 0.5],'MarkerSize',20,'LineWidth',2,'Parent',subplot2);
plot([yWindow2 yWindow2],[xWindow1 xWindow2],'--','Color',[0.5 0.5 0.5],'MarkerSize',20,'LineWidth',2,'Parent',subplot2);
plot([yWindow1 yWindow2],[xWindow2 xWindow2],'--','Color',[0.5 0.5 0.5],'MarkerSize',20,'LineWidth',2,'Parent',subplot2);
plot([yWindow1 yWindow1],[xWindow1 xWindow2],'--','Color',[0.5 0.5 0.5],'MarkerSize',20,'LineWidth',2,'Parent',subplot2);
% Plot features
axis equal
ylim(subplot2,[xWindow1 xWindow2]);
xlim(subplot2,[yWindow1 yWindow2]);
% Box for next plot
xWindow1 = Interface.xc-2*abs(Interface.r); xWindow2 = Interface.xc+2*abs(Interface.r); 
yWindow1 = Interface.yc-2*abs(Interface.r); yWindow2 = Interface.yc+2*abs(Interface.r);
plot([yWindow1 yWindow2],[xWindow1 xWindow1],'-','Color',[0.75 0.75 0.75],'MarkerSize',20,'LineWidth',2,'Parent',subplot2);
plot([yWindow2 yWindow2],[xWindow1 xWindow2],'-','Color',[0.75 0.75 0.75],'MarkerSize',20,'LineWidth',2,'Parent',subplot2);
plot([yWindow1 yWindow2],[xWindow2 xWindow2],'-','Color',[0.75 0.75 0.75],'MarkerSize',20,'LineWidth',2,'Parent',subplot2);
plot([yWindow1 yWindow1],[xWindow1 xWindow2],'-','Color',[0.75 0.75 0.75],'MarkerSize',20,'LineWidth',2,'Parent',subplot2);
axis off

axes1 = axes('Position',[0.05 0.05 0.25 0.25]);
hold(axes1,'on');
% Interface from MATLAB
plot(Interface.Bas(:,2),Interface.Bas(:,1),'c-','MarkerSize',10,'LineWidth',3,'Parent',axes1);
% Interface from Basilisk
plot(Bas.Yf,Bas.Xf,'k-','MarkerSize',20,'LineWidth',3,'Parent',axes1);
% Contour from Basilisk
% contourf(Y,X,f,[0.5 0.5],'facecolor', 'r','Parent',axes1)
% Box
plot([yWindow1 yWindow2],[xWindow1 xWindow1],'-','Color',[0.75 0.75 0.75],'MarkerSize',20,'LineWidth',2,'Parent',axes1);
plot([yWindow2 yWindow2],[xWindow1 xWindow2],'-','Color',[0.75 0.75 0.75],'MarkerSize',20,'LineWidth',2,'Parent',axes1);
plot([yWindow1 yWindow2],[xWindow2 xWindow2],'-','Color',[0.75 0.75 0.75],'MarkerSize',20,'LineWidth',2,'Parent',axes1);
plot([yWindow1 yWindow1],[xWindow1 xWindow2],'-','Color',[0.75 0.75 0.75],'MarkerSize',20,'LineWidth',2,'Parent',axes1);
% Plot features
axis equal
ylim(axes1,[xWindow1 xWindow2]);
xlim(axes1,[yWindow1 yWindow2]);
axis off
annotation(figure1,'arrow',[0.323809523809524 0.250595238095238],...
    [0.433433433433433 0.304304304304304],'LineWidth',4,'HeadWidth',25,...
    'HeadLength',25);

% Saving files
name = [folder '/' sprintf('Bo%5.4f.mat',Bond)];
save(name,'Interface');
name = [folder '/' sprintf('Bo%5.4f.png',Bond)];
set(figure1,'pos',[1 1 1120 630]);
print(name,'-dpng','-r300')
toc
%% Necessary Functions
function Interface = shape(Bond, nParams)
warning('off','all'); % Turn on during debugging
% Bond - actual bond number based on equivalent radius of the drop/bubble. 
% RhoDS - the density ratio between the Drop/Bubble & surrounding
%% Initial Conditions
x0 = 1e-8; y0 = 0;
odeOptions = odeset('AbsTol',1e-10,'RelTol',1e-8);
%% Solution of the Bubble inside the Pool
phi = [linspace(0,(pi/2-pi/1e8),nParams.n11) linspace(pi/2,pi,nParams.n1)];
counter1 = 0; Err = 1.0; 
RbMax = nParams.RbGuess; RbMin = 1.0;
while (Err > nParams.VolTol && counter1 < nParams.VolMaxCounter)
    counter1 = counter1+1;
    Rb = (RbMax+RbMin)/2.0;
    [phi,y] = ode45(@eqns1,phi,[x0; y0], odeOptions, Bond, Rb);
    %% Solving the tail part of the interface
    counter2 = 0; TailErr = 1.0; phiMin = max(pi/2.0, 0.75*(pi - 0.5*sqrt(Bond/3))); phiMax = pi;
    while (TailErr > nParams.TailTol && counter2 < nParams.TailMaxCounter)
        phic = (phiMin+phiMax)/2.0; 
        % Choosing the initial condition
        [~,ind] = min(abs(phi-phic)); ratio = (phic - phi(ind))/(phi(ind+1)+phi(ind)); 
        % Choosing the critical point: Interpolation might be required to get
        % xc and yc: if the ratio defined above is non-zero
        xc = y(ind,1) + ratio*(y(ind+1,1)-y(ind,1)); yc = y(ind,2) + ratio*(y(ind+1,2)-y(ind,2));
        Rc = xc/sin(phic); % Radius of curvature at the critical point
        hinf = (2/Bond)*(2/Rc - 1/Rb); % Height of the meniscus when x -> \inf
        xTail = linspace(xc,nParams.TailxMax,nParams.nTail);
        IC = [yc, tan(phic)]; % IC for [y, y']    
        [xTail,yTail] = ode45(@eqnsTAIL,xTail,IC, odeOptions, Bond, hinf);
        TailErr = abs(yTail(end,1) - hinf);
        if (yTail(end,1) > hinf)
            phiMax = phic;
        else
            phiMin = phic;
        end
        if (rem(counter2,10) == 0)
            fprintf('Solving Tail: iteration = %d & Err = %4.3e\n', counter2, TailErr);
        end
        counter2 = counter2 + 1;
    end
    x2 = linspace(y(ind,1), 0, nParams.n2); y2 = yc + xc/tan(phic) + sqrt((xc/sin(phic))^2 - x2.^2);
    x1 = y(1:ind-1,1); y1 = y(1:ind-1,2);
    Vol1 = 0.75*(trapz(y1,x1.^2));
    Vol2 = 0.75*(trapz(y2,x2.^2));
    Vol = Vol1 + Vol2;
%     fprintf('Volume 1 = %f\n',Vol1);
%     fprintf('Volume 2 = %f\n',Vol2);
%     fprintf('Total Vol. = %f\n',Vol);
    Err = abs(1.0 - Vol);
    if (Vol > 1.0)
        RbMax = Rb;
    else
        RbMin = Rb;
    end
    fprintf('Radius of Curvature at the bottom = %5.4e at Iteration = %d and Error = %4.3e\n',Rb,counter1,Err);
end
Interface.Rb = Rb;
if (TailErr > 1e1*nParams.TailTol)
    plot(x1,y1,'r.','MarkerSize',20,'LineWidth',2);
    hold on;
    plot(x2,y2,'b.','MarkerSize',20,'LineWidth',2);
    plot(xTail,yTail(:,1),'g.','MarkerSize',20,'LineWidth',2);
    axis equal
    error('Tail did not converge. Error is %4.3e Check!\n',TailErr);
end
if (Err > 1e2*nParams.VolTol)
    plot(x1,y1,'r.','MarkerSize',20,'LineWidth',2);
    hold on;
    plot(x2,y2,'b.','MarkerSize',20,'LineWidth',2);
    plot(xTail,yTail(:,1),'g.','MarkerSize',20,'LineWidth',2);
    axis equal
    error('Vol did not converge. Error is %4.3e Check!\n', Err);
end
if (xTail(end) < 32)
    xTail = [xTail; linspace(xTail(end),32,500)'];
    yTail = [yTail(:,1); yTail(end,1)*ones(500,1)];
else
    yTail = yTail(:,1); 
end
%% Filleting the interfaces
XfilletTail = xTail(nParams.nfillet); YfilletTail = yTail(nParams.nfillet); 
[~,ind] = min(sqrt((x1-XfilletTail).^2+(y1-YfilletTail).^2));
XfilletBubble = x1(ind); YfilletBubble = y1(ind); PhiFillet = phi(ind);
r = 0.5*(sec(PhiFillet))*((XfilletTail-XfilletBubble)^2+(YfilletTail-YfilletBubble)^2)/((YfilletTail-YfilletBubble) - (XfilletTail-XfilletBubble)*tan(PhiFillet));
xc = XfilletBubble - r*sin(PhiFillet); yc = YfilletBubble + r*(cos(PhiFillet)); 
xfillet1 = linspace(XfilletBubble, XfilletTail, 2e3)'; yfillet1 = yc + sqrt(r^2 - (xfillet1 - xc).^2);
yfillet2 = linspace(YfilletBubble, YfilletTail, 2e3)'; xfillet2 = xc - sqrt(r^2 - (yfillet2 - yc).^2);
%% Concatenating all the interfaces together
Interface.Bas = [[y1(1:ind)-yTail(end) x1(1:ind)]; [yfillet2-yTail(end) xfillet2]; [yfillet1-yTail(end) xfillet1]; [yTail(nParams.nfillet:end)-yTail(end) xTail(nParams.nfillet:end)]]; 
%% Transforming the coordinates for Basilisk
Interface.x1 = y1-yTail(end); Interface.y1 = x1;
Interface.x2 = y2-yTail(end); Interface.y2 = x2;
Interface.xTail = yTail-yTail(end); Interface.yTail = xTail;
Interface.xc = yc-yTail(end); Interface.yc = xc;
Interface.r = r;
Interface.xfillet = [yfillet1-yTail(end);yfillet2-yTail(end)]; 
Interface.yfillet = [xfillet2;xfillet1];
end
%% Function defining equation handles for solution of the part of bubble inside the pool
function dydt = eqns1(phi,y,Bo,Rb)
dydt = [y(1)*cos(phi)/( y(1)*( 2/Rb + Bo*y(2) ) - sin(phi));...
    (y(1)*cos(phi)/( y(1)*( 2/Rb + Bo*y(2) ) - sin(phi)))*tan(phi)];
end
%% Function defining equation handles for solution of the tail part of the interface
function dydt = eqnsTAIL(x,y,Bo,hinf)
dydt = [y(2);...
    Bo*(y(1)-hinf)*(1+y(2)^2)^1.5 - (y(2)/x)*(1+y(2)^2)];
end

