function [x,y] = makeFoil(shape,camber,fileSaveName)

%% DEFINE FOIL SHAPE
load('foilLibrary.mat','A_mean','U_a','S_a','V_a','sampleX')

[~,xloc] = min(abs(sampleX-camber.loc)); %nearest integer where camber pieces meet
camber.y1 = (camber.max/camber.loc^2) * (2*camber.loc*sampleX(1:xloc) - sampleX(1:xloc).^2);                        %camber piece 1 (front)
camber.y2 = camber.max/(1-camber.loc)^2*((1-2*camber.loc)+2*camber.loc*sampleX(xloc+1:end)-sampleX(xloc+1:end).^2); %camber piece 2 (back)
camber.y = [camber.y1 camber.y2];

coefficients = [shape.A shape.B shape.C shape.D].*[0.0008 0.0023 -0.0126 -0.0038]; 

%% ENSURE FOIL IS PHYSICAL (SURFACE CAN'T CROSS CHORDLINE)
U_b = U_a*S_a;
U = U_b(:,1:length(coefficients));
S = S_a(1:length(coefficients),1:length(coefficients));
tmin = 0.000;
options = optimset('Display','off');
% coefficientsN = quadprog(S'*S,-coefficients*S'*S,-U*S,A_mean,[],[],[],[],[],options)';
coefficientsN = quadprog(U'*U,-U'*U*coefficients',-U,A_mean - tmin);

modeNum = length(coefficientsN);
    
%% MAKE FOIL
foil_SVD = zeros(length(sampleX),1);
for k = 1:modeNum
    mode = U_a(:,k)*S_a(k,k)*coefficientsN(k)';
    mode(1) = 0;
    mode(end) = 0;
    foil_SVD = foil_SVD + mode;
end
A_mean(1) = 0;
A_mean(end) = 0;
foil_SVD = foil_SVD + A_mean;

%% PLOT FOIL
plot(sampleX,foil_SVD'+camber.y,'.-r')
hold on
plot(sampleX,-foil_SVD'+camber.y,'.-r')
set(gca,'fontname','Source Sans Pro Light','fontsize',12,'ygrid','on','xgrid','on')
xlim([0 1])
set(gcf,'position',[520   600   560   198])
plot(sampleX,camber.y,'--k')
axis equal

%% SAVE FOIL
x = [sampleX(end:-1:1) sampleX];
top = foil_SVD'+camber.y;
bottom = -foil_SVD'+camber.y;
y = [top(end:-1:1) bottom];
mat = [x' y'];
save(fileSaveName,'mat','-ascii')