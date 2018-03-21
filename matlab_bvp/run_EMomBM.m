%  Script to set up and run a steady-state version of the EMomBM
%  Brian Rose
%
%  loop over an array of ice edges, solve for the constant and also
%  the wind stress.

%  draws a plot of wind stresses versus latitude and associated ice extent

Bout = 1.7;
Aout = 212;
epsilon = 0.04;
Omega = 2*pi/3600/24;
f0 = 1E-4;
R = 287;
a = 6.373E6;
SD = 40;
Ld = 5E5;
Tf = -10;
p0=1E5;
Ca=1E7;
Ka=2.2E6;
g=9.8;

param.s2 = -0.48;
param.ao = [0.7, -0.078];
param.ai = 0.38;
param.ltheta = Ka*Ca/Bout/a^2;
param.A1 = Aout/Bout/SD;
param.Tf = Tf/SD;
param.s = R*SD/f0/2/Omega/Ld^2;
param.ltau = Ka*p0/g/epsilon/a^2;
param.ld = (Ld/a)^2/2;
param.solvetau=true;
param.qinitial = 341.5/Bout/SD;
param.kappaguess = 0.2;

icelat = [90;80;70;60;50;40];
xsarray = sin(icelat/180*pi);

for n=1:length(xsarray)
    solution(n) = YetAnotherEMomBM(xsarray(n),param);
    tau(n,:) = solution(n).taustar.*sqrt(1-solution(n).x.^2)*2*Omega*a*epsilon;
    lat(n,:) = solution(n).lat;
    ices(n,:) = [icelat(n) 90];
end
ices(1,1) = 89;
offsets = repmat((-0.1:-0.025:-0.225)',[1 2]);


colors=[1 0 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1; 0 0 0];
figure
hold on
for n=length(xsarray):-1:1
    plot(lat(n,:),tau(n,:),'LineWidth',2,'Color',colors(n,:))
    plot(ices(n,:),offsets(n,:),'LineWidth',2,'Color',colors(n,:))
end
hold off
grid on
set(gca,'XLim',[0 90],'XTick',0:15:90,'YLim',[-0.25,0.25],'FontSize',14)
xlabel('Latitude')
ylabel('Wind stress (N/m^2)')

