tRange=[0 0.1];
Y0=[0 0 6.5*10^5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
[tSol,YSol] = ode45(@odefun,tRange,Y0);

y1=YSol(:,17);
plot(tSol,y1);
hold on;
y2=YSol(:,16);
plot(tSol,y2);
hold on;
y3=YSol(:,12);
plot(tSol,y3);
hold on;
y4=YSol(:,15);
plot(tSol,y4);


xlabel('Time');
ylabel('Cellulose');
legend('c1','c2','c6','c15');

function dYdt = odefun(t,Y)
%Ci-soluble cello-oligosaccharide chain of length i (i < 7) 
eg=Y(1);%eg-endoglucanase 
cbh=Y(2);%cbh-cellobiohydrolase 
c15=Y(3);
c14=Y(4);
c13=Y(5);
c12=Y(6);
c11=Y(7);
c10=Y(8);
c9=Y(9);
c8=Y(10);
c7=Y(11);
c6=Y(12);
c5=Y(13);
c4=Y(14);
c3=Y(15);
c2=Y(16);
c1=Y(17);
r=Y(18);
a=19403;
b=6.65*10^5;%*max-free sites
c=2452;
d=2937;
e=7.8*10^7;
ri = [1,2,3,4,5,6];
Sum = 0;
for i = 1:6
%ri-rate of formation of a length i soluble sugar from solid
substrate
Sum = Sum + ri(i) * i;
end
%kEG2-ads=8640
%kx-ads=adsorption constant for a cellulase of type x onto a cellulose
surface
%kCBH1-ads=8640
%SA_cellulose surface area
%MW-cellulose monomer molecular weight 
SA = 0.162; MW = 8;
drdt = -r*(SA*MW/3)* Sum;
%EEG2-total=1.6(enzyme loadings set) 
%ECBH2-total=1.6(Enzyme loadings set) 
degdt=8640*(1.6-eg)*(6.6*10^-5)-19.3*eg;
%degdt=kcat-EG-cell/kM-EG-cell*(EG-total-eg)*(*max)-kEG-des*eg; 
dcbhdt=8640*(1.6-cbh)*(6.6*10^-5)-164*cbh;
%dcbhdt=kcat-CBH-cell/kM-CBH-cell*(CBH-total-cbh)*(*max)-kCBH-des*cbh;
dc15dt=a*eg*(-14*c15)-b*cbh*c15;
dc14dt=a*eg*(c15-13*c14)-b*cbh*c14;
dc13dt=a*eg*(c14+c15-12*c13)+b*cbh*(c15-c13);
dc12dt=a*eg*(c13+c14+c15-11*c12)+b*cbh*(c14-c12);
dc11dt=a*eg*(c12+c13+c14+c15-10*c11)+b*cbh*(c13-c11);
dc10dt=a*eg*(c11+c12+c13+c14+c15-9*c10)+b*cbh*(c12-c10);
dc9dt=a*eg*(c10+c11+c12+c13+c14+c15-8*c9)+b*cbh*(c11-c9);
dc8dt=a*eg*(c9+c10+c11+c12+c13+c14+c15-7*c8)+b*cbh*(c10-c8);
dc7dt=a*eg*(c8+c9+c10+c11+c12+c13+c14+c15-6*c7)+b*cbh*(c9-c7);
dc6dt=(a*eg*(c7+c8+c9+c10+c11+c12+c13+c14+c15)+(b*cbh*c8)+c*eg*(c6+c5+c4+c3-5*c6)-d*cbh*c6)*r*r+c*(1.6-eg)*(c6+c5+c4+c3-5*c6)-d*(1.6-cbh)*c6;
dc5dt=(a*eg*(c7+c8+c9+c10+c11+c12+c13+c14+c15)+(b*cbh*c7)+c*eg*(c6+c5+c4+c3-4*c5)-d*cbh*c5)*r*r+c*(1.6-eg)*(c6+c5+c4+c3-4*c5)-d*(1.6-cbh)*c5;
dc4dt=(a*eg*(c7+c8+c9+c10+c11+c12+c13+c14+c15)+c*eg*(c6+c5+c4-3*c4)+d*cbh*(c6-c4))*r*r+c*(1.6-eg)*(c6+c5+c4-3*c4)+d*(1.6-cbh)*(c6-c4);
dc3dt=(a*eg*(c7+c8+c9+c10+c11+c12+c13+c14+c15)+c*eg*(c6+c5+c4-2*c3)+d*cbh*(c5-c3))*r*r+c*(1.6-eg)*(c6+c5+c4-2*c3)+d*(1.6-cbh)*(c5-c3);
dc2dt=(a*eg*(c7+c8+c9+c10+c11+c12+c13+c14+c15)+b*cbh*(c7+c8+c9+c10+c11+c12+c13+c14+c15)+c*eg*(c6+c5+c4+c3)+d*cbh*(c6+c5+2*c4+c3))*r*r+c*(1.6-eg)*(c6+c5+c4+c3)+d*(1.6-cbh)*(c6+c5+2*c4+c3);
dc1dt=(a*eg*(c7+c8+c9+c10+c11+c12+c13+c14+c15)+d*cbh*c3+eg*(c6+c5+c4+c3))*r*r+d*(1.6-cbh)*c3+(1.6-eg)*(c6+c5+c4+c3);

dYdt=[degdt; dcbhdt; dc15dt; dc14dt; dc13dt; dc12dt; dc11dt; dc10dt; dc9dt; dc8dt; dc7dt; dc6dt; dc5dt; dc4dt; dc3dt; dc2dt; dc1dt; drdt];
end
