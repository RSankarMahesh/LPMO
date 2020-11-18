[t,r] = ode45(@odefun,[0 50], 1);
plot(t,r);
%time is arbitrary for simulation purpose I've fixed it as 100
function dr = odefun(t,r)
ri = [1,2,3,4,5,6];% this value is random and represents rate of formation of a length i soluble sugar from solid substrate
Sum = 0;
for i = 1:6
Sum = Sum + ri(i) * i;
end
SA = 0.162; MW = 8;
dr = r*(SA*MW/3)* Sum; %r = radius of the particle ; SA = the surface area per mass of the cellulose; MW = Molecular weight of anhydrous glucose
end
