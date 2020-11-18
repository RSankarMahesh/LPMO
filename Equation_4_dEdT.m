[t,Eads] = ode45(@odefun,[0 10], 1);% time value can be arbitrary
plot(Eads,t);
function Eads = odefun(t,Eads)
Kxads = 8640 ;Efree = 10;f = 4;Kxdes = 19.3; % Let x be EG2; Efree and 'f' - free sites here is an arbitrary value, 
dEads = Kxads*Efree*f - Kxdes*Eads;
end
