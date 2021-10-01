% set up reference system
b2 = 1;
a0 = 1;
a1 = 2;
a2 = 3;
system = tf([b2],[a0 a1 a2]); 
step(system) % validate against this solution


% numerical sim
dt = 0.00001;
t = 0:dt:6;
y = zeros(length(t),1);

for i=3:length(t)
   y(i) = second_order_filter(a0,a1,a2,b2,dt,1,1,1,y(i-1),y(i-2));
end

hold on
plot(t,y,'r--');