% L1 adaptive control motivational example
% Saving the time delay margin

s = tf('s');
gam = [10 100 1000];
mrac = [];
L1 = [];
w = 1;
C = w/(s+w);
for i=1:length(gam)
    mrac = [mrac gam(i)/(s*(s+1))]; % integral action (perfect mrac)
    L1 =  [L1 gam(i)*C/(s*(s+1) + gam(i)*(1-C))];
end


figure()
hold on
for i=1:length(gam)
    nyquist(mrac(i));
end
title('MRAC');
circle(0,0,1)
xlim([-1.1 0.1]);
ylim([-1.1 0.1]);


figure()
hold on
for i=1:length(gam)
    nyquist(L1(i));
end
title('L1 Adaptive Control')
circle(0,0,1)
xlim([-1.1 0.1]);
ylim([-1.1 0.1]);

