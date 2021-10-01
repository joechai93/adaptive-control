% Computing the inverse Hm in laplace domain
% Computing the State space realisation of inv(Hm)*Hum
% SISO only at the moment

clearvars;
sym s;
s = tf('s');

Am = [-10 -40 -80;
       1 0 0;
       0 1 0]; 
Bm = [1;
      0;
      0]; 
Bum = [0 0;
       1 0;
       0 1];
c = [0 0 1];

n = length(Am);
assert(norm(Bm'*Bum) == 0);
assert(rank([Bm Bum]) == n);

I = eye(n);

Hm = c * inv(s*I - Am) * Bm;
Hm_i = inv(Hm);

Hum = c * inv(s*I - Am) * Bum;

combined = Hm_i * Hum;
num = combined.Numerator{1};
den = combined.Denominator{1};


% State space realisation of Hm_i * Hum
[A,B,C,D] = tf2ss(num,den);

%% Same but for 2 state simple system
clearvars;
sym s;
s = tf('s');

Am = [-10 -40;
       1 0]; 
Bm = [1;
      0]; 
c = [0 1];
kg = -inv(c*inv(Am)*Bm);
Bum = [0; 
       2];
   
n = length(Am);
assert(norm(Bm'*Bum) == 0);
assert(rank([Bm Bum]) == n);

I = eye(n);

Hm = c * inv(s*I - Am) * Bm;
Hm_i = inv(Hm);

Hum = c * inv(s*I - Am) * Bum;

combined = Hm_i * Hum;
num = combined.Numerator{1};
den = combined.Denominator{1};


% State space realisation of Hm_i * Hum
[A,B,C,D] = tf2ss(num,den);