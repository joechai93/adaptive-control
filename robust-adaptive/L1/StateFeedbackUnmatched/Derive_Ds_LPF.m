clearvars
sym s;
s = tf('s');
tf_num = 1;
tf_den = (s*((s/25)+1)*((s/70)+1)*((s*s/(40*40)) + (1.8*s/40) + 1));
lpf_tf = tf_num/tf_den;
num = lpf_tf.Numerator{1};
den = lpf_tf.Denominator{1};

[A,B,C,D] = tf2ss(num,den);

%%
clearvars
sym s;
s = tf('s');
w = 50;
K = 1;
C = w*K*inv(s + w*K)