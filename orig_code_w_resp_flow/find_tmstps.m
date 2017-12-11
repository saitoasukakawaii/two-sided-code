nhb = 8; % number of heartbeats in a breathing cycle
Tper = 0.7;
Tper = Tper*nhb;

tmstps = 4096;
tmstps = tmstps*nhb

n=1;
while (2^n)<tmstps
n=n+1;
end
n;
tmstps;
tmstps = 2^n + 1

A = (importdata('gvPA_4096.dat'));
A = [A; A; A; A; A; A; A; A]; plot(A); shg
A = A(1:(tmstps-1));
plot(A)
save gvPA_3_16384.dat A -ascii
figure(2);plot(importdata('gvPA_3_16384.dat')); shg
zero = zeros([(length(A)-1)*2, 1]); 
save zero_3.dat zero -ascii
length(importdata('zero_3.dat'))