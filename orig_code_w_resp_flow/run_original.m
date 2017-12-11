unix(sprintf('make clean'));
unix(sprintf('make -f Makefile-jam'));

rm = 0.005;  f1 = 0; f2 = 1; f3 = 260000;  fa1 = 0; fa2 = 1; fa3 = 50000;
asym = 0.41; expo = 2.76;  lrrA = 15.76; lrrV = 14.54;
verb = 0; nhb = 6;
id = 1;
str = sprintf('./sor06 %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g, %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g',...
                       rm,   f1,   f2,   f3,   fa1,  fa2,  fa3,  asym,  expo, lrrA, lrrV, verb, nhb,  id);
unix(str);

load -ascii pu1_1.2d;
[t,x1,p1,q1,~,~] = gnuplot(pu1_1);       % Sorting out data from 3D to 2D arrays
t =  t(:, 1);
p1 = p1(:,1);                  % P at x = 0
q1 = q1(:,1);                  % Q at x = 0
figure(1)
subplot(2, 1, 1);
plot(t, p1, 'LineWidth', 2)
set(gcf, 'Units', 'centimeters', 'Position', [5 5 30 20]);
xlabel('t [s]'); ylabel('p [mmHg]'); grid on;
subplot(2, 1, 2)
plot(t,q1,'LineWidth', 2);
grid on; xlabel('t [s]'); ylabel('q [cm^3/s]');
shg;