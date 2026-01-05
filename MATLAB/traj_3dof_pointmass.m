%% RAVEN – 3DOF Point Mass Test
clear; clc; close all;

t = linspace(0,10,100);
x = t.^2;
y = t;
z = sin(t);

figure;
plot3(x,y,z,'LineWidth',2);
grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('RAVEN sanity check');

