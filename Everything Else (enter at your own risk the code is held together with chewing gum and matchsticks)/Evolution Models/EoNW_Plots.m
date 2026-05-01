%% Plots for the evolution of niche width chapter in my dissertation

%% Parameters

%% Functions

function r = GrowthRate(x,y,c,R_0,w,s,T)
    r = R_0*w*s/(sqrt(s^2+y^2))*exp(-(c*y)-x.^2/(2*(s^2+y^2)))-T;
end

%% Plots

x = -5:0.1:5;
y = 0:0.1:10;