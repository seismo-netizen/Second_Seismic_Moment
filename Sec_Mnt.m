clear; clc; close all;

strike = 0; dip = 90;

Lc = 1; % km along strike
Wc = Lc/2; % km along dip
Vr = 3.5613*0.9; % rupture velocity


nn = 1000;
x = linspace(-Lc,Lc,nn); % along strike
y = linspace(-Wc,Wc,round(nn/Lc*Wc)); % along dip
dx = x(2)- x(1);
dy = y(2)- y(1);

t = linspace(0,2*Lc/Vr,nn);
dt = t(2) - t(1);


%% starting point
xc = -Lc; yc = 0;
%xc = 0; yc = 0;

% coef matrix
W = zeros(length(x),length(y)); % scaled slip on each grid
T = W; % rupture timing
for i = 1:length(x)
    for j = 1:length(y)
        ra2 = (x(i)/Lc)^2 + (y(j)/Wc)^2;
        if ra2 <= 1
            W(i,j) = 1;
            dist = sqrt((x(i) - xc)^2 + (y(j) - yc)^2);
            T(i,j) = dist/Vr;
        else
            T(i,j) = nan;
        end
    end
end

%%
x_elps1 = x; x_elps2= flip(x);
for i = 1:length(x)
    y_elps1(i) = Wc*sqrt(1-(x(i)/Lc)^2); 
end
y_elps2 = -flip(y_elps1);
x_elps = [x_elps1 x_elps2];
y_elps = [y_elps1 y_elps2];

figure
plot(x_elps,y_elps,'k','linewidth',2); hold on; colorbar; caxis([0 0.6])
contour(x,y,T',12,'linewidth',2); shading interp
axis equal; colormap(jet)
xlabel('Along Strike (km)');
ylabel('Along Dip (km)');
set(gca,'Fontsize',20)


%% intigrate spatiotemporal function
% assuming STF at each grid is a Drac delta
sum_f = 0;
for i = 1:length(x)
    for j = 1:length(y)
        if W(i,j) == 1
            sum_f = sum_f + dx*dy;
        end
    end
end


%% centroid
x0 = 0; y0 = 0; t0 = 0;
for i = 1:length(x)
    for j = 1:length(y)
        if W(i,j) == 1
            x0 = x0 + x(i)*dx*dy/sum_f;
            y0 = y0 + y(j)*dx*dy/sum_f;
            dist = sqrt((x(i) - xc)^2 + (y(j) - yc)^2);
            t0 = t0 + (dist/Vr)*dx*dy/sum_f;
        end
    end
end

%% 2nd moment
xx = 0; xy = 0; yy = 0;
xt = 0; yt = 0;
tt = 0;
for i = 1:length(x)
    for j = 1:length(y)
        if W(i,j) == 1
            
            xx = xx + (x(i)-x0)^2*dx*dy/sum_f;
            yy = yy + (y(j)-y0)^2*dx*dy/sum_f;
            xy = xy + (x(i)-x0)*(y(j)-y0)*dx*dy/sum_f;
            
            dist = sqrt((x(i) - xc)^2 + (y(j) - yc)^2);
            
            xt = xt + (x(i)-x0)*(dist/Vr-t0)*dx*dy/sum_f;
            yt = yt + (y(j)-y0)*(dist/Vr-t0)*dx*dy/sum_f;
            
            tt = tt + (dist/Vr-t0)^2*dx*dy/sum_f;
        end
    end
end

m2 = [tt xt yt xx xy yy]; 

X = [m2(4), m2(5); m2(5), m2(6);];
[U,S,V]=svd(X);
L_c=2*sqrt(max(max(S))); W_c=2*sqrt(S(2,2));
v0=m2(2:3)/m2(1); mv0=sqrt(sum(v0.^2));
v0_up = -v0(2)*sin(dip/180*pi);
v0_east = v0(1)*sin(strike/180*pi) + v0(2)*cos(dip/180*pi)*cos(strike/180*pi);
v0_north = v0(1)*cos(strike/180*pi) - v0(2)*cos(dip/180*pi)*sin(strike/180*pi);
tauc=2*sqrt(m2(1)); L0=tauc*mv0; dir=L0/L_c;
display(L_c); display(W_c); display(tauc); display(mv0); display(dir);
display(v0); display(v0_up); display(v0_east); display(v0_north);




