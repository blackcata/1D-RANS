clc
clear
close all 

filename = 'CHANNEL_0180_'

data1 = 'mean_prof'
data2 = 'vel_fluc_prof'

for k=1:30
   
    date = sprintf('%02i',k)
    
    fn1  = strcat(filename,data1,'.plt')
    fn2  = strcat(filename,data2,'.plt')

    data_p1 = fopen(fn1);
    data_p2 = fopen(fn2);
    
    A = fscanf(data_p1,'%f %f %f',[6 inf]); 
    B = fscanf(data_p2,'%f %f %f',[9 inf]); 
 
end
%% Data Transfer %%
Y    = A(1,:);
Yp   = A(2,:);
U    = A(3,:);
DUDY = A(4,:);
W    = A(5,:);
P    = A(6,:);

uu   = B(3,:);
vv   = B(4,:);
ww   = B(5,:);
uv   = B(6,:);
uw   = B(7,:);
vw   = B(8,:);
k    = B(9,:);


%% Plotting %%
semilogx(Yp,U,'ro')
xlabel('y+')
ylabel('u+')
axis([1 max(Yp)+20 0 max(U)+2])

%% %%

Re_tau = 180
nu = 3.5e-4
del = 1.0

u_tau = Re_tau * nu /del

k_max = max(k*u_tau*u_tau) 

%% %%
u_f = sqrt(uu) * u_tau
v_f = sqrt(vv) * u_tau
w_f = sqrt(ww) * u_tau

for it = 2:length(uu)-1
    dy = (Y(it+1) -Y(it-1))/2;
    diss_u = (u_f(it+1) - u_f(it-1))/(2*dy);
    diss_v = (v_f(it+1) - v_f(it-1))/(2*dy);
    diss_w = (w_f(it+1) - w_f(it-1))/(2*dy);
    diss(it)= nu*(diss_u + diss_v + diss_w);
    dissp(it) = diss(it) /(nu *(u_tau^2/nu)^2);
end

semilogx(Yp(1:95),dissp)
xlabel('y+')
ylabel('diss')
axis([1 max(Yp)+20 0 max(dissp)])

%% Writing Y,Yp,U,k %%
fileID = fopen('CHANNEL_0180_profile.plt','w');
fprintf(fileID,'VARIABLES="Y","Y+","U","k","k+","diss","diss+" \n');
fprintf(fileID,'%f %f %f %f \n',[Y' Yp' U' k' (k'*u_tau^2), diss', dissp']');
fclose(fileID);

