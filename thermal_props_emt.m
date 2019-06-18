k0_p = 0.4795;
k1_p = 0.001923
a0_p = 0.1329;
a1_p = 0.00011;
k0_m = 0.4994;
k1_m = 0.001102;
a0_m = 0.1278;
a1_m = 0.00055;
k0_c = 0.4989;
k1_c = 0.001288;
a0_c = 0.1266;
a1_c = 0.00055;
T = linspace(3,45,7);
k_p = k0_p + (k1_p)*T;
k_m =  k0_m + (k1_m)*T;
k_c =  k0_c + (k1_c)*T;
a_p =  a0_p + (a1_p)*T;
a_p = a_p/(1000)^2;
a_c =  (a0_c + (a1_c)*T)/(1000)^2;
a_m =  (a0_m + (a1_m)*T)/(1000)^2;

C_k=3924;
C_w=4208;


plot(T,k_p,'-o','LineWidth',0.5,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b', 'MarkerSize',7);
hold on
plot(T,k_m,'-s','LineWidth',0.5,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r', 'MarkerSize',7);
hold on
plot(T,k_c,'-^','LineWidth',0.5,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','g', 'MarkerSize',7);
figure
plot(T,a_p,'-ob');
hold on
plot(T,a_m,'-sr');
plot(T,a_c,'-^g');

%thermal conductivity of water from 0C to 100C
tw = [0 10 20 30 40 50 60 70 80 90 100];
kw = [0.561 0.580 0.5984 0.6154 0.6305 0.6435 0.6543 0.6631 0.67 0.6753 0.6791];
Km = k0_m + (k1_m)*tw;
Kc =  k0_c + (k1_c)*tw;
Kp = k0_p + (k1_p)*tw;
%assume equal densities of water and dry material
%volume fraction of dry material 
v2 = 0.173

%volume fraction of water
v1 =0.827

k2m = ((-(v1/v2)*((kw-Km)/(kw+2*Km))*2*Km)+Km)/(1+(v1/v2)*((kw-Km)/(kw+2*Km)));
k2c = ((-(v1/v2)*((kw-Kc)/(kw+2*Kc))*2*Kc)+Kc)/(1+(v1/v2)*((kw-Kc)/(kw+2*Kc)));
k2p = ((-(v1/v2)*((kw-Kp)/(kw+2*Kp))*2*Kp)+Kp)/(1+(v1/v2)*((kw-Kp)/(kw+2*Kp)));

c2m = ((-(v1/v2)*((C_w-C_k)/(C_w+2*C_k))*2*C_k)+C_k)/(1+(v1/v2)*((C_w-C_k)/(C_w+2*C_k)))


figure(1)
hold on
plot(tw,k2p,'-b')
plot(tw,k2m,'-r')
plot(tw,k2c,'-g')

hold on
plot(tw,kw,'k')

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, 6.83, 6]);
set(gcf, 'PaperSize', [6.83 6]);
set(gca,'Fontsize',12,'linewidth',0.5);
xlim([0 55]);
ylim([0 1]);
box on
colordef white
hold on
ylabel('Thermal Conductivity');
xlabel('Temperature, {\circ}C');
h_legend = legend('k_{measured}, pelvis','k_{measured}, medulla', 'k_{measured}, cortex','k_{calculated}, pelvis dry material', 'k_{calculated}, medulla dry material','k_{calculated}, cortex dry material','k_{measured}, water')
legend('boxoff')
set(h_legend,'FontSize',12);

load ('C:\Users\Student\Desktop\3omega\Transient Hot Wire\M22\TH_M22_s44_06_28_16_T4.mat')

Gs = vertcat(G(1:48,:),G(66:end,:));
km22f = polyfit(Gs(:,4),Gs(:,1),4)
km22 = polyval(km22f,G(:,4))
figure
plot(G(:,4),G(:,1),'ob')
hold on

plot(G(:,4),km22,'-r')
hold on

xlim([-140 55]);
ylim([0 1]);

plot(tw,k2p,'-b')
plot(tw,k2m,'-r')
plot(tw,k2c,'-g')

