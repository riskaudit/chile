clear, clc, close

%% load data
GSR_10m = readtable("GridShare_Rel_ExpM.xlsx");
GSR_2m = readtable("GridShare_Rel_ExpM_HighRes.xlsx");
GSR_90m = readtable("GridShare_Rel_WSF_ExpM.xlsx");

%% HBET_1_2

% CR_1_3
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.CR_1_3(i,1) = GSR_2m.CR_1_3(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.CR_1_3(i,2) = GSR_10m.CR_1_3(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.CR_1_3(i,3) = GSR_90m.CR_1_3(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
[R.CR_1_3, P.CR_1_3, RL.CR_1_3, RU.CR_1_3] = corrcoef(GSR_combined.CR_1_3);

[f,xi] = ksdensity(GSR_combined.CR_1_3(:,1)); 
figure
plot(xi,f,'LineWidth',2);
hold on
[f,xi] = ksdensity(GSR_combined.CR_1_3(:,2)); 
plot(xi,f,'LineWidth',2);
[f,xi] = ksdensity(GSR_combined.CR_1_3(:,3)); 
plot(xi,f,'LineWidth',2);
hold off
grid on
legend('2m','10m','90m','Location','best');
xlabel('Number of Buildings')
ylabel('PDF')
title('CR-1-3')

% MATO_1_2
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MATO_1_2(i,1) = GSR_2m.MATO_1_(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MATO_1_2(i,2) = GSR_10m.MATO_1_(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MATO_1_2(i,3) = GSR_90m.MATO_1_(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
[R.MATO_1_2, P.MATO_1_2, RL.MATO_1_2, RU.MATO_1_2] = corrcoef(GSR_combined.MATO_1_2);

[f,xi] = ksdensity(GSR_combined.MATO_1_2(:,1)); 
figure
plot(xi,f,'LineWidth',2);
hold on
[f,xi] = ksdensity(GSR_combined.MATO_1_2(:,2)); 
plot(xi,f,'LineWidth',2);
[f,xi] = ksdensity(GSR_combined.MATO_1_2(:,3)); 
plot(xi,f,'LineWidth',2);
hold off
grid on
legend('2m','10m','90m','Location','best');
xlabel('Number of Buildings')
ylabel('PDF')
title('MATO-1-2')

% MCF_MR_1_2
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MCF_MR_1_2(i,1) = GSR_2m.MCF_MR_1(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MR_1_2(i,2) = GSR_10m.MCF_MR_1(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MR_1_2(i,3) = GSR_90m.MCF_MR_1(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
[R.MCF_MR_1_2, P.MCF_MR_1_2, RL.MCF_MR_1_2, RU.MCF_MR_1_2] = corrcoef(GSR_combined.MCF_MR_1_2);

[f,xi] = ksdensity(GSR_combined.MCF_MR_1_2(:,1)); 
figure
plot(xi,f,'LineWidth',2);
hold on
[f,xi] = ksdensity(GSR_combined.MCF_MR_1_2(:,2)); 
plot(xi,f,'LineWidth',2);
[f,xi] = ksdensity(GSR_combined.MCF_MR_1_2(:,3)); 
plot(xi,f,'LineWidth',2);
hold off
grid on
legend('2m','10m','90m','Location','best');
xlabel('Number of Buildings')
ylabel('PDF')
title('MCF-MR-1-2')

% MCF_MCF_1_2
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MCF_MCF_1_2(i,1) = GSR_2m.MCF_MCF_1(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MCF_1_2(i,2) = GSR_10m.MCF_MCF_1(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MCF_1_2(i,3) = GSR_90m.MCF_MCF_1(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
[R.MCF_MCF_1_2, P.MCF_MCF_1_2, RL.MCF_MCF_1_2, RU.MCF_MCF_1_2] = corrcoef(GSR_combined.MCF_MCF_1_2);

[f,xi] = ksdensity(GSR_combined.MCF_MCF_1_2(:,1)); 
figure
plot(xi,f,'LineWidth',2);
hold on
[f,xi] = ksdensity(GSR_combined.MCF_MCF_1_2(:,2)); 
plot(xi,f,'LineWidth',2);
[f,xi] = ksdensity(GSR_combined.MCF_MCF_1_2(:,3)); 
plot(xi,f,'LineWidth',2);
hold off
grid on
legend('2m','10m','90m','Location','best');
xlabel('Number of Buildings')
ylabel('PDF')
title('MCF-MCF-1-2')

% MR_1_2
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MR_1_2(i,1) = GSR_2m.MR_1_2(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MR_1_2(i,2) = GSR_10m.MR_1_2(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MR_1_2(i,3) = GSR_90m.MR_1_2(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
[R.MR_1_2, P.MR_1_2, RL.MR_1_2, RU.MR_1_2] = corrcoef(GSR_combined.MR_1_2);

[f,xi] = ksdensity(GSR_combined.MR_1_2(:,1)); 
figure
plot(xi,f,'LineWidth',2);
hold on
[f,xi] = ksdensity(GSR_combined.MR_1_2(:,2)); 
plot(xi,f,'LineWidth',2);
[f,xi] = ksdensity(GSR_combined.MR_1_2(:,3)); 
plot(xi,f,'LineWidth',2);
hold off
grid on
legend('2m','10m','90m','Location','best');
xlabel('Number of Buildings')
ylabel('PDF')
title('MR-1-2')

% MURADO_1_2
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MURADO_1_2(i,1) = GSR_2m.MURADO_(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MURADO_1_2(i,2) = GSR_10m.MURADO_(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MURADO_1_2(i,3) = GSR_90m.MURADO_(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
[R.MURADO_1_2, P.MURADO_1_2, RL.MURADO_1_2, RU.MURADO_1_2] = corrcoef(GSR_combined.MURADO_1_2);

[f,xi] = ksdensity(GSR_combined.MURADO_1_2(:,1)); 
figure
plot(xi,f,'LineWidth',2);
hold on
[f,xi] = ksdensity(GSR_combined.MURADO_1_2(:,2)); 
plot(xi,f,'LineWidth',2);
[f,xi] = ksdensity(GSR_combined.MURADO_1_2(:,3)); 
plot(xi,f,'LineWidth',2);
hold off
grid on
legend('2m','10m','90m','Location','best');
xlabel('Number of Buildings')
ylabel('PDF')
title('MURADO-1-2')

% MUR_1_2
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MUR_1_2(i,1) = GSR_2m.MUR_1_2(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MUR_1_2(i,2) = GSR_10m.MUR_1_2(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MUR_1_2(i,3) = GSR_90m.MUR_1_2(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
[R.MUR_1_2, P.MUR_1_2, RL.MUR_1_2, RU.MUR_1_2] = corrcoef(GSR_combined.MUR_1_2);

[f,xi] = ksdensity(GSR_combined.MUR_1_2(:,1)); 
figure
plot(xi,f,'LineWidth',2);
hold on
[f,xi] = ksdensity(GSR_combined.MUR_1_2(:,2)); 
plot(xi,f,'LineWidth',2);
[f,xi] = ksdensity(GSR_combined.MUR_1_2(:,3)); 
plot(xi,f,'LineWidth',2);
hold off
grid on
legend('2m','10m','90m','Location','best');
xlabel('Number of Buildings')
ylabel('PDF')
title('MUR-1-2')

% W_1_3
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.W_1_3(i,1) = GSR_2m.W_1_3(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.W_1_3(i,2) = GSR_10m.W_1_3(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.W_1_3(i,3) = GSR_90m.W_1_3(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
[R.W_1_3, P.W_1_3, RL.W_1_3, RU.W_1_3] = corrcoef(GSR_combined.W_1_3);

[f,xi] = ksdensity(GSR_combined.W_1_3(:,1)); 
figure
plot(xi,f,'LineWidth',2);
hold on
[f,xi] = ksdensity(GSR_combined.W_1_3(:,2)); 
plot(xi,f,'LineWidth',2);
[f,xi] = ksdensity(GSR_combined.W_1_3(:,3)); 
plot(xi,f,'LineWidth',2);
hold off
grid on
legend('2m','10m','90m','Location','best');
xlabel('Number of Buildings')
ylabel('PDF')
title('W-1-3')

% WDNO_1_2
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.WDNO_1_2(i,1) = GSR_2m.WDNO_1_(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.WDNO_1_2(i,2) = GSR_10m.WDNO_1_(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.WDNO_1_2(i,3) = GSR_90m.WDNO_1_(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
[R.WDNO_1_2, P.WDNO_1_2, RL.WDNO_1_2, RU.WDNO_1_2] = corrcoef(GSR_combined.WDNO_1_2);

[f,xi] = ksdensity(GSR_combined.WDNO_1_2(:,1)); 
figure
plot(xi,f,'LineWidth',2);
hold on
[f,xi] = ksdensity(GSR_combined.WDNO_1_2(:,2)); 
plot(xi,f,'LineWidth',2);
[f,xi] = ksdensity(GSR_combined.WDNO_1_2(:,3)); 
plot(xi,f,'LineWidth',2);
hold off
grid on
legend('2m','10m','90m','Location','best');
xlabel('Number of Buildings')
ylabel('PDF')
title('WDNO-1-2')

%% HBET_3

% MCF_MR_3
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MCF_MR_3(i,1) = GSR_2m.MCF_MR_3(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MR_3(i,2) = GSR_10m.MCF_MR_3(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MR_3(i,3) = GSR_90m.MCF_MR_3(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
[R.MCF_MR_3, P.MCF_MR_3, RL.MCF_MR_3, RU.MCF_MR_3] = corrcoef(GSR_combined.MCF_MR_3);

[f,xi] = ksdensity(GSR_combined.MCF_MR_3(:,1)); 
figure
plot(xi,f,'LineWidth',2);
hold on
[f,xi] = ksdensity(GSR_combined.MCF_MR_3(:,2)); 
plot(xi,f,'LineWidth',2);
[f,xi] = ksdensity(GSR_combined.MCF_MR_3(:,3)); 
plot(xi,f,'LineWidth',2);
hold off
grid on
legend('2m','10m','90m','Location','best');
xlabel('Number of Buildings')
ylabel('PDF')
title('MCF-MR-3')

% MCF_MCF_3
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MCF_MCF_3(i,1) = GSR_2m.MCF_MCF_3(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MCF_3(i,2) = GSR_10m.MCF_MCF_3(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MCF_3(i,3) = GSR_90m.MCF_MCF_3(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
[R.MCF_MCF_3, P.MCF_MCF_3, RL.MCF_MCF_3, RU.MCF_MCF_3] = corrcoef(GSR_combined.MCF_MCF_3);

[f,xi] = ksdensity(GSR_combined.MCF_MCF_3(:,1)); 
figure
plot(xi,f,'LineWidth',2);
hold on
[f,xi] = ksdensity(GSR_combined.MCF_MCF_3(:,2)); 
plot(xi,f,'LineWidth',2);
[f,xi] = ksdensity(GSR_combined.MCF_MCF_3(:,3)); 
plot(xi,f,'LineWidth',2);
hold off
grid on
legend('2m','10m','90m','Location','best');
xlabel('Number of Buildings')
ylabel('PDF')
title('MCF-MCF-3')


% MR_3
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MR_3(i,1) = GSR_2m.MR_3(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MR_3(i,2) = GSR_10m.MR_3(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MR_3(i,3) = GSR_90m.MR_3(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
[R.MR_3, P.MR_3, RL.MR_3, RU.MR_3] = corrcoef(GSR_combined.MR_3);

[f,xi] = ksdensity(GSR_combined.MR_3(:,1)); 
figure
plot(xi,f,'LineWidth',2);
hold on
[f,xi] = ksdensity(GSR_combined.MR_3(:,2)); 
plot(xi,f,'LineWidth',2);
[f,xi] = ksdensity(GSR_combined.MR_3(:,3)); 
plot(xi,f,'LineWidth',2);
hold off
grid on
legend('2m','10m','90m','Location','best');
xlabel('Number of Buildings')
ylabel('PDF')
title('MR-3')

%% HBET_4_5

% MCF_MR_4_5
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MCF_MR_4_5(i,1) = GSR_2m.MCF_MR_4(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MR_4_5(i,2) = GSR_10m.MCF_MR_4(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MR_4_5(i,3) = GSR_90m.MCF_MR_4(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
[R.MCF_MR_4_5, P.MCF_MR_4_5, RL.MCF_MR_4_5, RU.MCF_MR_4_5] = corrcoef(GSR_combined.MCF_MR_4_5);

[f,xi] = ksdensity(GSR_combined.MCF_MR_4_5(:,1)); 
figure
plot(xi,f,'LineWidth',2);
hold on
[f,xi] = ksdensity(GSR_combined.MCF_MR_4_5(:,2)); 
plot(xi,f,'LineWidth',2);
[f,xi] = ksdensity(GSR_combined.MCF_MR_4_5(:,3)); 
plot(xi,f,'LineWidth',2);
hold off
grid on
legend('2m','10m','90m','Location','best');
xlabel('Number of Buildings')
ylabel('PDF')
title('MCF-MR-4-5')

% MR_4_5
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MR_4_5(i,1) = GSR_2m.MR_4_5(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MR_4_5(i,2) = GSR_10m.MR_4_5(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MR_4_5(i,3) = GSR_90m.MR_4_5(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
[R.MR_4_5, P.MR_4_5, RL.MR_4_5, RU.MR_4_5] = corrcoef(GSR_combined.MR_4_5);

[f,xi] = ksdensity(GSR_combined.MR_4_5(:,1)); 
figure
plot(xi,f,'LineWidth',2);
hold on
[f,xi] = ksdensity(GSR_combined.MR_4_5(:,2)); 
plot(xi,f,'LineWidth',2);
[f,xi] = ksdensity(GSR_combined.MR_4_5(:,3)); 
plot(xi,f,'LineWidth',2);
hold off
grid on
legend('2m','10m','90m','Location','best');
xlabel('Number of Buildings')
ylabel('PDF')
title('MR-4-5')

% MCF_MCF_4_5
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MCF_MCF_4_5(i,1) = GSR_2m.MCF_MCF_4(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MCF_4_5(i,2) = GSR_10m.MCF_MCF_4(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MCF_4_5(i,3) = GSR_90m.MCF_MCF_4(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
[R.MCF_MCF_4_5, P.MCF_MCF_4_5, RL.MCF_MCF_4_5, RU.MCF_MCF_4_5] = corrcoef(GSR_combined.MCF_MCF_4_5);

[f,xi] = ksdensity(GSR_combined.MCF_MCF_4_5(:,1)); 
figure
plot(xi,f,'LineWidth',2);
hold on
[f,xi] = ksdensity(GSR_combined.MCF_MCF_4_5(:,2)); 
plot(xi,f,'LineWidth',2);
[f,xi] = ksdensity(GSR_combined.MCF_MCF_4_5(:,3)); 
plot(xi,f,'LineWidth',2);
hold off
grid on
legend('2m','10m','90m','Location','best');
xlabel('Number of Buildings')
ylabel('PDF')
title('MCF-MCF-4-5')

%% HBET_6_9

% CR_4_9
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.CR_4_9(i,1) = GSR_2m.CR_4_9(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.CR_4_9(i,2) = GSR_10m.CR_4_9(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.CR_4_9(i,3) = GSR_90m.CR_4_9(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
[R.CR_4_9, P.CR_4_9, RL.CR_4_9, RU.CR_4_9] = corrcoef(GSR_combined.CR_4_9);

[f,xi] = ksdensity(GSR_combined.CR_4_9(:,1)); 
figure
plot(xi,f,'LineWidth',2);
hold on
[f,xi] = ksdensity(GSR_combined.CR_4_9(:,2)); 
plot(xi,f,'LineWidth',2);
[f,xi] = ksdensity(GSR_combined.CR_4_9(:,3)); 
plot(xi,f,'LineWidth',2);
hold off
grid on
legend('2m','10m','90m','Location','best');
xlabel('Number of Buildings')
ylabel('PDF')
title('CR-4-9')

%% HBET_10_24

% CR_10_24
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.CR_10_24(i,1) = GSR_2m.CR_10_2(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.CR_10_24(i,2) = GSR_10m.CR_10_2(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.CR_10_24(i,3) = GSR_90m.CR_10_2(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
[R.CR_10_24, P.CR_10_24, RL.CR_10_24, RU.CR_10_24] = corrcoef(GSR_combined.CR_10_24);

[f,xi] = ksdensity(GSR_combined.CR_10_24(:,1)); 
figure
plot(xi,f,'LineWidth',2);
hold on
[f,xi] = ksdensity(GSR_combined.CR_10_24(:,2)); 
plot(xi,f,'LineWidth',2);
[f,xi] = ksdensity(GSR_combined.CR_10_24(:,3)); 
plot(xi,f,'LineWidth',2);
hold off
grid on
legend('2m','10m','90m','Location','best');
xlabel('Number of Buildings')
ylabel('PDF')
title('CR-10-24')

%% HBET_25_40

% CR_25_40
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.CR_25_40(i,1) = GSR_2m.CR_25_4(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.CR_25_40(i,2) = GSR_10m.CR_25_4(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.CR_25_40(i,3) = GSR_90m.CR_25_4(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
[R.CR_25_40, P.CR_25_40, RL.CR_25_40, RU.CR_25_40] = corrcoef(GSR_combined.CR_25_40);

[f,xi] = ksdensity(GSR_combined.CR_25_40(:,1)); 
figure
plot(xi,f,'LineWidth',2);
hold on
[f,xi] = ksdensity(GSR_combined.CR_25_40(:,2)); 
plot(xi,f,'LineWidth',2);
[f,xi] = ksdensity(GSR_combined.CR_25_40(:,3)); 
plot(xi,f,'LineWidth',2);
hold off
grid on
legend('2m','10m','90m','Location','best');
xlabel('Number of Buildings')
ylabel('PDF')
title('CR-25-40')
