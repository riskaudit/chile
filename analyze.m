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
R_2m_10m_90m.CR_1_3 = corrcoef(GSR_combined.CR_1_3);

% MATO_1_2
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MATO_1_2(i,1) = GSR_2m.MATO_1_(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MATO_1_2(i,2) = GSR_10m.MATO_1_(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MATO_1_2(i,3) = GSR_90m.MATO_1_(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
R_2m_10m_90m.MATO_1_2 = corrcoef(GSR_combined.MATO_1_2);

% MCF_MR_1_2
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MCF_MR_1_2(i,1) = GSR_2m.MCF_MR_1(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MR_1_2(i,2) = GSR_10m.MCF_MR_1(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MR_1_2(i,3) = GSR_90m.MCF_MR_1(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
R_2m_10m_90m.MCF_MR_1_2 = corrcoef(GSR_combined.MCF_MR_1_2);

% MCF_MCF_1_2
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MCF_MCF_1_2(i,1) = GSR_2m.MCF_MCF_1(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MCF_1_2(i,2) = GSR_10m.MCF_MCF_1(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MCF_1_2(i,3) = GSR_90m.MCF_MCF_1(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
R_2m_10m_90m.MCF_MCF_1_2 = corrcoef(GSR_combined.MCF_MCF_1_2);

% MR_1_2
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MR_1_2(i,1) = GSR_2m.MR_1_2(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MR_1_2(i,2) = GSR_10m.MR_1_2(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MR_1_2(i,3) = GSR_90m.MR_1_2(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
R_2m_10m_90m.MR_1_2 = corrcoef(GSR_combined.MR_1_2);

% MURADO_1_2
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MURADO_1_2(i,1) = GSR_2m.MURADO_(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MURADO_1_2(i,2) = GSR_10m.MURADO_(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MURADO_1_2(i,3) = GSR_90m.MURADO_(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
R_2m_10m_90m.MURADO_1_2 = corrcoef(GSR_combined.MURADO_1_2);

% MUR_1_2
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MUR_1_2(i,1) = GSR_2m.MUR_1_2(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MUR_1_2(i,2) = GSR_10m.MUR_1_2(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MUR_1_2(i,3) = GSR_90m.MUR_1_2(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
R_2m_10m_90m.MUR_1_2 = corrcoef(GSR_combined.MUR_1_2);

% W_1_3
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.W_1_3(i,1) = GSR_2m.W_1_3(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.W_1_3(i,2) = GSR_10m.W_1_3(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.W_1_3(i,3) = GSR_90m.W_1_3(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
R_2m_10m_90m.W_1_3 = corrcoef(GSR_combined.W_1_3);

% WDNO_1_2
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.WDNO_1_2(i,1) = GSR_2m.WDNO_1_(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.WDNO_1_2(i,2) = GSR_10m.WDNO_1_(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.WDNO_1_2(i,3) = GSR_90m.WDNO_1_(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
R_2m_10m_90m.WDNO_1_2 = corrcoef(GSR_combined.WDNO_1_2);

%% HBET_3

% MCF_MR_3
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MCF_MR_3(i,1) = GSR_2m.MCF_MR_3(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MR_3(i,2) = GSR_10m.MCF_MR_3(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MR_3(i,3) = GSR_90m.MCF_MR_3(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
R_2m_10m_90m.MCF_MR_3 = corrcoef(GSR_combined.MCF_MR_3);

% MCF_MCF_3
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MCF_MCF_3(i,1) = GSR_2m.MCF_MCF_3(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MCF_3(i,2) = GSR_10m.MCF_MCF_3(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MCF_3(i,3) = GSR_90m.MCF_MCF_3(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
R_2m_10m_90m.MCF_MCF_3 = corrcoef(GSR_combined.MCF_MCF_3);

% MR_3
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MR_3(i,1) = GSR_2m.MR_3(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MR_3(i,2) = GSR_10m.MR_3(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MR_3(i,3) = GSR_90m.MR_3(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
R_2m_10m_90m.MR_3 = corrcoef(GSR_combined.MR_3);

%% HBET_4_5

% MCF_MR_4_5
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MCF_MR_4_5(i,1) = GSR_2m.MCF_MR_4(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MR_4_5(i,2) = GSR_10m.MCF_MR_4(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MR_4_5(i,3) = GSR_90m.MCF_MR_4(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
R_2m_10m_90m.MCF_MR_4_5 = corrcoef(GSR_combined.MCF_MR_4_5);

% MR_4_5
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MR_4_5(i,1) = GSR_2m.MR_4_5(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MR_4_5(i,2) = GSR_10m.MR_4_5(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MR_4_5(i,3) = GSR_90m.MR_4_5(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
R_2m_10m_90m.MR_4_5 = corrcoef(GSR_combined.MR_4_5);

% MCF_MCF_4_5
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.MCF_MCF_4_5(i,1) = GSR_2m.MCF_MCF_4(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MCF_4_5(i,2) = GSR_10m.MCF_MCF_4(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.MCF_MCF_4_5(i,3) = GSR_90m.MCF_MCF_4(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
R_2m_10m_90m.MCF_MCF_4_5 = corrcoef(GSR_combined.MCF_MCF_4_5);

%% HBET_6_9

% CR_4_9
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.CR_4_9(i,1) = GSR_2m.CR_4_9(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.CR_4_9(i,2) = GSR_10m.CR_4_9(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.CR_4_9(i,3) = GSR_90m.CR_4_9(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
R_2m_10m_90m.CR_4_9 = corrcoef(GSR_combined.CR_4_9);

%% HBET_10_24

% CR_10_24
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.CR_10_24(i,1) = GSR_2m.CR_10_2(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.CR_10_24(i,2) = GSR_10m.CR_10_2(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.CR_10_24(i,3) = GSR_90m.CR_10_2(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
R_2m_10m_90m.CR_10_24 = corrcoef(GSR_combined.CR_10_24);

%% HBET_25_40

% CR_10_24
GSR_combined.GridNr = GSR_90m.GridNr;
for i = 1:numel(GSR_combined.GridNr) 
    GSR_combined.CR_25_40(i,1) = GSR_2m.CR_25_4(GSR_2m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.CR_25_40(i,2) = GSR_10m.CR_25_4(GSR_10m.GridNr==GSR_combined.GridNr(i,1));
    GSR_combined.CR_25_40(i,3) = GSR_90m.CR_25_4(GSR_90m.GridNr==GSR_combined.GridNr(i,1));
end
R_2m_10m_90m.CR_25_40 = corrcoef(GSR_combined.CR_25_40);

%% PLOTTING CDF
cdfplot(GSR_combined.CR_1_3(:,1))
hold on
cdfplot(GSR_combined.CR_1_3(:,2))
cdfplot(GSR_combined.CR_1_3(:,3))
legend('2m','10m','90m','Location','best')
hold off
