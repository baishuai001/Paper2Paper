# CRC K=10 逐程序联合稳定性与技术筛查人工审查报告

## 一句话结论

有资格进入 Che 确认性外部复核的主程序为：P1, P7, P8。
这里的‘有资格’只表示 Liu 内部证据足以接受下一步外部检验，不表示已经被 Che 验证。

## 本轮做了什么

每个 Liu 主分析 K=10 程序都接受了跨 seed、细胞重抽样、患者 A/B 独立留出、两项恶性细胞集合敏感性、患者贡献和技术成分检查。
所有阈值在生成资格结论前冻结；没有运行 Che，没有做程序生物学命名。

## 技术指标可用性说明

冻结 cNMF 输入不能得到有变化的逐细胞线粒体比例，因此该相关性记为不可用，而不是记为通过。
但所有程序 top-50 中线粒体基因最多为 0 个；预注册相关性淘汰规则还要求至少 5 个线粒体 top gene，
所以这个缺失指标不会改变本轮任何程序的资格。总 counts 和检测基因数相关性均已计算。

## 逐程序结论

### P1：可进入 Che 确认性复核

- top genes：SELENOP,BEX4,CEACAM6,PLCB4,CST3,TSPAN8,SPINK1,FGGY,CEACAM5,B2M,RRBP1,SMIM26,APOC2,FZD10,IL2RG,CYSTM1,SEC62,CD24,ENSG00000236449,ST6GAL1
- 五项稳定性：5/5 通过。
- 患者贡献：通过；三套完整输入中最小 effective patient=7.21，最少 ≥5% 贡献患者数=5。
- 技术主导：否。
- 最终依据：五项稳定性、患者贡献及技术筛查均通过。

### P2：不可进入 Che 确认性复核

- top genes：CLDN4,CD55,EFNB2,CD9,PCCA,MIR22HG,EMP1,TUBB2A,ALDOA,MISP,NR4A1,COMMD6,TACSTD2,TMEM219,MIDN,ITPKC,RGCC,CAMK2N1,KLF4,PDLIM4
- 五项稳定性：4/5 通过。
- 患者贡献：通过；三套完整输入中最小 effective patient=5.89，最少 ≥5% 贡献患者数=4。
- 技术主导：否。
- 最终依据：未通过：patient_A_B_joint。

### P3：不可进入 Che 确认性复核

- top genes：PFN1,GNAS,C1QBP,TAF10,ECH1,TMEM238,HSPB6,MZT2B,PPDPF,MDH2,NME4,TRAPPC1,SNRPN,ASCL2,FCGRT,CD81,PSMB6,MZT2A,SCAND1,RNASEK
- 五项稳定性：2/5 通过。
- 患者贡献：通过；三套完整输入中最小 effective patient=5.59，最少 ≥5% 贡献患者数=3。
- 技术主导：否。
- 最终依据：未通过：cell_resample；exclude_high_unresolved；include_one_method。

### P4：不可进入 Che 确认性复核

- top genes：JUN,BTG2,TOX,EGR1,HES6,PCSK1,TUBA1A,ATF3,RHOB,ITPR2,UBC,MATN2,DUSP1,CA8,FOSB,JUND,ROBO2,KLK12,HEPACAM2,SCG3
- 五项稳定性：4/5 通过。
- 患者贡献：通过；三套完整输入中最小 effective patient=7.46，最少 ≥5% 贡献患者数=6。
- 技术主导：是（stress_dissociation_top50_ge_10）。
- 最终依据：未通过：patient_A_B_joint；technical_dominant。

### P5：不可进入 Che 确认性复核

- top genes：ANPEP,SLC25A5,PRAC1,NDUFS5,AGR2,TMSB4X,CYCS,TPM2,ENO1,MESP1,NOS3,PCK1,EIF3I,SSBP1,DDT,MAOA,ATP5MF,CYP2W1,TPM1,SERPINA6
- 五项稳定性：4/5 通过。
- 患者贡献：通过；三套完整输入中最小 effective patient=6.91，最少 ≥5% 贡献患者数=6。
- 技术主导：否。
- 最终依据：未通过：patient_A_B_joint。

### P6：不可进入 Che 确认性复核

- top genes：CYC1,PSMG1,G3BP1,LARS1,L1TD1,CSNK1A1,DNAJA2,PDGFRA,CDX1,SMIM3,SRP72,PAICS,TCOF1,RBCK1,ADISSP,EXOSC4,NDUFA6,XRN2,COX5A,FERMT1
- 五项稳定性：3/5 通过。
- 患者贡献：通过；三套完整输入中最小 effective patient=5.91，最少 ≥5% 贡献患者数=7。
- 技术主导：否。
- 最终依据：未通过：patient_A_B_joint；exclude_high_unresolved。

### P7：可进入 Che 确认性复核

- top genes：S100A11,KLK6,RBP1,ARL4C,TM4SF1,PHLDA3,HLA-A,DCBLD2,TESC,S100A10,PHLDA1,S100A6,FSCN1,ANXA2,TMSB10,SERPINB5,SERPINB1,FLNA,BMP4,MSN
- 五项稳定性：5/5 通过。
- 患者贡献：通过；三套完整输入中最小 effective patient=5.23，最少 ≥5% 贡献患者数=5。
- 技术主导：否。
- 最终依据：五项稳定性、患者贡献及技术筛查均通过。

### P8：可进入 Che 确认性复核

- top genes：TPX2,PLK1,NUSAP1,UBE2C,TOP2A,TUBA1B,CCNB1,CCNA2,BIRC5,CDK1,HMGB2,CDKN3,HMMR,CKAP2,CDCA3,PTTG1,CDC20,KIF20B,CKS1B,CENPA
- 五项稳定性：5/5 通过。
- 患者贡献：通过；三套完整输入中最小 effective patient=10.24，最少 ≥5% 贡献患者数=10。
- 技术主导：否。
- 仍需注意的非淘汰标记：cell_cycle_dominant。
- 最终依据：五项稳定性、患者贡献及技术筛查均通过。

### P9：不可进入 Che 确认性复核

- top genes：AGR3,SNX3,PCSK1N,LINC01819,ST8SIA6-AS1,FAM3B,TSPAN6,CSRP2,PTGES3,AKR7A3,TCEA1,HSPD1,HSPE1,RFPL2,ST13,CCT6A,MOCS2,CACYBP,OTUD6B-AS1,ZNF706
- 五项稳定性：4/5 通过。
- 患者贡献：通过；三套完整输入中最小 effective patient=5.06，最少 ≥5% 贡献患者数=4。
- 技术主导：否。
- 最终依据：未通过：patient_A_B_joint。

### P10：不可进入 Che 确认性复核

- top genes：KLK7,KLK10,TSPAN1,PLAAT4,SPRR3,LCN2,UCA1,APOL1,KRT7,CTSE,IGFL2,CRIP1,CALB1,PSORS1C2,TIMP2,CLIC3,SPNS2,NXN,ITGB4,C19orf33
- 五项稳定性：4/5 通过。
- 患者贡献：通过；三套完整输入中最小 effective patient=4.22，最少 ≥5% 贡献患者数=3。
- 技术主导：否。
- 最终依据：未通过：patient_A_B_joint。

## 如何进入下一步

Che 投射时建议盲法计算全部 10 个程序，以保留阴性对照并避免选择性报告；
但只有本表预先标为 eligible_for_confirmatory_che 的程序，才能承担确认性外部复核结论。
其余程序即使 Che 中偶然出现较高得分，也不能反向抹去 Liu 内部稳定性失败。

## 仍需人工审查

- 程序是否具有明确生物学含义尚未判断；本轮没有做命名或富集。
- 细胞周期程序可以是真实肿瘤状态，因此只作控制/复核标记，不自动删除。
- Che 只有 5 位配对患者，后续必须报告效应方向和不确定性，不能把阴性简单写成‘未验证’。
