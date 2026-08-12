# CRC TF一级闸门

该目录实现 `analysis/gate1-protocol.md`。代码可以在本地接受合成输入测试，但真实CRC数据的
下载和运行只允许在项目云服务器进行。

云端固定布局：

```text
tmp/tnbc-chromatin-tf-nc-2026/
├── raw/          # CRC H5AD只读链接、CollecTRI和TCGA下载
├── env/          # 独立Python环境及freeze
├── work/         # 患者pseudobulk和中间矩阵
├── outputs/      # 完整结果
├── logs/         # 命令日志
└── manifests/    # 输入、环境、运行和判定收据
```

入口 `run_gate1.sh` 依次执行审计、网络冻结、pseudobulk、TF统计、验证和判定。任一阶段失败即
退出，不会继续到锚点的ATAC或药物层。

