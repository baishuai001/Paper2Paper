#!/usr/bin/env python3
"""Patient-level TF models, random-effects meta-analysis and no-leak LODO validation."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

from gate1_common import bh_fdr, sha256_file, write_json


def fit_ols_many(y: np.ndarray, metadata: pd.DataFrame) -> pd.DataFrame:
    label = metadata["immune_label"].eq("M").astype(float).to_numpy()
    dataset = pd.get_dummies(metadata["dataset"], drop_first=True, dtype=float)
    log_cells = np.log10(metadata["cancer_cells_aggregated"].to_numpy(float) + 1.0)
    X = np.column_stack([np.ones(len(metadata)), label, dataset.to_numpy(), log_cells])
    xtx_inv = np.linalg.pinv(X.T @ X)
    beta = xtx_inv @ X.T @ y
    residual = y - X @ beta
    dof = len(metadata) - np.linalg.matrix_rank(X)
    sigma2 = (residual * residual).sum(axis=0) / dof
    se = np.sqrt(np.maximum(xtx_inv[1, 1] * sigma2, np.finfo(float).tiny))
    z = beta[1] / se
    p = 2 * norm.sf(np.abs(z))
    return pd.DataFrame({"adjusted_effect": beta[1], "SE": se, "z": z, "p": p, "FDR": bh_fdr(p)})


def hedges_g(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    n1, n0 = len(x), len(y)
    pooled = np.sqrt(((n1 - 1) * x.var(ddof=1) + (n0 - 1) * y.var(ddof=1)) / (n1 + n0 - 2))
    if pooled == 0 or not np.isfinite(pooled):
        return 0.0, np.inf
    correction = 1.0 - 3.0 / (4.0 * (n1 + n0) - 9.0)
    g = correction * (x.mean() - y.mean()) / pooled
    variance = (n1 + n0) / (n1 * n0) + g * g / (2 * (n1 + n0 - 2))
    return float(g), float(variance)


def random_effects(effects: np.ndarray, variances: np.ndarray) -> tuple[float, float, float, float, float]:
    weights = 1.0 / variances
    fixed = np.sum(weights * effects) / np.sum(weights)
    Q = np.sum(weights * (effects - fixed) ** 2)
    df = len(effects) - 1
    c = np.sum(weights) - np.sum(weights * weights) / np.sum(weights)
    tau2 = max(0.0, (Q - df) / c) if c > 0 else 0.0
    random_w = 1.0 / (variances + tau2)
    pooled = np.sum(random_w * effects) / np.sum(random_w)
    se = np.sqrt(1.0 / np.sum(random_w))
    p = float(2 * norm.sf(abs(pooled / se)))
    i2 = max(0.0, (Q - df) / Q) * 100.0 if Q > 0 else 0.0
    return float(pooled), float(se), p, float(i2), float(tau2)


def meta_many(activity: np.ndarray, tfs: np.ndarray, metadata: pd.DataFrame, primary: set[str]) -> pd.DataFrame:
    rows=[]
    for j, tf in enumerate(tfs):
        if tf not in primary: continue
        effects=[]; variances=[]; directions=[]
        for dataset, idx in metadata.groupby("dataset").groups.items():
            idx=np.asarray(list(idx),dtype=int); label=metadata.iloc[idx]["immune_label"].eq("M").to_numpy()
            if label.sum()<3 or (~label).sum()<3: continue
            g,v=hedges_g(activity[idx,j][label],activity[idx,j][~label]); effects.append(g);variances.append(v);directions.append(np.sign(g))
        if len(effects)<3: continue
        pooled,se,p,i2,tau2=random_effects(np.asarray(effects),np.asarray(variances))
        same=np.mean(np.asarray(directions)==np.sign(pooled)) if pooled!=0 else 0.0
        rows.append({"TF":tf,"meta_effect":pooled,"meta_SE":se,"meta_p":p,"I2":i2,"tau2":tau2,"informative_datasets":len(effects),"same_direction_fraction":same})
    result=pd.DataFrame(rows)
    if not result.empty: result["meta_FDR"]=bh_fdr(result["meta_p"].to_numpy())
    return result


def select_features(activity: np.ndarray, label: np.ndarray, train: np.ndarray, top_n: int=10) -> np.ndarray:
    scale=StandardScaler().fit(activity[train]); x=scale.transform(activity[train]); y=label[train]
    effects=x[y==1].mean(axis=0)-x[y==0].mean(axis=0)
    return np.argsort(np.abs(effects))[-top_n:]


def lodo_predict(activity: np.ndarray, label: np.ndarray, datasets: np.ndarray) -> np.ndarray:
    pred=np.full(len(label),np.nan)
    for held in sorted(set(datasets)):
        test=datasets==held;train=~test
        if len(np.unique(label[train]))<2: continue
        selected=select_features(activity,label,train)
        scale=StandardScaler().fit(activity[train][:,selected])
        model=LogisticRegression(class_weight="balanced",penalty="l2",C=1.0,max_iter=5000,random_state=1729)
        model.fit(scale.transform(activity[train][:,selected]),label[train])
        pred[test]=model.predict_proba(scale.transform(activity[test][:,selected]))[:,1]
    return pred


def validation(activity: np.ndarray, metadata: pd.DataFrame, permutations:int, bootstraps:int, seed:int=1729) -> tuple[dict[str,object],pd.DataFrame,pd.DataFrame]:
    label=metadata["immune_label"].eq("M").astype(int).to_numpy(); datasets=metadata["dataset"].to_numpy(str)
    observed=lodo_predict(activity,label,datasets); valid=np.isfinite(observed); auc=float(roc_auc_score(label[valid],observed[valid]))
    rng=np.random.default_rng(seed); boot=[]
    strata={(d,y):np.flatnonzero((datasets==d)&(label==y)) for d in sorted(set(datasets)) for y in (0,1)}
    for _ in range(bootstraps):
        idx=np.concatenate([rng.choice(v,size=len(v),replace=True) for v in strata.values() if len(v)])
        if len(np.unique(label[idx]))==2: boot.append(roc_auc_score(label[idx],observed[idx]))
    perm_rows=[]
    for p in range(permutations):
        perm=label.copy()
        for d in sorted(set(datasets)):
            idx=np.flatnonzero(datasets==d); perm[idx]=rng.permutation(perm[idx])
        pp=lodo_predict(activity,perm,datasets); pv=np.isfinite(pp)
        perm_rows.append({"permutation":p,"AUROC":float(roc_auc_score(perm[pv],pp[pv]))})
    perm=pd.DataFrame(perm_rows); empirical=float((1+(perm["AUROC"]>=auc).sum())/(1+len(perm)))
    pred=pd.DataFrame({"donor_id":metadata["donor_id"],"dataset":datasets,"M":label,"prediction":observed})
    receipt={"LODO_AUROC":auc,"bootstrap_95_CI":[float(np.quantile(boot,.025)),float(np.quantile(boot,.975))],"bootstrap_replicates":len(boot),"permutations":len(perm),"permutation_empirical_p":empirical}
    return receipt,pred,perm


def run(metadata_path:Path,activity_path:Path,output_dir:Path,threshold:int=50,permutations:int=500,bootstraps:int=1000)->dict[str,object]:
    metadata=pd.read_csv(metadata_path,sep="\t"); loaded=np.load(activity_path,allow_pickle=False); activity=loaded["activity"].astype(float); tfs=loaded["TF"].astype(str)
    keep=metadata["cancer_cells_aggregated"].ge(threshold).to_numpy(); metadata=metadata.loc[keep].reset_index(drop=True); activity=activity[keep]
    zactivity=(activity-activity.mean(axis=0))/np.where(activity.std(axis=0,ddof=1)>0,activity.std(axis=0,ddof=1),1)
    main=fit_ols_many(zactivity,metadata); main.insert(0,"TF",tfs); primary=set(main.loc[(main.FDR<=.05)&(main.adjusted_effect.abs()>=.40),"TF"])
    meta=meta_many(zactivity,tfs,metadata,primary); validation_receipt,pred,perm=validation(zactivity,metadata,permutations,bootstraps)
    meta_pass=meta[(meta.meta_FDR<=.10)&(meta.same_direction_fraction>=.70)&(meta.I2<=75)] if not meta.empty else meta
    conditions={"primary_TFs_ge_10":len(primary)>=10,"meta_TFs_ge_5":len(meta_pass)>=5,"LODO_AUROC_ge_0_65":validation_receipt["LODO_AUROC"]>=.65,"bootstrap_lower_gt_0_55":validation_receipt["bootstrap_95_CI"][0]>.55,"permutation_p_le_0_05":validation_receipt["permutation_empirical_p"]<=.05}
    output_dir.mkdir(parents=True,exist_ok=True); main.to_csv(output_dir/"G5_primary_TF_models.tsv",sep="\t",index=False,lineterminator="\n");meta.to_csv(output_dir/"G5_random_effects_meta.tsv",sep="\t",index=False,lineterminator="\n");pred.to_csv(output_dir/"G5_lodo_predictions.tsv",sep="\t",index=False,lineterminator="\n");perm.to_csv(output_dir/"G5_lodo_permutations.tsv",sep="\t",index=False,lineterminator="\n")
    receipt={"status":"passed" if all(conditions.values()) else "failed","generated_at":datetime.now(timezone.utc).isoformat(),"threshold":threshold,"patients":len(metadata),"M_patients":int(metadata.immune_label.eq("M").sum()),"primary_TFs":len(primary),"meta_reproduced_TFs":len(meta_pass),**validation_receipt,"conditions":conditions,"outputs":{}}
    for path in sorted(output_dir.glob("G5_*.tsv")):receipt["outputs"][path.name]={"bytes":path.stat().st_size,"sha256":sha256_file(path)}
    write_json(output_dir/"G5_tf_statistics_receipt.json",receipt);print(receipt);return receipt


def main()->int:
    p=argparse.ArgumentParser();p.add_argument("--metadata",required=True,type=Path);p.add_argument("--activity",required=True,type=Path);p.add_argument("--output-dir",required=True,type=Path);p.add_argument("--threshold",type=int,default=50);p.add_argument("--permutations",type=int,default=500);p.add_argument("--bootstraps",type=int,default=1000);a=p.parse_args();run(a.metadata,a.activity,a.output_dir,a.threshold,a.permutations,a.bootstraps);return 0
if __name__=="__main__":raise SystemExit(main())
