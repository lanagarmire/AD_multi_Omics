"""Extract baseline education without imputation."""
from pathlib import Path
import pandas as pd
HERE=Path(__file__).resolve().parent
(HERE/'results').mkdir(exist_ok=True)
clinical=pd.read_csv(HERE/'inputs/cn_emci_lmci_ad_clin_pseudotime.csv')
ids=clinical.loc[clinical.SNF.isin(['emcisubtype1','emcisubtype2','lmcisubtype1','lmcisubtype2']),'PID']
v=pd.read_csv(HERE/'inputs/ADNIMERGE.csv',low_memory=False)
a=v[v.PTID.isin(ids)&v.VISCODE.eq('bl')].copy()
assert not a.PTID.duplicated().any() and set(a.PTID)==set(ids)
a[['PTID','RID','EXAMDATE','COLPROT','ORIGPROT','PTEDUCAT']].to_csv(HERE/'results/additional_covariates.csv',index=False)
