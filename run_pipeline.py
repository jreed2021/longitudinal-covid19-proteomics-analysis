"""Reproducible secondary analysis of published longitudinal COVID-19 Olink NPX data.
All findings are exploratory. No additional NPX normalization, imputation, or automatic outlier removal.
"""
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
OUT = ROOT / 'outputs'
FIG = ROOT / 'figures'
OUT.mkdir(exist_ok=True); FIG.mkdir(exist_ok=True)

def save(df, name):
    df.to_csv(OUT/name, index=False)

p = pd.read_csv(DATA/'plasma_npx_level.csv')
m = pd.read_csv(DATA/'plasma_sample_level.csv')
key = ['SampleID','Panel','Assay']
for name, df, cols in [
    ('proteomics',p,key+['Individual_ID','NPX']),
    ('metadata',m,['SampleID','Individual_ID','Plate_ID','Case_Control','WHO_Severity_Peak','Time_From_First_Swab'])
]:
    absent = set(cols) - set(df.columns)
    if absent: raise ValueError(f'{name}: missing {absent}')
if p.duplicated(key).any(): raise ValueError('Nonunique sample-panel-assay keys; reconcile before continuing')
if m.SampleID.duplicated().any(): raise ValueError('Nonunique metadata SampleID')
if p[key].isna().any().any() or m.SampleID.isna().any(): raise ValueError('Missing essential IDs')
p['NPX'] = pd.to_numeric(p.NPX, errors='raise')
if not np.isfinite(p.NPX).all(): raise ValueError('Nonfinite NPX values')
link = p[['SampleID','Individual_ID']].drop_duplicates().merge(m[['SampleID','Individual_ID']],on='SampleID',how='outer',indicator=True,validate='one_to_one',suffixes=('_protein','_metadata'))
if not link._merge.eq('both').all() or not link.Individual_ID_protein.eq(link.Individual_ID_metadata).all(): raise ValueError('Sample to patient metadata reconciliation failed')

# Assay completeness reference comes from assays observed in this particular file,
# not from an independent instrument manifest.
ref = p[['Panel','Assay']].drop_duplicates()
expected = m[['SampleID']].merge(ref,how='cross')
observed = p[key].drop_duplicates()
absence = expected.merge(observed,on=key,how='left',indicator=True,validate='one_to_one')
absence = absence.loc[absence._merge.eq('left_only'),key].merge(m[['SampleID','Plate_ID','Case_Control']],on='SampleID',validate='many_to_one')
save(absence,'structurally_absent_assays.csv')
save(p.isna().sum().rename_axis('column').reset_index(name='missing'), 'proteomics_column_missingness.csv')
save(m.isna().sum().rename_axis('column').reset_index(name='missing'), 'metadata_column_missingness.csv')

# Sample-panel QC (do not compare NPX across individual assays as absolute concentration).
joined = p.merge(m[['SampleID','Plate_ID','Case_Control','Individual_ID','WHO_Severity_Peak','Time_From_First_Swab']],on='SampleID',how='left',validate='many_to_one',suffixes=('','_meta'))
sample_panel = joined.groupby(['SampleID','Panel','Plate_ID','Case_Control'],as_index=False).agg(median_npx=('NPX','median'),n_assays=('Assay','nunique'))
save(sample_panel,'sample_panel_qc.csv')

# Long-format assay-level QC, no imputation.
assay_qc = joined.groupby(['Panel','Assay'],as_index=False).agg(n_samples=('SampleID','nunique'),median_npx=('NPX','median'),sd_npx=('NPX','std'))
assay_qc['coverage_pct'] = assay_qc.n_samples / m.SampleID.nunique()*100
save(assay_qc,'assay_qc.csv')

# PCA restricted to assay columns measured for all samples: complete-case assays.
# This avoids constructing imputed values and is a QC visualization only.
p['feature'] = p['Panel'].astype(str)+'::'+p['Assay'].astype(str)
mat = p.pivot(index='SampleID',columns='feature',values='NPX')
complete = mat.loc[:,mat.notna().all(axis=0)]
if complete.shape[1] < 3: raise ValueError('Not enough complete assays for PCA')
z = StandardScaler().fit_transform(complete)
pca = PCA(n_components=2)
pcs = pca.fit_transform(z)
pca_df = pd.DataFrame({'SampleID':complete.index,'PC1':pcs[:,0],'PC2':pcs[:,1]}).merge(m[['SampleID','Plate_ID','Case_Control']],on='SampleID',validate='one_to_one')
save(pca_df,'pca_coordinates.csv')
fig,ax=plt.subplots(figsize=(9,6)); sns.scatterplot(data=pca_df,x='PC1',y='PC2',hue='Plate_ID',style='Case_Control',alpha=.8,ax=ax)
ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)'); ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)')
ax.set_title('PCA of complete-case Olink NPX assays (descriptive QC)'); fig.tight_layout(); fig.savefig(FIG/'pca_plate_covid.png',dpi=170); plt.close(fig)

# Longitudinal association: separate exploratory OLS for every assay, restricted
# to COVID-POSITIVE observations. Inference uses PATIENT-CLUSTERED ROBUST SEs
# to address repeated measures; NOT equivalent to a mixed-effects model.
# A fixed plate term is descriptive, NOT proof of removing technical batch effects.
positive = joined.loc[(joined.Case_Control=='POSITIVE') & joined.Time_From_First_Swab.notna()].copy()
positive['time_days']=pd.to_numeric(positive.Time_From_First_Swab,errors='coerce')
rows=[]
for (panel,assay),g in positive.groupby(['Panel','Assay'],sort=True):
    g=g.dropna(subset=['time_days','NPX','Plate_ID','Individual_ID'])
    n=g.Individual_ID.nunique()
    if n < 20 or len(g) < 40 or g.time_days.nunique()<5:
        rows.append({'Panel':panel,'Assay':assay,'status':'insufficient_data','n_observations':len(g),'n_patients':n})
        continue
    try:
        fit=smf.ols('NPX ~ time_days + C(Plate_ID)',data=g).fit(cov_type='cluster',cov_kwds={'groups':g.Individual_ID,'use_correction':True})
        b=float(fit.params['time_days']); ci=fit.conf_int().loc['time_days']
        rows.append({'Panel':panel,'Assay':assay,'status':'ok','n_observations':len(g),'n_patients':n,'beta_npx_per_day':b,'ci_low':float(ci.iloc[0]),'ci_high':float(ci.iloc[1]),'p_value':float(fit.pvalues['time_days'])})
    except Exception as exc:
        rows.append({'Panel':panel,'Assay':assay,'status':'model_error','n_observations':len(g),'n_patients':n,'error':str(exc)[:180]})
result=pd.DataFrame(rows)
result['fdr_bh']=np.nan
valid=result.status.eq('ok') & result.p_value.notna()
if valid.any(): result.loc[valid,'fdr_bh']=multipletests(result.loc[valid,'p_value'],method='fdr_bh')[1]
result=result.sort_values(['status','fdr_bh','p_value'],na_position='last')
save(result,'all_assay_longitudinal_cluster_robust.csv')

# Sample-level severity metadata is patient-constant? Examine before modeling severity.
severity_var = m.loc[m.Case_Control.eq('POSITIVE')].groupby('Individual_ID')['WHO_Severity_Peak'].nunique()
severity_consistent = bool(severity_var.le(1).all())

# Compact QA summary as actual computed findings, no invented biological interpretations.
summary={
 'proteomics_rows':len(p),'metadata_samples':len(m),'patients':m.Individual_ID.nunique(),
 'assay_definitions':len(ref),'sample_assay_absences':len(absence),
 'fully_absent_sample_panels':int((absence.groupby(['SampleID','Panel']).size().reset_index(name='absent').merge(ref.groupby('Panel').size().reset_index(name='expected'),on='Panel').eval('absent == expected')).sum()),
 'pca_complete_assays':complete.shape[1],'pca_pc1_variance_pct':round(pca.explained_variance_ratio_[0]*100,2),
 'n_modeled_assays':int(valid.sum()),'n_fdr_lt_005':int((result.fdr_bh<.05).sum()),
 'positive_patients':positive.Individual_ID.nunique(),'positive_measurements':int(((p.Assay=='IL6') & (p.Panel=='In')).sum()),
 'peak_severity_constant_within_patient':severity_consistent
}
pd.Series(summary).rename_axis('metric').reset_index(name='value').to_csv(OUT/'run_summary.csv',index=False)
print(pd.Series(summary).to_string())
print('Output files:',len(list(OUT.glob('*.csv'))),'Figure files:',len(list(FIG.glob('*.png'))))
