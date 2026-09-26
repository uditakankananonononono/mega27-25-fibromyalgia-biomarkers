"""P16: frozen GSE67311 eight-gene panel in GSE221921 blood; no clinical inference."""
import csv,gzip,hashlib,io,json,re,tarfile
from pathlib import Path
import numpy as np,pandas as pd,requests
from scipy.stats import ttest_ind
from ubiomark.stats import hedges_g
G='GSE221921'; URL='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221921/'
TAR=Path('data/geo/GSE221921/raw.tar')
assert hashlib.sha256(TAR.read_bytes()).hexdigest()=='79ff6dcbebcbbb50202b86996f5c5f95ab30b483a62044e0cca85b5402c9b795'
r=requests.get(URL+'matrix/GSE221921_series_matrix.txt.gz',timeout=30);r.raise_for_status();meta_text=gzip.decompress(r.content).decode('utf8','replace').split('!series_matrix_table_begin')[0]
d={}
for line in meta_text.splitlines():
 if line.startswith('!Sample_'):
  row=next(csv.reader(io.StringIO(line),delimiter='\t'));d.setdefault(row[0],[]).append(row[1:])
titles=d['!Sample_title'][0];gsm=d['!Sample_geo_accession'][0];chars=d['!Sample_characteristics_ch1'];assert len(titles)==len(gsm)==189
assert all(len(z)==189 for a in d.values() for z in a)
labels=[s.split(':',1)[1].strip() for s in chars[3]]
assert set(labels)=={'Fibromylagia patient (FMA)','Healthy patient (CONTROL)'} and labels.count('Fibromylagia patient (FMA)')==96 and labels.count('Healthy patient (CONTROL)')==93
bios=[next((z[i].split('/')[-1] for z in d['!Sample_relation'] if 'BioSample:' in z[i]),'') for i in range(189)]
sra=[next((z[i].split('=')[-1] for z in d['!Sample_relation'] if 'SRA:' in z[i]),'') for i in range(189)]
assert len(set(titles))==len(set(gsm))==len(set(bios))==len(set(sra))==189 and all(bios) and all(sra)
meta=pd.DataFrame({'title':titles,'gsm':gsm,'biosample':bios,'sra':sra,'label':labels,'sex':[s.split(':',1)[1].strip() for s in chars[2]]})
files={};tables=[];dup=None
with tarfile.open(TAR) as archive:
 for member in archive:
  match=re.fullmatch(r'(GSM\d+)_(Sample_\d+)\.txt\.gz',member.name)
  assert match and member.name not in files
  accession,title=match.groups();files[title]=accession
  content=gzip.decompress(archive.extractfile(member).read())
  x=pd.read_csv(io.BytesIO(content),sep=r'\s+',engine='python',keep_default_na=False)
  assert x.columns.tolist()==['gene',title+'.txt'] and len(x)==21915
  assert x.iloc[:,1].notna().all() and np.isfinite(x.iloc[:,1]).all() and (x.iloc[:,1]>=0).all()
  if dup is None:dup=x.gene[x.gene.duplicated(keep=False)].unique().tolist();first=x.gene.tolist()
  assert x.gene.tolist()==first
  tables.append(x.iloc[:,1].to_numpy())
assert len(files)==189 and all(files[t]==g for t,g in zip(titles,gsm))
A=np.stack([tables[list(files).index(t)] for t in titles],axis=1);assert A.shape==(21915,189)
counts=pd.Series(first).value_counts();keep=np.array([counts[g]==1 for g in first]);genes=np.array(first,dtype=str)[keep];X=A[keep];assert len(set(genes))==len(genes)
labels=np.array(labels);case=labels=='Fibromylagia patient (FMA)';control=~case
frozen=['CPA3','GCSAML','MS4A2','FCER1A','ITGB8','UQCC3','HDC','GATA2']
assert len(set(frozen))==8
results=[]
for g in frozen:
 if g not in genes:results.append({'gene':g,'measured':False,'duplicate_in_deposit':g in dup});continue
 row=X[np.where(genes==g)[0][0]];gval,var=hedges_g(row[case][None,:],row[control][None,:]);gval=float(gval[0]);var=float(var[0]);p=float(ttest_ind(row[case],row[control],equal_var=False).pvalue)
 results.append({'gene':g,'measured':True,'duplicate_in_deposit':False,'g':gval,'se':var**.5,'normal_95_CI':[gval-1.96*var**.5,gval+1.96*var**.5],'welch_p':p,'direction':'down' if gval<0 else 'up','mean_fm':float(row[case].mean()),'mean_control':float(row[control].mean())})
# Same-direction/coverage reference: all frozen genes are discovery-DOWN and unique-row measured.
disc=pd.read_csv('results/meta_discovery/fibromyalgia.csv.gz');D=disc.set_index('gene');eligible=np.array(sorted(set(genes)&set(D.index)-set(frozen)));eligible=np.array([g for g in eligible if D.loc[g,'mu']<0]);
idx={g:i for i,g in enumerate(genes)};M=X[[idx[g] for g in eligible]]
negative=(M[:,case].mean(axis=1)-M[:,control].mean(axis=1))<0
k=sum(r.get('measured',False) for r in results);observed=sum(r.get('direction')=='down' for r in results)
rng=np.random.default_rng(20260925);B=10000
assert len(eligible)>=k and k>=7
null=np.array([negative[rng.choice(len(eligible),size=k,replace=False)].sum() for _ in range(B)])
p_emp=float((1+np.sum(null>=observed))/(B+1))
# Post-outcome sex sensitivity, never a replacement for the registered all-sample outcome.
import statsmodels.api as sm
female=meta.sex.eq('Female').to_numpy(); sensitivity=[]
for g in frozen:
 if g not in idx:continue
 row=X[idx[g]]; a=row[case&female];b=row[control&female];v,_=hedges_g(a[None,:],b[None,:])
 design=pd.DataFrame({'const':np.ones(len(meta)),'FM':case.astype(float),'female':female.astype(float)})
 fit=sm.OLS(row,design).fit(cov_type='HC3')
 sensitivity.append({'gene':g,'female_n_fm':len(a),'female_n_control':len(b),'female_g':float(v[0]),'female_welch_p':float(ttest_ind(a,b,equal_var=False).pvalue),'sex_adjusted_beta_FM':float(fit.params['FM']),'sex_adjusted_robust_p':float(fit.pvalues['FM']),'sex_adjusted_95_CI':[float(z) for z in fit.conf_int().loc['FM']]})
meta.to_csv('results/fibromyalgia_p16_sample_map.csv',index=False)
out={'source':f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={G}','raw_tar_url':URL+'suppl/GSE221921_RAW.tar','raw_tar_sha256':hashlib.sha256(TAR.read_bytes()).hexdigest(),'deposited_scale':'nonnegative apparently logged abundance despite processing note about TPM/TMM; no second transform','n_fm':int(case.sum()),'n_control':int(control.sum()),'rows_raw':len(first),'duplicate_symbol_rows':int((~keep).sum()),'duplicate_symbols':len(dup),'unique_gene_rows':int(keep.sum()),'frozen':results,'null':{'seed':20260925,'draws':B,'universe':'unique-row measured genes with GSE67311 discovery mu<0, excluding eight targets','universe_size':len(eligible),'observed_down':observed,'measured_frozen':k,'null_mean_down':float(null.mean()),'p_ge_observed':p_emp,'limitations':'Same discovery direction and coverage but not expression/variance matched; gene-gene correlation and cell composition not controlled; post-selection exploratory diagnostics.'},'descriptive_threshold_pass':bool(k>=7 and observed>=7 and p_emp<.01),'post_outcome_sex_sensitivity':sensitivity,'independence_caveat':'Distinct GEO GSM, title, BioSample, SRA. No direct patient crosswalk with GSE67311; identifiers alone cannot guarantee participant independence. Published prior FM transcriptomic signatures; no novel clinical claim.'}
Path('results/fibromyalgia_p16_result.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'results':results,'null':out['null'],'pass':out['descriptive_threshold_pass']},indent=2))
