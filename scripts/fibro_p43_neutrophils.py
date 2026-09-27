"""P43 prespecified eight whole-blood DOWN genes in untreated neutrophils."""
import csv,hashlib,json,sys
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import ttest_ind
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from ubiomark import geo,stats
root=Path(__file__).resolve().parents[1];folder=root/'data/geo/p43';p=folder/'GSE229750_FM_HC.xlsx';expected='84e9475d9794bc89d544a441eba159f6a372e7d594300bc1eb7d19452f179214';assert hashlib.sha256(p.read_bytes()).hexdigest()==expected
rows=list(csv.DictReader((root/'results/fibro_p43_samples.csv').open()));assert len(rows)==12 and len({z['gsm'] for z in rows})==12
for z in rows:assert hashlib.sha256((folder/'gsm'/(z['gsm']+'.txt')).read_bytes()).hexdigest()==z['sha256']
baseline=[z for z in rows if z['timepoint']=='week 0' and z['treatment']=='none'];assert len(baseline)==10 and len({z['donor'] for z in baseline})==10 and {g:sum(z['group']==g for z in baseline) for g in ('fibromyalgia','healthy')}=={'fibromyalgia':5,'healthy':5}
x=pd.read_excel(p,dtype={'Geneid':str,'Symbol':str});assert x.Geneid.str.fullmatch(r'ENSG\d+').all() and x.Geneid.is_unique;keys={z['gsm']:'Counts_'+z['token']+('(AB)' if z['group']=='fibromyalgia' else '(C)') for z in baseline};assert set(keys.values())=={z for z in x.columns if z.startswith('Counts_')}
y=x.set_index('Geneid')[list(keys.values())].copy();assert np.isfinite(y.to_numpy()).all() and (y.to_numpy()>=0).all() and np.equal(y.to_numpy(),np.floor(y.to_numpy())).all()
h=pd.read_csv(geo.HGNC_PATH,sep='\t',dtype=str,usecols=['symbol','ensembl_gene_id']).dropna();amb=set(h.loc[h.ensembl_gene_id.duplicated(keep=False),'ensembl_gene_id']);h=h[~h.ensembl_gene_id.isin(amb)];mapping=dict(zip(h.ensembl_gene_id,h.symbol));y['gene']=y.index.map(mapping);expr=y.dropna(subset=['gene']).groupby('gene').sum();assert expr.index.is_unique;lib=expr.sum(axis=0);assert (lib>0).all();cpm=expr.div(lib,axis=1)*1e6;log=np.log2(cpm.loc[(cpm>1).mean(axis=1)>=.2]+1)
case=[keys[z['gsm']] for z in baseline if z['group']=='fibromyalgia'];ctrl=[keys[z['gsm']] for z in baseline if z['group']=='healthy'];g,v=stats.hedges_g(log[case].to_numpy(),log[ctrl].to_numpy());E=pd.DataFrame({'g':g,'v':v},index=log.index).replace([np.inf,-np.inf],np.nan).dropna();panel=['CPA3','GCSAML','MS4A2','FCER1A','ITGB8','UQCC3','HDC','GATA2'];rec=[]
for gene in panel:
 if gene in E.index:rec.append(dict(gene=gene,status='measured',g=float(E.loc[gene,'g']),v=float(E.loc[gene,'v']),welch_p=float(ttest_ind(log.loc[gene,case],log.loc[gene,ctrl],equal_var=False).pvalue),down=bool(E.loc[gene,'g']<0)))
 else:rec.append(dict(gene=gene,status='not uniquely mapped or expression filtered'))
count=sum(z.get('down',False) for z in rec);corrected=sum(z.get('down',False) and z.get('welch_p',1)<.05/8 for z in rec)
out=dict(source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229750',source_sha256=expected,n_case=5,n_control=5,n_repeat_treated_excluded=2,n_raw_ensg_rows=len(x),n_unique_mapped_genes=len(expr),n_effect_genes=len(E),n_fixed_measured=sum(z['status']=='measured' for z in rec),n_down=count,n_corrected_down=corrected,registered_descriptive_support=bool(sum(z['status']=='measured' for z in rec)>=7 and count>=7 and corrected>=1),ambiguous_ensembl_ids=len(amb),genes=rec,limitation='Five donors per group, neutrophils versus old whole blood, historical selected genes; not an independent biomarker or benchmark.')
(root/'results/fibro_p43_result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
