# Fibromyalgia - proposed separate project

Status: project scaffold and gate audit only, not a finished paper or a positive benchmark. The parent repository is shared core; this directory must hold a disease-specific protocol, accession and external-service evidence ledger, reproducible results and a 50-page substantive paper before its gates can be claimed.

Current disease-tagged manifest records: 205 = 3 GSE studies + 201 nested GSM samples + 1 other. These are record units, not independent datasets or patients. The shared 40-service and 49-page PDF do not transfer as automatic per-project passes. Benchmark/discovery endpoint is open.

Disease-specific documents and source logs are not yet split from shared core; use the shared source code and result filenames by disease as leads, then verify original record attribution before copying.

## Discovery-cohort comorbidity boundary

[GSE67311](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67311), the original whole-blood discovery source, has 142 array columns with 67 GEO-labeled FM and 75 healthy controls, despite the series summary rounding its analysis as 70 and 70. Its per-GSM matrix metadata show ten FM cases also marked chronic fatigue syndrome "Yes", fifty-one marked "No", and six missing/"-"; 72 controls are CFS "No", two "-", and one blank. These are 142 source-labeled sample records, not evidence that an FM-only signature is specific against ME/CFS. The historical effects did not stratify CFS comorbidity. The source's unique titles are a lead, not a separate donor crosswalk; before clinical differential-diagnosis claims, audit individual records, lock comorbidity handling and compare FM to a real symptomatic comparator, not just healthy controls. No new endpoint or accession record is claimed here.
