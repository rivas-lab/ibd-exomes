import pandas as pd
import numpy as np

manifest = pd.read_table('manifest/IBD_WES_MANIFEST_AUGUST2019.txt')

genomes_manifest = pd.read_table('manifest/IMMUNE_CCDG_WGS_MANIFEST_MAY2018.txt')

pop_dict = {'1KG': 'CONTROL',
            'ADSP': 'CONTROL',
            'AHMAD': 'NFE',
            'ALLEZ': 'NFE',
            'ALM': '',
            'CCFA': 'US- Other Cohorts',
            'CHO': 'AJ',
            'DEUUKB': 'CONTROL',
            'DEUUKL': 'CONTROL',
            'DEUULG': 'CONTROL',
            'DEUUTB': 'CONTROL',
            'ESP': 'CONTROL',
            'FARKKILA': 'FIN',
            'FARMER': '',
            'FINNISH': 'FIN',
            'FRANCHIMONT': 'NFE',
            'FRANKE': 'NFE',
            'IRLRCI': 'CONTROL',
            'KELSEN': 'US- Other Cohorts',
            'KUGATHASAN': 'US- Other Cohorts',
            'KUPCINSKAS': 'LIT',
            'LOUIS': 'NFE',
            'MCGOVERN': 'AJ',
            'MIGEN-LEICESTER': 'CONTROL',
            'MIGEN-OTTAWA': 'CONTROL',
            'NEC': '',
            'NIDDK-BRANT': 'US- NIDDK',
            'NIDDK-CHO': 'US- NIDDK',
            'NIDDK-DUERR': 'US- NIDDK',
            'NIDDK-MCGOVERN': 'US- NIDDK',
            'NIDDK-RIOUX': 'US- NIDDK',
            'NIDDK-SILVERBERG': 'US- NIDDK',
            'NIMH': 'CONTROL',
            'OSTRER': 'CONTROL',
            'PALOTIE': 'FIN',
            'PLAGNOL': 'NFE',
            'PULVER': 'CONTROL',
            'RIOUX-GENIZON': 'NFE',
            'SHARE-CEDARS': 'US- Other Cohorts',
            'SHARE-MAYO': 'US- Other Cohorts',
            'SHARE-MGH': 'US- Other Cohorts',
            'SHARE-MSSM': 'US- Other Cohorts',
            'SHARE-UC': 'US- Other Cohorts',
            'SHARE-UNC': 'US- Other Cohorts',
            'SHARE-WASHU': 'US- Other Cohorts',
            'SOKOL': 'NFE',
            'SONG': 'EAS',
            'T2D': 'CONTROL',
            'TURNER': 'AJ',
            'UK-IRELAND': 'CONTROL',
            'USACHP': 'CONTROL',
            'USAEGP': 'CONTROL',
            'USAHEP': 'CONTROL',
            'USAUPN': 'CONTROL',
            'WEERSMA': 'NFE',
            'WINTER': 'US- Other Cohorts',
            'XAVIER': 'US- MGH PRISM'}

manifest['RACE/ETHNICITY'] = 'NA'

for i, row in manifest.iterrows():
	manifest.at[i, 'RACE/ETHNICITY'] = pop_dict[manifest.at[i, 'COHORT']]

manifest.to_csv('manifest/IBD_WES_MANIFEST_AUGUST2019_POP.txt', sep='\t', index=False)
genomes_manifest = genomes_manifest[['SAMPID', 'RACE/ETHNICITY', 'DIAGNOSIS', 'COHORT']]
genomes_manifest['SAMPLE_ID'] = genomes_manifest['SAMPID']
genomes_manifest = genomes_manifest[['SAMPLE_ID','RACE/ETHNICITY', 'DIAGNOSIS', 'COHORT']]

pd.concat([manifest[['SAMPLE_ID', 'RACE/ETHNICITY', 'DIAGNOSIS', 'COHORT']], genomes_manifest]).to_csv('v36+ccdg_pop+diagnosis.tsv', sep='\t', index=False)


