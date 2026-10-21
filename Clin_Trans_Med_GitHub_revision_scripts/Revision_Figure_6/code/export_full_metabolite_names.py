from pathlib import Path
import pandas as pd,shutil,re
r=Path(__file__).resolve().parents[1];m=pd.read_csv(r/'inputs/metabolite_mapping.csv')
labels={'HMDB0000201':'Acetylcarnitine','HMDB0000564':'Dipalmitoylphosphatidylcholine (32:0)','HMDB0000651':'Decanoylcarnitine','HMDB0000791':'Octanoylcarnitine','HMDB0000870':'Histamine','HMDB0002250':'Lauroylcarnitine','HMDB0010384':'Lysophosphatidylcholine (18:0)','HMDB07869':'Phosphatidylcholine (30:0) / phosphatidylethanolamine (33:0) / N-methylphosphatidylethanolamine (32:0)','HMDB07890':'Phosphatidylcholine (36:5)','HMDB08029':'Ether-linked phosphatidylcholine (34:2 or 34:3)','HMDB08055':'Phosphatidylcholine (40:5)','HMDB10382':'Lysophosphatidylcholine (16:0)','HMDB10383':'Lysophosphatidylcholine (16:1)','HMDB10393':'Lysophosphatidylcholine (20:3)','HMDB11481':'Lysophosphatidylethanolamine (20:0) / lysophosphatidylcholine (17:0)','HMDB13331':'3,5-Tetradecadienoylcarnitine','HMDB13464':'Sphingomyelin (d18:0/16:1)','HMDB07921':'Phosphatidylcholine (36:3)','HMDB10386':'Lysophosphatidylcholine (18:2)','HMDB10395':'Lysophosphatidylcholine (20:4)','HMDB12100':'Sphingomyelin (d18:1/18:1)'}
prefixes={'PC':'Phosphatidylcholine','PE':'Phosphatidylethanolamine','PE-NMe':'N-methylphosphatidylethanolamine','LysoPC':'Lysophosphatidylcholine','LysoPE':'Lysophosphatidylethanolamine','SM':'Sphingomyelin'}
def expand(s):
 for old,new in prefixes.items():s=re.sub(r'(?<![A-Za-z-])'+re.escape(old)+r'(?=\()',new+' ',s)
 return s.rstrip(';')
m['figure_label']=[labels.get(h,expand(str(n)))for h,n in zip(m.target_hmdb,m.Name)];m['expanded_source_annotation']=m.Name.map(expand);m.to_csv(r/'results/full_metabolite_labels.csv',index=False)
(r/'docs/metabolite_names.md').write_text('# Full metabolite names\n\n'+ '\n\n'.join('**'+v.figure_label+'**\n\nComplete source alternatives: '+v.expanded_source_annotation for _,v in m.iterrows())+'\n\nLipid carbon:double-bond compositions are retained. Multiple source candidates are presented as alternatives; they are not asserted to be independently identified compounds.\n')
print('Exported full metabolite names:',len(m))
