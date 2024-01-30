import os
import scanpy as sc
from scrna_tools import r, helpers

from scrna_tools._core import _vdata as mvdata

output_folder = helpers.mkdir('/Users/kkruse/project-results/stewen/efn-eph-huvecs-scran-strict')
adata_dir = helpers.mkdir(os.path.join(output_folder, 'adata'))
adata = sc.read(os.path.join(adata_dir, 'adata.h5ad'))

#
# DE
#
de_folder = helpers.mkdir(output_folder, 'de')
for name, celltype in [
    ('mitotic', ['Mitotic']),
    ('tip_cell_like', ['Tip cell-like']),
    ('venous_like', ['Venous-like']),
    ('tip_cell_like_and_venous_like', ['Tip cell-like', 'Venous-like']),
    ('tip_cell_like_and_venous_like_and_mitotic', ['Tip cell-like', 'Venous-like', 'Mitotic']),
]:
    for sample1, sample2 in [
        ('siEFNB2', 'siCTRL'),
        ('siEPHB4', 'siCTRL'),
    ]:
        print(f'{name}: {sample1} vs {sample2}')
        de_sub_folder = helpers.mkdir(de_folder, 'main_celltype', f'{sample1}_vs_{sample2}')
        vdata_sub = mvdata.VData(adata)
        vdata_sub.add_categorical_obs_constraint('main_celltype', celltype)
        
        de = r.de_pseudobulk(vdata_sub, 'label', sample1, sample2)
        de.to_csv(os.path.join(de_sub_folder, f'de_{sample1}_vs_{sample2}_in_{name}.txt'), sep="\t")
