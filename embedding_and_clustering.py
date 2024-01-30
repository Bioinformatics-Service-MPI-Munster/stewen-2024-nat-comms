import os
import scanpy as sc
from scrna_tools import helpers, process
import logging

from scrna_tools._core import _vdata as mvdata

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


output_folder = helpers.mkdir('/Users/kkruse/project-results/stewen/efn-eph-huvecs-scran-strict')
adata_dir = helpers.mkdir(os.path.join(output_folder, 'adata'))
adata = sc.read(os.path.join(adata_dir, 'adata.h5ad'))

vdata = mvdata.VData(adata)
vdata.rename_groups('sample', {
    'stewen_huvec_si_ctrl': 'siCTRL', 
    'stewen_huvec_si_efnb2': 'siEFNB2',
    'stewen_huvec_si_ephb4': 'siEPHB4',
}, key_added='label')


vdata = mvdata.VData(adata)
process.relayout(adata, n_pcs=50, n_top_genes=2000, batch_algorithm='harmony')
process.recluster(adata)

vdata = mvdata.VData(adata)
vdata.rename_groups(
    'leiden_0.2', 
    {
        '0': 'Venous-like',
        '1': 'Tip cell-like',
        '2': 'Mitotic',
    }, 
    'main_celltype'
)
adata.write(os.path.join(adata_dir, 'adata.h5ad'))

vdata.rename_groups(
    'leiden_0.8', 
    {
        '0': 'Tip cell-like 1',
        '1': 'Venous-like 1',
        '2': 'Venous-like 2',
        '3': 'Venous-like 3',
        '4': 'Mitotic',
        '5': 'Tip cell-like 2',
        '6': 'Mitotic',
    }, 
    'sub_celltype'
)
adata.write(os.path.join(adata_dir, 'adata.h5ad'))

vdata_sub = mvdata.VData(adata)
vdata_sub.add_categorical_obs_constraint('main_celltype', ['Tip cell-like', 'Venous-like'])
non_mitotic = vdata_sub.adata_view
process.relayout(non_mitotic, n_pcs=50, n_top_genes=2000, batch_algorithm='harmony')
process.recluster(non_mitotic)

v_non_mitotic = mvdata.VData(non_mitotic)
v_non_mitotic.rename_groups(
    'leiden_0.2', 
    {
        '0': 'Venous-like',
        '1': 'Tip cell-like',
    }, 
    'main_celltype'
)

non_mitotic.write(os.path.join(adata_dir, 'adata_non_mitotic.h5ad'))
