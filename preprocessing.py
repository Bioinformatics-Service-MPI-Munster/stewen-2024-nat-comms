import os
import scanpy as sc
import matplotlib.pyplot as plt
from scrna_tools import helpers, preprocessing
import logging


logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


output_folder = helpers.mkdir('/Users/kkruse/project-results/stewen/efn-eph-huvecs-scran-strict')
adata_dir = helpers.mkdir(os.path.join(output_folder, 'adata'))
preprocessing_dir = helpers.mkdir(output_folder, 'preprocessing')


counts_folders = {
    'stewen_huvec_si_ctrl': '/home/kkruse/project-results/stewen/huvec/'
                            'star/stewen_huvec_si_ctrl/'
                            'Solo.out/Gene/raw/',
    'stewen_huvec_si_ephb4': '/home/kkruse/project-results/stewen/huvec/'
                            'star/stewen_huvec_si_ephb4/'
                            'Solo.out/Gene/raw/',
    'stewen_huvec_si_efnb2': '/home/kkruse/project-results/stewen/huvec/'
                            'star/stewen_huvec_si_efnb2/'
                            'Solo.out/Gene/raw/',
}

for sample_name, folder in counts_folders.items():
    print(folder)
    assert os.path.exists(folder)
    assert os.path.exists(os.path.join(folder, 'barcodes.tsv'))
    assert os.path.exists(os.path.join(folder, 'features.tsv'))
    assert os.path.exists(os.path.join(folder, 'matrix.mtx'))


adatas = dict()
for sample_name, folder in counts_folders.items():
    print(sample_name)
    
    adatas[sample_name] = preprocessing.read_star_data(
        os.path.join(folder, 'matrix.mtx'),
        os.path.join(folder, 'barcodes.tsv'),
        os.path.join(folder, 'features.tsv'),
        barcodes_suffix=sample_name
    )





#
# Remove low complexity barcodes
#
for name, adata in adatas.items():
    print(name)
    fig, ax = plt.subplots(dpi=150)
    preprocessing.combined_knee_and_cumulative_count_depth_plot(adata, ax=ax)
    fig.savefig(os.path.join(preprocessing_dir, f'knee_{name}.pdf'))
    plt.show()


pairs = [
    ('stewen_huvec_si_ctrl', 1100),
    ('stewen_huvec_si_ephb4', 1300),
    ('stewen_huvec_si_efnb2', 1100),
]

processed_adatas = dict()
for name, cutoff in pairs:
    sc.pp.filter_cells(adatas[name], min_counts=cutoff)

#
# All other annotations
#
for name, adata in adatas.items():
    sub_folder = helpers.mkdir(preprocessing_dir, name)
    preprocessing.annotate_base_stats(adata, genome='hg', output_folder=sub_folder)

#
# Filter mitochondrial excess and doublets
#
for name, adata in adatas.items():
    adata = adata[adata.obs['pct_counts_mt'] < 20]
    adata = adata[adata.obs['predicted_doublets'] == False, :]
    adatas[name] = adata

for name, adata in adatas.items():
    adata.write(os.path.join(output_folder, 'adata', f'{name}.adata.h5ad'))


#
# scran norm
#
for name, adata in adatas.items():
    preprocessing.scran_norm(adata)
    

# sample info
for name, adata in adatas.items():
    adata.obs['sample'] = name

for name, adata in adatas.items():
    adata.write(os.path.join(output_folder, 'adata', f'{name}.adata.h5ad'))


#
# Merge
#
adata = adatas['stewen_huvec_si_ctrl'].concatenate(
    [adatas['stewen_huvec_si_ephb4'], adatas['stewen_huvec_si_efnb2']],
    join="outer",
    index_unique=None
)
sc.pp.filter_genes(adata, min_cells=10)
adata.write(os.path.join(output_folder, 'adata', 'adata.h5ad'))
