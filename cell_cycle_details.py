import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import scanpy as sc
import scrna_tools
from scrna_tools import helpers, plotting
import numpy as np


output_folder = helpers.mkdir('/Users/kkruse/project-results/stewen/efn-eph-huvecs-scran-strict')
adata_file = os.path.join(output_folder, 'adata', 'adata.h5ad')

adata = sc.read(adata_file)


def label_plot(ax, label, x=1., y=0.5, **kwargs):
    kwargs.setdefault('ha', 'right')
    kwargs.setdefault('va', 'center')

    ax.text(x, y, label, transform=ax.transAxes, **kwargs)
    ax.set_axis_off()


#
# Cell cycle figure
#
fig = plt.figure(figsize=(10, 10))
gs = GridSpec(3, 5, width_ratios=[4, 10, 10, 10, 8])

ax_phase_label = plt.subplot(gs[0, 0])
label_plot(ax_phase_label, label='Phase', x=0.5)
ax_clusters_label = plt.subplot(gs[1, 0])
label_plot(ax_clusters_label, label='Clusters', x=0.5)
ax_bar_label = plt.subplot(gs[2, 0])
label_plot(ax_bar_label, label='Phase\npercentages', x=0.5)

samples = [('stewen_huvec_si_ctrl', 'siCTRL'),
           ('stewen_huvec_si_ephb4', 'siEPHB4'),
           ('stewen_huvec_si_efnb2', 'siEFNB2')]

for i, (sample_id, sample_name) in enumerate(samples):
    if i == 0:
        legend = True
        lax_phase = plt.subplot(gs[0, 4])
        lax_phase.axis('off')
        lax_clusters = plt.subplot(gs[1, 4])
        lax_clusters.axis('off')
    else:
        legend = False
        lax_phase = None
        lax_clusters = None

    vdata = scrna_tools.VData(adata)
    vdata.add_categorical_obs_constraint('sample', [sample_id])

    ax_phase = plt.subplot(gs[0, i + 1])
    plotting.embedding_plot(vdata, 'umap', groupby='phase', ax=ax_phase,
                                colors='Set2', alpha=0.7, groups_rename={sample_id: sample_name},
                                legend=legend, lax=lax_phase)
    ax_phase.set_title(sample_name)

    ax_clusters = plt.subplot(gs[1, i + 1])
    plotting.embedding_plot(vdata, 'umap',
                                groupby='main_celltype', ax=ax_clusters,
                                groups=['Venous-like', 'Tip cell-like', 'Mitotic'],
                                alpha=0.7, groups_rename={sample_id: sample_name},
                                legend=legend, lax=lax_clusters)



for i, (phase, ylim) in enumerate([('G1', 100), ('G2M', 30), ('S', 10)]):
    ax_bar = plt.subplot(gs[2, i + 1])
    ixs_phase = adata.obs['phase'] == phase

    percentages = []
    colors = []
    for sample_id, sample_name in samples:
        ixs_sample = adata.obs['sample'] == sample_id
        ixs_both = np.logical_and(ixs_phase, ixs_sample)
        p = np.sum(ixs_both) / np.sum(ixs_sample) * 100
        percentages.append(p)
        print(f'{phase}, {sample_id[13:]}: {int(p*100)/100}')
        #colors.append(sample_colors[sample_id])

    x = np.arange(0, len(samples))
    ax_bar.bar(x=x, height=percentages)#, color=colors)
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels([s[1] for s in samples], rotation=90, ha='right')
    ax_bar.set_title(phase)
    ax_bar.set_ylim((0, ylim))

fig.tight_layout()
fig.savefig(os.path.join(output_folder, 'cell_cycle_details.png'))
fig.savefig(os.path.join(output_folder, 'cell_cycle_details.pdf'))
plt.close(fig)
