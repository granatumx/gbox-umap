import matplotlib.pyplot as plt

import scanpy as sc
import numpy as np
import umap

from granatum_sdk import Granatum
import time

def main():
    tic = time.perf_counter()

    gn = Granatum()

    df = gn.pandas_from_assay(gn.get_import('assay'))
    n_neighbors = gn.get_arg('n_neighbors')
    min_dist = gn.get_arg('min_dist')
    n_pcs = gn.get_arg('n_pcs')
    if n_pcs == 0:
       n_pcs = None
    metric = gn.get_arg('metric')
    random_seed = gn.get_arg('random_seed')

    #adata = sc.AnnData(df.copy(), dtype=np.float64)  # Do not want to destroy the data variable
    #sc.tl.pca(adata, svd_solver='arpack', random_state=random_seed)
    #sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, metric=metric)
    #sc.tl.umap(adata, min_dist=min_dist, random_state=random_seed)

    embedding = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, metric=metric, random_state=random_seed).fit_transform(df.values.T)
    embedding = adata.obsm["X_umap"]

    plt.figure()
    plt.scatter(embedding[:, 0], embedding[:, 1], min(5000 / df.shape[0], 36.0))
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    plt.tight_layout()

    gn.add_current_figure_to_results('UMAP plot: each dot represents a cell', dpi=75)

    umap_export = {
        'dimNames': ['UMAP 1', 'UMAP 2'],
        'coords': {sample_id: embedding[i, :].tolist() for i, sample_id in enumerate(df.columns)},
    }
    gn.export_statically(umap_export, 'UMAP coordinates')

    toc = time.perf_counter()
    time_passed = round(toc - tic, 2)

    timing = "* Finished UMAP step in {} seconds*".format(time_passed)
    gn.add_result(timing, "markdown")

    gn.commit()


if __name__ == '__main__':
    main()
