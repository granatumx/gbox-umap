import matplotlib.pyplot as plt

import umap
from granatum_sdk import Granatum
import time

def main():
    tic = time.perf_counter()

    gn = Granatum()

    df = gn.pandas_from_assay(gn.get_import('assay'))
    n_neighbors = gn.get_arg('n_neighbors')
    min_dist = gn.get_arg('min_dist')
    metric = gn.get_arg('metric')

    embedding = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, metric=metric).fit_transform(df.values.T)

    plt.figure()
    plt.scatter(embedding[:, 0], embedding[:, 1], min(5000 / df.shape[0], 36.0))
    plt.xlabel('UMAP dim. 1')
    plt.ylabel('UMAP dim. 2')
    plt.tight_layout()

    gn.add_current_figure_to_results('UMAP plot: each dot represents a cell', dpi=75)

    pca_export = {
        'dimNames': ['UMAP dim. 1', 'UMAP dim. 2'],
        'coords': {sample_id: embedding[i, :].tolist() for i, sample_id in enumerate(df.columns)},
    }
    gn.export_statically(pca_export, 'UMAP coordinates')

    toc = time.perf_counter()
    time_passed = round(toc - tic, 2)

    timing = "* Finished UMAP step in {} seconds*".format(time_passed)
    gn.add_result(timing, "markdown")

    gn.commit()


if __name__ == '__main__':
    main()
