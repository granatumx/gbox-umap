import matplotlib.pyplot as plt

import numpy as np
import umap
import ast
import pandas as pd

from granatum_sdk import Granatum
import time

def main():
    tic = time.perf_counter()

    gn = Granatum()

    df = gn.pandas_from_assay(gn.get_import('assay')).T
    n_neighbors = gn.get_arg('n_neighbors')
    min_dist = gn.get_arg('min_dist')
    dens_lambda = gn.get_arg('dens_lambda')
    metric = gn.get_arg('metric')
    eval_along_line = gn.get_arg('eval_along_line')
    random_seed = gn.get_arg('random_seed')
    whether_parametric = gn.get_arg('whether_parametric')
    n_epochs = gn.get_arg('n_epochs')

    mapper = None

    if whether_parametric:
        from umap.parametric_umap import ParametricUMAP
        mapper = ParametricUMAP(n_epochs=n_epochs, random_state=random_seed).fit(df.values)
    else:
        if dens_lambda > 0.0:
            mapper = umap.UMAP(densmap=True, dens_lambda=dens_lambda, n_neighbors=n_neighbors, min_dist=min_dist, metric=metric, random_state=random_seed).fit(df.values)
        else:
            mapper = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, metric=metric, random_state=random_seed).fit(df.values)
    embedding = mapper.transform(df.values)

    plt.figure()
    plt.scatter(embedding[:, 0], embedding[:, 1], min(5000 / df.shape[0], 36.0))
        
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    plt.tight_layout()


    if eval_along_line != "":
        myline = ast.literal_eval(eval_along_line)
        start_pt = np.array([myline[0][0], myline[0][1]])
        end_pt = np.array([myline[1][0], myline[1][1]])
        plt.plot([myline[0][0], myline[1][0]], [myline[0][1], myline[1][1]], 'k-', color = 'r')
        test_pts = [start_pt*x + end_pt*(1.0-x) for x in np.linspace(0, 1, 100)]
        inv_xform = mapper.inverse_transform(test_pts)
        inverse_mapped_points = pd.DataFrame(inv_xform, columns=df.columns)
        print(inverse_mapped_points.head(), flush=True)
        inverse_mapped_points["xlocs"] = np.linspace(0, 1, 100)
        variable_genes = inverse_mapped_points.std().sort_values(ascending=False).iloc[:10].index

    gn.add_current_figure_to_results('UMAP plot: each dot represents a cell', dpi=75)

    if eval_along_line != "":
        plt.figure()
        inverse_mapped_points.plot(y=variable_genes)
        plt.legend()
        gn.add_current_figure_to_results('UMAP top gene variation along cut line (red)', dpi=75)
        

    umap_export = {
        'dimNames': ['UMAP 1', 'UMAP 2'],
        'coords': {sample_id: embedding[i, :].tolist() for i, sample_id in enumerate(df.columns)},
    }
    gn.export(umap_export, "{}".format(gn.get_arg("coord_name")), kind='sampleCoords', meta=None)
    #gn.export_statically(umap_export, 'UMAP coordinates')

    toc = time.perf_counter()
    time_passed = round(toc - tic, 2)

    timing = "* Finished UMAP step in {} seconds*".format(time_passed)
    gn.add_result(timing, "markdown")

    gn.commit()


if __name__ == '__main__':
    main()
