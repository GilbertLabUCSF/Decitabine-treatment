import pandas as pd
from CanDI import candi as can
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from adjustText import adjust_text 
# https://stackoverflow.com/questions/34693991/repel-annotations-in-matplotlib


def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))


def pretty_print_attr(obj):
    attr = []
    ls_attr = []
    meth = []
    for i in dir(obj):
        if "_" != i[0]:
            if type(getattr(obj, i)) == str or type(getattr(obj, i)) == int:
                attr.append(i)
            elif type(getattr(obj, i)) == list:
                ls_attr.append(i)
            else:
                meth.append(i)

    print("Attributes:\n")
    for i in attr: print(i+":", getattr(obj, i))
    for i in ls_attr: print(i+" list first item:", getattr(obj, i)[0])
    for i in ls_attr: print(i+" length:", len(getattr(obj, i)))
    print("\nMethods:\n")
    for i in meth: print(i)


def mt_wt_objs(cancer, genes):
    # Cell lines with mutation
    mt_list = cancer.mutated(genes)
    mt = can.CellLineCluster(mt_list) 

    # CellLineCluster ojbect must be instantiated with a mutable sequence
    # I use set operations to get wild type cell line ids and convert to a list
    wt_list = list(set(cancer.depmap_ids) - set(mt_list)) 

    wt = can.CellLineCluster(wt_list) 
    print (f'#of mutated cell lines:\n\t{len(mt.depmap_ids)}')
    print (f'#of wildtype cell lines:\n\t{len(wt.depmap_ids)}')
    
    return mt, wt


def gene_effect_heatmap(obj1, obj2, genes, name = None):
    #Make Figure appropriate size, dpi, and font
    plt.rcParams.update({"figure.figsize": (16, 6),
                        "savefig.dpi": 300,
                        "font.size": 12
                        })

    #One figure with one subplot
    fig, ax = plt.subplots(1,1)

    #Construcat matrix to make heatmap and cell line labels
    data = pd.concat([obj1.effect_of(genes), obj2.effect_of(genes)], axis=1)
    names = can.data.cell_lines.loc[data.columns, "cell_line_name"]

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(names)))
    ax.set_yticks(np.arange(len(genes)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(names)
    ax.set_yticklabels(genes)

    #make heatmap
    im = ax.imshow(data, cmap="RdBu")

    #Make colorbar scale to axis
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = ax.figure.colorbar(im, ax = ax, cax = cax)
    cbar.ax.set_ylabel("Gene Effect", rotation=-90, va="bottom")

    #Draw Dividing line btween mutant and
    ax.axvline(x=obj1.gene_effect.shape[1] - 0.5, c = "black", linewidth = 3)
    plt.setp(ax.get_xticklabels(), rotation=-90, ha="left", va="center",
         rotation_mode="anchor")
    plt.tight_layout()
    plt.show()

    if name:
        fig.savefig(name, dpi=300)

        
def gene_effect_scatter(mt, wt, genes=None, name_scatter='', return_effect=False, plot=True, name=None):

    # Average Gene Effect 
    mt_effect = mt.gene_dependency.mean(1)
    wt_effect = wt.gene_dependency.mean(1)

    if return_effect:
        return mt_effect, wt_effect
    
    if plot:
        #Make Figure appropriate size, dpi, and font
        plt.rcParams.update({"figure.figsize": (8, 8),
                            "savefig.dpi": 300,
                            "font.family": "sans-serif",
                            "font.size": 12
                            })

        #Generate Figure and Axis objects
        fig, ax = plt.subplots(1,1)

        #Label Axes
        ax.set_xlabel(f"{name_scatter} MT Average Gene Effect (CERES Score)")
        ax.set_ylabel(f"{name_scatter} WT Average Gene Effect (CERES Score)")

        #Draw Line at median common essential value
        ax.axhline(y = 0.50,
                   c = "black",
                   linewidth=0.5,
                   label = "Minimun Gene Dependencey Probability"
                  )

        ax.axvline(x = 0.50,
                   c= "black",
                   linewidth=0.5)

        #Plot all genes
        ax.scatter(mt_effect,
                   wt_effect,
                   c = "#2166ac",
                   alpha = 0.1,
                   s = 50
                  )

        ax.legend()

        if genes:
            #For Labeling
            if type(genes) != list:
                mt_lab = mt_effect.loc[[genes]]
                wt_lab = wt_effect.loc[[genes]]
            else:
                mt_lab = mt_effect.loc[genes]
                wt_lab = wt_effect.loc[genes]
            # Outline Genes To label
            ax.scatter(mt_lab,
                       wt_lab,
                       c = "#2166ac",
                       s = 50,
                       edgecolor = "black",
                       linewidth = 2,
                       alpha = 0.7
                      )
        #     for i in range(mt_lab.shape[0]):
        #         text = list(mt_lab.index)
        #         ax.annotate(text[i],
        #                     xy = (mt_lab[i], wt_lab[i]),
        #                     xytext = label[i],
        #                     xycoords = "data",
        #                     arrowprops = {"arrowstyle": "-"}
        #                    )
            texts = []
            for x, y, s in zip(mt_lab, wt_lab, mt_lab.index.tolist()):
                texts.append(plt.text(x, y, s))

            adjust_text(
                texts, force_points=0.2, force_text=0.2,
                expand_points=(2, 2), 
                expand_text=(2, 2),
                arrowprops=dict(arrowstyle="-", color='black', lw=0.2)
            )

        plt.show()

        if name:
            fig.savefig(name, dpi=300)
