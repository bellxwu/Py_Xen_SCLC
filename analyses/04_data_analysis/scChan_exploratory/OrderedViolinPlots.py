"""
Created on 2026.01.16

Description: Created an ordered violin plot of different hits.

@author: bellwu
"""
# %%
def OrderedViolinPlots(adata, GeneList, ColumnOfInterest, AddEmpty=True):
    def OrderedViolinOfHit(adata, GeneOfInterest, ColumnOfInterest, AddEmpty):
        '''
        Generate a descending ordered ScanPy violin plot of gene hit from anndata
        '''
        import scanpy as sc
        import warnings
        # ensure list
        if not isinstance(GeneOfInterest, list):
            GeneOfInterest = [GeneOfInterest]
            warnings.warn(
                "Input Gene of Interest as a list",
                FutureWarning,
                stacklevel=2)
        # ensure string
        if isinstance(ColumnOfInterest, list):
            ColumnOfInterest = ColumnOfInterest[0]
            warnings.warn(
                "Input Column of Interest as a string",
                FutureWarning,
                stacklevel=2)
            
        # subset out gene row from anndata
        mask = adata.var_names.isin(GeneOfInterest)
        adata_GoI = adata[:, mask].copy()
        # collect all nonzero values and subset
        mask = adata_GoI.X > 0 
        adata_GoI = adata_GoI[mask].copy()
        # create a counts table and sort
        ValueCountsOutput = adata_GoI.obs[ColumnOfInterest].value_counts()
        SortedOutput = ValueCountsOutput.sort_values(ascending=False)
        OrderedList = list(SortedOutput.index)
        # take all zero categorical from original anndata and append to list
        if AddEmpty == True:
            TotalList = adata.obs[ColumnOfInterest].unique()
            OrderedList = OrderedList + list(TotalList[~TotalList.isin(OrderedList)])
        # plot
        OrderedViolinPlot = sc.pl.violin(
            adata=adata,
            keys=GeneOfInterest,
            groupby=ColumnOfInterest,
            jitter=0.4,
            rotation=45,
            order=OrderedList,
            use_raw=False,
            show=False)
        return(OrderedViolinPlot)
    # run for loop to output plot of gene hits as dictionary key-value
    plots = {}
    for Gene in GeneList:
        print(f"Creating violin plot of {Gene}")
        ViolinPlot = OrderedViolinOfHit(
            adata=adata,
            GeneOfInterest=Gene,
            ColumnOfInterest=ColumnOfInterest,
            AddEmpty=AddEmpty)
        plots[Gene] = ViolinPlot
    return(plots)
# %%
