"""Analysis tools for clustering and dimensionality reduction."""

import logging
from typing import Optional, List
import numpy as np
import anndata
import scanpy as sc
from sklearn.cluster import KMeans

logger = logging.getLogger(__name__)


class ClusteringTools:
    """Clustering algorithms."""

    @staticmethod
    def leiden(
        adata: anndata.AnnData,
        resolution: float = 1.0,
        random_state: int = 0,
        key_added: str = "leiden",
    ) -> anndata.AnnData:
        """Apply Leiden clustering.
        
        Args:
            adata: Input data with neighborhood graph
            resolution: Resolution parameter (higher = more clusters)
            random_state: Random seed
            key_added: Key to store clusters in obs
            
        Returns:
            AnnData with clusters in obs[key_added]
        """
        logger.info(f"Running Leiden clustering (resolution={resolution})...")
        adata = adata.copy()
        
        # Ensure neighbor graph exists
        if "neighbors" not in adata.obsp:
            logger.info("Computing neighbors first...")
            sc.pp.neighbors(adata, n_neighbors=15)
        
        sc.tl.leiden(adata, resolution=resolution, random_state=random_state, key_added=key_added)
        
        n_clusters = len(adata.obs[key_added].unique())
        logger.info(f"Found {n_clusters} clusters")
        
        return adata

    @staticmethod
    def louvain(
        adata: anndata.AnnData,
        resolution: float = 1.0,
        random_state: int = 0,
        key_added: str = "louvain",
    ) -> anndata.AnnData:
        """Apply Louvain clustering.
        
        Args:
            adata: Input data with neighborhood graph
            resolution: Resolution parameter
            random_state: Random seed
            key_added: Key to store clusters in obs
            
        Returns:
            AnnData with clusters in obs[key_added]
        """
        logger.info(f"Running Louvain clustering (resolution={resolution})...")
        adata = adata.copy()
        
        # Ensure neighbor graph exists
        if "neighbors" not in adata.obsp:
            logger.info("Computing neighbors first...")
            sc.pp.neighbors(adata, n_neighbors=15)
        
        sc.tl.louvain(adata, resolution=resolution, random_state=random_state, key_added=key_added)
        
        n_clusters = len(adata.obs[key_added].unique())
        logger.info(f"Found {n_clusters} clusters")
        
        return adata

    @staticmethod
    def kmeans(
        adata: anndata.AnnData,
        n_clusters: int = 10,
        random_state: int = 0,
        key_added: str = "kmeans",
    ) -> anndata.AnnData:
        """Apply K-means clustering.
        
        Args:
            adata: Input data
            n_clusters: Number of clusters
            random_state: Random seed
            key_added: Key to store clusters in obs
            
        Returns:
            AnnData with clusters in obs[key_added]
        """
        logger.info(f"Running K-means clustering (n_clusters={n_clusters})...")
        adata = adata.copy()
        
        kmeans = KMeans(n_clusters=n_clusters, random_state=random_state, n_init=10)
        adata.obs[key_added] = kmeans.fit_predict(adata.X).astype(str)
        
        logger.info(f"Found {n_clusters} clusters")
        
        return adata


class DimensionReductionTools:
    """Dimensionality reduction methods."""

    @staticmethod
    def pca(
        adata: anndata.AnnData,
        n_components: int = 50,
        use_highly_variable: bool = True,
    ) -> anndata.AnnData:
        """Apply PCA.
        
        Args:
            adata: Input data
            n_components: Number of principal components
            use_highly_variable: Use only highly variable genes
            
        Returns:
            AnnData with PCA embeddings
        """
        logger.info(f"Running PCA (n_components={n_components})...")
        adata = adata.copy()
        
        if use_highly_variable and "highly_variable" in adata.var:
            adata = adata[:, adata.var["highly_variable"]]
        
        sc.tl.pca(adata, n_comps=n_components)
        
        logger.info(f"Explained variance: {adata.uns['pca']['variance_ratio'][:5].sum():.2%}")
        
        return adata

    @staticmethod
    def umap(
        adata: anndata.AnnData,
        n_components: int = 2,
        n_neighbors: int = 15,
        min_dist: float = 0.1,
        metric: str = "euclidean",
        random_state: int = 0,
    ) -> anndata.AnnData:
        """Apply UMAP.
        
        Args:
            adata: Input data (must have PCA or neighbors)
            n_components: Number of dimensions
            n_neighbors: Number of neighbors
            min_dist: Minimum distance between points
            metric: Distance metric
            random_state: Random seed
            
        Returns:
            AnnData with UMAP embeddings
        """
        logger.info(f"Running UMAP (n_components={n_components})...")
        adata = adata.copy()
        
        # Compute neighbors if not present
        if "neighbors" not in adata.obsp:
            logger.info("Computing neighbors first...")
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, metric=metric)
        
        sc.tl.umap(adata, n_components=n_components, min_dist=min_dist, random_state=random_state)
        
        logger.info("UMAP complete")
        
        return adata

    @staticmethod
    def tsne(
        adata: anndata.AnnData,
        n_components: int = 2,
        perplexity: float = 30.0,
        random_state: int = 0,
        use_rep: str = "X_pca",
    ) -> anndata.AnnData:
        """Apply t-SNE.
        
        Args:
            adata: Input data
            n_components: Number of dimensions
            perplexity: Perplexity parameter
            random_state: Random seed
            use_rep: Representation to use (X_pca or X)
            
        Returns:
            AnnData with t-SNE embeddings
        """
        logger.info(f"Running t-SNE (n_components={n_components})...")
        adata = adata.copy()
        
        sc.tl.tsne(
            adata,
            n_pcs=50 if use_rep == "X_pca" else None,
            perplexity=perplexity,
            random_state=random_state
        )
        
        logger.info("t-SNE complete")
        
        return adata


class AnalysisTools:
    """General analysis tools combining clustering and dimensionality reduction."""

    @staticmethod
    def standard_workflow(
        adata: anndata.AnnData,
        clustering_method: str = "leiden",
        resolution: float = 1.0,
        dim_reduction: str = "umap",
        n_neighbors: int = 15,
        random_state: int = 0,
    ) -> anndata.AnnData:
        """Run standard clustering + dimensionality reduction workflow.
        
        Args:
            adata: Preprocessed data
            clustering_method: 'leiden', 'louvain', or 'kmeans'
            resolution: Resolution for leiden/louvain
            dim_reduction: 'umap' or 'tsne'
            n_neighbors: Number of neighbors
            random_state: Random seed
            
        Returns:
            Analyzed AnnData object
        """
        logger.info("Running standard analysis workflow...")
        adata = adata.copy()
        
        # PCA
        logger.info("Step 1: PCA")
        adata = DimensionReductionTools.pca(adata, n_components=50)
        
        # Compute neighbors
        logger.info("Step 2: Computing neighbors")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors)
        
        # Clustering
        logger.info(f"Step 3: Clustering ({clustering_method})")
        if clustering_method == "leiden":
            adata = ClusteringTools.leiden(adata, resolution=resolution, random_state=random_state)
        elif clustering_method == "louvain":
            adata = ClusteringTools.louvain(adata, resolution=resolution, random_state=random_state)
        elif clustering_method == "kmeans":
            adata = ClusteringTools.kmeans(adata, random_state=random_state)
        
        # Dimensionality reduction for visualization
        logger.info(f"Step 4: Dimensionality reduction ({dim_reduction})")
        if dim_reduction == "umap":
            adata = DimensionReductionTools.umap(adata, n_neighbors=n_neighbors, random_state=random_state)
        elif dim_reduction == "tsne":
            adata = DimensionReductionTools.tsne(adata, random_state=random_state)
        
        logger.info("Workflow complete")
        
        return adata
