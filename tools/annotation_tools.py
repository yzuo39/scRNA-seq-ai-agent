"""Cell type annotation tools."""

import logging
from typing import Dict, Any, List
import numpy as np
import pandas as pd
import anndata
import scanpy as sc

logger = logging.getLogger(__name__)


class AnnotationTools:
    """Tools for cell type annotation and marker gene identification."""

    @staticmethod
    def find_markers(
        data_path: str,
        method: str = "wilcoxon",
        min_log2fc: float = 0.25,
        min_pct: float = 0.1,
        n_top_genes: int = 10,
    ) -> Dict[str, Any]:
        """Find marker genes for each cluster.
        
        Args:
            data_path: Path to clustered data
            method: 'wilcoxon', 't-test', or 'logreg'
            min_log2fc: Minimum log2 fold-change
            min_pct: Minimum percentage of cells expressing marker
            n_top_genes: Top N genes per cluster
            
        Returns:
            Dictionary with marker gene results
        """
        logger.info(f"Loading data from {data_path}...")
        from data import DataLoader
        
        # Load data
        adata = DataLoader.load(data_path)
        logger.info(f"Loaded data: {adata.shape}")
        
        # Determine cluster column
        cluster_key = None
        for key in ["leiden", "louvain", "kmeans", "cluster"]:
            if key in adata.obs:
                cluster_key = key
                break
        
        if cluster_key is None:
            logger.warning("No cluster column found, using first categorical obs")
            cluster_key = adata.obs.select_dtypes(include=["category"]).columns[0]
        
        logger.info(f"Finding markers using {method} method...")
        
        # Find markers
        sc.tl.rank_genes_groups(
            adata,
            groupby=cluster_key,
            method=method,
            key_added="rank_genes_groups",
        )
        
        # Extract top markers
        markers = {}
        for cluster in adata.obs[cluster_key].unique():
            cluster_markers = sc.get.rank_genes_groups_df(
                adata, group=cluster, key="rank_genes_groups"
            ).head(n_top_genes)
            markers[str(cluster)] = cluster_markers[["names", "logfoldchanges", "pvals_adj"]].to_dict(orient="list")
        
        # Save results
        output_file = data_path.replace(".h5ad", f"_markers_{method}.h5ad")
        adata.write_h5ad(output_file)
        
        results = {
            "status": "markers_found",
            "method": method,
            "n_clusters": len(markers),
            "n_genes_per_cluster": n_top_genes,
            "markers": markers,
            "output_file": output_file,
        }
        
        logger.info(f"Found markers for {len(markers)} clusters")
        return results

    @staticmethod
    def annotate_cells(
        data_path: str,
        method: str = "celltypist",
        reference: str = "Immune_All_Low.pkl",
        confidence_threshold: float = 0.5,
    ) -> Dict[str, Any]:
        """Annotate cell types using marker genes or CellTypist.
        
        Args:
            data_path: Path to clustered data
            method: 'celltypist', 'marker_based', or 'reference'
            reference: Reference dataset or marker file
            confidence_threshold: Confidence threshold for annotations
            
        Returns:
            Dictionary with annotation results
        """
        logger.info(f"Loading data from {data_path}...")
        from data import DataLoader
        
        # Load data
        adata = DataLoader.load(data_path)
        logger.info(f"Loaded data: {adata.shape}")
        
        if method == "celltypist":
            logger.info("Running CellTypist annotation...")
            try:
                import celltypist
                
                # Predict cell types
                predictions = celltypist.predict(
                    adata,
                    model=celltypist.models.Model.fetch(reference),
                    majority_voting=True,
                )
                
                # Extract annotations
                adata.obs["cell_type"] = predictions.pred_label
                adata.obs["cell_type_confidence"] = predictions.confidence_score
                
                # Filter by confidence
                adata.obs["cell_type_filtered"] = adata.obs["cell_type"].where(
                    adata.obs["cell_type_confidence"] >= confidence_threshold,
                    "Unknown"
                )
                
                logger.info(f"Annotated {(adata.obs['cell_type_confidence'] >= confidence_threshold).sum()} cells")
                
            except ImportError:
                logger.warning("CellTypist not installed, using marker-based annotation instead")
                return AnnotationTools._annotate_marker_based(adata)
        
        elif method == "marker_based":
            return AnnotationTools._annotate_marker_based(adata)
        
        else:
            raise ValueError(f"Unknown annotation method: {method}")
        
        # Save annotated data
        output_file = data_path.replace(".h5ad", f"_annotated_{method}.h5ad")
        adata.write_h5ad(output_file)
        
        # Get cell type distribution
        cell_types = adata.obs["cell_type"].value_counts()
        
        results = {
            "status": "annotated",
            "method": method,
            "n_cell_types": len(cell_types),
            "cell_type_distribution": cell_types.to_dict(),
            "output_file": output_file,
        }
        
        return results

    @staticmethod
    def _annotate_marker_based(adata: anndata.AnnData) -> Dict[str, Any]:
        """Marker-based cell type annotation.
        
        Args:
            adata: AnnData object
            
        Returns:
            Dictionary with annotation results
        """
        logger.info("Using marker-based annotation...")
        
        # Define marker genes for common cell types
        marker_genes = {
            "T_cells": ["CD3D", "CD3E", "CD3G"],
            "B_cells": ["CD19", "CD79A", "CD79B"],
            "Monocytes": ["CD14", "CD16", "FCGR3A"],
            "Macrophages": ["CD14", "CD68", "MARCO"],
            "NK_cells": ["NCAM1", "NKG7", "GZMB"],
            "Neutrophils": ["FCGR3B", "CSF3R", "ELANE"],
            "Dendritic_cells": ["FCER1A", "DC-SIGN", "CD83"],
        }
        
        # Score cell types
        for cell_type, genes in marker_genes.items():
            available_genes = [g for g in genes if g in adata.var_names]
            if available_genes:
                sc.tl.score_genes(adata, available_genes, score_name=f"{cell_type}_score")
        
        # Assign cell types based on highest score
        score_cols = [col for col in adata.obs.columns if col.endswith("_score")]
        if score_cols:
            adata.obs["cell_type"] = adata.obs[score_cols].idxmax(axis=1).str.replace("_score", "")
        else:
            logger.warning("No marker genes found in data")
            adata.obs["cell_type"] = "Unknown"
        
        return {
            "status": "annotated",
            "method": "marker_based",
            "n_cell_types": len(adata.obs["cell_type"].unique()),
            "cell_type_distribution": adata.obs["cell_type"].value_counts().to_dict(),
        }

    @staticmethod
    def get_annotation_summary(adata: anndata.AnnData, annotation_key: str = "cell_type") -> Dict[str, Any]:
        """Get summary of cell type annotations.
        
        Args:
            adata: AnnData object with annotations
            annotation_key: Key in obs with annotations
            
        Returns:
            Dictionary with annotation summary
        """
        if annotation_key not in adata.obs:
            return {"error": f"Annotation key '{annotation_key}' not found"}
        
        summary = {
            "n_cell_types": len(adata.obs[annotation_key].unique()),
            "cell_type_distribution": adata.obs[annotation_key].value_counts().to_dict(),
            "total_cells": adata.n_obs,
        }
        
        return summary
