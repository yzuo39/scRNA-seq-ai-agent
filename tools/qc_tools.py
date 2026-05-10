"""Quality control tools for scRNA-seq analysis."""

import logging
from typing import Dict, Any
import numpy as np
import anndata
import scanpy as sc

logger = logging.getLogger(__name__)


class QCTools:
    """Quality control and filtering tools."""

    @staticmethod
    def run_qc(
        data_path: str,
        min_genes: int = 200,
        min_cells: int = 3,
        max_genes: int = None,
        mt_threshold: float = 0.2,
    ) -> Dict[str, Any]:
        """Perform quality control on scRNA-seq data.
        
        Args:
            data_path: Path to input data
            min_genes: Minimum genes per cell
            min_cells: Minimum cells per gene
            max_genes: Maximum genes per cell (for doublet detection)
            mt_threshold: Mitochondrial percentage threshold
            
        Returns:
            Dictionary with QC statistics
        """
        logger.info(f"Loading data from {data_path}...")
        from data import DataLoader, Preprocessor
        
        # Load data
        adata = DataLoader.load(data_path)
        logger.info(f"Loaded data: {adata.shape}")
        
        # Calculate QC metrics
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
        
        initial_cells = adata.n_obs
        initial_genes = adata.n_vars
        
        # Apply QC
        adata = Preprocessor.quality_control(
            adata,
            min_genes=min_genes,
            min_cells=min_cells,
            max_genes=max_genes,
            mt_threshold=mt_threshold,
        )
        
        # Prepare results
        results = {
            "status": "qc_complete",
            "initial_cells": int(initial_cells),
            "final_cells": int(adata.n_obs),
            "cells_removed": int(initial_cells - adata.n_obs),
            "initial_genes": int(initial_genes),
            "final_genes": int(adata.n_vars),
            "genes_removed": int(initial_genes - adata.n_vars),
            "mean_genes_per_cell": float(adata.obs["n_genes_by_counts"].mean()),
            "mean_counts_per_cell": float(adata.obs["total_counts"].mean()),
            "mean_mt_percentage": float(adata.obs["pct_counts_mt"].mean()),
            "output_file": data_path.replace(".h5ad", "_qc.h5ad"),
        }
        
        # Save QC'd data
        adata.write_h5ad(results["output_file"])
        logger.info(f"Saved QC'd data to {results['output_file']}")
        
        return results

    @staticmethod
    def get_qc_summary(adata: anndata.AnnData) -> Dict[str, Any]:
        """Get summary statistics from QC-filtered data.
        
        Args:
            adata: AnnData object with QC metrics
            
        Returns:
            Dictionary with QC statistics
        """
        summary = {
            "n_cells": adata.n_obs,
            "n_genes": adata.n_vars,
            "mean_genes_per_cell": float(adata.obs["n_genes_by_counts"].mean()),
            "median_genes_per_cell": float(adata.obs["n_genes_by_counts"].median()),
            "mean_counts_per_cell": float(adata.obs["total_counts"].mean()),
            "median_counts_per_cell": float(adata.obs["total_counts"].median()),
        }
        
        if "pct_counts_mt" in adata.obs:
            summary["mean_mt_percentage"] = float(adata.obs["pct_counts_mt"].mean())
            summary["median_mt_percentage"] = float(adata.obs["pct_counts_mt"].median())
        
        return summary
