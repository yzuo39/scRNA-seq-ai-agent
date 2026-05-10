"""Preprocessing tools for normalization and feature selection."""

import logging
from typing import Dict, Any
import numpy as np
import anndata
import scanpy as sc

logger = logging.getLogger(__name__)


class PreprocessingTools:
    """Tools for data normalization and scaling."""

    @staticmethod
    def normalize(
        data_path: str,
        method: str = "log",
        target_sum: float = 10000,
    ) -> Dict[str, Any]:
        """Normalize scRNA-seq expression data.
        
        Args:
            data_path: Path to input data
            method: 'log', 'sqrt', or 'none'
            target_sum: Target sum for normalization
            
        Returns:
            Dictionary with normalization results
        """
        logger.info(f"Loading data from {data_path}...")
        from data import DataLoader, Preprocessor
        
        # Load data
        adata = DataLoader.load(data_path)
        logger.info(f"Loaded data: {adata.shape}")
        
        # Apply normalization
        adata = Preprocessor.normalize(adata, method=method, target_sum=target_sum)
        
        # Save normalized data
        output_file = data_path.replace(".h5ad", f"_normalized_{method}.h5ad")
        adata.write_h5ad(output_file)
        
        results = {
            "status": "normalized",
            "method": method,
            "target_sum": target_sum,
            "n_cells": adata.n_obs,
            "n_genes": adata.n_vars,
            "output_file": output_file,
        }
        
        logger.info(f"Saved normalized data to {output_file}")
        return results

    @staticmethod
    def scale(
        data_path: str,
        zero_center: bool = True,
        max_value: float = 10.0,
    ) -> Dict[str, Any]:
        """Scale normalized expression data.
        
        Args:
            data_path: Path to input data
            zero_center: Whether to center data
            max_value: Maximum value for clipping
            
        Returns:
            Dictionary with scaling results
        """
        logger.info(f"Loading data from {data_path}...")
        from data import DataLoader
        
        # Load data
        adata = DataLoader.load(data_path)
        
        # Scale
        logger.info("Scaling data...")
        sc.pp.scale(adata, zero_center=zero_center, max_value=max_value)
        
        # Save scaled data
        output_file = data_path.replace(".h5ad", "_scaled.h5ad")
        adata.write_h5ad(output_file)
        
        results = {
            "status": "scaled",
            "zero_center": zero_center,
            "max_value": max_value,
            "output_file": output_file,
        }
        
        logger.info(f"Saved scaled data to {output_file}")
        return results

    @staticmethod
    def select_hvg(
        data_path: str,
        n_genes: int = 2000,
        flavor: str = "seurat_v3",
    ) -> Dict[str, Any]:
        """Select highly variable genes.
        
        Args:
            data_path: Path to input data
            n_genes: Number of genes to select
            flavor: HVG selection method
            
        Returns:
            Dictionary with HVG results
        """
        logger.info(f"Loading data from {data_path}...")
        from data import DataLoader, Preprocessor
        
        # Load data
        adata = DataLoader.load(data_path)
        
        # Select HVG
        adata = Preprocessor.select_hvg(adata, n_genes=n_genes, flavor=flavor)
        
        # Save HVG data
        output_file = data_path.replace(".h5ad", f"_hvg_{n_genes}.h5ad")
        adata.write_h5ad(output_file)
        
        results = {
            "status": "hvg_selected",
            "n_genes_selected": adata.n_vars,
            "n_genes_total": adata.n_vars,
            "flavor": flavor,
            "output_file": output_file,
        }
        
        logger.info(f"Saved HVG-selected data to {output_file}")
        return results
