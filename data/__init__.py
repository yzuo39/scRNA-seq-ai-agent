"""Data loading and preprocessing utilities."""

import logging
from pathlib import Path
from typing import Union, Optional
import numpy as np
import pandas as pd
import anndata
from scipy.io import mmread

logger = logging.getLogger(__name__)


class DataLoader:
    """Load scRNA-seq data from various formats."""

    SUPPORTED_FORMATS = [".h5ad", ".mtx", ".csv", ".npz", ".xlsx"]

    @staticmethod
    def load(filepath: Union[str, Path], first_column_names: bool = False) -> anndata.AnnData:
        """Load data from file.

        Args:
            filepath: Path to data file
            first_column_names: Whether first column contains cell/gene names

        Returns:
            AnnData object
        """
        filepath = Path(filepath)
        suffix = filepath.suffix.lower()

        if suffix == ".h5ad":
            return DataLoader._load_h5ad(filepath)
        elif suffix == ".mtx":
            return DataLoader._load_mtx(filepath)
        elif suffix == ".csv":
            return DataLoader._load_csv(filepath, first_column_names)
        elif suffix == ".npz":
            return DataLoader._load_npz(filepath)
        elif suffix == ".xlsx":
            return DataLoader._load_xlsx(filepath, first_column_names)
        else:
            raise ValueError(
                f"Unsupported format: {suffix}. Supported: {DataLoader.SUPPORTED_FORMATS}"
            )

    @staticmethod
    def _load_h5ad(filepath: Path) -> anndata.AnnData:
        """Load H5AD format."""
        logger.info(f"Loading H5AD file: {filepath}")
        return anndata.read_h5ad(filepath)

    @staticmethod
    def _load_mtx(filepath: Path) -> anndata.AnnData:
        """Load Matrix Market format."""
        logger.info(f"Loading MTX file: {filepath}")
        X = mmread(filepath).T.tocsr()

        # Look for barcodes and features files
        parent = filepath.parent
        barcodes_file = parent / "barcodes.tsv.gz"
        features_file = parent / "features.tsv.gz"

        obs_names = None
        var_names = None

        if barcodes_file.exists():
            obs_names = pd.read_csv(barcodes_file, header=None, sep="\t")[0].values
        if features_file.exists():
            features_df = pd.read_csv(features_file, header=None, sep="\t")
            var_names = features_df[1].values if features_df.shape[1] > 1 else features_df[0].values

        adata = anndata.AnnData(X)
        if obs_names is not None:
            adata.obs_names = obs_names
        if var_names is not None:
            adata.var_names = var_names

        return adata

    @staticmethod
    def _load_csv(filepath: Path, first_column_names: bool = False) -> anndata.AnnData:
        """Load CSV format."""
        logger.info(f"Loading CSV file: {filepath}")
        df = pd.read_csv(filepath, index_col=0 if first_column_names else None)
        return anndata.AnnData(
            df.values,
            obs=pd.DataFrame(index=df.index),
            var=pd.DataFrame(index=df.columns)
        )

    @staticmethod
    def _load_npz(filepath: Path) -> anndata.AnnData:
        """Load NPZ format."""
        logger.info(f"Loading NPZ file: {filepath}")
        data = np.load(filepath, allow_pickle=True)
        X = data["X"]
        obs_names = data.get("obs_names", None)
        var_names = data.get("var_names", None)

        adata = anndata.AnnData(X)
        if obs_names is not None:
            adata.obs_names = obs_names
        if var_names is not None:
            adata.var_names = var_names

        return adata

    @staticmethod
    def _load_xlsx(filepath: Path, first_column_names: bool = False) -> anndata.AnnData:
        """Load Excel format."""
        logger.info(f"Loading Excel file: {filepath}")
        df = pd.read_excel(filepath, index_col=0 if first_column_names else None)
        return anndata.AnnData(
            df.values,
            obs=pd.DataFrame(index=df.index),
            var=pd.DataFrame(index=df.columns)
        )


class Preprocessor:
    """Preprocessing utilities for scRNA-seq data."""

    @staticmethod
    def quality_control(
        adata: anndata.AnnData,
        min_genes: int = 200,
        min_cells: int = 3,
        max_genes: Optional[int] = None,
        mt_threshold: float = 0.2,
    ) -> anndata.AnnData:
        """Apply quality control filtering.
        
        Args:
            adata: Input AnnData object
            min_genes: Minimum genes per cell
            min_cells: Minimum cells per gene
            max_genes: Maximum genes per cell
            mt_threshold: Mitochondrial gene percentage threshold
            
        Returns:
            Filtered AnnData object
        """
        logger.info("Running quality control...")
        adata = adata.copy()

        # Calculate QC metrics
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        import scanpy as sc
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

        # Filter cells
        initial_cells = adata.n_obs
        adata = adata[adata.obs["n_genes_by_counts"] >= min_genes]
        if max_genes is not None:
            adata = adata[adata.obs["n_genes_by_counts"] <= max_genes]
        adata = adata[adata.obs["pct_counts_mt"] < mt_threshold * 100]

        # Filter genes
        adata = adata[:, adata.var["n_cells_by_counts"] >= min_cells]

        logger.info(f"Cells: {initial_cells} → {adata.n_obs}")
        logger.info(f"Genes: {adata.n_vars}")

        return adata

    @staticmethod
    def normalize(
        adata: anndata.AnnData,
        method: str = "log",
        target_sum: float = 10000,
    ) -> anndata.AnnData:
        """Normalize expression data.
        
        Args:
            adata: Input data
            method: 'log', 'sqrt', or 'none'
            target_sum: Target sum for normalization
            
        Returns:
            Normalized AnnData object
        """
        logger.info(f"Normalizing with {method} method...")
        adata = adata.copy()

        import scanpy as sc

        if method == "log":
            sc.pp.normalize_total(adata, target_sum=target_sum)
            sc.pp.log1p(adata)
        elif method == "sqrt":
            sc.pp.normalize_total(adata, target_sum=target_sum)
            adata.X = np.sqrt(adata.X.data) if hasattr(adata.X, 'data') else np.sqrt(adata.X)
        elif method != "none":
            raise ValueError(f"Unknown normalization method: {method}")

        return adata

    @staticmethod
    def select_hvg(
        adata: anndata.AnnData,
        n_genes: int = 2000,
        flavor: str = "seurat_v3",
    ) -> anndata.AnnData:
        """Select highly variable genes.
        
        Args:
            adata: Input data
            n_genes: Number of genes to select
            flavor: HVG selection method
            
        Returns:
            Filtered AnnData object with HVG only
        """
        logger.info(f"Selecting {n_genes} highly variable genes...")
        adata = adata.copy()

        import scanpy as sc
        sc.pp.highly_variable_genes(adata, n_top_genes=n_genes, flavor=flavor)
        adata = adata[:, adata.var["highly_variable"]]

        return adata
