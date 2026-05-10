"""Tool definitions for Claude to use."""

def get_all_tools():
    """Return all available tools formatted for Claude."""
    return [
        {
            "name": "run_qc",
            "description": "Perform quality control filtering on scRNA-seq data. Removes low-quality cells and genes based on gene/UMI counts.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "data_path": {
                        "type": "string",
                        "description": "Path to the input data file (H5AD, CSV, or MTX)"
                    },
                    "min_genes": {
                        "type": "integer",
                        "description": "Minimum number of genes per cell (default: 200)",
                        "default": 200
                    },
                    "min_cells": {
                        "type": "integer",
                        "description": "Minimum number of cells per gene (default: 3)",
                        "default": 3
                    },
                    "max_genes": {
                        "type": "integer",
                        "description": "Maximum number of genes per cell (optional, for removing doublets)"
                    },
                    "mt_threshold": {
                        "type": "number",
                        "description": "Mitochondrial gene percentage threshold (default: 0.2)",
                        "default": 0.2
                    }
                },
                "required": ["data_path"]
            }
        },
        {
            "name": "normalize",
            "description": "Normalize and scale scRNA-seq expression data.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "data_path": {
                        "type": "string",
                        "description": "Path to the preprocessed data"
                    },
                    "method": {
                        "type": "string",
                        "enum": ["log", "sqrt", "none"],
                        "description": "Normalization method (log: log-normalization, sqrt: square root, none: no normalization)",
                        "default": "log"
                    },
                    "target_sum": {
                        "type": "number",
                        "description": "Target sum for normalization (default: 10000)",
                        "default": 10000
                    }
                },
                "required": ["data_path"]
            }
        },
        {
            "name": "select_hvg",
            "description": "Select highly variable genes for downstream analysis.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "data_path": {
                        "type": "string",
                        "description": "Path to the normalized data"
                    },
                    "n_genes": {
                        "type": "integer",
                        "description": "Number of highly variable genes to select (default: 2000)",
                        "default": 2000
                    },
                    "flavor": {
                        "type": "string",
                        "enum": ["seurat", "cell_ranger", "seurat_v3"],
                        "description": "Method for HVG selection",
                        "default": "seurat_v3"
                    }
                },
                "required": ["data_path"]
            }
        },
        {
            "name": "reduce_dimensions",
            "description": "Apply dimensionality reduction techniques (PCA, UMAP, t-SNE).",
            "input_schema": {
                "type": "object",
                "properties": {
                    "data_path": {
                        "type": "string",
                        "description": "Path to the preprocessed data"
                    },
                    "method": {
                        "type": "string",
                        "enum": ["pca", "umap", "tsne"],
                        "description": "Dimensionality reduction method",
                        "default": "umap"
                    },
                    "n_components": {
                        "type": "integer",
                        "description": "Number of components (default: 2 for visualization, 50 for PCA)",
                        "default": 2
                    },
                    "n_neighbors": {
                        "type": "integer",
                        "description": "Number of neighbors for UMAP/t-SNE (default: 15)",
                        "default": 15
                    }
                },
                "required": ["data_path", "method"]
            }
        },
        {
            "name": "cluster",
            "description": "Identify cell clusters using various clustering algorithms.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "data_path": {
                        "type": "string",
                        "description": "Path to the dimensionality-reduced data"
                    },
                    "method": {
                        "type": "string",
                        "enum": ["leiden", "louvain", "kmeans"],
                        "description": "Clustering algorithm",
                        "default": "leiden"
                    },
                    "resolution": {
                        "type": "number",
                        "description": "Resolution parameter for Leiden/Louvain (higher = more clusters, default: 1.0)",
                        "default": 1.0
                    },
                    "n_clusters": {
                        "type": "integer",
                        "description": "Number of clusters for K-means (default: 10)"
                    },
                    "random_state": {
                        "type": "integer",
                        "description": "Random seed for reproducibility (default: 0)",
                        "default": 0
                    }
                },
                "required": ["data_path", "method"]
            }
        },
        {
            "name": "find_markers",
            "description": "Identify marker genes that define each cluster.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "data_path": {
                        "type": "string",
                        "description": "Path to the clustered data"
                    },
                    "method": {
                        "type": "string",
                        "enum": ["wilcoxon", "t-test", "logreg"],
                        "description": "Statistical test for marker detection",
                        "default": "wilcoxon"
                    },
                    "min_log2fc": {
                        "type": "number",
                        "description": "Minimum log2 fold-change for markers (default: 0.25)",
                        "default": 0.25
                    },
                    "min_pct": {
                        "type": "number",
                        "description": "Minimum percentage of cells expressing marker (default: 0.1)",
                        "default": 0.1
                    },
                    "n_top_genes": {
                        "type": "integer",
                        "description": "Top N genes to return per cluster (default: 10)",
                        "default": 10
                    }
                },
                "required": ["data_path"]
            }
        },
        {
            "name": "annotate_cells",
            "description": "Automatically annotate cell types based on marker gene expression.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "data_path": {
                        "type": "string",
                        "description": "Path to the processed data with clusters"
                    },
                    "method": {
                        "type": "string",
                        "enum": ["celltypist", "marker_based", "reference"],
                        "description": "Cell type annotation method",
                        "default": "celltypist"
                    },
                    "reference": {
                        "type": "string",
                        "description": "Reference dataset or marker file (for marker-based annotation)",
                        "default": "Immune_All_Low.pkl"
                    },
                    "confidence_threshold": {
                        "type": "number",
                        "description": "Confidence threshold for annotations (default: 0.5)",
                        "default": 0.5
                    }
                },
                "required": ["data_path"]
            }
        }
    ]
