"""Example: Basic scRNA-seq analysis with the AI agent."""

import os
import sys
from pathlib import Path

# Add project to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from agent import ScRNAAgent


def main():
    """Run basic analysis example."""
    
    # Initialize agent
    print("Initializing ScRNA Agent...")
    agent = ScRNAAgent(api_key=os.getenv("ANTHROPIC_API_KEY"))
    
    # Example data path (you should replace with your actual data)
    data_path = "data/example_data.h5ad"
    
    # Check if example data exists
    if not Path(data_path).exists():
        print(f"Note: Example data not found at {data_path}")
        print("To use this script, provide your scRNA-seq data in H5AD format.")
        print("\nExample usage:")
        print("  python examples/basic_analysis.py")
        print("\nThe agent will ask Claude to:")
        print("  1. Run quality control filtering")
        print("  2. Normalize expression data")
        print("  3. Select highly variable genes")
        print("  4. Perform dimensionality reduction (UMAP)")
        print("  5. Cluster cells (Leiden algorithm)")
        print("  6. Find marker genes for each cluster")
        print("  7. Annotate cell types")
        return
    
    # Define analysis query
    query = """
    Please perform a complete single-cell RNA-seq analysis:
    1. Run quality control with default parameters
    2. Normalize the data using log normalization
    3. Select the top 2000 highly variable genes
    4. Reduce dimensions with UMAP
    5. Cluster cells using the Leiden algorithm (resolution=1.0)
    6. Find marker genes for each cluster
    7. Annotate cell types based on marker genes
    8. Provide a summary of the findings
    """
    
    print("\n" + "="*80)
    print("STARTING SCRNA-SEQ ANALYSIS")
    print("="*80 + "\n")
    
    # Run analysis
    results = agent.analyze(
        data_path=data_path,
        query=query,
        verify=True
    )
    
    # Display results
    print("\n" + "="*80)
    print("ANALYSIS PLAN")
    print("="*80)
    print(results["plan"])
    
    if results["verification"]:
        print("\n" + "="*80)
        print("VERIFICATION RESULTS")
        print("="*80)
        print(results["verification"])
    
    print("\n" + "="*80)
    print("FINAL REPORT")
    print("="*80)
    print(results["report"])
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)


if __name__ == "__main__":
    main()
