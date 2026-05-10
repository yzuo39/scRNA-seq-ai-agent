# scRNA-seq AI Agent

An intelligent AI agent for analyzing single-cell RNA sequencing (scRNA-seq) data. Uses Claude as the backbone with specialized agents for planning, observation, verification, and reporting.

## 🏗️ Architecture

```
User Query
    ↓
Planning Agent (Claude)
    ↓ decides which tools to call
Tool Layer
├── run_qc()
├── normalize()
├── cluster()
├── find_markers()
└── annotate()
    ↓
Observation Agent (Claude)
    ↓ interprets results
Verification Agent (Claude)
    ↓ checks for errors
Report Agent (Claude)
    ↓
Final Report
```

## 🚀 Quick Start

### Installation

```bash
git clone https://github.com/yzuo39/scRNA-seq-ai-agent.git
cd scRNA-seq-ai-agent

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Basic Usage

```python
from agent import ScRNAAgent

# Initialize agent
agent = ScRNAAgent(api_key="your-anthropic-api-key")

# Analyze data
results = agent.analyze(
    data_path="path/to/data.h5ad",
    query="Perform quality control, normalize, and cluster the data"
)

# Get report
report = results.get_report()
print(report)
```

## 📋 Features

- ✅ **Quality Control**: Filter low-quality cells and genes
- ✅ **Normalization**: Multiple normalization methods (log, sqrt, etc.)
- ✅ **Feature Selection**: Highly variable gene selection
- ✅ **Dimensionality Reduction**: PCA, t-SNE, UMAP
- ✅ **Clustering**: Leiden, Louvain, K-means
- ✅ **Cell Type Annotation**: Marker gene identification
- ✅ **Intelligent Planning**: Claude decides analysis workflow
- ✅ **Result Interpretation**: Claude interprets findings
- ✅ **Error Verification**: Automatic validation
- ✅ **Report Generation**: Comprehensive analysis reports

## 📁 Project Structure

```
scRNA-seq-ai-agent/
├── agent/
│   ├── __init__.py
│   ├── main.py                  # Main ScRNAAgent class
│   ├── planner.py               # Planning agent
│   ├── observer.py              # Observation agent
│   ├── verifier.py              # Verification agent
│   └── reporter.py              # Report generation agent
├── tools/
│   ├── __init__.py
│   ├── data_tools.py            # Data loading/saving
│   ├── qc_tools.py              # Quality control
│   ├── preprocessing_tools.py   # Normalization, scaling
│   ├── analysis_tools.py        # Clustering, reduction
│   └── annotation_tools.py      # Cell type annotation
├── data/
│   ├── __init__.py
│   ├── loader.py                # Load H5AD, CSV, MTX
│   ├── preprocessor.py          # QC & preprocessing
│   └── normalizer.py            # Normalization methods
├── visualization/
│   ├── __init__.py
│   └── plots.py                 # Plotting utilities
├── examples/
│   ├── basic_analysis.py
│   ├── custom_workflow.py
│   └── example_data/
├── tests/
│   ├── test_agent.py
│   ├── test_tools.py
│   └── test_data.py
├── requirements.txt
├── setup.py
└── config.yaml
```

## 🧠 Agent Components

### 1. Planning Agent
- Receives user query
- Analyzes data characteristics
- Plans optimal analysis workflow
- Selects appropriate tools

### 2. Tool Layer
Core analysis tools:
- `run_qc()` - Quality control filtering
- `normalize()` - Normalization & scaling
- `select_hvg()` - Highly variable gene selection
- `reduce_dimensions()` - PCA, t-SNE, UMAP
- `cluster()` - Leiden, Louvain, K-means
- `find_markers()` - Marker gene identification
- `annotate_cells()` - Cell type annotation

### 3. Observation Agent
- Interprets tool execution results
- Extracts biological insights
- Identifies cell populations
- Detects patterns and correlations

### 4. Verification Agent
- Validates results
- Checks for errors
- Ensures quality metrics are met
- Suggests corrections if needed

### 5. Report Agent
- Compiles findings
- Generates visualizations
- Creates summary statistics
- Produces final report

## 🔧 Configuration

Edit `config.yaml` to customize:

```yaml
preprocessing:
  min_genes: 200
  min_cells: 3
  n_top_genes: 2000
  
analysis:
  clustering_method: leiden
  resolution: 1.0
  n_neighbors: 15
  
visualization:
  figsize: [12, 8]
  dpi: 100
```

## 📚 Supported Data Formats

- **H5AD** (.h5ad) - AnnData format
- **Matrix Market** (.mtx) - Sparse matrix
- **CSV** (.csv) - Dense expression matrix
- **NPZ** (.npz) - NumPy compressed

## 🧪 Testing

```bash
pytest tests/ -v
pytest tests/ --cov=agent --cov-report=html
```

## 📦 Dependencies

- Python 3.8+
- anthropic >= 0.7.0
- torch >= 1.9.0
- scanpy >= 1.8.0
- anndata >= 0.8.0
- numpy, pandas, scipy
- matplotlib, seaborn
- scikit-learn

## 📄 License

MIT License

## 🚀 Roadmap

- [ ] Web UI for interactive analysis
- [ ] Multi-sample batch analysis
- [ ] Trajectory inference
- [ ] Gene ontology enrichment
- [ ] Integration with public databases
- [ ] Custom tool registration

---

**Last Updated**: 2026-05-10
