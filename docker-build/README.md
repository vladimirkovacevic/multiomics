# eQTL Analysis Tutorial - Python/Jupyter Notebook

[![Docker Pulls](https://img.shields.io/docker/pulls/vladimirkovacevic/eqtl_intro)](https://hub.docker.com/r/vladimirkovacevic/eqtl_intro)
[![Docker Image Size](https://img.shields.io/docker/image-size/vladimirkovacevic/eqtl_intro/latest)](https://hub.docker.com/r/vladimirkovacevic/eqtl_intro)

A comprehensive Docker image for learning **expression Quantitative Trait Locus (eQTL) analysis** using Python, Jupyter Notebook, and modern statistical libraries.

This is a Python/Jupyter port of the original R-based eQTL tutorial from [https://jknightlab.github.io/eqtl-intro/](https://jknightlab.github.io/eqtl-intro/).

## 📚 What's Included

### Tutorial Content
- **Section 1**: Simple SNP-Gene Expression Associations
- **Section 2**: Handling Confounding Variation
- **Section 3**: Principal Components Analysis
- **Section 4**: Genome-Wide eQTL Analysis

### Software Stack
- Python 3.11
- Jupyter Notebook & JupyterLab
- pandas, numpy, matplotlib, seaborn
- scipy, scikit-learn, statsmodels

### Data
- Simulated genotype data (300 individuals, 10 SNPs)
- Gene expression data (10 genes)
- Covariate data for confounding analysis

---

## 🚀 Quick Start

### Option 1: Run with Data Included (Recommended)

Pull and run the Docker image with data directory mounted:

```bash
# Pull the Docker image
docker pull vladimirkovacevic/eqtl_intro:latest

# Run the container (mount your data directory)
docker run -d -p 8887:8887 \
  -v /path/to/your/data:/home/eqtl/data \
  --name eqtl_tutorial \
  vladimirkovacevic/eqtl_intro:latest
```

**If you don't have data**, you can extract it from the original R Docker image:

```bash
# Create data directory
mkdir -p eqtl_data/simulated

# Extract data from the original R image
docker run --rm humburg/eqtl-intro cat /data/simulated/sim_genotypes.tab > eqtl_data/simulated/sim_genotypes.tab
docker run --rm humburg/eqtl-intro cat /data/simulated/sim_expression1.tab > eqtl_data/simulated/sim_expression1.tab
docker run --rm humburg/eqtl-intro cat /data/simulated/sim_expression2.tab > eqtl_data/simulated/sim_expression2.tab
docker run --rm humburg/eqtl-intro cat /data/simulated/sim_covariates.tab > eqtl_data/simulated/sim_covariates.tab

# Now run with the data
docker run -d -p 8887:8887 \
  -v $(pwd)/eqtl_data:/home/eqtl/data \
  --name eqtl_tutorial \
  vladimirkovacevic/eqtl_intro:latest
```

### Option 2: One-Line Quick Start

```bash
docker run -d -p 8887:8887 --name eqtl_tutorial vladimirkovacevic/eqtl_intro:latest
```

⚠️ **Note**: This option doesn't mount data, so you'll need to add data paths manually in the notebook or use dummy data.

---

## 🌐 Access Jupyter Notebook

Once the container is running:

1. **Open your browser** and navigate to: **`http://localhost:8887`**

2. **Open the tutorial notebook**:
   - In the Jupyter file browser, navigate to: `notebooks/eQTL_tutorial.ipynb`

3. **Start learning!**
   - Run cells sequentially with `Shift + Enter`
   - Read the explanations and interpret the visualizations
   - Modify code to explore different analyses

---

## 📖 Tutorial Structure

### Section 1: SNP-Gene Expression Associations
- Load and explore genotype and expression data
- Calculate Minor Allele Frequencies (MAF)
- Visualize expression distributions by genotype
- Perform linear regression for eQTL mapping
- Create volcano plots and effect size distributions

### Section 2: Confounding Variation
- Load covariate data
- Compare models with and without covariates
- Understand the impact of confounders
- Visualize effect size changes

### Section 3: Principal Components Analysis
- Compute PCs from covariate data
- Select optimal number of components
- Perform eQTL analysis with PC adjustment
- Compare three model types

### Section 4: Genome-Wide Analysis
- Apply multiple testing corrections (Bonferroni, FDR)
- Create Manhattan plots
- Generate effect size heatmaps
- Export results for further analysis

---

## 🛠️ Container Management

### View logs
```bash
docker logs eqtl_tutorial
```

### Stop the container
```bash
docker stop eqtl_tutorial
```

### Start the container again
```bash
docker start eqtl_tutorial
```

### Remove the container
```bash
docker stop eqtl_tutorial
docker rm eqtl_tutorial
```

### Access container shell (for debugging)
```bash
docker exec -it eqtl_tutorial /bin/bash
```

---

## 🧪 Testing the Installation

### Verify Jupyter is running:
```bash
curl http://localhost:8887
```

### Check data availability (if mounted):
```bash
docker exec eqtl_tutorial ls -la /home/eqtl/data/simulated/
```

Expected output:
```
sim_genotypes.tab
sim_expression1.tab
sim_expression2.tab
sim_covariates.tab
```

### Run a test cell:
Open the notebook and run the first cell to verify all imports work correctly.

---

## 📊 Data Format

### Genotype Data (`sim_genotypes.tab`)
- Tab-separated file
- Rows: individuals (300)
- Columns: sample ID + SNPs (10)
- Values: 0, 1, 2 (copies of alternative allele)

### Expression Data (`sim_expression1.tab`)
- Tab-separated file
- Rows: individuals (300)
- Columns: sample ID + genes (10)
- Values: log2-transformed normalized expression

### Covariate Data (`sim_covariates.tab`)
- Tab-separated file
- Contains technical and biological confounders
- Used for adjusted eQTL analysis

---

## 🐍 Python vs R: Key Differences

| Task | R (Original) | Python (This Tutorial) |
|------|--------------|------------------------|
| Data Loading | `readr::read_tsv()` | `pandas.read_csv()` |
| Data Reshaping | `tidyr::gather()` | `pandas.melt()` |
| Visualization | `ggplot2` | `matplotlib + seaborn` |
| Linear Regression | base R `lm()` | `statsmodels.OLS()` |
| PCA | `prcomp()` | `scikit-learn.PCA()` |
| Multiple Testing | custom | `statsmodels.multipletests()` |

---

## 🔧 Advanced Usage

### Use JupyterLab instead of Jupyter Notebook
The image includes JupyterLab. Access it at: `http://localhost:8887/lab`

### Custom Port Mapping
```bash
docker run -d -p 9999:8887 \
  -v $(pwd)/eqtl_data:/home/eqtl/data \
  --name eqtl_tutorial \
  vladimirkovacevic/eqtl_intro:latest
```
Access at: `http://localhost:9999`

### Mount Custom Notebooks
```bash
docker run -d -p 8887:8887 \
  -v $(pwd)/eqtl_data:/home/eqtl/data \
  -v $(pwd)/my_notebooks:/home/eqtl/notebooks \
  --name eqtl_tutorial \
  vladimirkovacevic/eqtl_intro:latest
```

### Save Output Files
Mount a directory for saving results:
```bash
docker run -d -p 8887:8887 \
  -v $(pwd)/eqtl_data:/home/eqtl/data \
  -v $(pwd)/results:/home/eqtl/results \
  --name eqtl_tutorial \
  vladimirkovacevic/eqtl_intro:latest
```

---

## 📚 Learning Resources

### eQTL Analysis
- [GTEx Consortium](https://gtexportal.org) - Genotype-Tissue Expression project
- [eQTL Catalogue](https://www.ebi.ac.uk/eqtl) - Public eQTL datasets

### Python Libraries
- [pandas](https://pandas.pydata.org/docs/) - Data manipulation
- [statsmodels](https://www.statsmodels.org/) - Statistical modeling
- [scikit-learn](https://scikit-learn.org/) - Machine learning
- [seaborn](https://seaborn.pydata.org/) - Statistical visualization

### Original Tutorial
- [eQTL Intro (R version)](https://jknightlab.github.io/eqtl-intro/)

---

## 🐛 Troubleshooting

### Issue: Port 8887 already in use
```bash
# Check what's using the port
lsof -i :8887

# Use a different port
docker run -d -p 9999:8887 --name eqtl_tutorial vladimirkovacevic/eqtl_intro:latest
```

### Issue: Data not found in notebook
- Ensure you mounted the data directory correctly with `-v`
- Check data path in notebook matches `/home/eqtl/data/simulated/`
- Verify files exist: `docker exec eqtl_tutorial ls /home/eqtl/data/simulated/`

### Issue: Out of memory errors
```bash
# Increase Docker memory limit (Docker Desktop settings)
# Or run with memory limit:
docker run -d -p 8887:8887 -m 4g --name eqtl_tutorial vladimirkovacevic/eqtl_intro:latest
```

### Issue: Jupyter doesn't load
- Wait 10-15 seconds for Jupyter to fully start
- Check logs: `docker logs eqtl_tutorial`
- Verify container is running: `docker ps`

---

## 📦 Building from Source

If you want to customize the image:

```bash
git clone <repository>
cd docker-build

# Build the image
docker build -t my-eqtl-tutorial .

# Run your custom image
docker run -d -p 8887:8887 my-eqtl-tutorial
```

---

## 📝 License

This tutorial is based on the original eQTL introduction by Julian Knight's Lab at the University of Oxford, converted to Python for educational purposes.

## 🙏 Acknowledgments

- Original R tutorial: [Julian Knight Lab, University of Oxford](https://jknightlab.github.io/eqtl-intro/)
- Python conversion and Docker packaging: Vladimir Kovacevic

---

## 💬 Support & Feedback

- **Issues**: Report bugs or request features on the repository issue tracker
- **Questions**: Open a discussion or issue for help
- **Contributions**: Pull requests welcome!

---

**Happy Learning! 🎓📊🧬**
