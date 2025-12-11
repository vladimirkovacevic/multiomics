# eQTL Analysis Tutorial - Python Edition

[![Docker Pulls](https://img.shields.io/docker/pulls/vladimirkovacevic/eqtl_intro)](https://hub.docker.com/r/vladimirkovacevic/eqtl_intro)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/)

A comprehensive tutorial for learning **expression Quantitative Trait Locus (eQTL) analysis** using Python, Jupyter Notebook, and modern statistical libraries. This is a Python port of the excellent [R-based eQTL tutorial](https://jknightlab.github.io/eqtl-intro/) from Julian Knight's Lab at the University of Oxford, specifically based on the [exercises and solutions](https://jknightlab.github.io/eqtl-intro/exercises/exercises_and_solutions.html) section.

## 📚 Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Tutorial Content](#tutorial-content)
- [Installation](#installation)
  - [Option 1: Docker (Recommended)](#option-1-docker-recommended)
  - [Option 2: Local Installation](#option-2-local-installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Data](#data)
- [Learning Path](#learning-path)
- [Repository Structure](#repository-structure)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)
- [Acknowledgments](#acknowledgments)

---

## 🎯 Overview

Expression Quantitative Trait Loci (eQTL) analysis identifies genetic variants that influence gene expression levels. This tutorial provides a hands-on introduction to eQTL mapping, covering fundamental concepts through advanced genome-wide analysis.

### What You'll Learn

- **SNP-Gene Association Testing**: Basic linear regression for eQTL discovery
- **Confounding Factor Analysis**: Controlling for technical and biological variation
- **Principal Components Analysis**: Dimensionality reduction for improved power
- **Genome-Wide Analysis**: Multiple testing correction and result interpretation
- **Professional Visualization**: Creating publication-quality figures

---

## ✨ Features

- 📓 **Complete Jupyter Notebook** with detailed explanations and code
- 🐍 **Modern Python Stack** (pandas, numpy, matplotlib, seaborn, scikit-learn, statsmodels)
- 🐳 **Docker Support** for reproducible environments
- 📊 **Professional Visualizations** (volcano plots, Manhattan plots, heatmaps)
- 🎓 **Educational Focus** with step-by-step guidance
- 🔬 **Real Data Examples** using simulated genotype and expression data
- 📈 **Statistical Rigor** with proper multiple testing corrections

---

## 📖 Tutorial Content

### Section 1: Simple SNP-Gene Expression Associations
- Data loading and quality control
- Minor Allele Frequency (MAF) calculation
- Expression visualization by genotype
- Linear regression for effect size estimation
- Confidence intervals and p-values

### Section 2: Handling Confounding Variation
- Understanding technical and biological confounders
- Incorporating covariates in regression models
- Comparing simple vs. adjusted models
- Visualizing effect size changes

### Section 3: Principal Components Analysis (PCA)
- Computing principal components from covariates
- Selecting optimal number of components
- eQTL analysis with PC-adjusted models
- Model comparison and interpretation

### Section 4: Genome-Wide eQTL Analysis
- Multiple testing corrections (Bonferroni, FDR)
- Manhattan plots for visualization
- Effect size heatmaps
- Result summarization and export

---

## 🚀 Installation

### Option 1: Docker (Recommended)

Docker provides a complete, reproducible environment with all dependencies pre-installed.

#### Prerequisites
- [Docker](https://docs.docker.com/get-docker/) installed on your system
- At least 4GB RAM available
- Internet connection for pulling the image

#### Pull and Run

```bash
# Pull the Docker image
docker pull vladimirkovacevic/eqtl_intro:latest

# Extract data from the original R tutorial
mkdir -p eqtl_data/simulated
docker run --rm humburg/eqtl-intro cat /data/simulated/sim_genotypes.tab > eqtl_data/simulated/sim_genotypes.tab
docker run --rm humburg/eqtl-intro cat /data/simulated/sim_expression1.tab > eqtl_data/simulated/sim_expression1.tab
docker run --rm humburg/eqtl-intro cat /data/simulated/sim_expression2.tab > eqtl_data/simulated/sim_expression2.tab
docker run --rm humburg/eqtl-intro cat /data/simulated/sim_covariates.tab > eqtl_data/simulated/sim_covariates.tab

# Run the container
docker run -d -p 8887:8887 \
  -v $(pwd)/eqtl_data:/home/eqtl/data \
  --name eqtl_tutorial \
  vladimirkovacevic/eqtl_intro:latest
```

#### Access Jupyter
Open your browser and navigate to: **http://localhost:8887**

Open the notebook: `notebooks/eQTL_tutorial.ipynb`

#### Managing the Container

```bash
# View logs
docker logs eqtl_tutorial

# Stop the container
docker stop eqtl_tutorial

# Start the container again
docker start eqtl_tutorial

# Remove the container
docker stop eqtl_tutorial
docker rm eqtl_tutorial
```

---

### Option 2: Local Installation

For local installation without Docker, you'll need Python 3.11+ and the required packages.

#### Prerequisites

- Python 3.11 or higher
- pip (Python package manager)
- Git

#### Clone the Repository

```bash
git clone https://github.com/yourusername/eqtl-tutorial-python.git
cd eqtl-tutorial-python
```

#### Install Dependencies

```bash
# Create a virtual environment (recommended)
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install required packages
pip install -r docker-build/requirements.txt
```

Required packages:
- jupyter==1.0.0
- jupyterlab==4.0.9
- pandas==2.2.3
- numpy==1.26.4
- matplotlib==3.8.2
- seaborn==0.13.0
- scipy==1.11.4
- scikit-learn==1.3.2
- statsmodels==0.14.1

#### Extract Data

```bash
# Create data directory
mkdir -p data/simulated

# Extract from Docker image
docker run --rm humburg/eqtl-intro cat /data/simulated/sim_genotypes.tab > data/simulated/sim_genotypes.tab
docker run --rm humburg/eqtl-intro cat /data/simulated/sim_expression1.tab > data/simulated/sim_expression1.tab
docker run --rm humburg/eqtl-intro cat /data/simulated/sim_expression2.tab > data/simulated/sim_expression2.tab
docker run --rm humburg/eqtl-intro cat /data/simulated/sim_covariates.tab > data/simulated/sim_covariates.tab
```

Or download directly from the [original tutorial](https://jknightlab.github.io/eqtl-intro/).

#### Start Jupyter

```bash
# Start Jupyter Notebook
jupyter notebook notebooks/eQTL_tutorial.ipynb

# Or use JupyterLab
jupyter lab
```

---

## 🎓 Quick Start

### For Docker Users

```bash
# One-command setup
docker run -d -p 8887:8887 --name eqtl_tutorial vladimirkovacevic/eqtl_intro:latest

# Access at http://localhost:8887
```

⚠️ **Note**: Without data mounting, you'll need to adjust file paths in the notebook.

### For Local Users

```bash
# Activate environment
source venv/bin/activate

# Launch Jupyter
jupyter notebook notebooks/eQTL_tutorial.ipynb
```

---

## 💻 Usage

### Running the Tutorial

1. **Open the notebook**: Navigate to `notebooks/eQTL_tutorial.ipynb`
2. **Run cells sequentially**: Press `Shift + Enter` to execute each cell
3. **Read explanations**: Each section includes detailed markdown explanations
4. **Explore visualizations**: Examine plots to understand the data
5. **Modify and experiment**: Try changing parameters to see different results

### Tips for Learning

- **Start from the beginning**: Sections build on each other
- **Take your time**: Read the explanations carefully
- **Run all cells**: Ensure each cell executes successfully
- **Check outputs**: Compare your results with expected outputs
- **Experiment**: Modify code to deepen understanding

### Common Tasks

#### Change data paths
If your data is located elsewhere, update the paths in the notebook:

```python
geno = pd.read_csv('/path/to/your/sim_genotypes.tab', sep='\t')
expr = pd.read_csv('/path/to/your/sim_expression1.tab', sep='\t')
```

#### Save results
Export analysis results:

```python
# Save eQTL results
results_final.to_csv('eqtl_results.tsv', sep='\t', index=False)

# Save figures
plt.savefig('manhattan_plot.pdf', dpi=300, bbox_inches='tight')
```

#### Use your own data
Replace the simulated data with your own:

```python
# Load your genotype data (samples × SNPs)
geno = pd.read_csv('your_genotypes.tab', sep='\t')

# Load your expression data (samples × genes)
expr = pd.read_csv('your_expression.tab', sep='\t')
```

---

## 📊 Data

### Simulated Dataset

The tutorial uses simulated data with:
- **300 individuals**
- **10 SNPs** (Single Nucleotide Polymorphisms)
- **10 genes** with expression measurements
- **Covariates** for confounding analysis

### Data Format

#### Genotype Data (`sim_genotypes.tab`)
```
sample    snp_1  snp_2  snp_3  ...
sample_1    0      1      2    ...
sample_2    1      0      1    ...
...
```
- Rows: Individuals
- Columns: Sample ID + SNPs
- Values: 0, 1, 2 (copies of alternative allele)

#### Expression Data (`sim_expression1.tab`)
```
sample    gene_1  gene_2  gene_3  ...
sample_1   7.45    6.82    8.13  ...
sample_2   6.91    7.23    7.56  ...
...
```
- Rows: Individuals
- Columns: Sample ID + Genes
- Values: log2-transformed normalized expression

#### Covariate Data (`sim_covariates.tab`)
- Technical and biological confounders
- Used for adjusted eQTL analysis

---

## 🗺️ Learning Path

### Beginner (2-3 hours)
1. Complete Section 1: Basic associations
2. Understand MAF and linear regression
3. Interpret p-values and effect sizes

### Intermediate (3-4 hours)
1. Work through Section 2: Confounders
2. Learn about covariate adjustment
3. Complete Section 3: PCA

### Advanced (4-5 hours)
1. Master Section 4: Genome-wide analysis
2. Apply multiple testing corrections
3. Create publication-quality figures
4. Analyze your own data

---

## 📁 Repository Structure

```
eqtl-tutorial-python/
├── README.md                          # This file
├── LICENSE                            # GNU GPL v3 license
├── DOCKER_QUICKSTART.md              # Quick Docker guide
├── notebooks/
│   └── eQTL_tutorial.ipynb           # Main tutorial notebook
├── data/
│   └── simulated/                    # Simulated data files
│       ├── sim_genotypes.tab
│       ├── sim_expression1.tab
│       ├── sim_expression2.tab
│       └── sim_covariates.tab
├── docker-build/
│   ├── Dockerfile                    # Docker image definition
│   ├── requirements.txt              # Python dependencies
│   └── README.md                     # Docker documentation
└── results/                          # Output files (generated)
```

---

## 🤝 Contributing

Contributions are welcome! Here's how you can help:

### Reporting Issues
- Use the [GitHub issue tracker](https://github.com/yourusername/eqtl-tutorial-python/issues)
- Describe the problem clearly
- Include error messages and system info
- Provide steps to reproduce

### Suggesting Enhancements
- Open an issue with your suggestion
- Explain the use case
- Discuss implementation approaches

### Pull Requests
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

### Development Setup
```bash
git clone https://github.com/yourusername/eqtl-tutorial-python.git
cd eqtl-tutorial-python
python3 -m venv venv
source venv/bin/activate
pip install -r docker-build/requirements.txt
jupyter notebook
```

---

## 📖 Citation

If you use this tutorial in your research or teaching, please cite:

### This Tutorial
```bibtex
@software{eqtl_tutorial_python,
  author = {Kovacevic, Vladimir},
  title = {eQTL Analysis Tutorial - Python Edition},
  year = {2025},
  url = {https://github.com/yourusername/eqtl-tutorial-python}
}
```

### Original R Tutorial
```bibtex
@misc{knight_eqtl_tutorial,
  author = {Knight Lab, University of Oxford},
  title = {Introduction to eQTL Analysis},
  year = {2016},
  url = {https://jknightlab.github.io/eqtl-intro/},
  note = {Python conversion based on exercises: https://jknightlab.github.io/eqtl-intro/exercises/exercises_and_solutions.html}
}
```

### Key References

**GTEx Consortium** (2020). The GTEx Consortium atlas of genetic regulatory effects across human tissues. *Science*, 369(6509), 1318-1330.

**Shabalin, A.A.** (2012). Matrix eQTL: ultra fast eQTL analysis via large matrix operations. *Bioinformatics*, 28(10), 1353-1358.

**Ongen, H. et al.** (2016). Fast and efficient QTL mapper for thousands of molecular phenotypes. *Bioinformatics*, 32(10), 1479-1485.

---

## 📜 License

This project is licensed under the **GNU General Public License v3.0** - see the [LICENSE](LICENSE) file for details.

### What this means:
- ✅ You can use this tutorial freely for education and research
- ✅ You can modify and distribute modified versions
- ✅ You must make source code available if you distribute
- ✅ Modified versions must also be GPL-licensed
- ❌ No warranty is provided

### Third-Party Components

This tutorial uses several open-source libraries:
- **Python**: PSF License
- **pandas**: BSD 3-Clause License
- **NumPy**: BSD 3-Clause License
- **matplotlib**: PSF-based License
- **seaborn**: BSD 3-Clause License
- **scikit-learn**: BSD 3-Clause License
- **statsmodels**: BSD 3-Clause License
- **Jupyter**: BSD 3-Clause License

All dependencies maintain their original licenses.

---

## 🙏 Acknowledgments

### Original Tutorial
This work is based on the excellent [eQTL Introduction tutorial](https://jknightlab.github.io/eqtl-intro/) created by:
- **Julian Knight's Lab**, University of Oxford
- **Peter Humburg**, primary developer

The Python conversion specifically follows the structure and exercises from:
- [Exercises and Solutions](https://jknightlab.github.io/eqtl-intro/exercises/exercises_and_solutions.html)

### Python Conversion
- **Vladimir Kovacevic** - Python conversion, Docker packaging, additional documentation

### Resources
- **GTEx Consortium** - For advancing eQTL research
- **Bioconductor Community** - For R/Bioconductor packages that inspired Python equivalents
- **Python Scientific Community** - For excellent statistical libraries

### Inspiration
This tutorial was created to make eQTL analysis more accessible to Python users and to provide a modern, containerized learning environment.

---

## 📞 Support & Contact

### Getting Help
- 📚 **Documentation**: Read the [full README](docker-build/README.md)
- 💬 **Issues**: Open an [issue on GitHub](https://github.com/yourusername/eqtl-tutorial-python/issues)
- 📧 **Email**: contact@example.com (for private inquiries)

### Community
- **Discussions**: Use GitHub Discussions for questions
- **Bug Reports**: Use GitHub Issues
- **Feature Requests**: Use GitHub Issues with `enhancement` label

### Related Resources
- [Original R Tutorial](https://jknightlab.github.io/eqtl-intro/)
- [Original R Exercises (basis for this tutorial)](https://jknightlab.github.io/eqtl-intro/exercises/exercises_and_solutions.html)
- [GTEx Portal](https://gtexportal.org)
- [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/)
- [statsmodels Documentation](https://www.statsmodels.org/)
- [pandas Documentation](https://pandas.pydata.org/)

---

## 🌟 Star History

If you find this tutorial helpful, please consider giving it a star ⭐ on GitHub!

---

## 🔄 Updates & Changelog

### Version 1.0.0 (2025-12-11)
- Initial release
- Complete Python conversion of R tutorial
- Docker support with Jupyter Notebook
- Comprehensive documentation
- All 4 sections implemented
- Professional visualizations

### Planned Features
- [ ] Integration with real GTEx data
- [ ] Additional visualization options
- [ ] Support for trans-eQTL analysis
- [ ] Batch processing scripts
- [ ] Video tutorials
- [ ] Interactive widgets

---

**Made with ❤️ for the genomics community**

**Happy Learning! 🎓📊🧬**
