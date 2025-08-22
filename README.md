# TPP Plotter - Thermal Proteome Profiling Analysis Tool

## Download Ready-to-Use Application

**For End Users (No Python Required):**

**Windows:** Download `TPP-Plotter-Windows.zip` from [Releases](../../releases)
- Extract the zip file
- Double-click `TPP Plotter.exe` to run

**macOS:** Download `TPP-Plotter-macOS.zip` from [Releases](../../releases)
- Extract the zip file
- Drag `TPP Plotter.app` to your Applications folder
- Double-click to run (you may need to right-click → Open the first time)

Both versions include all dependencies and provide a complete GUI experience.

---

## Overview

TPP Plotter is a desktop application for analyzing thermal proteome profiling (TPP) data. TPP is a systems-wide method to study protein thermal stability and drug-protein interactions by analyzing protein aggregation at different temperatures and compound concentrations.

## Application Features

### Data Input
- **Data File**: TSV/CSV files containing proteomics data with thermal stability measurements
- **Output Folder**: Directory where generated plots and analysis results will be saved
- **Gene List File**: Optional text file containing specific genes to analyze (one gene per line)

### Plot Configuration

#### Plot Modes
- **Plot One Line**: Analyze a single compound variant
  - **Compound Name Field**: When this mode is selected, you MUST specify which compound variant to plot
  - **Purpose**: The compound name field becomes editable and required
  - **Example**: If your data contains variants "DMSO", "Compound_A", "Compound_B", enter "Compound_A" to plot only that variant

- **Plot Two Lines**: Compare two compound variants automatically
  - **Compound Name Field**: Becomes disabled (grayed out) as both variants are plotted automatically
  - **Purpose**: Automatically detects and plots both available variants for comparison

> **Important**: The compound name field is only editable in "Plot One Line" mode because that's when you need to specify which single variant to analyze. In "Plot Two Lines" mode, the application automatically uses all available variants.

### Analysis Modes
- **Plot Single Gene**: Analyze one specific gene (enter gene symbol)
- **Plot from Gene List**: Batch analysis of genes from uploaded file
- **Plot All Proteins**: Comprehensive analysis of entire dataset

### Advanced Filtering Options

The application includes filtering capabilities for data quality control and biological insight extraction:

#### 1. No Filtering
- Processes all data without quality filters
- Useful for initial data exploration

#### 2. Area Between Curves Filtering
- **Purpose**: Identifies proteins with significant thermal stability differences
- **Method**: Calculates area between thermal curves using trapezoidal integration
- **Algorithm**: Computes |control(T) - treatment(T)| at each temperature and integrates
- **Parameters**:
  - **Area Threshold**: Minimum area for significance (typical: 1.0-10.0)
  - **Reference Variant**: Which compound to use as control
  - **Require Flat Control**: Ensures control curve is stable
  - **Flatness Threshold**: Maximum variability for "flat" designation
- **Use Cases**: Drug target identification, hit compound screening

#### 3. Filter for Stable Curves (Quality Control)
- **Purpose**: Select proteins with stable thermal profiles
- **Method**: Simple standard deviation threshold
- **Algorithm**: σ < flatness_threshold
- **Parameters**:
  - **Flatness Threshold**: Maximum standard deviation (typical: 0.35)
  - **Reference Variant**: Which variant to analyze
- **Use Cases**: Housekeeping protein identification, quality control

#### 4. Filter Against Stable Curves (Dose-Response)
- **Purpose**: Select proteins showing thermal response variation
- **Method**: Simple standard deviation threshold (opposite of flatness filter)
- **Algorithm**: σ ≥ flatness_threshold
- **Parameters**:
  - **Variability Threshold**: Minimum standard deviation (typical: 0.35)
  - **Reference Variant**: Which variant to analyze
- **Use Cases**: Dose-response analysis, drug target validation

### Execution Options
- **Add Error Bars**: Include standard error bars in plots (if available in data)
- **Save Detailed Filter Metrics**: Export filtering statistics to CSV files
- **Real-time Progress**: Live status updates during processing
- **Stop Functionality**: Safely interrupt long-running analyses

## Expected Data Format

### Input TSV/CSV Structure
```
Genes	Comparison (group1/group2)	AVG Log2 Ratio	Standard Error	# Unique Total Peptides
GAPDH	CellLine_R_1uM conc	0.12	0.05	5
GAPDH	CellLine_R_3uM conc	0.25	0.08	5
GAPDH	CellLine_S_1uM conc	0.15	0.06	4
GAPDH	CellLine_S_3uM conc	0.45	0.12	4
ACTB	CellLine_R_1uM conc	-0.05	0.03	8
...
```

### Required Columns
- **Genes**: Protein/gene identifiers
- **Comparison (group1/group2)**: Format: `CellLine_Variant_Concentration`
- **AVG Log2 Ratio**: Thermal stability measurements
- **# Unique Total Peptides**: Quality metric (minimum 2 required)
- **Standard Error**: For error bars (optional)

### Data Quality Requirements
- At least 2 concentration points
- Valid numeric values for thermal measurements
- Consistent variant naming across concentrations

## Usage Workflow

### Basic Analysis
1. **Load Data**: Select your TSV/CSV data file
2. **Set Output**: Choose destination folder for results
3. **Choose Plot Mode**:
   - Single line: Enter specific compound name
   - Two lines: Automatic comparison
4. **Select Analysis Scope**: Single gene, gene list, or all proteins
5. **Configure Filtering**: Choose appropriate quality control method
6. **Execute**: Click "Run" and monitor progress

### Advanced Analysis
1. **Quality Control**: Start with "Filter for Stable Curves" to identify reliable data
2. **Target Discovery**: Use "Area Between Curves" with area threshold 2.0-5.0
3. **Dose-Response**: Apply "Filter Against Stable Curves" for responsive proteins
4. **Validation**: Review saved metrics and plots for biological relevance

## Output Files

### Generated Plots
- **Format**: PNG files (300 DPI)
- **Naming**: `{GeneName}_{CellLine}.png`
- **Features**:
  - Consistent color scheme and styling
  - Optional error bars
  - Clear axis labels and legends
  - Horizontal reference line at y=0

### Filtering Metrics (Optional)
- **Format**: CSV files saved to output directory
- **Content**:
  - Filter method and parameters used
  - Statistical measures
  - Processing time and data quality metrics
  - Pass/fail decisions with rationale
- **Files Generated**:
  - Single gene: `{gene_name}_filter_metrics.csv`
  - Gene list: `gene_list_filter_metrics.csv`
  - All genes: `all_genes_filter_metrics.csv`

#### Example CSV Output
```csv
gene_name,filter_method,passes_filter,processing_time_ms,max_area,threshold,control_std,control_mean,is_flat
MNT,area_between_curves,false,6.3,1.184,2.0,0.294,0.145,true
GAPDH,area_between_curves,true,5.8,3.456,2.0,0.121,0.234,true
ACTB,filter_for_flatness,true,4.2,0.0,0.35,0.089,0.156,true
HSP90,filter_against_flatness,true,3.1,0.0,0.35,0.512,0.089,false
```

#### Working with CSV Metrics
The CSV files can be easily opened in Excel, imported into R/Python for further analysis, or processed with command-line tools:

```bash
# View metrics in terminal
cat gene_list_filter_metrics.csv

# Count genes that passed filtering
awk -F',' '$3=="true" {count++} END {print count}' gene_list_filter_metrics.csv

# Filter for high-area genes
awk -F',' 'NR==1 || $5>3.0' all_genes_filter_metrics.csv
```

## Developer Installation and Setup

### Prerequisites
- Python 3.8 or higher
- [uv](https://docs.astral.sh/uv/getting-started/installation/) for fast dependency management

### Install uv
```bash
# Linux/macOS
curl -LsSf https://astral.sh/uv/install.sh | sh

# Windows
powershell -c "irm https://astral.sh/uv/install.ps1 | iex"

# Or with pip
pip install uv
```

### Clone and Setup
```bash
# Clone the repository
git clone https://github.com/Jack-Coutts/tpp.git
cd thermal_proteome_profiling

# Install dependencies and create virtual environment
uv sync
```

### Running the Application
```bash
# GUI Mode (Recommended)
./run.sh

# Or manually:
uv run python src/main.py
```

### Command Line Interface
```bash
# Single gene analysis
uv run python src/main.py -d data/file.tsv -o outputs/ -g GENE_NAME

# Gene list analysis
uv run python src/main.py -d data/file.tsv -o outputs/ -gl gene_list.txt

# Full dataset with filtering
uv run python src/main.py -d data/file.tsv -o outputs/ -f

# With error bars
uv run python src/main.py -d data/file.tsv -o outputs/ -e
```

## Creating Standalone Executables

For non-technical end users, you can create standalone executables that don't require Python installation.

### Prerequisites
- PyInstaller: `uv add pyinstaller --dev`

### Building Executables

#### For macOS
```bash
python build_mac.py
```
**Output:** `dist/TPP Plotter.app` - Users can drag this to their Applications folder and double-click to run.

#### For Windows
```bash
python build_windows.py
```
**Output:** `dist/TPP Plotter.exe` - Users can double-click this file to run the application.

### Distribution
- **macOS:** Share the entire `TPP Plotter.app` bundle
- **Windows:** Share the single `TPP Plotter.exe` file
- Both versions include all dependencies and don't show a console window

### Cleanup (Optional)
You can remove development files before distribution:
- `build/` directory (build artifacts)
- `*.spec` files (PyInstaller configuration)
- Keep `dist/` for the final executables
