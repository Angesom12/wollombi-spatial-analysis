Wollombi Spatial Analysis
Overview
This project performs spatial analysis on the Wollombi area, including processing of Digital Elevation Models (DEM), channel lines, and gauging station data.

Project Structure
wollombi-spatial-analysis/
├── data/              # Storage for processed data
├── input_data/        # Raw data storage
├── outputs/           # Analysis outputs and visualizations
├── src/               # Source code
└── requirements.txt   # Python dependencies

Features

Processing of 10m and 30m resolution DEMs
Analysis of channel lines and gauging stations
Elevation comparison between different resolution DEMs
Spatial relationship analysis between channels and gauging stations

Installation

1.Clone the repository:

bash
git clone https://github.com/Angesom12/wollombi-spatial-analysis.git
cd wollombi-spatial-analysis
2.Install dependencies:
pip install -r requirements.txt

Data Requirements
The following data files are required:

wollombi_channelline.rar
wollombi_gauging.rar
wollombi_DEM_10m.rar
wollombi_DEM_30m.rar

Place these files in the input_data directory.
Usage
Run the spatial analysis script:

python src/spatial_analysis.py

Outputs
The script generates:

Visualizations of DEMs
Channel line and gauging station plots
Elevation analysis results
GeoPackage files with processed data
