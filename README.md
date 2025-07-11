# Description

This repository contains the code required to reproduce analyses and figures from Taff et al. "Organismal robustness and resilience influence the ability to cope with cold snaps". Repository and code curated by Conor Taff (cct63@cornell.edu).

The data files required for these analyses are stored in separate repository here: XXX. In order to reproduce analyses, the data should first be downloaded and unzipped, then the folders should all be added to the same local folder as the code repository. The code is structured to access data based on relative file paths.

Code repository permanent archive: [![DOI](https://zenodo.org/badge/827961817.svg)](https://doi.org/10.5281/zenodo.15864802)

Data repository permanent archive: [![DOI](https://zenodo.org/badge/827961817.svg)](https://doi.org/10.5281/zenodo.15864830)

# Code structure

The code repository includes three scripts.

- **incubation_probe_processing**: This script processes the raw data from fake egg thermocouples used to assess incubation behavior. It outputs bout level data which is used in downstream analyses.

- **raw_rfid_processing**: This script processes the raw RFID data from RFID readers at each box. It outputs a single file with all reads and also summarizes feeding rate per hour. The summarized data is used in downstream analyses.

- **main_analysis_script**: This script performs all the analyses and generates all figures used in the manuscript. It is heavily annotated and organized into sections corresponding to the manuscript sections. It relies on the output of the raw data processing scripts above, but those intermediate data products have already been saved so that those scripts should not need to be run again in order to start using the main analysis script.

# Data structure

The data repository includes several folders that hold various kinds of raw or processed data.

- **hobo_csvs**: Contains the raw data files from temperature loggers at each nest before any processing. Each file is from one nest.

- **hobo_plots**: Contains diagnostic plots used to identify excluded sections from temperature loggers. Each plot is from one nest.

- **rfid_csvs**: Contains the raw data files from RFID readers at each box before any processing. Each file is from one nest.

- **cleaned_rfids**: Contains intermediate processed RFID data after cleaning and before summary. Each file is from one nest.

- **reference_files**: These files include the reference data that is added to raw RFID and thermocouple data, such as nest dates, band numbers, and sections to be excluded. These are all integrated to produce the summary information.

- **output_plots**: Contains plots produced by the various sections of analysis script used either for intermediate checks or for final figures.

- **saved_models**: Contains saved model objects for some sections of the analysis script. These can be loaded and inspected without having to re-run the model fitting.

- **processed_data_output**: Contains processed data that was produced from the raw data processing scripts. These are saved so that the full raw scripts do not need to be run in order to work with the data. Files are described in detailed annotation within the processing scripts.

- **raw_data**: This folder contains the main raw data used in the primary manuscript analysis. There are five data files that are described in more detail below.

# Raw data folder files

These are the primary files used in the analyses.

- **GFR_weather** and **ithaca_airport_weather** include the temperature records used to match hourly ambient temperature to all records from the long term database and in all behavioral analyses related to current ambient temperature.

- **excluded_treatments** is used to exclude records from the long term database in which manipulative treatments were applied.

- **data_by_capture** includes one row for each capture of adults or nestlings from the full database. Column names are self explanatory and include all pertinent information from that capture.

- **data_by_nest** includes one row for each nest from the full database. Column names are self explanatory and include all pertinent information from that nest. Nest records are joined to capture records by the nest_key.
