# SPAC Technical Details & Project Report

This document provides comprehensive technical information about SPAC development, performance optimizations, and project contributions.

## Project Overview and Contributions

This project involved a collaborative effort between the Frederick National Laboratory for Cancer Research (FNL), Purdue University, and The Data Mine project. Our team worked closely with FNL throughout the project, coordinating efforts to significantly enhance the functionality and performance of the SPAC library.

## Architecture Overview

SPAC employs a four-tier architecture that includes:

1. **Modular Python-based analysis engine** - Built on scverse ecosystem (AnnData, SCANPY, Squidpy, SCIMAP)
2. **Seamless integration with HPC and GPU acceleration** - Automated workload distribution for large datasets
3. **Interactive browser interface** - No-code workflow configuration powered by enterprise platforms
4. **Real-time visualization layer** - Shiny for Python dashboards with dynamic exploration capabilities

This design empowers distinct user roles: data scientists can extend and customize analysis modules, while bench scientists can execute complete workflows and interactively explore results without coding.

## Performance Optimizations & Benchmarks

### Front-End Interface Expansion
We dramatically expanded the functionality of the front-end interface, adding nearly **50 new UI features**. These enhancements include:
*   **Customization of parameter options:** Allowing users to fine-tune visualization settings
*   **UI elements for adjusting styling:** Providing control over the aesthetic elements of generated plots
*   **Novel features for flexible subsetting and data wrangling:** Enabling more dynamic data exploration and display within the UI

### Visualization Performance Improvements
A major focus was on characterizing and optimizing the run times of SPAC visualizations. By refactoring existing functions to utilize optimal computational methods, we achieved substantial performance gains:

**Optimization Results:**
We observed significant improvements in the performance of key visualization functions:

| Visualization         | Dataset Size           | Original Implementation | New Implementation | Efficiency Increase |
| :-------------------- | :--------------------- | :---------------------- | :----------------- | :------------------ |
| Boxplot               | 10 Million cells with annotations | 20.91 seconds           | 2.28 seconds       | **917%**            |
| Histogram             | 10 Million cells with annotations | 1.21 seconds            | 0.17 seconds       | **712%**            |

These optimizations resulted in over  **5x speedup** for existing boxplot and histogram visualization functions.

## Technical Implementation Details

### Foundational Dependencies
SPAC builds on established open-source packages for single-cell and spatial proteomics analysis:
- **AnnData** - Standardized data structure for single-cell analysis
- **SCANPY** - Core single-cell analysis toolkit
- **Squidpy** - Spatial single-cell analysis methods
- **SCIMAP** - Spatial analysis and visualization tools

These tools are integral to the "scverse" ecosystem, ensuring compatibility with community standards and robust core functionalities.

## Acknowledgements and Credits

### Institutional Partners
*   **Frederick National Laboratory for Cancer Research (FNL)** - Guidance, expertise, and project coordination
*   **The Data Mine Project at Purdue University** - Framework and resources for collaborative student-led research

### Development Team
*   **Project Team Members:**
    *   Aileen Chow
    *   Ahmad Abdallah
    *   Arshnoor Randhawa
    *   Ella Delvecchio
    *   Jiazhen Li
    *   Liza Shchehlik
    *   Megan Lawson
    *   Mustapha Braimoh
    *   Ruhi Sharmin
    *   Sam Ying
    *   Shamita Yediapalli
    *   Anthony Cusimano
    *   Suriya Selvarajan
    *   Qianyue Wang
    *   Andree Kolliegbo
    *   Boqiang Zhang
    *   Sungmin Lee
    *   Saran Nagubandi
    *   Noah Lee
    *   Heaven Golladay-Watkins
    *   Arjun Chhabra
    *   Nila Kumar

*   **Teaching Assistants from Purdue's Data Mine:**
    *   Alex Liu
    *   Omar Eldaghar
    *   Thomas Sheeley
    *   Ramya Rajaram

*   **Additional Support:**
    *   Fang Liu and Rui He
    *   George Zaki (FNL mentor)
    *   The entire Data Mine staff


## Contact and Support

For technical questions and development inquiries:
- **Primary Contact**: george.zaki@nih.gov
- **GitHub Issues**: [SPAC_Shiny Issues](https://github.com/FNLCR-DMAP/SPAC_Shiny/issues)

---

*This document provides detailed technical information about SPAC development and optimization. For user-facing information and quick start instructions, see the main [README.md](README.md).*
