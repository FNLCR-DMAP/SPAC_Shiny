# SPAC_Shiny

SPAC_Shiny is the The Shiny Interactive Realtime Dashboard for SPAC, designed to empower researchers with intuitive, dynamic, and customizable visualization tools for spatial single-cell datasets derived from cancerous tumors. Built using the Shiny platform, SPAC_SHINY streamlines the exploration and analysis of complex biological data through a user-friendly dashboard interface.

## Project Overview and Contributions

This project involved a collaborative effort between the Frederick National Laboratory for Cancer Research (FNL), Purdue University, and The Data Mine project. Our team worked closely with FNL throughout the project, coordinating efforts to significantly enhance the functionality and performance of the SPAC library.

## User Experience

A look at the redesigned SPAC_SHINY dashboard, highlighting new features and enhanced interactivity:
![SPAC_SHINY New UI](https://github.com/user-attachments/assets/e5b42da3-0af6-4439-bd27-4a0bba71d04d)

## Key Features & Enhancements

### Front-End Interface Expansion
We dramatically expanded the functionality of the front-end interface, adding nearly **50 new UI features**. These enhancements include:
*   **Customization of parameter options:** Allowing users to fine-tune visualization settings.
*   **UI elements for adjusting styling:** Providing control over the aesthetic elements of generated plots.
*   **Novel features for flexible subsetting and data wrangling:** Enabling more dynamic data exploration and display within the UI.

### Optimization of SPAC Visualizations
A major focus was on characterizing and optimizing the run times of SPAC visualizations. By refactoring existing functions to utilize optimal computational methods, we achieved substantial performance gains:
*   **Characterizing Visualization Run Times:** Thorough analysis was conducted to pinpoint performance bottlenecks.
*   **Refactoring for Efficiency:** Existing visualization functions were re-engineered for improved speed.

**Optimization Results:**
We observed significant improvements in the performance of key visualization functions:

| Visualization         | Dataset Size           | Original Implementation | New Implementation | Efficiency Increase |
| :-------------------- | :--------------------- | :---------------------- | :----------------- | :------------------ |
| Boxplot               | 10 Million cells with annotations | 20.91 seconds           | 2.28 seconds       | **917%**            |
| Histogram             | 10 Million cells with annotations | 1.21 seconds            | 0.17 seconds       | **712%**            |

These optimizations resulted in a nearly **10x speedup** for existing boxplot and histogram visualization functions.

### Implementation of Clustering Algorithms
To enhance the analytical capabilities of SPAC, we integrated new clustering algorithms into the Python package:
*   **K-Nearest Neighbors**
*   **K-Means**

## Quickstart Tutorials

1. [Quickstart for Development on Windows WSL](./tutorials/Quickstart_WSL.md)

## Acknowledgements and Credits

This project was a success thanks to the invaluable collaboration and support from several key institutions and individuals:

*   **Frederick National Laboratory for Cancer Research (FNL):** For their guidance, expertise, and close coordination throughout the project. Special thanks to our FNL mentor, George Zaki.
*   **Purdue University:** For facilitating this research and development effort.
*   **The Data Mine Project at Purdue University:** For providing the framework and resources for this collaborative student-led initiative.
*   **Purdue Institute for Cancer Research:** For facilitating collaboration between Purdue University and the Frederick National Laboratory for Cancer Research
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
*   **Teaching Assistants (TAs) from Purdue's Data Mine:**
    *   Alex Liu
    *   Omar Eldaghar
    *   Thomas Sheeley
*   **Additional Support:** Kang Liu and Rui He, and the entire Data Mine staff.
