## Pre-requisites

- [x] `miniconda` for Macbook M1/M2/M3 (`arm64`)
- [x] `homebrew`, `git`, `XCode`
- [x] Github account and authentication tokens
- [x] VS Code

## Steps

1.  Open a Terminal application and `git clone` this repository inside a folder, e.g., `~/summer2025/`:

    ```sh
    mkdir -p ~/summer2025;
    git clone [https://github.com/Summer2025-SPAC/SPAC_Shiny](https://github.com/Summer2025-SPAC/SPAC_Shiny);
    cd ~/summer2025/SPAC_Shiny
    ```

2.  `git checkout` your development branch, e.g., `gh_iss2_rk`:

    ```sh
    git checkout branch_name;
    ```

3.  Open the folder in VS Code:

    ```sh
    code ./
    ```

4.  Inside the VS Code `terminal`, create a `conda` environment using the `environment.yml` file. This file specifies all necessary packages and ensures compatibility.

    
    Ensure you are in the SPAC_Shiny directory

    ```sh
    cd ~/summer2025/SPAC_Shiny;
    ```

    Make edits these edits to the `environment.yml` file:

        Under dependencies, add these lines:

        ```sh
        - tables>=3.8.0
        - c-blosc2
        - libtiff
        ```

        Under pip, remove this line:

        ```sh
        - tables==3.8.0
        ```
        
        Installing tables and c-blosc2 with conda guarantees a smooth setup on MacOS by using pre-built, compatible binaries, while avoiding build errors that occur with pip installations.

        By adding libtiff directly to your conda dependencies, you instruct conda to download this library and place it correctly within your shiny environment's path. This ensures that when Pillow is imported, its dependencies are available and properly linked, resolving the error.

        # Create the environment from environment.yml
        # This will install Python 3.9 and all specified packages.
        conda env create -f environment.yml

        # Activate the newly created environment using the environment name

        ```sh
        conda activate shiny
        ```

5.  Assuming all installation works, run the shiny app in the terminal.

    ```sh
    shiny run app.py
    ```

6. If you encounter ImportErrors with libtiff, carry out these commands:

    ```sh
    conda activate shiny
    pip uninstall pillow
    conda install -c conda-forge --force-reinstall pillow libtiff
    ```
    Then run the shiny app again in the terminal:

    ```sh
    shiny run app.py
    ```