Installation
===============

This section details the steps to installing Cedar and its dependencies.

.. note::
    If you wish to run any models that utilize OpenMC, then:

    - Cedar must be executed on Linux (Windows users can use WSL).
    - OpenMC must be included in the Conda installation.

(Optional) Windows Users
------------------------
Any Cedar model that utilizes OpenMC can only be used if you install and run
Cedar on Linux. If you're a user on a Windows Operating System, the easiest way
for you to use Linux is to take advantage of Windows Subsystem for Linux (WSL).

Check here for instructions:
`How to Install Windows Subsystem for Linux <https://learn.microsoft.com/en-us/windows/wsl/install>`_.

Clone Cedar from GitHub
-----------------------
From the Ubuntu command line, use these commands:

.. code-block:: bash

    git clone https://github.com/jcstonehill/cedar.git

Create Conda Environment
------------------------
It's recommended that Cedar users use Conda and Miniforge to manage their Python
environment. 

The following commands will download Miniforge and initialize Conda, then exit
the terminal. Upon restarting the terminal, the Conda environment prefix should
show before the command line.

.. code-block:: bash

    curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
    
    bash Miniforge3-Linux-x86_64.sh -b -p ~/miniforge
    
    export PATH=$HOME/miniforge/bin:$PATH
    
    conda init --all
    exit

Using the ``exit`` command will exit your terminal. Once you re-open your terminal, ``conda`` will be active.

Then, update the Conda base envionment and add the conda-forge channel.

.. code-block:: bash

    conda update --all --yes
    conda config --add channels conda-forge

Finally, create a new Conda environment and using the included environment.yml
file.

.. code-block:: bash

    cd cedar
    conda env create -f environment.yml

Activate Environment
--------------------

.. code-block:: bash

    conda activate cedar

Install Python API
------------------
Navigate to the root Cedar directory (where ``pyproject.toml`` is located) and
install the API using pip. Make sure that the Cedar Conda environment is active
before installing the API.

.. code-block:: bash

    python -m pip install .