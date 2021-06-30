# Development
The project is structured to have 
- reusable functions in S2S/
- scripts in the scripts/
- datasets in data/
- to install the environment: conda env create -f conda_environment.yml 
- to update the environment: conda env update --file conda_environment.yml --prune 

- Run '''pip install -e .''' inside the project root folder and with the relevant python environment active. Then script files will find the function hierarchy below S2S/

All file paths are stored in S2S/local_configuration.py. Change the paths there for local use. The file should contain a dictionary '''config = {'S2S_DIR': ..., ...} ''' with paths on local machine.
