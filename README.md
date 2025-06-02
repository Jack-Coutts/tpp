# Thermal Proteome Profiling Script

The code is contained in the `thermal_proteome_profiling/` directory.

Ensure you create a `data/` directory (at the same level as `thermal_proteome_profiling/`) containing the data file(s).

Also, create empty `outputs/` directory also at the same level as `thermal_proteome_profiling/`.

Directory structure should look like:

```
|-- data/
|   |-- data.csv
|-- outputs/
|-- thermal_proteome_profiling/
|   |-- __init__.py
|   |-- main.py
|-- .gitignore
|-- poetry.lock
|-- pyproject.toml
|-- README.md
|-- run.sh
|-- test.sh
```



* Install poetry
* `poetry install`
* Run script with `./run.sh`
