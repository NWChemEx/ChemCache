Building the Documentation
==========================

Documentation for ChemCache has three parts:

- The Sphinx front-end which contains things like tutorials
- The Sphinx autodoc-generated Python documentation
- The Doxygen-generated C++ API documentation

**Note:** If the C++ or Python documentation does not exist, building the Sphinx 
front-end will still work, but the links to the C++ and Python API won't work.

Building the Sphinx Front-End
-----------------------------

The Sphinx front-end of the documentation requires the following Python packages
(all can be installed with `pip`):

- sphinx_rtd_theme (The Read-The-Docs theme for sphinx)
- sphinx (The thing that makes the documenation)

To build the Sphinx front-end of the documentation it is recommended that you
use Python's virtual environments. Assuming you want to call your virtual
environment `chemcache_venv` the steps are (run from this directory):

```.bash
# First two lines are optional to generate virtual environment
$ python3 -m venv chemcache_venv
$ . chemcache_venv/bin/activate

# Install dependencies and generate documentation
(venv) $ pip install -r requirements.txt
(venv) $ make html
```

You can then point your web browser at `build/html/index.html` to view the
documentation.

Building the Python API Documentation
-------------------------------------

The Python documentation is currently built manually. It can be built by 
running (from the `docs/` directory, assuming the virtual environment above is 
active):

```bash
$ sphinx-apidoc --force --separate --no-toc \
    -o ../docs/source/python_api/reference_data ../reference_data
```

Building the CXX API Documentation
----------------------------------

The C++ API documentation is built as part of the build process and then moved
to `docs/source/chemcache_cxx_api`. When the Sphinx front-end documentation is 
built it will link to the C++ documentation. 