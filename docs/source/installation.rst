.. Copyright 2023 NWChemEx-Project
..
.. Licensed under the Apache License, Version 2.0 (the "License");
.. you may not use this file except in compliance with the License.
.. You may obtain a copy of the License at
..
.. http://www.apache.org/licenses/LICENSE-2.0
..
.. Unless required by applicable law or agreed to in writing, software
.. distributed under the License is distributed on an "AS IS" BASIS,
.. WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
.. See the License for the specific language governing permissions and
.. limitations under the License.

############
Installation
############

***************
Install Command
***************

ChemCache uses `CMaize <https://cmakepp.github.io/CMaize/index.html>`__ as
its build system, more detailed build instructions can be found
`here <https://cmakepp.github.io/CMaize/getting_started/building/index.html>`__.
ChemCache can be installed using the usual CMake commands:

.. code-block:: bash

   cmake -H. \
         -B<build_dir> \
         -DCMAKE_INSTALL_PREFIX:PATH=<where/to/install/libraries> \
         -DCMAKE_TOOLCHAIN_FILE:PATH=<path/to/toolchain.cmake>
   cmake --build build --target install --parallel 2

Here ``<build_dir>`` is the name of the build directory CMake should use (most
users just set this to ``-Bbuild``), ``<where/to/install/libraries>`` should
be set to where you want to install the dependencies ChemCache builds for you,
and ``<path/to/toolchain.cmake>`` should point to your ``toolchain.cmake`` file.
Of particular note, make sure that in your toolchain file you set
``NWX_MODULE_PATH`` to where you want ChemCache installed and you may want to
set both ``Python_EXECUTABLE`` and ``Python3_EXECUTABLE`` to the Python
interpreter from your virtual environment (with the environment activated
run ``which python3`` to get it's path).

*********************
Configuration Options
*********************

This is a list of configuration options recognized by ChemCache's build
system.

``BUILD_TESTING``.
   Off by default. Set to a truth-y value to enable testing.
``BUILD_DOCS``.
   Off by default. Set to a truth-y value to build the C++ API documentation.
   This variable is defined by the ``nwx_cxx_api_docs`` CMake module.
``ONLY_BUILD_DOCS``.
   Off by default. If ``BUILD_DOCS`` is set to a truth-y value and
   ``ONLY_BUILD_DOCS`` is also set to a truth-y value, then the configure
   process will skip all other aspects of the configuration aside from creating
   the ``chemcache_cxx_api`` target. This variable is defined by the
   ``nwx_cxx_api_docs`` CMake module.
``BUILD_PYBIND11_PYBINDINGS``.
  On by default. When enabled the optional Python API is built.
``ENABLE_EXPERIMENTAL_FEATURES``.
  Off by default. When set to a truth-y value classes and functions which are in
  a pre-release state will be built.

************
Dependencies
************

Required Dependencies
=====================

These are dependencies which must be pre-installed and can not be built by
ChemCache's build system.

CMake
-----

CMake is the basis of CMaize, and minimum version of 3.14 is required to
properly build ChemCache.


C++ Compiler
------------

ChemCache relies on the C++17 standard and should work with any C++17
compliant compiler (GCC 9.x or newer).

Optional Dependencies
=====================

These are dependencies that ChemCache's build system can not build; however,
they are only required if certain features are enabled.

Doxygen
-------

`Home Page <https://www.doxygen.nl/>`__

Used to generate the C++ API documentation. Only needed if ``BUILD_DOCS`` is
set to a truth-y value.

Python
------

Needed if ``BUILD_PYBIND11_PYTHONBINDINGS`` is enabled. You will need the
developer headers and libraries for Python.

Other Dependencies
==================

The dependencies in this section can be built by ChemCache's build system
when they are not located. Under normal circumstances users can ignore them.
They are listed here primarily for completeness.

SimDE
-----

URL: `<https://github.com/NWChemEx/SimDE>`__

ChemCache is build off of the Simulation Development Environment (SimDE). See
`here <https://nwchemex.github.io/SimDE/install.html#simde-dependencies>`__ for
the list of dependencies inherited from SimDE.


Catch2
------

URL: `<https://github.com/catchorg/Catch2>`__

Used for unit testing. Only needed if unit testing is enabled (controlled by
the CMake variable ``BUILD_TESTING``, which is ``OFF`` by default).
