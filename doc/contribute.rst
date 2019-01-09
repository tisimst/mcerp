Contribute
==========

This is an open source project on Github: https://github.com/tisimst/mcerp

Contributions are very welcome, feel free to file an issue with a suggestion,
or send a pull request with an improvement any time.

We don't have a mailing list of other means of communication, so if you have
any question, please just open a Github issue and ask there.

Develop
-------

To work on ``mcerp``, first get the latest version::

    git clone https://github.com/tisimst/mcerp.git
    cd mcerp

The files to edit are here:

- Code in ``mcerp``
- Tests in ``mcerp/tests``
- Documentation in ``doc``

Install
-------

To hack on ``mcerp``, you need to have a development environment
with all packages and tools installed.

If you're using ``conda``, it's easy to create a development environment::

    conda env create -f environment.yml
    conda activate mcerp

With the conda environment active, run this command::

    pip install -e .

This installs ``mcerp`` in editable mode, meaning a pointer
is put in your site-packages to the current source folder, so
that after editing the code you only have to re-start python
and re-run to get this new version, and not run an install command again.

Tests
-----

Run all tests::

    pytest -v

Run tests and check coverage::

    pytest -v --cov=mcerp --cov-report=html
    open htmlcov/index.html

Code style
----------

We use the `black`_ code style. To apply it in-place to all files::

    black mcerp readme_code.py setup.py

Docs
----

To build the docs::

    cd docs
    make clean && make html
    open _build/html/index.html

Then for any other tasks go back to the top level of the package::

    cd ..

Release
-------

To make a release for this package, follow the following steps

#. check that the tests and docs build are OK
#. check via ``git tag`` or on Github or PyPI what the next version should be
#. ``git clean -fdx``
#. ``git tag vX.Y`` (substitute actual version number here and in the following steps)
#. ``python setup.py build sdist``
#. check the package in ``dist`` (should automate somehow)
#.  ``twine upload dist/*``
#. ``git push --tags``

.. _black: https://black.readthedocs.io
