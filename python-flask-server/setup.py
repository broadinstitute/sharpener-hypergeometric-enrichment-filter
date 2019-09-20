# coding: utf-8

import sys
from setuptools import setup, find_packages

NAME = "hypergeometric-enrichment-filter"
VERSION = "1.2.0"

# To install the library, run the following
#
# python setup.py install
#
# prerequisite: setuptools
# http://pypi.python.org/pypi/setuptools

REQUIRES = ["connexion"]

setup(
    name=NAME,
    version=VERSION,
    description="Hypergeometric enrichment filter",
    author_email="",
    url="",
    keywords=["Swagger", "Hypergeometric enrichment filter"],
    install_requires=REQUIRES,
    packages=find_packages(),
    package_data={'': ['swagger/swagger.yaml']},
    include_package_data=True,
    entry_points={
        'console_scripts': ['swagger_server=swagger_server.__main__:main']},
    long_description="""\
    Gene set filter based on hypergeometric enrichment test.
    """
)

