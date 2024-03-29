# coding: utf-8

import sys
from setuptools import setup, find_packages

NAME = "hypergeometric-enrichment-filter"
VERSION = "1.3.0"

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
    description="MSigDB hypergeometric enrichment",
    author_email="",
    url="",
    keywords=["Swagger", "MSigDB hypergeometric enrichment"],
    install_requires=REQUIRES,
    packages=find_packages(),
    package_data={'': ['swagger/swagger.yaml']},
    include_package_data=True,
    entry_points={
        'console_scripts': ['swagger_server=swagger_server.__main__:main']},
    long_description="""\
    Gene-list filter based on hypergeometric enrichment in MSigDB gene sets (http://software.broadinstitute.org/gsea/index.jsp).
    """
)

