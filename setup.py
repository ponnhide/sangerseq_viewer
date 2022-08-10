#!/usr/bin/env python3
#
# Copyright (C) Hideto Mori


DESCRIPTION      = "Visualization tool for Sanger sequencing data"
DISTNAME         = 'sangerseq_viewer'
MAINTAINER       = 'Hideto Mori'
MAINTAINER_EMAIL = 'hidto7592@gmail.com'
URL              = 'https://github.com/ponnhide/patchworklib'
LICENSE          = 'GNU General Public License v3.0'
DOWNLOAD_URL     = 'https://github.com/ponnhide/patchworklib'
VERSION          = '0.1.3'
PYTHON_REQUIRES  = ">=3.7"

INSTALL_REQUIRES = [
    'biopython>=1.78',
    'matplotlib>=3.4.3',
    'patchworklib>=0.4.5',
    'python-queen>=1.1.0',
]

PACKAGES = [
    'sangerseq_viewer'
]

SCRIPTS=['bin/sangerseq_viewer']

CLASSIFIERS = [
    'Intended Audience :: Science/Research',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
]

with open('README.md', 'r', encoding='utf-8') as fp:
    readme = fp.read()
LONG_DESCRIPTION = readme
LONG_DESCRIPTION_CONTENT_TYPE = 'text/markdown'

if __name__ == "__main__":
    from distutils.core import setup
    try:
        # Setuptools only needed for building the package
        import setuptools  # noqa
    except ImportError:
        pass
    
    import sys
    if sys.version_info[:2] < (3, 7):
        raise RuntimeError("sangerseq_viewer requires python >= 3.7.")

    setup(
        name=DISTNAME,
        author=MAINTAINER,
        author_email=MAINTAINER_EMAIL,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        long_description_content_type=LONG_DESCRIPTION_CONTENT_TYPE,
        license=LICENSE,
        url=URL,
        version=VERSION,
        download_url=DOWNLOAD_URL,
        python_requires=PYTHON_REQUIRES,
        install_requires=INSTALL_REQUIRES,
        scripts=SCRIPTS,
        packages=PACKAGES,
        classifiers=CLASSIFIERS
    )
