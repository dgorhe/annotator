#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='annotator',
    version='0.2.0',
    url='github.com/byee4/annotator',
    license='',
    author='brianyee',
    author_email='',
    description='Annotate BED files with genic information given a priority list',
    packages=['annotator'],
    install_requires=[
        'tqdm',
        'gffutils',
    ],
    extras_require={
        'regions': ['pybedtools>=0.9', 'pandas>=1.0'],
        'utils': ['pandas>=1.0'],
        'dev': [
            'pytest>=7.0',
            'autopep8>=2.0',
            'flake8>=7.0',
            'flake8-pyproject>=1.2.0',
            'mkdocs>=1.6',
            'mkdocs-material>=9.5',
            'mkdocstrings[python]>=0.27',
            'mkdocs-llmstxt>=0.5.0',
        ],
    },
    package_dir={
        'annotator': 'annotator',
    },
    entry_points={
        'console_scripts': [
            'annotator = annotator.annotator:main',
            'create_region_bedfiles = annotator.create_region_bedfiles:main',
            'miRNA_name2id = annotator.miRNA_name2id:main',
            'gene_name2id = annotator.gene_name2id:main',
            'build_gffutils_db = annotator.build_gffutils_db:main',
            'get_region_lengths = annotator.get_region_lengths:main',
            'create_as_structure = annotator.create_AS_STRUCTURE:main',
        ]
    }
)
