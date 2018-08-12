#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

__version__ = '0.11'


requirements = [
    'numpy',
    'scipy',
    'setuptools-scm'
]

test_requirements = [
    # TODO: put package test requirements here
]

setup_requirements = [
    'setuptools_scm'
]


setup(
    name = 'focus',
    use_scm_version = True,
    description = 'FOCUS: An Agile Profiler for Metagenomic Data',
    author = 'Genivaldo Silva',
    author_email = 'genivaldo.gueiros@gmail.com',
    url = 'https://github.com/metageni/focus',
    packages = [
        'focus_app',
    ],
    package_dir = {'focus_app': 'focus_app'},
    include_package_data = True,
    install_requires = requirements,
    setup_requires = setup_requirements,
    zip_safe = False,
    keywords = 'focus_app',
    classifiers = [
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Bioinformatics',
    ],
    test_suite = 'tests',
    tests_require = test_requirements,
    entry_points={
        'console_scripts': [
            'focus = focus_app.focus:main',
        ]
    },
)
