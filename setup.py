# -*- coding: utf-8 -*-
from distutils.core import setup
from setuptools import find_packages

setup(
    name='Sally_KDE',
    author='sallyda',
    version='0.0',
    url='https://github.com/sallyda/Sally_KDE',
    download_url='https://github.com/sallyda/Sally_KDE/tarball/v0.0',
    packages=find_packages(),
    license='MIT',
    description='',
    classifiers=[
        # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
    ],
    include_package_data=True,
    install_requires=[],
    setup_requires=['pytest-runner'],
)
