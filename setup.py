from setuptools import setup, find_packages

import os 

setup(
    name="lmfilter",
    version="1.0",
    packages="lmfilter",
    author="James Boocock",
    author_email="james.boocock@otago.ac.nz",
    description="JoinMap Locus File Filter",
    license="Mit",
    zip_safe=False,
     entry_points={
        'console_scripts': [
            'lmfilter = jmfilter.locus_filter:main',
        ]
        },
    url="github.com/theboocock/jmfilter",
    use_2to3=True,
)
