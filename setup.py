from setuptools import setup

setup(
    name='reseqtb-uvp',
    version='2.5.1',
    description='Mycobacterium tuberculosis next generation sequence analysis',
    url='http://github.com/CPTR-ReSeqTB/UVP',
    author='Matthew Ezewudo',
    author_email='',
    license='MIT',
    packages=['uvp'],
    scripts=['bin/uvp'],
    zip_safe=False,
    install_requires = [
        "pyyaml",
    ],
)
