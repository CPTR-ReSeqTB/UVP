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
    scripts=[
        'bin/uvp',
        'bin/UVP',
        'scripts/coverage_estimator.py',
        'scripts/del_parse.py',
        'scripts/lineage_parser.py',
        'scripts/parse_annotation.py',
        'scripts/raw_vcf_parse.py',
        'scripts/resis_parser.py',
    ],
    zip_safe=False,
    install_requires = [
        "pyyaml",
    ],
)
