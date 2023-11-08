from setuptools import setup

setup(
    name='haplot',
    version='0.0.4',
    description='A package for genomic data visualization',
    url='https://github.com/swu1019lab/haplot',
    author='XiaoDong Li',
    author_email='lxd1997xy@163.com',
    license='BSD 3-clause',
    packages=['haplot'],
    install_requires=[
        'numpy>=1.22.3',
        'pandas>=1.4.2',
        'matplotlib>=3.6.2',
        'scipy>=1.9.3',
        'geopandas>=0.13.2',
        'networkx>=3.1',
    ],
    zip_safe=False,
    python_requires='>=3.10'
)
