from setuptools import setup

setup(name='recycler',
      version='0.62',
      description='Recycler: an algorithm for detecting plasmids from de novo assembly graphs',
      url='https://github.com/Shamir-Lab/Recycler',
      author='Roye Rozov',
      author_email = '',
      license='BSD-3-Clause',
      scripts = ['bin/recycle.py', 'bin/make_fasta_from_fastg.py', 'bin/get_simple_cycs.py'],
      packages = ['recyclelib'],
      requires=['python (<3.0)'],
      install_requires=[
        'networkx',
        'pysam',
        'nose',
        'numpy']
      )
