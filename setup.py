from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()


setup(name='EpigeneticPacemaker',
      version='0.0.3',
      description='Epigenetic State Modeling Utility',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/NuttyLogic/EpigeneticPacemaker',
      project_urls={'Documentation': 'https://epigeneticpacemaker.readthedocs.io'},
      author='Colin P. Farrell, Sagi Snir',
      author_email='colinpfarrell@gmail.com',
      packages=['EpigeneticPacemaker', 'EpigeneticPacemaker.ExampleData', ],
      requires=['scipy', 'numpy', 'tqdm'],
      install_requires=['numpy>=1.16.3', 'tqdm>=4.31.1', 'scipy>=1.3.0'],
      python_requires='>=3.6',
      test_suite='tests',
      license='MIT',
      zip_safe=False,
      include_package_data=True,
      )
