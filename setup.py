from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()


setup(name='EpigeneticPacemaker',
      version='2.0.0',
      description='Epigenetic State Modeling Utility',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/NuttyLogic/EpigeneticPacemaker',
      project_urls={'Documentation': 'https://epigeneticpacemaker.readthedocs.io'},
      author='Colin P. Farrell, Sagi Snir',
      author_email='colinpfarrell@gmail.com',
      packages=['EpigeneticPacemaker', 'EpigeneticPacemaker.ExampleData', ],
      requires=['numpy', 'tqdm'],
      install_requires=['numpy>=1.20.3', 'tqdm>=4.61.0',
                        'joblib>=1.0.1'],
      python_requires='>=3.6',
      test_suite='tests',
      license='MIT',
      zip_safe=False,
      include_package_data=True,
      )
