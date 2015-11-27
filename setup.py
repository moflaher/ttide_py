from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()

setup(name='ttide',
      version='0.3_exp',
      description='Python distribution of the MatLab package TTide.',
      long_description=readme(),
      url='https://github.com/moflaher/ttide_py',
      author='Mitchell O\'Flaherty-Sproul',
      author_email='073208o@acadiau.ca',
      license='MIT',
      packages=['ttide'],
      package_data={'ttide': ['data/*.nc']},
      zip_safe=False)
