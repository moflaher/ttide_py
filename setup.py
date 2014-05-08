from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()

setup(name='ttide',
      version='0.1_exp',
      description='Python distribution of the MatLab package TTide.',
      long_description=readme(),
      url='https://github.com/moflaher/ttide_py',
      author='Mitchell O\'Flaherty-Sproul',
      author_email='073208o@acadiau.ca',
      license='MIT',
      packages=['ttide'],
      zip_safe=False)
