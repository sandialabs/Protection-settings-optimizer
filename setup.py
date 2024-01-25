from setuptools import setup, find_packages

with open('README.md',encoding="utf8") as f:
        long_description = f.read()

reqs = open('requirements.txt').readlines()

setup(
        name='RSO_pack',
        version='1.0.0',
        description='Sandia Protection Settings Optimizer Algorithms',
        long_description_content_type='text/markdown',
        long_description=long_description,
        #author='',
        #author_email='',
        url='https://github.com/sandialabs/Protection-settings-optimizer',
        packages = find_packages(),
        include_package_data=True,
        setup_requires=reqs,
        install_requires=reqs,
        license = 'GPLv3'
)
