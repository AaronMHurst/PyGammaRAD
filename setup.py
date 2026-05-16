import setuptools

requirements = ["numpy", "pandas"]

setuptools.setup(
    name="PyGammaRAD",
    version="0.2.1",
    url="https://github.com/AaronMHurst/PyGammaRAD",
    author="Aaron M. Hurst",
    author_email="amhurst@berkeley.edu",
    description="Calculations of angular distribution coefficients in aligned and partially-aligned nuclei and general vector-coupling calculator for quantum mechanical applications based on the theory of angular momentum",
    long_description=open('README.md').read(),
    license_files=('LICENSE'),
    #packages=setuptools.find_packages(),
    packages=setuptools.find_namespace_packages(),
    #install_requires=[],
    install_requires=requirements,
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    include_package_data=True,
    package_data={'': ['data/*.json', 'data/*.csv']},
)
