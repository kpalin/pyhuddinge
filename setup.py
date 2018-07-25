from setuptools import setup,Extension
import versioneer

requirements = [
    "numpy","scipy","pandas","six"
    # package requirements go here
]



setup(
    name='pyhuddinge',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Python interface for computing huddinge distance from python code",
    author="Kimmo Palin",
    author_email='kimmo.palin at helsinki.fi',
    url='https://github.com/kpalin/pyhuddinge',
    packages=['pyhuddinge'],
    #ext_modules = [moder_ext],
    install_requires=requirements,
    keywords='pyhuddinge',
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
    ]
)
