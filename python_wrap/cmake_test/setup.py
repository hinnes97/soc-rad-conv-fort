from skbuild import setup

setup(
    name="add_cython",
    version="0.1",
    description="a description",
    author='Hamish Innes',
    license="MIT",
    packages=['add_cython'],
    install_requires=['cython','numpy'],
)
