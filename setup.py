from setuptools import setup, find_packages

setup(
    name="Framework_package",  
    version="0.0.1",  
    author="HsinYu-Hsu",  
    author_email="shinyu229@gmail.com",  
    description="A single-sample network-based framework",  
    long_description=open("README.md").read(),  
    long_description_content_type="text/markdown",  
    url="https://github.com/HsinYu-Hsu/Framework_package", 
    packages=find_packages(),  
    include_package_data=True, 
    install_requires=[ 
        "numpy",
        "pandas",
        "scipy",
        "statsmodels",
        "matplotlib",
        "argparse",
        "networkx"
    ],
    classifiers=[  
        "Programming Language :: Python :: 3.9.17",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",  # 支援的 Python 版本
)