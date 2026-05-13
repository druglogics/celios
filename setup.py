from pathlib import Path

from setuptools import find_packages, setup


ROOT = Path(__file__).parent
README = ROOT / "README.md"

setup(
	name="celios",
	version="0.1.0",
	description="CELIOS: a DrugLogics pipeline module that processes CELl LIne OmicS (CELIOS) for calibration of Boolean models in the DrugLogics and TRAFIKK pipeline.",
	long_description=README.read_text(encoding="utf-8"),
	long_description_content_type="text/markdown",
	author='Viviam Solangeli Bermudez',
	author_email='viviamsb@ntnu.no',
	packages=find_packages(where="src"),
	package_dir={"": "src"},
	package_data={
		"celios.features": ["Model.csv"],
	},
	include_package_data=True,
	python_requires=">=3.8",
	install_requires=[
		"pandas>=1.0.0",
	],
	extras_require={
		"dev": [
			"pytest>=6.0.0",
		],
	},
	entry_points={
		"console_scripts": [
			"celios=celios.cli:main",
		]
	},
	project_urls={
		"Homepage": "https://github.com/druglogics/celios",
		"Source": "https://github.com/druglogics/celios",
		"Issues": "https://github.com/druglogics/celios/issues",
	},
)
